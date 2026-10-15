#pragma once

#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/shared_memory_object.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace bi = boost::interprocess;
using ShmemVecAllocator =
    bi::allocator<double, bi::managed_shared_memory::segment_manager>;
using BoostIpVector = bi::vector<double, ShmemVecAllocator>;
using ShmemIdAllocator =
    bi::allocator<uint64_t, bi::managed_shared_memory::segment_manager>;
using BoostIpIdVector = bi::vector<uint64_t, ShmemIdAllocator>;

// Version the object name as well as the payload: an old consumer must not
// interpret this layout as the former single-inventory PhotonSharingData.
class PhotonSharingData {
public:
  static constexpr const char *objectName() { return "photon_inventories_v2"; }

  explicit PhotonSharingData(bi::managed_shared_memory &segment)
      : _photon_fluxes(segment.get_segment_manager()),
        _elem_strength(segment.get_segment_manager()),
        _photon_bins(segment.get_segment_manager()),
        _global_element_ids(segment.get_segment_manager()) {}

  static size_t checkedProduct(size_t a, size_t b) {
    if (b && a > std::numeric_limits<size_t>::max() / b)
      throw std::overflow_error("Photon inventory storage size overflow");
    return a * b;
  }

  // Fixed arrays: never resize after views have been handed to the producer.
  void allocate(size_t inventories, size_t elements, size_t bins) {
    if (_num_inventories)
      throw std::logic_error("Photon inventory arrays are already allocated");
    if (!inventories || !bins)
      throw std::invalid_argument("Photon inventories and bins must be nonzero");
    _num_inventories = inventories;
    _num_local_elems = elements;
    _num_photon_bins = bins;
    const size_t rows = checkedProduct(inventories, elements);
    _photon_fluxes.resize(checkedProduct(rows, bins), 0.0);
    _elem_strength.resize(rows, 0.0);
    _global_element_ids.resize(elements);
  }

  size_t strengthIndex(size_t inventory, size_t element) const {
    if (inventory >= _num_inventories || element >= _num_local_elems)
      throw std::out_of_range("Photon inventory/element index out of range");
    return inventory * _num_local_elems + element;
  }

  const double *spectrum(size_t inventory, size_t element) const {
    return _photon_fluxes.data() +
           strengthIndex(inventory, element) * _num_photon_bins;
  }

  double elementStrength(size_t inventory, size_t element) const {
    return _elem_strength[strengthIndex(inventory, element)];
  }

  // Called only between transport solves. Arrays are immutable while consumers
  // sample. The generation is publication metadata, NOT an IPC synchronization
  // primitive for independently executing producer/consumer jobs.
  void publish(size_t inventory, double local_strength, double total_strength,
               bool uniform) {
    if (!_populated || inventory >= _num_inventories)
      throw std::logic_error("Cannot publish unavailable photon inventory");
    _selected_inventory = inventory;
    _local_domain_strength = local_strength;
    _total_domain_strength = total_strength;
    _uniform = uniform;
    ++_generation;
  }

  size_t _num_inventories = 0;
  size_t _num_photon_bins = 0;
  size_t _num_local_elems = 0;
  size_t _selected_inventory = 0;
  uint64_t _generation = 0;
  double _local_domain_strength = 0;
  double _total_domain_strength = 0;
  bool _populated = false;
  bool _uniform = false;

  // [inventory][local element][group], with group-integrated emission rates.
  BoostIpVector _photon_fluxes;
  BoostIpVector _elem_strength;
  BoostIpVector _photon_bins;
  // Explicitly in spectrum-row order, independent of global ID ordering.
  BoostIpIdVector _global_element_ids;
};

// Sample histogram data without allocating PDF/CDF arrays. Both random inputs
// must be in [0,1). total must be the sum of the nonnegative group rates.
inline double samplePhotonHistogram(const double *rates, const double *bounds,
                                    size_t bins, double total, double group_u,
                                    double energy_u) {
  if (!bins || !(total > 0) || !std::isfinite(total))
    throw std::invalid_argument("Cannot sample an empty/zero photon spectrum");
  const double target = group_u * total;
  double cumulative = 0;
  size_t last_positive = bins;
  for (size_t group = 0; group < bins; ++group) {
    if (rates[group] > 0) {
      last_positive = group;
      cumulative += rates[group];
      if (target < cumulative)
        return bounds[group] + energy_u * (bounds[group + 1] - bounds[group]);
    }
  }
  // Accumulation/product rounding can put the target at the upper endpoint.
  if (last_positive == bins)
    throw std::invalid_argument("Photon spectrum has no positive groups");
  return bounds[last_positive] +
         energy_u * (bounds[last_positive + 1] - bounds[last_positive]);
}
