#pragma once

#include <stdexcept>
#include <vector>

// Same indexing for locally owned storage and non-owning shared-memory views.
// The owner of an external buffer must outlive this view and never resize it.
class PhotonSpectra {
public:
  PhotonSpectra(size_t nt, size_t ne, size_t nb)
      : _nt(nt), _ne(ne), _nb(nb), _owned(nt * ne * nb, 0),
        _data(_owned.empty() ? &_empty : _owned.data()) {}

  PhotonSpectra(size_t nt, size_t ne, size_t nb, double *buffer)
      : _nt(nt), _ne(ne), _nb(nb), _data(buffer ? buffer : &_empty) {
    if (nt && ne && nb && !buffer)
      throw std::invalid_argument("Missing photon storage buffer");
  }

  PhotonSpectra(const PhotonSpectra &) = delete;
  PhotonSpectra &operator=(const PhotonSpectra &) = delete;

  double &at(size_t t, size_t e, size_t b) { return _data[idx(t, e, b)]; }
  const double &at(size_t t, size_t e, size_t b) const { return _data[idx(t, e, b)]; }

  double *spectrum_begin(size_t t, size_t e) { return _data + idx(t, e, 0); }
  double *spectrum_end(size_t t, size_t e) { return spectrum_begin(t, e) + _nb; }

  double *time_begin(size_t t) {
    if (t >= _nt)
      throw std::out_of_range("Photon inventory index out of range");
    return _data + t * _ne * _nb;
  }
  double *time_end(size_t t) { return time_begin(t) + _ne * _nb; }
  double *data_buffer() { return _data; }

private:
  size_t idx(size_t t, size_t e, size_t b) const {
    if (t >= _nt || e >= _ne || b >= _nb)
      throw std::out_of_range("Photon spectrum index out of range");
    return (t * _ne + e) * _nb + b;
  }

  size_t _nt, _ne, _nb;
  std::vector<double> _owned;
  double _empty = 0;
  double *_data;
};
