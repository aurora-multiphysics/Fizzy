#include "FispactFluxInput.h"

class DummyFluxInput : public FispactFluxInput {

public:
  static InputParameters validParams() {

    InputParameters params = FispactFluxInput::validParams();
    return params;
  };

  DummyFluxInput(const InputParameters &params) : FispactFluxInput(params) {

    _n_bins = 10;
    _wall_loading = 1;
    _flux_energy_groups = std::vector<double>(_n_bins);

    for (int i = 0; i < _n_bins; i++) {

      _flux_energy_groups[i] = i + 1;
    }
  }

  virtual std::vector<double> getElemFlux(dof_id_type elem_id) {

    return std::vector<double>(10, 1);
  };

protected:
  size_t _flux_tally_id;

  size_t _energy_filter_id;

  std::string _statepoint_filename;
};
