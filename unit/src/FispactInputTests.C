#include "FispactProblemTest.h"
#include "gmock/gmock.h"

class FispactProblemInputTests : public FispactProblemTestUserObjects {
public:
  FispactProblemInputTests() : FispactProblemTestUserObjects() {
    _density = 100;
  }

  ~FispactProblemInputTests() = default;

protected:
  double _rtol, _atol, _density, _wall_loading, _atoms_threshold;

  std::vector<double> _gamma_energy_bounds, _time, _input_flux_amplitude, _zais,
      _atomic_numbers, _atoms, _mass_fractions, _input_flux_spectra,
      _gamma_bins;

  std::string _flux_name;

  bool _exclude_xray;
};

TEST_F(FispactProblemInputTests, setInputDataMassMat) {
  std::vector<double> input_flux = {1, 2, 3, 4};
  std::vector<double> input_flux_energy_groups = {1, 2, 3, 4, 5};
  std::vector<double> gamma_bins = {1, 2, 3, 4, 5};
  double volume = 1.0;
  double rtol = 1e-05;
  double atol = 1e-06;
  double wall_loading = 1.0;
  std::string flux_name = "test_flux";
  bool exclude_xray = false;

  auto problemParams = problemParameters();
  problemParams.setParameters<double>("rtol", rtol);
  problemParams.setParameters<double>("atol", atol);
  problemParams.setParameters<double>("atol", atol);
  problemParams.setParameters<bool>("exclude_xrays", exclude_xray);
  problemParams.setParameters<std::vector<double>>("photon_bins", gamma_bins);

  buildProblem(problemParams);
  buildObjectsMass();
  _fe_problem->resolveFispactUserObjects();

  ON_CALL(mockUtils(), GetAtomicNumberFromElementName("H"))
      .WillByDefault(::testing::Return(1));

  auto [atomic_numbers, nuclide_fractions] =
      _fe_problem->calculateMassInput(*_mat);

  setTargetEnergyGroups(input_flux_energy_groups);

  EXPECT_CALL(mockInput(), setGammaEnergyBounds(gamma_bins))
      .WillOnce(::testing::SaveArg<0>(&_gamma_bins));

  EXPECT_CALL(mockInput(), setFlux(input_flux_energy_groups, input_flux))
      .WillOnce(::testing::SaveArg<1>(&_input_flux_spectra));

  EXPECT_CALL(mockInput(), setFluxWallLoading(wall_loading))
      .WillOnce(::testing::SaveArg<0>(&_wall_loading));

  EXPECT_CALL(mockInput(), setSolverTolerance(rtol, atol))
      .WillOnce(::testing::DoAll(::testing::SaveArg<0>(&_rtol),
                                 ::testing::SaveArg<1>(&_atol)));
  EXPECT_CALL(mockInput(), setExcludeXrays(exclude_xray))
      .WillOnce(::testing::SaveArg<0>(&_exclude_xray));

  EXPECT_CALL(mockUtils(), GetAtomicNumberFromElementName("H")).Times(1);

  EXPECT_CALL(mockInput(), setDensity(_mat->getDensity()))
      .WillOnce(::testing::SaveArg<0>(&_density));

  EXPECT_CALL(mockInput(), setMassTotal(::testing::_)).Times(1);
  EXPECT_CALL(mockInput(), setMass(atomic_numbers, nuclide_fractions)).Times(1);

  // Should not get to either of these, if we have we took the wrong branch
  EXPECT_CALL(mockInput(), setFuel(::testing::_, ::testing::_)).Times(0);

  _fe_problem->setFispactInputData(*_mat, input_flux, volume, mockInput());

  ASSERT_DOUBLE_EQ(_density, 1);

  ASSERT_EQ(input_flux.size(), _input_flux_spectra.size());
  for (int i = 0; i < input_flux.size(); i++) {
    ASSERT_DOUBLE_EQ(_input_flux_spectra[i], input_flux[i]);
  }
}

TEST_F(FispactProblemInputTests, setInputDataFuelMat) {
  std::vector<double> input_flux = {1, 2, 3, 4};
  std::vector<double> input_flux_energy_groups = {1, 2, 3, 4, 5};
  std::vector<double> gamma_bins = {1, 2, 3, 4, 5};
  double volume = 1.0;
  double rtol = 1e-05;
  double atol = 1e-06;
  double wall_loading = 1.0;
  std::string flux_name = "test_flux";
  bool exclude_xray = false;

  auto problemParams = problemParameters();
  problemParams.setParameters<double>("rtol", rtol);
  problemParams.setParameters<double>("atol", atol);
  problemParams.setParameters<double>("atol", atol);
  problemParams.setParameters<bool>("exclude_xrays", exclude_xray);
  problemParams.setParameters<std::vector<double>>("photon_bins", gamma_bins);

  buildProblem(problemParams);
  buildObjectsFuel();

  _fe_problem->resolveFispactUserObjects();

  double total_mass_grams = volume * _mat->getDensity();

  auto [zais, atoms] = _fe_problem->calculateFuelInput(*_mat, total_mass_grams);

  setTargetEnergyGroups(input_flux_energy_groups);

  EXPECT_CALL(mockInput(), setGammaEnergyBounds(gamma_bins))
      .WillOnce(::testing::SaveArg<0>(&_gamma_bins));

  EXPECT_CALL(mockInput(), setFlux(input_flux_energy_groups, input_flux))
      .WillOnce(::testing::SaveArg<1>(&_input_flux_spectra));

  EXPECT_CALL(mockInput(), setFluxWallLoading(wall_loading))
      .WillOnce(::testing::SaveArg<0>(&_wall_loading));

  EXPECT_CALL(mockInput(), setSolverTolerance(rtol, atol))
      .WillOnce(::testing::DoAll(::testing::SaveArg<0>(&_rtol),
                                 ::testing::SaveArg<1>(&_atol)));

  EXPECT_CALL(mockInput(), setDensity(_mat->getDensity()))
      .WillOnce(::testing::SaveArg<0>(&_density));

  EXPECT_CALL(mockInput(), setExcludeXrays(exclude_xray))
      .WillOnce(::testing::SaveArg<0>(&_exclude_xray));

  EXPECT_CALL(mockInput(), setFuel(zais, atoms)).Times(1);

  _fe_problem->setFispactInputData(*_mat, input_flux, volume, mockInput());

  ASSERT_DOUBLE_EQ(_density, 1);

  ASSERT_EQ(input_flux.size(), _input_flux_spectra.size());
  for (int i = 0; i < input_flux.size(); i++) {
    ASSERT_DOUBLE_EQ(_input_flux_spectra[i], input_flux[i]);
  }
}
