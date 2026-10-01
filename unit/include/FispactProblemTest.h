#pragma once

#include "FizzyObjectUnitTest.h"
#include "OpenMCFluxInput.h"

class FispactProblemTest : public FizzyObjectUnitTest {
public:
  FispactProblemTest() : FizzyObjectUnitTest("FizzyApp") {}
};

class FispactProblemTestUserObjects : public FizzyObjectUnitTest {
public:
  FispactProblemTestUserObjects()
      : FizzyObjectUnitTest("FizzyApp", ProblemSetup::Deferred) {}

protected:
  void buildObjectsMass() {
    buildMaterialMass();
    buildFluxInput();
    buildSchedule();
    buildNuclearDataPaths();
  }

  void buildObjectsFuel() {
    buildMaterialFuel();
    buildFluxInput();
    buildSchedule();
    buildNuclearDataPaths();
  }

  void buildMaterialMass() {
    InputParameters mat_params = _factory.getValidParams("FispactMaterial");

    mat_params.set<MooseEnum>("material_type") = "MASS";
    mat_params.set<std::vector<std::string>>("nuclides") = {"H"};
    mat_params.set<std::vector<double>>("nuclide_fraction") = {1};
    mat_params.set<double>("density") = 1;
    mat_params.set<MooseEnum>("density_units") = "g/cm3";
    mat_params.set<MooseEnum>("fraction_type") = "wo";
    mat_params.set<std::vector<SubdomainName>>("block") = {"1"};
    _fe_problem->addObject<FispactMaterial>("FispactMaterial", "test_material",
                                            mat_params);
    _mat = &_fe_problem->getUserObject<FispactMaterial>("test_material");
  }

  void buildMaterialFuel() {
    InputParameters mat_params = _factory.getValidParams("FispactMaterial");

    mat_params.set<MooseEnum>("material_type") = "FUEL";
    mat_params.set<std::vector<std::string>>("nuclides") = {"H1"};
    mat_params.set<std::vector<double>>("nuclide_fraction") = {1};
    mat_params.set<double>("density") = 1;
    mat_params.set<MooseEnum>("density_units") = "g/cm3";
    mat_params.set<MooseEnum>("fraction_type") = "wo";
    mat_params.set<std::vector<SubdomainName>>("block") = {"1"};
    _fe_problem->addObject<FispactMaterial>("FispactMaterial", "test_material",
                                            mat_params);
    _mat = &_fe_problem->getUserObject<FispactMaterial>("test_material");
  }

  void buildSchedule() {

    InputParameters schedule_params =
        _factory.getValidParams("FispactSchedule");

    schedule_params.set<std::vector<double>>("times") = {10, 10, 10};
    schedule_params.set<std::vector<double>>("flux_amplitude") = {10, 10, 10};
    _fe_problem->addObject<FispactSchedule>("FispactSchedule", "test_schedule",
                                            schedule_params);
    _schedule = &_fe_problem->getUserObject<FispactSchedule>("test_schedule");
  }

  void buildFluxInput() {
    InputParameters flux_params = _factory.getValidParams("OpenMCFluxInput");

    flux_params.set<size_t>("flux_tally_id") = 1;
    flux_params.set<size_t>("energy_filter_id") = 1;
    flux_params.set<double>("wall_loading") = 0;
    flux_params.set<FileName>("statepoint_filename") = "statepoint.10.h5";
    _fe_problem->addObject<OpenMCFluxInput>("OpenMCFluxInput", "test_flux",
                                            flux_params);
    _flux_input = &_fe_problem->getUserObject<OpenMCFluxInput>("test_flux");
  }

  void buildNuclearDataPaths() {
    InputParameters nd_params =
        _factory.getValidParams("FispactNuclearDataPaths");

    _fe_problem->addObject<FispactNuclearDataPaths>("FispactNuclearDataPaths",
                                                    "test_nd", nd_params);
    _nuclear_data_paths =
        &_fe_problem->getUserObject<FispactNuclearDataPaths>("test_nd");
  }

  const FispactMaterial *_mat = nullptr;
  const FispactSchedule *_schedule = nullptr;
  const FispactFluxInput *_flux_input = nullptr;
  const FispactNuclearDataPaths *_nuclear_data_paths = nullptr;
};
