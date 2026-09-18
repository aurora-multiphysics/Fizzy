#include "AverageElementSize.h"
#include "FizzyObjectUnitTest.h"
#include "gtest/gtest.h"
#include <limits>
#include <stdexcept>

class FispactProblemTest : public FizzyObjectUnitTest {
public:
  FispactProblemTest() : FizzyObjectUnitTest("FizzyApp") {}
};

class FispactProblemPhotonBinsTest : public FizzyObjectUnitTest {
public:
  FispactProblemPhotonBinsTest()
      : FizzyObjectUnitTest("FizzyApp", ProblemSetup::Deferred) {}

protected:
  void buildFispactProblem(const std::vector<double> &bins) {
    auto params = problemParameters();
    params.set<std::vector<double>>("photon_bins") = bins;
    buildProblem(params);
  }
};

TEST_F(FispactProblemTest, isFlux) {
  std::vector<double> no_flux({0, 0, 0, 0});
  std::vector<double> flux({0, 0, 1, 0});
  EXPECT_TRUE(_fe_problem->isFlux(flux));
  EXPECT_FALSE(_fe_problem->isFlux(no_flux));
}

TEST_F(FispactProblemPhotonBinsTest, convertGammaEvToCount) {
  buildFispactProblem({1, 3, 5, 11});
  std::vector<double> input = {1, 2, 3};
  std::vector<double> output;
  _fe_problem->convertGammaEvToCount(input, output);
  EXPECT_EQ(output, std::vector<double>({5e5, 5e5, 3.75e5}));
}

TEST_F(FispactProblemPhotonBinsTest, acceptsPositiveIncreasingPhotonBins) {
  const std::vector<double> bins = {0.25, 0.5, 2.0, 10.0};
  ASSERT_FALSE(_fe_problem);
  buildFispactProblem(bins);
  const FispactProblem &const_problem = *_fe_problem;
  EXPECT_EQ(const_problem.getPhotonBins(), bins);
  EXPECT_EQ(_fe_problem->numPhotonBins(), 3u);
  EXPECT_EQ(_app->actionWarehouse().problemBase().get(), _fe_problem.get());
}

TEST_F(FispactProblemPhotonBinsTest, acceptsSinglePhotonBin) {
  buildFispactProblem({0.25, 0.5});
  EXPECT_EQ(_fe_problem->getPhotonBins(), std::vector<double>({0.25, 0.5}));
  EXPECT_EQ(_fe_problem->numPhotonBins(), 1u);
}

struct InvalidPhotonBins {
  std::string name;
  std::vector<double> bins;
  std::string message;
};

class FispactProblemInvalidPhotonBinsTest
    : public FispactProblemPhotonBinsTest,
      public ::testing::WithParamInterface<InvalidPhotonBins> {};

TEST_P(FispactProblemInvalidPhotonBinsTest, rejectsInvalidBoundaries) {
  const auto &test = GetParam();
  ASSERT_FALSE(_fe_problem);
  try {
    buildFispactProblem(test.bins);
    FAIL() << "Expected invalid photon boundaries to be rejected";
  } catch (const std::runtime_error &error) {
    const std::string message = error.what();
    EXPECT_NE(message.find(test.message), std::string::npos);
  }
}

INSTANTIATE_TEST_SUITE_P(
    PhotonBins, FispactProblemInvalidPhotonBinsTest,
    ::testing::Values(
        InvalidPhotonBins{"Empty", {}, "At least two"},
        InvalidPhotonBins{"OneBoundary", {1.0}, "At least two"},
        InvalidPhotonBins{"Zero", {0.0, 1.0}, "finite and strictly positive"},
        InvalidPhotonBins{
            "Negative", {-1.0, 1.0}, "finite and strictly positive"},
        InvalidPhotonBins{
            "NegativeLater", {1.0, -2.0}, "finite and strictly positive"},
        InvalidPhotonBins{"Duplicate", {1.0, 2.0, 2.0}, "strictly increasing"},
        InvalidPhotonBins{"Descending", {1.0, 3.0, 2.0}, "strictly increasing"},
        InvalidPhotonBins{"NaN",
                          {1.0, std::numeric_limits<double>::quiet_NaN()},
                          "finite and strictly positive"},
        InvalidPhotonBins{"Infinity",
                          {1.0, std::numeric_limits<double>::infinity()},
                          "finite and strictly positive"}),
    [](const ::testing::TestParamInfo<InvalidPhotonBins> &info) {
      return info.param.name;
    });

class FispactProblemStrengthTest : public FizzyObjectUnitTest {
public:
  FispactProblemStrengthTest() : FizzyObjectUnitTest("FizzyApp") {}

protected:
  void initializeStrengths() {
    // Normally allocated in initialSetup(), which also requires solve inputs.
    _fe_problem->_element_strengths.assign(_mesh->nActiveLocalElem(), -1.0);
  }

  const std::vector<double> &strengths() const {
    return _fe_problem->_element_strengths;
  }
};

TEST_F(FispactProblemStrengthTest, insertElementStrengthSumAndStore) {
  initializeStrengths();
  ASSERT_FALSE(strengths().empty());
  const auto *element = *_mesh->getActiveLocalElementRange()->begin();
  const auto local_index =
      _fe_problem->getLocalElemIndexMap().at(element->id());

  std::vector<double> expected(strengths().size(), -1.0);
  expected.at(local_index) = 7.5;

  _fe_problem->insertElementStrength(0, element, {1.25, 2.5, 3.75});

  EXPECT_EQ(strengths(), expected);
}

TEST_F(FispactProblemStrengthTest, insertElementStrengthUsesLocalElementIndex) {
  initializeStrengths();
  ASSERT_GE(strengths().size(), 2u);

  // Reverse the mapping so global IDs cannot stand in for local indices.
  auto &local_indices = _fe_problem->getLocalElemIndexMap();
  size_t local_index = strengths().size();
  for (const auto *element : *_mesh->getActiveLocalElementRange())
    local_indices.at(element->id()) = --local_index;

  auto element_it = _mesh->getActiveLocalElementRange()->begin();
  const auto *first = *element_it;
  const auto *second = *++element_it;
  ASSERT_NE(local_indices.at(first->id()), first->id());

  std::vector<double> expected(strengths().size(), -1.0);
  expected.back() = 6.0;
  expected.at(expected.size() - 2) = 9.0;

  _fe_problem->insertElementStrength(0, first, {1.0, 2.0, 3.0});
  _fe_problem->insertElementStrength(0, second, {4.0, 5.0});

  EXPECT_EQ(strengths(), expected);
}

TEST_F(FispactProblemTest, loadMolarMasses) {

  const std::unordered_map<std::string, double> &molar_masses =
      _fe_problem->getMolarMassMap();

  EXPECT_EQ(molar_masses.at("H1"), 1);
  EXPECT_EQ(molar_masses.at("H2"), 2);
  EXPECT_EQ(molar_masses.at("He3"), 3);

  EXPECT_EQ(molar_masses.count("H33"), 0u);
}

// class FispactProblemTestUserObjects : public FizzyObjectUnitTest {
// public:
//   FispactProblemTestUserObjects() : FizzyObjectUnitTest("FizzyApp") {
//     buildObjects();
//   }
//
// protected:
//   void buildObjects() {
//
//     InputParameters pars_mat1 = _factory.getValidParams("FispactMaterial");
//
//     pars_mat1.set<MooseEnum>("material_type") = "FUEL";
//     pars_mat1.set<std::vector<std::string>>("nuclides") = {"H1"};
//     pars_mat1.set<std::vector<double>>("nuclide_fraction") = {100};
//     pars_mat1.set<double>("density") = 1;
//     pars_mat1.set<MooseEnum>("density_units") = "g/cm3";
//     pars_mat1.set<MooseEnum>("fraction_type") = "wo";
//     pars_mat1.set<std::vector<SubdomainName>>("block") = {"1"};
//     _fe_problem->addObject<FispactMaterial>("FispactMaterial", "Hydrogen",
//                                             pars_mat1);
//     _mat1 = &_fe_problem->getUserObject<FispactMaterial>("Hydrogen");
//   }
//
//   const FispactMaterial *_mat1;
//   const FispactSchedule *_schedule;
//   // const FispactFluxInput *;
//   // const FispactNuclearDataPaths *;
// };

/// Test fixture to make sure getElementMaterial throws if two materials are
/// defined on the same block
class FispactProblemTestThrowMaterials : public FizzyObjectUnitTest {
public:
  FispactProblemTestThrowMaterials() : FizzyObjectUnitTest("FizzyApp") {
    buildObjects();
  }

protected:
  void buildObjects() {

    InputParameters pars_mat1 = _factory.getValidParams("FispactMaterial");
    InputParameters pars_mat2 = _factory.getValidParams("FispactMaterial");

    pars_mat1.set<MooseEnum>("material_type") = "FUEL";
    pars_mat1.set<std::vector<std::string>>("nuclides") = {"H1"};
    pars_mat1.set<std::vector<double>>("nuclide_fraction") = {1};
    pars_mat1.set<double>("density") = 1;
    pars_mat1.set<MooseEnum>("density_units") = "g/cm3";
    pars_mat1.set<MooseEnum>("fraction_type") = "wo";
    pars_mat1.set<std::vector<SubdomainName>>("block") = {"1"};
    _fe_problem->addObject<FispactMaterial>("FispactMaterial", "Hydrogen",
                                            pars_mat1);
    _mat1 = &_fe_problem->getUserObject<FispactMaterial>("Hydrogen");

    pars_mat2.set<MooseEnum>("material_type") = "FUEL";
    pars_mat2.set<std::vector<std::string>>("nuclides") = {"He3"};
    pars_mat2.set<std::vector<double>>("nuclide_fraction") = {1};
    pars_mat2.set<double>("density") = 3;
    pars_mat2.set<MooseEnum>("density_units") = "g/cm3";
    pars_mat2.set<MooseEnum>("fraction_type") = "wo";
    pars_mat2.set<std::vector<SubdomainName>>("block") = {"1"};
    _fe_problem->addUserObject("FispactMaterial", "Helium", pars_mat2);
    _mat2 = &_fe_problem->getUserObject<FispactMaterial>("Helium");
  }

  const FispactMaterial *_mat1;
  const FispactMaterial *_mat2;
};

TEST_F(FispactProblemTestThrowMaterials, getElementMaterial) {
  dof_id_type elem_id = _fe_problem->mesh().elemPtr(0)->id();
  EXPECT_THROW(_fe_problem->getElementMaterial(elem_id), std::runtime_error);
}

class FispactScheduleTest : public FizzyObjectUnitTest {
public:
  FispactScheduleTest() : FizzyObjectUnitTest("FizzyApp") {}

  void buildFispactSchedule(std::vector<double> times,
                            std::vector<double> flux_amplitude) {
    InputParameters fispact_schedule_params =
        _factory.getValidParams("FispactSchedule");
    fispact_schedule_params.set<std::vector<double>>("flux_amplitude") =
        flux_amplitude;
    fispact_schedule_params.set<std::vector<double>>("times") = times;
    _fe_problem->addObject<FispactSchedule>("FispactSchedule", "Schedule",
                                            fispact_schedule_params);
  }
};

struct InvalidScheduleParams {
  std::string name;
  std::vector<double> times;
  std::vector<double> flux_amplitude;
  std::string message;
};

class FispactScheduleInvalidParamsTest
    : public FispactScheduleTest,
      public ::testing::WithParamInterface<InvalidScheduleParams> {};

TEST_P(FispactScheduleInvalidParamsTest, rejectsInvalidParams) {
  const auto &test = GetParam();
  ASSERT_TRUE(_fe_problem);
  try {
    buildFispactSchedule(test.times, test.flux_amplitude);
    FAIL() << "Expected invalid schedule params to be rejected";
  } catch (const std::runtime_error &error) {
    const std::string message = error.what();
    std::cout << message << std::endl;
    EXPECT_NE(message.find(test.message), std::string::npos);
  }
}

INSTANTIATE_TEST_SUITE_P(
    ScheduleParams, FispactScheduleInvalidParamsTest,
    ::testing::Values(
        InvalidScheduleParams{"BothEmpty", {}, {}, "cannot be empty vectors."},
        InvalidScheduleParams{
            "DifferentSizes", {100, 200}, {1, 2, 3}, "must be equal."},

        InvalidScheduleParams{
            "Times0", {100, 0, 300}, {1, 2, 3}, "and strictly positive."},

        InvalidScheduleParams{"Timesnegative",
                              {100, -10, 300},
                              {1, 2, 3},
                              "and strictly positive."},

        InvalidScheduleParams{"flux_amplitudeNegative",
                              {100, 200, 300},
                              {-1, 2, 3},
                              "and strictly positive."},

        InvalidScheduleParams{"flux_amplitudeInf",
                              {100, 200, 300},
                              {std::numeric_limits<double>::infinity(), 2, 3},
                              "and strictly positive."},

        InvalidScheduleParams{
            "timesInf",
            {std::numeric_limits<double>::infinity(), 200, 300},
            {1, 2, 3},
            "and strictly positive."},

        InvalidScheduleParams{
            "flux_amplitudeAndtimesInf",
            {std::numeric_limits<double>::infinity(), 200, 300},
            {std::numeric_limits<double>::infinity(), 2, 3},
            "must be finite"}),
    [](const ::testing::TestParamInfo<InvalidScheduleParams> &info) {
      return info.param.name;
    });
