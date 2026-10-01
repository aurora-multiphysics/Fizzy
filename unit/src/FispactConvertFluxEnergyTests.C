#include "FispactProblemTest.h"

class FispactProblemConvertFluxEnergyTest : public FizzyObjectUnitTest {
public:
  FispactProblemConvertFluxEnergyTest()
      : FizzyObjectUnitTest("FizzyApp", ProblemSetup::Deferred) {}

protected:
  void buildConversionProblem(const std::string &conversion_type) {
    auto params = problemParameters();
    params.set<MooseEnum>("conversion_type") = conversion_type;
    buildProblem(params);
  }
};

TEST_F(FispactProblemConvertFluxEnergyTest, convertByLethargy) {
  buildConversionProblem("LETHARGY");

  const std::vector<double> input_bounds = {1, 2, 4, 8};
  const std::vector<double> original_flux = {10, 20, 30};
  const std::vector<double> target_bounds = {1, 4, 8};
  const std::vector<double> converted_flux = {11, 49};

  setTargetEnergyGroups(target_bounds);
  auto &utils = mockUtils();

  EXPECT_CALL(
      utils, GroupConvertByLethargy(input_bounds, original_flux, target_bounds))
      .Times(1)
      .WillOnce(::testing::Return(converted_flux));

  EXPECT_CALL(utils,
              GroupConvertByEnergy(::testing::_, ::testing::_, ::testing::_))
      .Times(0);

  auto flux = original_flux;
  _fe_problem->convertFluxEnergyGroups(flux, input_bounds);
  EXPECT_EQ(flux, converted_flux);
}

TEST_F(FispactProblemConvertFluxEnergyTest, convertByEnergy) {
  buildConversionProblem("ENERGY");

  const std::vector<double> input_bounds = {1, 2, 4, 8};
  const std::vector<double> original_flux = {10, 20, 30};
  const std::vector<double> target_bounds = {1, 4, 8};
  const std::vector<double> converted_flux = {22, 38};

  setTargetEnergyGroups(target_bounds);
  auto &utils = mockUtils();

  EXPECT_CALL(utils,
              GroupConvertByEnergy(input_bounds, original_flux, target_bounds))
      .Times(1)
      .WillOnce(::testing::Return(converted_flux));

  EXPECT_CALL(utils,
              GroupConvertByLethargy(::testing::_, ::testing::_, ::testing::_))
      .Times(0);

  auto flux = original_flux;
  _fe_problem->convertFluxEnergyGroups(flux, input_bounds);
  EXPECT_EQ(flux, converted_flux);
}
