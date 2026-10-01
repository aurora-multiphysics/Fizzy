#include "FispactMaterialTest.h"
#include "gtest/gtest.h"

TEST_F(FispactMaterialTests, get_average_molar_mass) {

  std::vector<std::string> nuclide_names = {"H1", "He3"};
  std::vector<double> nuclide_fractions = {0.5, 0.5};
  std::string fraction_type = "wo";

  EXPECT_DOUBLE_EQ(_mat1->getAverageMolarMass(nuclide_names, nuclide_fractions,
                                              fraction_type),
                   1.5);

  fraction_type = "ao";

  EXPECT_DOUBLE_EQ(_mat1->getAverageMolarMass(nuclide_names, nuclide_fractions,
                                              fraction_type),
                   2);

  nuclide_names = {"H1", "H2", "He3"};
  nuclide_fractions = {0.2, 0.3, 0.5};
  fraction_type = "wo";
  EXPECT_DOUBLE_EQ(_mat1->getAverageMolarMass(nuclide_names, nuclide_fractions,
                                              fraction_type),
                   60.0 / 31.0);

  fraction_type = "ao";
  EXPECT_DOUBLE_EQ(_mat1->getAverageMolarMass(nuclide_names, nuclide_fractions,
                                              fraction_type),
                   2.3);
}

TEST_F(FispactMaterialTests, convertAtomToMassFraction) {

  std::vector<std::string> nuclide_names = {"H1", "He3"};
  std::vector<double> nuclide_fractions = {0.5, 0.5};
  std::string fraction_type = "ao";

  const double average_molar_mass = 2.0;

  _mat1->convertFromAtomToMassFraction(nuclide_names, nuclide_fractions,
                                       average_molar_mass);

  std::vector<double> expected_mass_frac = {0.25, 0.75};

  EXPECT_EQ(nuclide_fractions, expected_mass_frac);
}

TEST_F(FispactMaterialConstructorTest, checkMaterialFractions) {
  InputParameters pars_mat1 = _factory.getValidParams("FispactMaterial");

  pars_mat1.set<MooseEnum>("material_type") = "FUEL";
  pars_mat1.set<std::vector<std::string>>("nuclides") = {"H1, H2"};
  pars_mat1.set<std::vector<double>>("nuclide_fraction") = {0.9, 0.4};
  pars_mat1.set<double>("density") = 1;
  pars_mat1.set<MooseEnum>("density_units") = "g/cm3";
  pars_mat1.set<MooseEnum>("fraction_type") = "wo";
  pars_mat1.set<std::vector<SubdomainName>>("block") = {"1"};

  EXPECT_THROW(_fe_problem->addObject<FispactMaterial>(
      "FispactMaterial", "test_mat_1", pars_mat1), std::runtime_error);
}
