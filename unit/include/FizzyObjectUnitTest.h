//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "gtest/gtest.h"
#include <memory>

#include "AppFactory.h"
#include "FEProblem.h"
#include "FileMesh.h"
#include "FileMeshGenerator.h"
#include "FispactProblem.h"
#include "MooseMain.h"
#include "MooseMesh.h"

class FizzyObjectUnitTest : public ::testing::Test {
public:
  enum class ProblemSetup { Automatic, Deferred };

  /**
   * @param app_name The name of client's application
   * @param setup Whether to construct the problem now or in the test body
   */
  FizzyObjectUnitTest(const std::string &app_name,
                      ProblemSetup setup = ProblemSetup::Automatic)
      : _app(Moose::createMooseApp(app_name, 0, nullptr)),
        _factory(_app->getFactory()) {
    buildMesh();
    if (setup == ProblemSetup::Automatic) {
      auto params = problemParameters();
      buildProblem(params);
    }
  }

protected:
  void buildMesh() {

    InputParameters mesh_params = _factory.getValidParams("FileMesh");
    mesh_params.set<MeshFileName>("file") = "../geometry/cube.e";

    _mesh = _factory.createUnique<FileMesh>("FileMesh", "name1", mesh_params);
    _mesh->setMeshBase(_mesh->buildMeshBaseObject());
    _mesh->buildMesh();
    _mesh->init();
    _mesh->prepare(&_mesh->getMesh());
  }

  InputParameters problemParameters() {
    InputParameters problem_params = _factory.getValidParams("FispactProblem");
    problem_params.set<MooseMesh *>("mesh") = _mesh.get();
    problem_params.set<std::string>(MooseBase::name_param) = "name2";
    problem_params.set<UserObjectName>("fispact_schedule_uo") = "name2";
    problem_params.set<UserObjectName>("fispact_nuclear_data_uo") = "name2";
    problem_params.set<UserObjectName>("fispact_input_flux_uo") = "name2";
    problem_params.set<FileName>("molar_mass_data") = "molar_mass_mock.h5";
    return problem_params;
  }

  void buildProblem(InputParameters &problem_params) {
    _fe_problem = _factory.create<FispactProblem>("FispactProblem", "problem",
                                                  problem_params);

    _fe_problem->createQRules(libMesh::QGAUSS, libMesh::FIRST, libMesh::FIRST,
                              libMesh::FIRST);

    _fe_problem->callFispactFactory();

    _app->actionWarehouse().problemBase() = _fe_problem;
  }

  template <typename T>
  T &addObject(const std::string &type, const std::string &name,
               InputParameters &params);

  std::unique_ptr<MooseMesh> _mesh;
  std::shared_ptr<MooseApp> _app;
  Factory &_factory;
  std::shared_ptr<FispactProblem> _fe_problem;
};

template <typename T>
T &FizzyObjectUnitTest::addObject(const std::string &type,
                                  const std::string &name,
                                  InputParameters &params) {
  auto objects = _fe_problem->addObject<T>(type, name, params);
  mooseAssert(objects.size() == 1, "Doesn't work with threading");
  return *objects[0];
}
