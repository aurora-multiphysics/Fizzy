#pragma once
#include "AuxKernel.h"
#include "FispactInventoryManager.h"
#include "FizzyEnums.h"
#include "libmesh/string_to_enum.h"

class FispactAuxKernel : public AuxKernel {

public:
  static InputParameters validParams() {
    InputParameters params = AuxKernel::validParams();
    return params;
  }

  FispactAuxKernel(const InputParameters &params) : AuxKernel(params) {
    const auto &fe_type = mooseVariableBase()->feType();
    if (fe_type.order != libMesh::CONSTANT ||
        fe_type.family != libMesh::MONOMIAL)
      paramError(
          "variable", "Must be of type CONSTANT MONOMIAL; variable '",
          mooseVariableBase()->name(), "' has order ",
          libMesh::Utility::enum_to_string<libMesh::Order>(fe_type.order),
          " and family ",
          libMesh::Utility::enum_to_string<libMesh::FEFamily>(fe_type.family));
  }

protected:
  virtual Real computeValue() override = 0;

  const FispactInventoryManager &getInventoryManager() {
    return getUserObject<FispactInventoryManager>("InvManager");
  }
};
