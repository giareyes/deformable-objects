#ifndef COMPOSITE_MATERIAL_H
#define COMPOSITE_MATERIAL_H

#include "MATERIAL.h"

///////////////////////////////////////////////////////////////////////
// Co-rotational length resistence,
///////////////////////////////////////////////////////////////////////
class COMPOSITE_MATERIAL : public MATERIAL
{
public:
  COMPOSITE_MATERIAL();
  ~COMPOSITE_MATERIAL(); 

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  MATRIX PK1(const MATRIX2& F);

  // S = second Piola-Kirchoff stress tensor
  MATRIX PK2(const MATRIX2& F);
  
  // derivative of PK1 w.r.t. F
  MATRIX DPDF(const MATRIX& F);

  // get the strain energy
  Real psi(const MATRIX2& F);

  // add a material to the composite
  void addMaterial(MATERIAL* material);

private:
  std::vector<MATERIAL*> _materials;
};

#endif
