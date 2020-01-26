#include "COMPOSITE_MATERIAL.h"
#include <iostream>

using namespace std;

COMPOSITE_MATERIAL::COMPOSITE_MATERIAL()
{
  _name = string("Composite");
}

COMPOSITE_MATERIAL::~COMPOSITE_MATERIAL()
{
  for (int x = 0; x < _materials.size(); x++)
    delete _materials[x];
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX COMPOSITE_MATERIAL::PK1(const MATRIX2& F)
{
  assert(_materials.size() > 0);
  MATRIX final = _materials[0]->PK1(F);

  for (int x = 1; x < _materials.size(); x++)
    final += _materials[x]->PK1(F);

  return final;
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor = 2 mu E
///////////////////////////////////////////////////////////////////////
MATRIX COMPOSITE_MATERIAL::PK2(const MATRIX2& F)
{
  assert(_materials.size() > 0);
  MATRIX final = _materials[0]->PK2(F);

  for (int x = 1; x < _materials.size(); x++)
    final += _materials[x]->PK2(F);

  return final;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX COMPOSITE_MATERIAL::DPDF(const MATRIX& F)
{
  assert(_materials.size() > 0);
  MATRIX final = _materials[0]->DPDF(F);

  for (int x = 1; x < _materials.size(); x++)
    final += _materials[x]->DPDF(F);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real COMPOSITE_MATERIAL::psi(const MATRIX2& F)
{
  assert(_materials.size() > 0);
  Real final = _materials[0]->psi(F);

  for (int x = 1; x < _materials.size(); x++)
    final += _materials[x]->psi(F);

  return final;
}

///////////////////////////////////////////////////////////////////////
// add a material to the composite
///////////////////////////////////////////////////////////////////////
void COMPOSITE_MATERIAL::addMaterial(MATERIAL* material)
{
  _materials.push_back(material);

  char buffer[1024];
  sprintf(buffer, "%s_%s", _name.c_str(), material->name().c_str());

  _name = string(buffer);

  cout << " Name: " << _name.c_str() << endl;
}
