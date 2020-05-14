#ifndef MATERIAL_H
#define MATERIAL_H

#include "SETTINGS.h"
#include <iostream>

class MATERIAL
{
public:
  MATERIAL() : _name("None") {};
  virtual ~MATERIAL() {};

  // P = first Piola-Kirchoff stress tensor
  // P = F * S
  virtual MATRIX PK1(const MATRIX2& F) = 0;

  // derivative of PK1 w.r.t. F
  virtual MATRIX DPDF(const MATRIX& F) = 0;

  // get the strain energy
  virtual Real psi(const MATRIX2& F) = 0;

  virtual Real getLambda() = 0;

  virtual Real getMu() = 0;

  const std::string& name() const { return _name; };

protected:
  // the name of the material, for display purposes
  std::string _name;
};

#endif
