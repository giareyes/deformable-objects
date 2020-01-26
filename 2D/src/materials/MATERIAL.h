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
  
  // S = first Piola-Kirchoff stress tensor
  virtual MATRIX PK2(const MATRIX2& F) = 0;

  // derivative of PK2 w.r.t E = 0.5 * (F^T * F - I)
  virtual MATRIX DSDE() {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " NOT IMPLEMENTED. " << std::endl;
    exit(0);
  };

  // needed by the HYPER_TAN interpolated material, gradient of psi
  // w.r.t. just one entry in F
  virtual Real DpsiDF(const MATRIX2& F, int index)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " NOT IMPLEMENTED. " << std::endl;
    exit(0);
  }
  
  // derivative of PK1 w.r.t. F
  virtual MATRIX DPDF(const MATRIX& F) = 0;

  // get the strain energy
  virtual Real psi(const MATRIX2& F) = 0;
 
  const std::string& name() const { return _name; };

  // draw PK1 to gnuplot, just for probing
  void plotPK1();

  // get the eigenvalues of DPDF, if we have expressions for them
  virtual VEC4 eigenvalues(const MATRIX2& F) { return VEC4::Zero(); };
  
  // get the eigensysem of DPDF, if we have expressions for them
  virtual void eigensystem(const MATRIX2& F, MATRIX4& vectors, VEC4& values) 
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " NOT IMPLEMENTED. " << std::endl;
    exit(0);
  }; 

  virtual MATRIX cauchy(const MATRIX2& F)
  {
    const Real J = F.determinant();
    return (1.0 / J) * PK1(F) * F.transpose();
  }

protected:
  MATRIX2 DFDF(const MATRIX2& F, int i, int j)
  {
    MATRIX2 zero;
    zero.setZero();
    zero(i,j) = 1;
    return zero;
  };

  // the name of the material, for display purposes
  std::string _name;
};

#endif
