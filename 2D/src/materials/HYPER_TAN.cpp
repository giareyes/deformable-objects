#include "HYPER_TAN.h"
#include <iostream>

using namespace std;

HYPER_TAN::HYPER_TAN(const Real mu0, const Real mu1, MATERIAL* psi1D) :
  _psi1D(psi1D), _mu0(mu0), _mu1(mu1)
{
  //_k = 50;
  //_k = 10;
  //_k = 100;
  //_k = 2.5;
  //_k = 0.1;
  //_k = 1;
  _k = 5;
  //_k = 0.01;

  Real hardeningLength = 1.5;
  //Real hardeningLength = 2.0;
  //Real hardeningLength = 1.0;

  // this is the psi of the unit test
  //_threshold = 0.34097;
  _threshold = pow(hardeningLength - 1.0, 2.0);
  //_threshold = 0.3;
  //_threshold = -10;
  _name = string("Tanh");
  std::cout << " Hypertan interpolation: " << std::endl;
  std::cout << " \t mu0: " << _mu0 << std::endl;
  std::cout << " \t mu1: " << _mu1 << std::endl;
  std::cout << " \t threshold: " << _threshold << std::endl;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX HYPER_TAN::PK1(const MATRIX2& F)
{
  return DmuDF(F) * _psi1D->psi(F) + mu(F) * _psi1D->PK1(F);
  //return _psi1D->PK1(F);
}

///////////////////////////////////////////////////////////////////////
// S = second Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX HYPER_TAN::PK2(const MATRIX2& F)
{
  return F.inverse() * PK1(F);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX HYPER_TAN::DPDF(const MATRIX& F)
{
  Real psi_1D = _psi1D->psi(F);
  MATRIX2 P_1D = _psi1D->PK1(F);

  Real tanhTerm = tanh(_k * (psi_1D - _threshold));
  Real sechTerm = 1.0 / cosh(_k * (psi_1D - _threshold));
  Real alpha = (_mu1 - _mu0) * _k * 0.5 * sechTerm * sechTerm;

  MATRIX final(4,4);
  final.setZero();
  for (int col = 0; col < 4; col++)
  {
    const int i = col % 2;
    const int j = col / 2;
    Real DalphaDF = -(_mu1 - _mu0) * _k * _k * tanhTerm * sechTerm * sechTerm * P_1D(i,j);
    Real DmuDF = (_mu1 - _mu0) * _k * 0.5 * sechTerm * sechTerm * P_1D(i,j);

    MATRIX2 DPDF_i = DalphaDF * psi_1D * P_1D;
    DPDF_i += alpha * P_1D(i,j) * P_1D;
    DPDF_i += DmuDF * P_1D;
    
    final.col(col) = flatten(DPDF_i);
  }

  final += (alpha * psi_1D + mu(F)) * _psi1D->DPDF(F);

  return final;
}
/*
{
  Real psi_1D = _psi1D->psi(F);
  // sech = 1 / cosh
  Real coeff = (0.5 * _k) * pow(1.0 / cosh(_k * (psi_1D - _threshold)), 2.0);
  Real omega = (_mu1 - _mu0) * coeff + mu(F);

  MATRIX2 P_1D = _psi1D->PK1(F);

  //MATRIX2 secondTerm = -(_mu1 - _mu0) * (_k * _k) * tanh(_k * (psi_1D - _threshold)) * pow(1.0 / cosh(_k * (psi_1D - _threshold)), 2.0) * P_1D;
  MATRIX2 secondTerm = -(_mu1 - _mu0) * (_k * _k) * tanh(_k * (psi_1D - _threshold)) * pow(1.0 / cosh(_k * (psi_1D - _threshold)), 2.0) * P_1D * psi_1D;
  MATRIX2 thirdTerm = (_mu1 - _mu0) * (0.5 * _k) * pow(1.0 / cosh(_k * (psi_1D - _threshold)), 2.0) * P_1D;
  MATRIX2 DmuDF_i = DmuDF(F);

  MATRIX final(4,4);
  final.setZero();
  for (int i = 0; i < 4; i++)
  {
    const int j = i % 2;
    const int k = i / 2;
    MATRIX2 DPDF_i = DmuDF_i(j,k) * P_1D;
    DPDF_i += secondTerm * P_1D(j,k);
    DPDF_i += thirdTerm * P_1D(j,k);
    final.col(i) = flatten(DPDF_i);
  }

  final += omega * _psi1D->DPDF(F);
  return final;
}
*/

///////////////////////////////////////////////////////////////////////
// flatten a matrix into a vector, stacking each of the columns
// on top of each other
///////////////////////////////////////////////////////////////////////
VECTOR HYPER_TAN::flatten(const MATRIX& A)
{
  VECTOR final(A.rows() * A.cols());

  int i = 0;
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++, i++)
      final[i] = A(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real HYPER_TAN::psi(const MATRIX2& F)
{
  return mu(F) * _psi1D->psi(F);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real HYPER_TAN::mu(const MATRIX2& F) const
{
  Real T = 0.5 + 0.5 * tanh(_k * (_psi1D->psi(F) - _threshold));

  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " Psi: " << _psi1D->psi(F) << endl;

  /*
  cout << " psi: " << _psi1D->psi(F) << endl;
  cout << " mu: " << _mu0 * (1.0 - T) + _mu1 * T << endl;
  exit(0);
  */

  return _mu0 * (1.0 - T) + _mu1 * T;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX HYPER_TAN::DmuDF(const MATRIX2& F)
{
  Real psi_1D = _psi1D->psi(F);
  //Real coeff = (0.5 * _k) * pow(sech(_k * (psi_1D - _threshold)), 2.0);
  
  // sech = 1 / cosh
  Real sech = 1.0 / cosh(_k * (psi_1D - _threshold));
  Real coeff = (0.5 * _k) * sech * sech;

  return (_mu1 - _mu0) * coeff * _psi1D->PK1(F);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real HYPER_TAN::DmuDF(const MATRIX2& F, int index)
{
  Real psi_1D = _psi1D->psi(F);
  Real coeff = (_mu1 - _mu0) * (0.5 * _k) * pow(1.0 / cosh(_k * (psi_1D - _threshold)), 2.0);

  const int i = index % 2;
  const int j = index / 2;

  return coeff * _psi1D->PK1(F)(i,j);
}
