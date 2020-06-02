#include "NEOHOOKEAN.h"
#include <iostream>

NEOHOOKEAN::NEOHOOKEAN(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _name = std::string("Neo-Hookean");
}

Real NEOHOOKEAN::getLambda() {
  return _lambda;
}

Real NEOHOOKEAN::getMu(){
  return _mu;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX NEOHOOKEAN::PK1(const MATRIX2& F)
{
  MATRIX2 ddetF; // derivative of det(F)
  ddetF(0,0) = F(1,1);
  ddetF(0,1) = -1*F(1,0);
  ddetF(1,0) = -1*F(0,1);
  ddetF(1,1) = F(0,0);

  // PK1 = (mu/2)*(d(trFTF)/dF) + (lambda*det(F) - lambda - mu)(d(detF)/dF)
  return _mu*F + (_lambda*F.determinant() - _lambda - _mu)*ddetF;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX NEOHOOKEAN::DPDF(const MATRIX& F)
{
  //derivative w.r.t. f0
  MATRIX2 df0;
  df0(0,0) = _lambda*pow(F(1,1), 2) + _mu;
  df0(1,0) = -1*_lambda*F(0,1)*F(1,1);
  df0(0,1) = -1*_lambda*F(1,0)*F(1,1);
  df0(1,1) = _lambda*(2*F(0,0)*F(1,1) - F(1,0)*F(0,1)) - _lambda - _mu;

  //derivative w.r.t. f1
  MATRIX2 df1;
  df1(0,0) = -1*_lambda*F(0,1)*F(1,1);
  df1(1,0) = _lambda*pow(F(0,1), 2) + _mu;
  df1(0,1) = _lambda*(2*F(1,0)*F(0,1) - F(0,0)*F(1,1)) + _lambda + _mu;
  df1(1,1) = -1*_lambda*F(0,0)*F(0,1);

  //derivative w.r.t. f2
  MATRIX2 df2;
  df2(0,0) = -1*_lambda*F(1,0)*F(1,1);
  df2(1,0) = _lambda*(2*F(1,0)*F(0,1) - F(0,0)*F(1,1)) + _lambda + _mu;
  df2(0,1) = _lambda*pow(F(1,0), 2) + _mu;
  df2(1,1) = -1*_lambda*F(0,0)*F(1,0);

  //derivative w.r.t. f3
  MATRIX2 df3;
  df3(0,0) = _lambda*(2*F(0,0)*F(1,1) - F(1,0)*F(0,1)) - _lambda - _mu;
  df3(1,0) = -1*_lambda*F(0,0)*F(0,1);
  df3(0,1) = -1*_lambda*F(0,0)*F(1,0);
  df3(1,1) = _lambda*pow(F(0,0), 2) + _mu;

  MATRIX4 hessian;
  hessian.col(0) = vectorize(df0);
  hessian.col(1) = vectorize(df1);
  hessian.col(2) = vectorize(df2);
  hessian.col(3) = vectorize(df3);

  return hessian;
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
Real NEOHOOKEAN::psi(const MATRIX2& F)
{
  MATRIX2 transpose = F.transpose();
  MATRIX2 mult = transpose*F;
  Real term1;
  Real term2;

  //calculate the first term, (mu/2)*(tr(FTF) - 3)
  term1 = (_mu/2) * (mult.trace() - 3);

  //calculate second term, (lambda/2)*(det(F) - 1 - mu/lambda)^2
  term2 = (_lambda/2)*pow(F.determinant() - 1 - _mu/_lambda, 2);

  return term1 + term2;
}

//finite difference method for the hessian matrix
void HessianDifference()
{
  MATRIX2 randomizedF;
  randomizedF(0,0) = rand()%10;
  randomizedF(0,1) = rand()%10;
  randomizedF(1,0) = rand()%10;
  randomizedF(1,1) = rand()%10;
  Real nu = 0.4;
  Real E = 1;

  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu = E / (2.0 * (1 + nu));

  NEOHOOKEAN* NEOHOOKEANAnalytic = new NEOHOOKEAN(lambda, mu);

  MATRIX hessianAnalytic = NEOHOOKEANAnalytic->DPDF(randomizedF);
  long double epsilon = 0.1;
  for(int x = 0; x < 6; x++)
  {
    Matrix4d hessianNumerical;
    Matrix4d matrixDiff;
    long double hessianDiff;

    MATRIX2 initial = randomizedF;
    MATRIX2 forceInitial = NEOHOOKEANAnalytic->PK1(initial);
    for(int i = 0; i < 2; i++)
    {
      MATRIX2 perturbed = randomizedF;
      perturbed(i,0) += epsilon;
      MATRIX2 forcePerturbed = NEOHOOKEANAnalytic->PK1(perturbed);
      MATRIX2 finiteDiff = (1/epsilon)*(forcePerturbed - forceInitial);
      hessianNumerical.col(i) = vectorize(finiteDiff);
    }
    for(int i = 0; i < 2; i++)
    {
      MATRIX2 perturbed = randomizedF;
      perturbed(i,1) += epsilon;
      MATRIX2 forcePerturbed = NEOHOOKEANAnalytic->PK1(perturbed);
      MATRIX2 finiteDiff = (1/epsilon)*(forcePerturbed - forceInitial);
      hessianNumerical.col(i+2) = vectorize(finiteDiff);
    }

    matrixDiff = hessianAnalytic - hessianNumerical;
    hessianDiff = matrixDiff.norm()/hessianAnalytic.norm();

    printf("hessianDiff: %Lf, epsilon: %Lf \n", hessianDiff, epsilon);
    epsilon = epsilon*0.1;
  }
}

//finite difference method for the PK1
void PK1Difference()
{
  MATRIX2 randomizedF;
  randomizedF(0,0) = rand()%10;
  randomizedF(0,1) = rand()%10;
  randomizedF(1,0) = rand()%10;
  randomizedF(1,1) = rand()%10;
  Real nu = 0.4;
  Real E = 1;

  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu = E / (2.0 * (1 + nu));

  NEOHOOKEAN* NEOHOOKEANAnalytic = new NEOHOOKEAN(lambda, mu);

  MATRIX PK1Analytic = NEOHOOKEANAnalytic->PK1(randomizedF);

  long double epsilon = 0.1;
  for(int x = 0; x < 6; x++)
  {
    MATRIX2 PK1Numerical;
    MATRIX2 matrixDiff;
    long double PK1Diff;

    MATRIX2 initial = randomizedF;
    Real forceInitial = NEOHOOKEANAnalytic->psi(initial);
    for(int i = 0; i < 2; i++)
    {
      MATRIX2 perturbed = randomizedF;
      perturbed(i,0) += epsilon;
      Real forcePerturbed = NEOHOOKEANAnalytic->psi(perturbed);
      Real finiteDiff = (1/epsilon)*(forcePerturbed - forceInitial);
      PK1Numerical(i,0) = finiteDiff;
    }

    for(int i = 0; i < 2; i++)
    {
      MATRIX2 perturbed = randomizedF;
      perturbed(i,1) += epsilon;
      Real forcePerturbed = NEOHOOKEANAnalytic->psi(perturbed);
      Real finiteDiff = (1/epsilon)*(forcePerturbed - forceInitial);
      PK1Numerical(i,1) = finiteDiff;
    }

    matrixDiff = PK1Analytic - PK1Numerical;
    PK1Diff = matrixDiff.norm()/PK1Analytic.norm();

    printf("PK1Diff: %Lf, epsilon: %Lf \n", PK1Diff, epsilon);
    epsilon = epsilon*0.1;
  }
}
