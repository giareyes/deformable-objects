#include "STVK.h"
#include <iostream>

STVK::STVK(const Real lambda, const Real mu) :
  _lambda(lambda), _mu(mu)
{
  _name = std::string("StVK");
}

Real STVK::getLambda() {
  return _lambda;
}

Real STVK::getMu(){
  return _mu;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// P = F * S
///////////////////////////////////////////////////////////////////////
MATRIX STVK::PK1(const MATRIX2& F)
{
  MATRIX2 transpose = F.transpose();
  //MATRIX2 identity;
  MATRIX2 mult = transpose*F;
  //identity.setIdentity();
  //mult = mult - identity;

  //first find the derivative of the first term
  //MATRIX2 dFirst = 4*_mu*(F*mult);

  //next find the derivative of the second term
  //chain rule gets rid of /2 so that we just have lambda * tr(blah) * d(blah)/dF
  MATRIX2 dSecond = _lambda*mult.trace()*2*F;

  return dSecond; // + dFirst;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX STVK::DPDF(const MATRIX& F)
{
  //derivative w.r.t. f0
  MATRIX2 df0;
  df0(0,0) = 4*_mu*(3*pow(F(0,0), 2) + pow(F(1,0), 2) + pow(F(0,1), 2))
              + 2*_lambda*(3*pow(F(0,0), 2) + pow(F(1,0), 2) + pow(F(0,1), 2) + pow(F(1,1), 2));
  df0(1,0) = 4*_mu*(2*F(1,0)*F(0,0) + F(0,1)*F(1,1)) + 4*_lambda*F(0,0)*F(1,0);
  df0(0,1) = 4*_mu*(2*F(0,1)*F(0,0) + F(1,0)*F(1,1)) + 4*_lambda*F(0,0)*F(0,1);
  df0(1,1) = 4*_mu*F(1,0)*F(0,1) + 4*_lambda*F(0,0)*F(1,1);

  //derivative w.r.t. f1
  MATRIX2 df1;
  df1(0,0) = 4*_mu*(2*F(0,0)*F(1,0) + F(0,1)*F(1,1)) + 4*_lambda*F(0,0)*F(1,0);
  df1(1,0) = 4*_mu*(pow(F(0,0), 2) + 3*pow(F(1,0), 2) + pow(F(1,1), 2))
              + 2*_lambda*(pow(F(0,0), 2) + 3*pow(F(1,0), 2) + pow(F(0,1), 2) + pow(F(1,1), 2));
  df1(0,1) = 4*_mu*F(0,0)*F(1,1) + 4*_lambda*F(1,0)*F(0,1);
  df1(1,1) = 4*_mu*(2*F(1,0)*F(1,1) + F(0,0)*F(0,1)) + 4*_lambda*F(1,0)*F(1,1);

  //derivative w.r.t. f2
  MATRIX2 df2;
  df2(0,0) = 4*_mu*(2*F(0,0)*F(0,1) + F(1,0)*F(1,1)) + 4*_lambda*F(0,0)*F(0,1);
  df2(1,0) = 4*_mu*F(0,0)*F(1,1) + 4*_lambda*F(1,0)*F(0,1);
  df2(0,1) = 4*_mu*(pow(F(0,0), 2) + 3*pow(F(0,1), 2) + pow(F(1,1), 2))
              + 2*_lambda*(pow(F(0,0), 2) + pow(F(1,0), 2) + 3*pow(F(0,1), 2) + pow(F(1,1), 2));
  df2(1,1) = 4*_mu*(2*F(0,1)*F(1,1) + F(0,0)*F(1,0)) + 4*_lambda*F(0,1)*F(1,1);

  //derivative w.r.t. f3
  MATRIX2 df3;
  df3(0,0) = 4*_mu*F(1,0)*F(0,1) + 4*_lambda*F(0,0)*F(1,1);
  df3(1,0) = 4*_mu*(F(0,1)*F(0,0) + 2*F(1,0)*F(1,1)) + 4*_lambda*F(1,0)*F(1,1);
  df3(0,1) = 4*_mu*(F(1,0)*F(0,0) + 2*F(0,1)*F(1,1)) + 4*_lambda*F(0,1)*F(1,1);
  df3(1,1) = 4*_mu*(3*pow(F(1,1), 2) + pow(F(1,0), 2) + pow(F(0,1), 2))
              + 2*_lambda*(pow(F(0,0), 2) + pow(F(1,0), 2) + pow(F(0,1), 2) + 3*pow(F(1,1), 2));

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
Real STVK::psi(const MATRIX2& F)
{
  MATRIX2 transpose = F.transpose();
  MATRIX2 mult = transpose*F;
  MATRIX2 identity;
  Real term1;
  Real term2;
  identity.setIdentity();

  //F transpose * F - I
  mult = mult - identity;

  //calculate the first term, mu*norm(f transpose f - i)^2
  term1 = _mu * pow(mult.norm(), 2);

  //calculate second term, (lambda/2)*trace(f transpose f - i)^2
  term2 = (_lambda/2)*pow(mult.trace(), 2);
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

  STVK* STVKAnalytic = new STVK(lambda, mu);

  MATRIX hessianAnalytic = STVKAnalytic->DPDF(randomizedF);
  long double epsilon = 0.1;
  for(int x = 0; x < 6; x++)
  {
    Matrix4d hessianNumerical;
    Matrix4d matrixDiff;
    long double hessianDiff;

    MATRIX2 initial = randomizedF;
    MATRIX2 forceInitial = STVKAnalytic->PK1(initial);
    for(int i = 0; i < 2; i++)
    {
      MATRIX2 perturbed = randomizedF;
      perturbed(i,0) += epsilon;
      MATRIX2 forcePerturbed = STVKAnalytic->PK1(perturbed);
      MATRIX2 finiteDiff = (1/epsilon)*(forcePerturbed - forceInitial);
      hessianNumerical.col(i) = vectorize(finiteDiff);
    }
    for(int i = 0; i < 2; i++)
    {
      MATRIX2 perturbed = randomizedF;
      perturbed(i,1) += epsilon;
      MATRIX2 forcePerturbed = STVKAnalytic->PK1(perturbed);
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

  STVK* STVKAnalytic = new STVK(lambda, mu);

  MATRIX PK1Analytic = STVKAnalytic->PK1(randomizedF);

  long double epsilon = 0.1;
  for(int x = 0; x < 6; x++)
  {
    MATRIX2 PK1Numerical;
    MATRIX2 matrixDiff;
    long double PK1Diff;

    MATRIX2 initial = randomizedF;
    Real forceInitial = STVKAnalytic->psi(initial);
    for(int i = 0; i < 2; i++)
    {
      MATRIX2 perturbed = randomizedF;
      perturbed(i,0) += epsilon;
      Real forcePerturbed = STVKAnalytic->psi(perturbed);
      Real finiteDiff = (1/epsilon)*(forcePerturbed - forceInitial);
      PK1Numerical(i,0) = finiteDiff;
    }

    for(int i = 0; i < 2; i++)
    {
      MATRIX2 perturbed = randomizedF;
      perturbed(i,1) += epsilon;
      Real forcePerturbed = STVKAnalytic->psi(perturbed);
      Real finiteDiff = (1/epsilon)*(forcePerturbed - forceInitial);
      PK1Numerical(i,1) = finiteDiff;
    }

    matrixDiff = PK1Analytic - PK1Numerical;
    PK1Diff = matrixDiff.norm()/PK1Analytic.norm();

    printf("PK1Diff: %Lf, epsilon: %Lf \n", PK1Diff, epsilon);
    epsilon = epsilon*0.1;
  }
}
