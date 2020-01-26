#include "MATERIAL.h"
#include "GNUPLOT.h"

///////////////////////////////////////////////////////////////////////
// draw PK1 to gnuplot, just for probing
///////////////////////////////////////////////////////////////////////
void MATERIAL::plotPK1()
{
  vector<Real> lambdas;
  vector<Real> PK1s;

  Real range = 1.5;
  Real delta = 0.01;
  Real lambda = -range;
  while (lambda < range)
  {
    MATRIX2 F = MATRIX2::Identity();
    F(0,0) = lambda;
    MATRIX currentPK1 = this->PK1(F);

    lambdas.push_back(lambda);
    PK1s.push_back(currentPK1(0,0));

    lambda += delta;
  }
  char buffer[512];
  sprintf(buffer, "%s", _name.c_str());

  GNUPLOT plot;
  plot.addPlot(lambdas, PK1s, buffer);
  sprintf(buffer, "%s.ps", _name.c_str());
  plot.outputPlots(buffer);
}
