// get sigma vs y for VideoSigma
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "henry_02.h"

double y_sigma, sigma_sigma;
// scalar PhiC[];
scalar c[], *stracers = {c};
// scalar f[];
char filename[80];

int main(int a, char const *arguments[])
{
  sprintf(filename, "%s", arguments[1]);
  restore(file = filename);
  boundary((scalar *){f});

  FILE *fp = fout;
  for (int i = 0; i < 1000; i++)
  {
    double csum = 0.0;
    int count = 0;
    for (int j = 0; j < 100; j++)
    {
      if (interpolate(f, j * 1.2e-2, i * 15e-3) > 1e-3)
      {
        csum += interpolate(c, j * 1.2e-2, i * 15e-3);
        count++;
      }
    }
    if (csum / count > 1e-6)
    {
      fprintf(ferr, "%g %g\n", i * 15e-3, csum / count);
    }
    else
    {
      fprintf(ferr, "%g %g\n", i * 15e-3, 0.0);
    }
  }
  fflush(fp);
  fclose(fp);
}