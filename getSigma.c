// get sigma vs y for VideoSigma
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "advdiff.h"
#include "curvature.h"


double y_sigma, sigma_sigma;
// scalar PhiC[];
char filename[80];
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  boundary((scalar *){f});
  
  scalar kappa[];
  curvature(f, kappa);

  FILE *fp = fout;

  foreach(){
    if (kappa[] != nodata){
        y_sigma = y;
        sigma_sigma = PhiC[];
        fprintf(ferr, "%g %g\n", y_sigma, sigma_sigma);
    }
  }
  fflush (fp);
  fclose (fp);
}