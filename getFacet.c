// get facets
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f[];
char filename[80];
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  #if TREE
    f.prolongation = fraction_refine;
  #endif
  boundary((scalar *){f});

  FILE * fp = ferr;
  output_facets(f,fp);
  fflush (fp);
  fclose (fp);
}