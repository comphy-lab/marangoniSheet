// get film thickness h

#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80], nameTrack[80];
scalar * list = NULL;
scalar f[];

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  restore (file = filename);
  boundary((scalar *){f, u.x, u.y});

  double xmax = -HUGE;
  double y_xmax = 0;
  
  foreach(){
    if (x > xmax && f[] > 1-1e-3 && y < 0.05)
      {
        xmax = x;
        y_xmax = y;
      }
  }

  FILE * fp = ferr;
  fprintf(ferr, "%f %7.6e %7.6e\n", t,  xmax, y_xmax);
  fflush (fp);
  fclose (fp);
}