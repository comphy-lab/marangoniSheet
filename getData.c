// get data

#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay, Oh = 0.;
scalar f[], D2c[], vel[];
scalar *list = NULL;

int main(int a, char const *arguments[])
{
    sprintf(filename, "%s", arguments[1]);
    xmin = atof(arguments[2]);
    ymin = atof(arguments[3]);
    xmax = atof(arguments[4]);
    ymax = atof(arguments[5]);
    ny = atoi(arguments[6]);
    Oh = atof(arguments[7]);

    list = list_add(list, D2c);
    list = list_add(list, vel);
    list = list_add(list, u.x);
    list = list_add(list, u.y);

    // boundary conditions
    u.n[right] = neumann(0.);
    p[right] = dirichlet(0.);

    restore(file = filename);
    f.prolongation = fraction_refine;
    boundary((scalar *){f, u.x, u.y});

    foreach ()
    {
        double D11 = (u.y[0, 1] - u.y[0, -1]) / (2 * Delta);
        double D22 = (u.y[] / y);
        double D33 = (u.x[1, 0] - u.x[-1, 0]) / (2 * Delta);
        double D13 = 0.5 * ((u.y[1, 0] - u.y[-1, 0] + u.x[0, 1] - u.x[0, -1]) / (2 * Delta));
        double D2 = (sq(D11) + sq(D22) + sq(D33) + 2.0 * sq(D13));
        D2c[] = 2 * (clamp(f[], 0., 1.) * (Oh - 2e-2 * Oh) + 2e-2 * Oh) * D2;

        D2c[] = (D2c[] > 0.) ? log(D2c[]) / log(10) : -10;
        vel[] = sqrt(sq(u.x[]) + sq(u.y[]));
    }
    boundary((scalar *){D2c, vel});

    FILE *fp = ferr;
    Deltay = (double)(ymax - ymin) / (ny);
    nx = (int)(xmax - xmin) / Deltay;
    Deltax = (double)(xmax - xmin) / (nx);
    len = list_len(list);

    double **field = (double **)matrix_new(nx, ny + 1, len * sizeof(double));
    for (int i = 0; i < nx; i++)
    {
        double x = Deltax * (i + 0.5) + xmin;
        for (int j = 0; j < ny; j++)
        {
            double y = Deltay * (j + 0.5) + ymin;
            int k = 0;
            for (scalar s in list)
            {
                field[i][len * j + k++] = interpolate(s, x, y);
            }
        }
    }

    for (int i = 0; i < nx; i++)
    {
        double x = Deltax * (i + 0.5) + xmin;
        for (int j = 0; j < ny; j++)
        {
            double y = Deltay * (j + 0.5) + ymin;
            fprintf(fp, "%g %g", x, y);
            int k = 0;
            for (scalar s in list)
            {
                fprintf(fp, " %g", field[i][len * j + k++]);
            }
            fputc('\n', fp);
        }
    }
    fflush(fp);
    fclose(fp);
    matrix_free(field);
}