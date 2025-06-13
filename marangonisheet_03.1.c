// Sheet breakup by surface tension gradient || Initialized with exponential shape disturbance || Use henry_02.h for finite Peclet number
// Correct boundary conditions
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "henry_02.h"

#define tsnap (5e-2) // 0.001 only for some cases.

// Error tolerancs
#define fErr (1e-3)   // error tolerance in f1 VOF
#define KErr (1e-6)   // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-3) // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J

int MAXlevel;
scalar c[], * stracers = {c};
double GammaR, Oh, D, D_air, alpha_inv, tmax, w0;
char nameOut[80], dumpFile[80];

// PhiC[top] = neumann(0.0);
uf.n[top] = dirichlet(0.);
u.n[top] = dirichlet(0.);
uf.t[top] = neumann(0.);
u.t[top] = neumann(0.);

uf.n[left] = dirichlet(0.);
uf.t[left] = neumann(0.);

u.t[left] = neumann(0.);
u.n[left] = dirichlet(0.);
d[top] = neumann(0.);

scalar sigmaf[];
int main(int argc, char const *argv[])
{
    dtmax = 1e-5;
    MAXlevel = atoi(argv[1]);
    Oh = atof(argv[2]);
    GammaR = atof(argv[3]);
    D = atof(argv[4]);
    D_air = atof(argv[5]);
    alpha_inv = atof(argv[6]);
    tmax = atof(argv[7]);
    if (argc < 8)
    {
        fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 6 - argc);
        return 1;
    }
    fprintf(ferr, "Level %d, Oh %2.1e, GammaR %4.3f, D %4.3f\n", MAXlevel, Oh, GammaR, D);

    init_grid(1 << MAXlevel);
    // Create intermediate for all snapshots.
    char comm[80];
    sprintf(comm, "mkdir -p intermediate");
    system(comm);
    // Name of the restart file, used in writingFiles event.
    sprintf(dumpFile, "dump");

    L0 = 15.0, w0 = 0.5;
    // Diff_C1 = D, Diff_C2 = D;
    c.D1 = D, c.D2 = D_air; 
    c.alpha_inv = alpha_inv;
    rho1 = 1., rho2 = 0.001;
    mu1 = Oh, mu2 = 2e-5;

    TOLERANCE = 1e-6;

    d.sigmaf = sigmaf;

    run();
}

event init(i = 0)
{
    if (!restore(file = dumpFile))
    {
        foreach ()
        {
             d[] = -(x  - 1.00001);
        }

        foreach ()
        {
            c[] = d[] > 0 ? (1/(w0*sqrt(2*pi)))*exp(-y*y/(2*sq(w0))) : 0 ; // Gaussian profile with 
            sigmaf[] = max(1 - GammaR * c[], 1e-3);
        }
    }
}

event properties(i++)
{
    foreach ()
        sigmaf[] = max(1 - GammaR * c[], 1e-3);
}

event logWriting(i++)
{
    double ke = 0.;
    foreach (reduction(+ : ke))
    {
        ke += (2 * pi * y) * (0.5 * rho(f[]) * (sq(u.x[]) + sq(u.y[]))) * sq(Delta);
    }
    static FILE *fp;
    if (pid() == 0)
    {
        if (i == 0)
        {
            fprintf(ferr, "i dt t ke\n");
            fp = fopen("log", "w");
            fprintf(fp, "Level %d, Oh %2.1e, GammaR %4.3f, D %4.3f\n", MAXlevel, Oh, GammaR, D);
            fprintf(fp, "i dt t ke\n");
        }
        else
        {
            fp = fopen("log", "a");
            fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
        }
        fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
        fclose(fp);
        fprintf(ferr, "%d %g %g %g\n", i, dt, t, ke);
    }

    assert(ke > -1e-10);
    assert(ke < 1e7);

    if ((ke > 1e7 || ke < -1e-10) && i > 1e1 && pid() == 0)
    {
        const char *message = ke > 1e7 ? "The kinetic energy blew up. Stopping simulation\n"
                                       : "kinetic energy too small now! Stopping!\n";
        fprintf(ferr, "%s", message);
        fp = fopen("log", "a");
        fprintf(fp, "%s", message);
        fclose(fp);
        dump(file = dumpFile);
        return 1;
    }
}

event adapt(i++)
{
    adapt_wavelet((scalar *){f, u.x, u.y},
                  (double[]){fErr, VelErr, VelErr}, MAXlevel, MAXlevel - 3);
}

event writingFiles(t = 0; t += tsnap; t <= tmax)
{
    dump(file = dumpFile);
    sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
    dump(file = nameOut);
}

event end(t = tmax)
{
}
