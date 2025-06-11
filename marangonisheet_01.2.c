// Sheet breakup by surface tension gradient || Initialized with gaussian surface tension profile
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
// #include "navier-stokes/conserving.h"
#include "integral.h"

#define tsnap (5e-2) // 0.001 only for some cases.

// Error tolerancs
#define fErr (1e-3)   // error tolerance in f1 VOF
#define KErr (1e-6)   // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-3) // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J

int MAXlevel;
double GammaR, Oh, tmax, w0;
char nameOut[80], dumpFile[80];

uf.n[bottom] = dirichlet(0.);
uf.n[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.n[top] = dirichlet(0.);

uf.t[bottom] = neumann(0.);
uf.t[top] = neumann(0.);

u.t[bottom] = neumann(0.);
u.t[top] = neumann(0.);

uf.n[left] = dirichlet(0.);
uf.t[left] = neumann(0.);

u.t[left] = neumann(0.);
u.n[left] = dirichlet(0.);

d[bottom] = neumann(0.);

scalar sigmaf[];
scalar PhiC[];
int main(int argc, char const *argv[])
{
    dtmax = 1e-5;
    MAXlevel = atoi(argv[1]);
    Oh = atof(argv[2]);
    GammaR = atof(argv[3]);
    tmax = atof(argv[4]);
    if (argc < 5)
    {
        fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 5 - argc);
        return 1;
    }
    fprintf(ferr, "Level %d, Oh %2.1e, GammaR %4.3f\n", MAXlevel, Oh, GammaR);

    init_grid(1 << MAXlevel);
    // Create intermediate for all snapshots.
    char comm[80];
    sprintf(comm, "mkdir -p intermediate");
    system(comm);
    // Name of the restart file, used in writingFiles event.
    sprintf(dumpFile, "dump");

    L0 = 15.0;
    rho1 = 1.;
    mu1 = Oh;
    rho2 = 0.001;
    mu2 = 2e-5;
    w0 = 0.5;
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
            d[] = -( x - 1 + ((double)L0) / ((double)N) / 2. + 0.1 * exp(-y * y / 0.1));
            // d[] = (y+((double)L0)/((double)N)/2.);
            // d[] = -x;
        }

        foreach ()
        {   
            PhiC[] = (1/(w0*sqrt(2*pi)))*exp(-y*y/(2*sq(w0)));
            sigmaf[] = max(1 - GammaR * PhiC[], 1e-3);
        }
    }
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
            fprintf(fp, "Level %d, Oh %2.1e, GammaR %4.3f\n", MAXlevel, Oh, GammaR);
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

    if ((ke > 1e7 || ke < 1e-6) && i > 1e1 && pid() == 0)
    {
        const char *message = ke > 1e3 ? "The kinetic energy blew up. Stopping simulation\n"
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
                  (double[]){fErr, VelErr, VelErr}, MAXlevel, MAXlevel - 5);
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
