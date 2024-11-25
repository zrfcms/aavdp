#ifndef _AAVDP_MC_H_
#define _AAVDP_MC_H_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
#include "../MODEL/MODEL.h"

#define SCATTERING_EVENT_NUMBER 300 //Number of scattering events along one trajectory

using namespace std;
extern void compute_Lambert_Projection(double xy[2], int &ierr, double xyz[3]);

class RNG 
{
public:
    RNG(int nseed=100);
    ~RNG();
    void   seed(int id);
    double random();
private:
    const char *rand_path="RandomSeeds.data";
    int nseed=0;
    int **default_seeds=nullptr;
    int ns=4;
    int state[4]={0};
};

class DED_MC
{
public:
    int    nump;
    int    numzbin;//int(depthmax/depthstep)+1
    int    numpz;
    int    **accum_E=nullptr;//Energy accumulator array in the modified Lambert projection
    int    ***accum_z=nullptr;//Depth accumulator array in the modified Lambert projection

    int    izmax=0;
    double dz, z;
    double *weight=nullptr;
    DED_MC(CELL *cell, double omega, double sigma, double Emax, double Emin, double depthmax, double depthstep, int nume, int nump);
    ~DED_MC();
    void   img(char *img_path, double dimension=6.0, int resolution=512);
private:
    RNG   *rng;
    double ave_M, ave_Z, density;
    void   compute(int &count_bse, int ne, double dir0[3], double E0, double Ex);
    void   update_free_path(double &step, double &alpha, double E, double E0);
    void   update_incident_energy(double &E, double step);
    void   update_incident_direction(double dir[3], double alpha);
    void   update_incident_coordinate(double xyz[3], double dir[3], double step);
};

#endif