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
#define CONST_IN_AREA_INV 1.5273987E19 //1.0/(5.21E-21*(4*PI))

using namespace std;

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

class MC
{
public:
    int    numpE=0;
    int    **accum_E=nullptr;//Energy accumulator array in the modified Lambert projection
    int    numzbin=0, numpz=0;//int(depthmax/depthstep)+1
    int    ***accum_z=nullptr;//Depth accumulator array in the modified Lambert projection

    int    izmax=0;
    double dz, z;
    double EkeV;
    double *weight=nullptr;
    MC(CELL *cell, double omega, double sigma, double Emax, double Emin, double depthmax, double depthstep, int nume, int nump);
    ~MC();
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