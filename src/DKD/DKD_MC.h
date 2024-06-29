#ifndef _AAVDP_DKD_MC_H_
#define _AAVDP_DKD_MC_H_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../MATH/MATH.h"
// #include "../HDF5/HDF5.h"
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

class DKD_MC
{
public:
    double omega=0.0;//sample rotate angle [in degree]
    double sigma=75.7;//sample tilt angle [in degree]
    double EkeV;//Energy of incident electron beam, i.e., accelerating voltage [in keV]
    double Emin;//Minimum energy to consider [in keV]
    double Ebin;
    double depthmax;//Maximum depth for which to keep track of exit energy statistics [in nm]
    double depthstep;//stepsize for depth-energy accumulator array [in nm]
    int    num_e;//Total number of electrons to try
    int    nump;
    int    numpz;

    int    numEbin;//int((EkeV-Emin)/Ebinsize), excluding Emin
    double *Ebins=nullptr;
    int    numzbin;//int(depthmax/depthstep)+1
    int    izmax;
    double *depths;
    int    ***accum_E=nullptr;//Energy accumulator array in the modified Lambert projection
    int    ****accum_z=nullptr;//Depth accumulator array in the modified Lambert projection
    DKD_MC(CELL *cell, double omega, double sigma, double Emax, double Emin, double Ebin, double zmax, double zstep, int num_e, int nump);
    // DKD_MC(const char* hdf5_path);
    ~DKD_MC();
    // void   hdf5(const char *hdf5_path);
    void   img(char *img_path, double dimension=6.0, int resolution=512);
private:
    RNG   *rng;
    const int nbatch=100;
    double ave_Z, ave_M, density;
    void   compute(int &count_bse, int &count_e, int ne, int id);
    void   update_free_path(double &step, double &alpha, double E, double E0);
    void   update_incident_energy(double &E, double step);
    void   update_incident_direction(double dir[3], double alpha);
    void   update_incident_coordinate(double xyz[3], double dir[3], double step);
    void   compute_depth_distribution();
};

#endif