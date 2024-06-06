#ifndef _AAVDP_DKD_MC_H_
#define _AAVDP_DKD_MC_H_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../MATH/MATH.h"
#include "../HDF5/HDF5.h"
#include "../MODEL/MODEL.h"

using namespace std;

class RNG 
{
public:
    int ns=4;
    int state[4];

    void seed(int seed){
        state[0]=seed;
        for(int i=1;i<ns;i++){
            state[i]=default_seed[i];
        }
    }

    double uniform(){
        int imz=state[0]-state[2];
        if(imz<0){
            imz+=2147483579;
        }
        state[0]=state[1];
        state[1]=state[2];
        state[2]=imz;
        state[3]=69069*state[3]+1013904243;
        imz+=state[3];
        return 0.5+0.23283064e-9*imz;
    }
private:
    int default_seed[4] = {521288629, 362436069, 16163801, 1131199299};
};

class DKD_MC
{
public:
    double omega;//sample rotate angle [in degree]
    double sigma;//sample tilt angle [in degree]
    double EkeV;//Energy of incident electron beam, i.e., accelerating voltage [in keV]
    double Ehistmin;//Minimum energy to consider in energy histogram [in keV]
    double Ebin;
    double depthmax;//Maximum depth for which to keep track of exit energy statistics [in nm]
    double depthstep;//stepsize for depth-energy accumulator array [in nm]
    int    num_e;//Total number of electrons to try
    int    nump;
    int    numpz;

    int    numEbin;//int((EkeV-Ehistmin)/Ebinsize)+1
    int    numzbin;//int(depthmax/depthstep)+1
    int    ***accum_E=nullptr;//Energy accumulator array in the modified Lambert projection
    int    ****accum_z=nullptr;//Depth accumulator array in the modified Lambert projection
    DKD_MC(const char *hdf5_path, double omega, double sigma, double Emax, double Emin, double Ebin, double zmax, double zstep, int num_e, int nump);
    DKD_MC(const char* hdf5_path);
    ~DKD_MC();
    void   hdf5(const char *hdf5_path);
    void   img(const char *img_path, double dimension=6.0, int resolution=256);

    int    izmax;
    double *Ebins, *depths;
    void   set_energies_and_depths();
private:
    RNG rng;
    int    prime_seed=932117;
    int    multiplier=8;
    double ave_Z, ave_M, density;
    void   compute(int &count, int seed=932117);
    void   update_free_path(double &step, double &alpha, double E);
    void   update_incident_energy(double &E, double step);
    void   update_incident_direction(double dir[3], double alpha);
    void   update_incident_coordinate(double xyz[3], double dir[3], double step);
};

#endif