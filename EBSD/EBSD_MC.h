#ifndef _EBSD_MC_H_
#define _EBSD_MC_H_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "EBSD_MATH.h"
#include "EBSD_CONST.h"
#include "EBSD_HDF5.h"

using namespace std;

class RNG 
{
public:
    int ns=4;
    int state[4];

    void seed(int seed){
        state[0]=seed;
        for(int i=1; i<ns; i++){
            state[i]=default_seed[i];
        }
    }

    double uniform(){
        int imz=state[0]-state[3];
        if(imz<0){
            imz+=2147483579;
        }
        state[0]=state[1];
        state[1]=state[2];
        state[2]=state[3];
        state[3]=imz;
        state[4]=69069*state[4]+1013904243;
        imz+=state[4];

        return 0.5+0.23283064e-9*imz;
    }

private:
    int default_seed[4] = {521288629, 362436069, 16163801, 1131199299};
};

class EBSD_MC
{
public:
    //parameters required
    double omega;//RD sample tilt angle [in degree]
    double sigma;//TD sample tilt angle [in degree]
    double EkeV;//Energy of incident electron beam [in keV]
    double Ehistmin;//Minimum energy in energy histogram [in keV]
    double Ebin;
    double depthmax;//Maximum depth for which to keep track of exit energy statistics [in nm]
    double depthstep;//stepsize for depth-energy accumulator array [in nm]
    int    num_e;//Total number of electrons to try
    int    multiplier;
    int    nump;

    int    numEbin;//int((EkeV-Ehistmin)/Ebinsize)+1
    int    numzbin;//int(depthmax/depthstep)+1
    size_t size_SP[3];
    double ***accum_SP=nullptr;
    size_t size_E[3];
    int    ***accum_E=nullptr;//Energy accumulator array in the modified Lambert projection
    size_t size_z[4];
    int    ****accum_z=nullptr;//Depth accumulator array in the modified Lambert projection
    EBSD_MC(const char* file_path);
    ~EBSD_MC();

    double *Ebins=nullptr;
    double *depths=nullptr;
    int    izmax;
    void   compute_energy_distribution();
    void   compute_depth_distribution();
    double Emin;         //Ehistmin-Ebin/2.0
    int    numpx;        //Number of scintillator points along x
    int    numpy;        //Number of scintillator points along y, equal to numsx
    int    num_E, num_z;
    int    count_E, count_z;
    //Material parameters
    double ave_Z=13.0;        //Average atomic number
    double ave_M=26.98154;        //Average atomic mass in g/mol
    double density=2.69781;      //Density in g/cm^3
private:
    //Monte Carlo related parameters
    RNG rng;
    double R;
    void   run(int seed);
    double compute_alphainfreepath(double E);
    double compute_freepath(double E, double alpha);
    double compute_energyloss(double E, double step);
    bool read_parameters_from_nml(const char *file_path);
    bool read_parameters_from_hdf5(const char *file_path);
    bool write_parameters_into_hdf5(const char *file_path);
};

#endif