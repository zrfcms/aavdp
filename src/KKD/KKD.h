#ifndef __AAVDP_KKD_H__
#define __AAVDP_KKD_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include "../MODEL/MODEL.h"
#include "../SED/SED.h"
#include "../MATH.h"
using namespace std;

//kinematic Kikuchi diffraction
class KKD
{
public:
    KKD(int zone[3], double threshold_ratio, double screenD, int screen_dpi, int screenW, int screenH);
    ~KKD();
private:
    double radiusK=0.0;
    int    numpx=0, numpy=0;
    double ***screenG=nullptr;
    double **screenI=nullptr;
    void compute_Kikuchi_sphere_projection(int screen_npx);
    void compute_Kikuchi_sphere_projection(int zone[3], int screen_npx, int screen_npy, double screenD, double screenW, double screenH);
    void compute_Kikuchi_intensity_projection(SED *sed, double threshold_ratio, double interceptK);
};

// class KKD
// {
// public:
//     KKD();
//     KKD(const char *file_path);
//     ~KKD();
//     void   compute_master_pattern(const char* file_path);
//     void   test_images(const char *file_path);
// private:
//     double lambda;
//     double **screenX=nullptr;
//     double **screenY=nullptr;
//     double **screenZ=nullptr;


//     double threshold;//intensity threshold
//     int    numKc=360*16;

//     double lambda=0.0;//excitation voltage [in keV], corresponding to a wavelength 0.0251 [in A]


//     double dmin;
//     int    nump=500;
//     int    npx, npy;
//     double **mLPNH=nullptr, **mLPSH=nullptr;
//     double **mSPNH=nullptr, **mSPSH=nullptr;
//     bool   read_parameters_from_nml(const char *file_path);
//     bool   write_parameters_into_hdf5(const char* file_path);
//     void   apply_anti_aliasing(double **master, double xy[2], int ix, int iy, double intensity);
//     double get_Lambert_interpolation(double xyz[3], double **mat, bool hexagonal_flag=false);
//     void   create_images(const char* png_path, double **mat);
// };

#endif