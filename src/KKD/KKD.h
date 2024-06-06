#ifndef __AAVDP_KKD_H__
#define __AAVDP_KKD_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <ctime>
#include "../SED/SED.h"
#include "../MATH/MATH.h"

using namespace std;

//kinematic Kikuchi diffraction
class KKD
{
public:
    KKD(const char *hdf5_path, double threshold, double screenD, int screenW, int screen_dpi);
    KKD(const char *hdf5_path, int zone[3], double threshold, double screenD, int screenW, int screenH, int screen_dpi);
    ~KKD();
    void img(const char* img_path, char mode='w');
private:
    int    numpx=0, numpy=0;//ncol, nrow
    double ***screenG=nullptr;
    double **screenI=nullptr;
    double intensity_min=1.0e8, intensity_max=0.0;
    void compute_Kikuchi_sphere_projection(SED *sed, int screenW, int screen_dpi);
    void compute_Kikuchi_sphere_projection(SED *sed, int zone[3], double screenD, double screenW, double screenH, int screen_dpi);
    void compute_Kikuchi_intensity_projection(SED *sed, double threshold, double screenD, int screen_dpi);
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