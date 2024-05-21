#ifndef __AAVDP_KKD_H__
#define __AAVDP_KKD_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <ctime>
#include "../MATH/MATH.h"

using namespace std;

//kinematic Kikuchi diffraction
class KKD
{
public:
    KKD(const char *kkd_path, double threshold_ratio, double screenD, int screen_dpi, int screenW, int screenH);
    KKD(const char *kkd_path, int zone[3], double threshold_ratio, double screenD, int screen_dpi, int screenW, int screenH);
    ~KKD();
    void img(const char* img_path, char mode='w');
private:
    int    numpx=0, numpy=0;//ncol, nrow
    double ***screenG=nullptr;
    double intensity_min=1.0e8, intensity_max=0.0;
    double **screenI=nullptr;

    double radiusK=39.8406374501992;//ELAMBDA
    int    numk=0;
    double **Kvectors=nullptr;
    double *Kintensity=nullptr;
    void read_diffraction_intensity_from_kkd(const char *kkd_path);
    void compute_Kikuchi_sphere_projection(int screen_npx);
    void compute_Kikuchi_sphere_projection(int zone[3], int screen_npx, int screen_npy, double screenD, double screenW, double screenH);
    void compute_Kikuchi_intensity_projection(double threshold_ratio, double interceptK);
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