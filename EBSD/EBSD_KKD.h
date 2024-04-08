#ifndef __EBSD_KKD_H__
#define __EBSD_KKD_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include "EBSD_CELL.h"
using namespace std;

struct QUATERNION
{
    double c1, c2, c3, c4;
};

extern QUATERNION quate_conjg(QUATERNION q);
extern QUATERNION quate_multi(QUATERNION q1, QUATERNION q2);
extern void quate_rotate(double r_v[3], double v[3], QUATERNION q);
extern QUATERNION quate_convert(double v[3], double angle);

//kinematic EBSD pattern
class EBSD_KKD
{
public:
    EBSD_KKD(const char *file_path);
    ~EBSD_KKD();
    void   compute_master_pattern(const char* file_path);
    void   test_images(const char *file_path);
private:
    double dmin;
    double ci;//intensity cutoff
    double voltage;
    int    nump=500;
    int    numphi=360*16;
    int    npx, npy;
    double **mLPNH=nullptr, **mLPSH=nullptr;
    double **mSPNH=nullptr, **mSPSH=nullptr;
    bool   read_parameters_from_nml(const char *file_path);
    bool   write_parameters_into_hdf5(const char* file_path);
    void   apply_anti_aliasing(double **master, double xy[2], int ix, int iy, double intensity);
    double get_Lambert_interpolation(double xyz[3], double **mat, bool hexagonal_flag=false);
    void   create_images(const char* png_path, double **mat);
};

#endif