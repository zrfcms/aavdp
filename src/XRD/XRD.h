#ifndef __AAVDP_XRD_H__
#define __AAVDP_XRD_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"
#define ZERO_LIMIT 0.000001
using namespace std;

//x-ray diffraction
class XRD
{
public:
    int    numk=0;
    int    **kvector=nullptr;
    double *Ktheta=nullptr;
    complex<double> *Kfactor=nullptr;
    double *Kintensity=nullptr;
    XRD(MODEL *model, double spacingK[3], double min2Theta, double max2Theta, bool is_lorentz, bool is_spacing_auto);
    XRD(MODEL *model, bool is_lorentz);
    XRD(const char *xrd_path);
    ~XRD();
    void   xrd(const char *xrd_path);
private:
    bool   is_lorentz_flag;
    int    kmin[3], kmax[3];
    double intensity_min=1.0e8, intensity_max=0.0;
    void   unique();
    void   count_diffraction_vector(MODEL *model, double spacingK[3], double minTheta, double maxTheta);
    void   compute_diffraction_intensity(MODEL *model, double spacingK[3], double minTheta, double maxTheta);
    void   count_diffraction_vector(MODEL *model);
    void   compute_diffraction_intensity(MODEL *model);
};

#endif