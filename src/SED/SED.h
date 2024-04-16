#ifndef __AAVDP_SED_H__
#define __AAVDP_SED_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"
using namespace std;

//electron diffraction
class SED
{
public:
    double radiusE=0.0;
    double spacingK[3]={1.0, 1.0, 1.0};
    int    kmax[3]={10000, 10000, 10000}, kmin[3]={-10000, -10000, -10000};
    double intensity_max=0.0, intensity_min=1.0e8;

    int    numk=0;
    int    **kvectors=nullptr;
    double **Kvectors=nullptr;
    double *Kintensity=nullptr;
    SED(EMODEL *model, double spacing[3], double Kmagnitude_max);
    SED(EMODEL *model, double spacing[3], int zone[3], double intercept_radius, double Kmagnitude_max);
    ~SED();
    void   result(const char *vtk_path, const char *text_path);
private:
    void   count_diffraction_vectors(double Kmagnitude_max);
    void   count_diffraction_vectors(int zone[3], double intercept_radiusE, double Kmagnitude_max);
    void   compute_diffraction_intensity(EMODEL *model, double Kmagnitude_max);
    void   compute_diffraction_intensity(EMODEL *model, int zone[3], double intercept_radiusE, double Kmagnitude_max);
    void   vtk(const char *vtk_path);
    void   text(const char *text_path);
};

#endif