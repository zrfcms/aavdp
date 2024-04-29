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
    SED(EMODEL *model, double spacing[3], double Kmagnitude_max);
    SED(EMODEL *model, double spacing[3], int zone[3], double thickness, double Kmagnitude_max);
    ~SED();
    void   vtk(const char *vtk_path);
    void   kkd(const char *kkd_path);
private:
    int    numk=0;
    int    **kvectors=nullptr;
    double **Kvectors=nullptr;
    double *Kintensity=nullptr;

    double radiusE=39.8406374501992;//ELAMBDA
    double spacingK[3]={0.1, 0.1, 0.1};
    int    kmax[3]={10000, 10000, 10000}, kmin[3]={-10000, -10000, -10000};
    void   count_diffraction_vectors(double Kmagnitude_max);
    void   count_diffraction_vectors(int zone[3], double thickness, double Kmagnitude_max);
    void   compute_diffraction_intensity(EMODEL *model, double Kmagnitude_max);
    void   compute_diffraction_intensity(EMODEL *model, int zone[3], double thickness, double Kmagnitude_max);
};

#endif