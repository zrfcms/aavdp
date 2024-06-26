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
// #include "../HDF5/HDF5.h"
using namespace std;

//electron diffraction
class SED
{
public:
    SED(EMODEL *model, double Kmagnitude_max, double spacing[3], bool is_spacing_auto);
    // SED(EMODEL *model, int zone[3], double thickness, double Kmagnitude_max, double spacing[3], bool is_spacing_auto);
    // SED(const char *hdf5_path);
    ~SED();
    double lambda;
    int    numk=0;
    int    **kvectors=nullptr;
    double **Kvectors=nullptr;
    double *Kintensity=nullptr;
    double intensity_min=1.0e8, intensity_max=0.0;
    // void   hdf5(const char *hdf5_path);
    void   vtk(const char *vtk_path, double threshold=1.0e-8);
    void   sed(const char *sed_path, double threshold);
    void   sed(const char *sed_path, int xaxis[3], int zaxis[3], double thickness, double threshold);
    friend class KKD;
private:
    double radiusE;
    double radiusK=1.70;
    double spacingK[3]={0.1, 0.1, 0.1};
    int    kmax[3]={-10000, -10000, -10000}, kmin[3]={10000, 10000, 10000};
    void   count_diffraction_vectors(EMODEL *model);
    // void   count_diffraction_vectors(EMODEL *model, int zone[3], double thickness);
    void   compute_diffraction_intensity(EMODEL *model);
    // void   compute_diffraction_intensity(EMODEL *model, int zone[3], double thickness);
};

#endif