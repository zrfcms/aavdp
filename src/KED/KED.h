#ifndef __AAVDP_KED_H__
#define __AAVDP_KED_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include "../MODEL/MODEL.h"
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
using namespace std;

struct KED_KNODE{
    int    k[3];
    double K[3];
    double Kmag;
    double intensity;
    KED_KNODE *next=nullptr;
};

//electron diffraction
class KED
{
public:
    KED(EMODEL *model, double Kmagnitude_max, double spacing[3], bool is_spacing_auto);
    KED(EMODEL *model, int zone[3], double thickness, double Kmagnitude_max, double spacing[3], bool is_spacing_auto);
    ~KED();
    int    numk=0;
    KED_KNODE  *khead=nullptr;
    KED_KNODE  *ktail=nullptr;
    void   vtk(const char *vtk_path, double threshold);
    void   ked(const char *ked_path, double threshold);
    void   ked(const char *ked_path, char *png_path, int xaxis[3], int yaxis[3], int zaxis[3], double threshold);
    void   restart(const char *restart_path);
private:
    double radiusE;
    double radiusK=1.70;
    double spacingK[3]={0.1, 0.1, 0.1};
    int    kmax[3]={-10000, -10000, -10000}, kmin[3]={10000, 10000, 10000};
    double intensity_min=1.0e8, intensity_max=0.0;
    void   add_k_node(int h, int k, int l, double K[3], double Kmag, double intensity);
    void   compute_diffraction_intensity(EMODEL *model);
    void   compute_diffraction_intensity(EMODEL *model, int zone[3], double thickness);
};

#endif