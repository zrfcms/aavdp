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

struct KNODE{
    int    k[3];
    double theta;
    double intensity;
    int    multiplicity=0;
    KNODE *next=nullptr;
};

//x-ray diffraction
class XRD
{
public:
    int    numk=0;
    KNODE  *khead=nullptr;
    KNODE  *ktail=nullptr;
    XRD(MODEL *model, double spacing[3], double min2Theta, double max2Theta, bool is_lorentz, bool is_spacing_auto);
    XRD(MODEL *model, bool is_lorentz);
    ~XRD();
    void   xrd(const char *xrd_path);
    void   xrd(const char *xrd_path, int nbin);
private:
    bool   is_lorentz_flag=true;
    int    kmin[3]={10000, 10000, 10000}, kmax[3]={-10000, -10000, -10000};
    double spacingK[3]={0.1, 0.1, 0.1};
    double minTheta=0.0, maxTheta=PI/2.0;
    double intensity_min=1.0e8, intensity_max=0.0;
    void   copy_knode_data(KNODE *knode1, KNODE *knode2);
    void   swap_knode_data(KNODE *knode1, KNODE *knode2);
    void   quick_sort(KNODE *kstart, KNODE *kend);
    void   quick_unique();
    void   add_k_node(int h, int k, int l, double theta, double intensity, int multiplicity);
    void   compute_diffraction_intensity(MODEL *model, double spacingK[3]);
    void   compute_diffraction_intensity(MODEL *model);
};

#endif