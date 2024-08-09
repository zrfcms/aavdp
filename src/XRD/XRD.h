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
#include "../MATH/GRAPH.h"

struct XRD_KNODE{
    int    k[3];
    double theta;
    double intensity;
    int    multiplicity=0;
    XRD_KNODE *next=nullptr;
};

//x-ray diffraction
class XRD
{
public:
    int    numk=0;
    XRD_KNODE  *khead=nullptr;
    XRD_KNODE  *ktail=nullptr;
    XRD(XMODEL *model, double min2Theta, double max2Theta, int lp_type, double spacing[3], bool is_spacing_auto);
    ~XRD();
    void   xrd(char *xrd_path, char *png_path);
    void   xrd(char *xrd_path, char *png_path, double peak_parameter, double scherrer_lambda, double scherrer_size, double dt);
private:
    double spacingK[3]={0.1, 0.1, 0.1};
    int    kmin[3]={10000, 10000, 10000}, kmax[3]={-10000, -10000, -10000};
    double minTheta=0.0, maxTheta=PI/2.0;
    double intensity_min=1.0e8, intensity_max=0.0;
    void   copy_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2);
    void   swap_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2);
    void   quick_sort(XRD_KNODE *kstart, XRD_KNODE *kend);
    void   quick_unique();
    void   add_k_node(int h, int k, int l, double theta, double intensity, int multiplicity);
    void   compute_diffraction_intensity(XMODEL *model, int lp_type);
};

#endif