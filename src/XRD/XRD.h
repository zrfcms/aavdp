#ifndef __AAVDP_XRD_H__
#define __AAVDP_XRD_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
#include "../MODEL/MODEL.h"

#define DEG_TO_RAD_HALF 0.008726646259971648
#define RAD_TO_DEG_TWO 114.59155902616465
#define XRD_INTENSITY_LIMIT 1.0e-6
#define XRD_THETA_LIMIT 1.0e-3
struct XRD_KNODE{
    int    hkl[3];
    double theta;//radian
    double intensity;
    int    multiplicity;
    XRD_KNODE *next=nullptr;
};

extern void copy_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2);
extern void swap_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2);
extern void quick_sort(XRD_KNODE *kstart, XRD_KNODE *kend);

#define SCHERRER_CONST 0.90
extern void pseudo_Voigt(double *y, double *x, int num, double x0, double y0, double eta, double w);

//x-ray diffraction
class XRD
{
public:
    int    numk=0;
    XRD_KNODE  *khead=nullptr;
    XRD_KNODE  *ktail=nullptr;
    double minTheta, maxTheta;
    double intensity_min=1.0e8, intensity_max=0.0;
    XRD(MODEL *model, double min2Theta, double max2Theta, double threshold, double spacing[3], bool is_spacing_auto);
    ~XRD();
    void   xrd(char *xrd_path);
    void   xrd(char *xrd_path, double mixing_param, double scherrer_lambda, double scherrer_diameter, double bin2Theta);
    void   xrd(char *xrd_path, double mixing_param, double FWHM, double bin2Theta);
private:
    void   add_k_node(int hkl[3], double theta, double intensity, int multiplicity);
    void   compute_diffraction_intensity(MODEL *model, double spacing[3], bool is_spacing_auto);
    void   filter_diffraction_intensity(double threshold);
    void   unique_diffraction_intensity();
    void   img(char *png_path, double *x, double *y, int num, double xmin, double xmax, char mode);
};

#endif