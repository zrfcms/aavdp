#ifndef _AAVDP_DED_H_
#define _AAVDP_DED_H_
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "../include/lapacke.h"
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
#include "../MODEL/MODEL.h"
#include "DED_MC.h"

#define DED_KMAG_LIMIT 1.0e-3
#define DED_INTENSITY_LIMIT 1.0e-6

using namespace std;

struct BETHE{
    double c1=4.0, c2=8.0, c3=50.0; // rg
    double c_sg=1.0;
};

struct DED_GNODE{
    double hkl[3];//Miller indices
    complex<double> Ug, qg;
    double sg;//excitation error
    bool   is_double_diffrac;
    bool   is_weak, is_strong;
    DED_GNODE *next=nullptr;
    DED_GNODE *nextw=nullptr;//weak
    DED_GNODE *nexts=nullptr;//strong
};

class DED_GVECTOR
{
public:
    DED_GVECTOR(CELL *cell, BETHE *bethe, double kk[3], double fn[3], double kn, double dmin);
    ~DED_GVECTOR();
    int    numg=0;
    int    nstrong=0, nweak=0;
    DED_GNODE *ghead=nullptr;
    DED_GNODE *heads=nullptr, *headw=nullptr;
private:
    DED_GNODE *gtail=nullptr;
    void   add_g_node(double hkl[3], complex<double> Ug, complex<double> qg, double sg, bool is_double_diffrac);
};

struct DED_KNODE{
    double hkl[3];
    double K[3];
    double Kmagnitude;
    double intensity;
    DED_KNODE *next=nullptr;
};

extern void copy_knode_data(DED_KNODE *knode1, DED_KNODE *knode2);
extern void swap_knode_data(DED_KNODE *knode1, DED_KNODE *knode2);
extern void quick_sort(DED_KNODE *kstart, DED_KNODE *kend);

class DED
{
public:
    int    numk=0;
    DED_KNODE  *khead=nullptr;
    DED_KNODE  *ktail=nullptr;
    double Kmagnitude_max=0.1;//nm
    double intensity_min=1.0e8, intensity_max=0.0;
    DED_KNODE *knearest_1=nullptr, *knearest_2=nullptr;
    DED(CELL *cell, BETHE *bethe, int zone[3], int fnorm[3], double fthick, double voltage, double Kmag_max, double threshold);
    ~DED();
    void   rotate(double x[3], double y[3]);
    void   ded(char *ded_path);
    void   ded(char *ded_path, double sigma, double dx);
private:
    double axes[3][3]={0.0};
    void   add_k_node(CELL *cell, double hkl[3], double intensity);
    void   compute_dynamic_matrix(complex<double> **dmat, CELL *cell, DED_GVECTOR *gvec);
    void   compute_diffraction_intensity(double *INTENS, complex<double> **DMAT, double Z, double KN, int NS);
    void   filter_diffraction_intensity(double threshold);
    void   find_first_and_second_knearests();
    void   rotate_by_first_knearest(CELL *cell, int zone[3]);
    void   img(char *png_path, double *x, double *y, double *value, int num, double limit);
    void   img(char *png_path, double **value, int nump);
};

struct DKD_KNODE{
    double k[3];
    double kn;
    int    i, j;
    int    hemisphere;//Northern = 1, Southern = -1
    double intensity=0.0;
    DKD_KNODE *next=nullptr;
};

class DKD_KVECTOR
{
public:
    DKD_KVECTOR(CELL *cell, int impx, int impy);//Rosca-Lambert
    ~DKD_KVECTOR();
    int    numk=0;
    DKD_KNODE *khead=nullptr;
    DKD_KNODE *ktail=nullptr;
    double intensity_min=1.0e8, intensity_max=0.0;
    void   add_k_intensity(DKD_KNODE *knode, complex<double> **Lgh, complex<double> ***Sgh, int nstrong, int napos, int npos);
private:
    double deltax, deltay;
    void   add_k_node(CELL *cell, int i, int j, double kn, bool southern_flag=false);
};

class DKD
{
public:
    int    numpx=0, numpy=0;
    int    impx, impy;
    double **mLPNH=nullptr, **mLPSH=nullptr; //The modified Lambert Projection Northern/Southern Hemisphere
    double **mSPNH=nullptr, **mSPSH=nullptr; //The master Stereographic Projection Northern/Southern Hemisphere
    DKD(CELL *cell, DED_MC *mc, BETHE *bethe, double voltage, double dmin, int nump);
    //DKD(CELL *cell, BETHE *bethe, double voltage, double fthick, double dmin);
    ~DKD();
    void   img(char *LPNH_path, char *LPSH_path, char *SPNH_path, char *SPSH_path, double dimension, int resolution);
private:
    void   compute_dynamic_matrix(complex<double> **dynmat, CELL *cell, DED_GVECTOR *gvec);
    void   compute_Sgh_matrices(complex<double> ***Sgh, CELL *cell, DED_GVECTOR *gvec);
    void   compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double *EWF, int IZMAX, double Z, double DZ, double KN, int NS);
    void   compute_Lambert_projection(CELL *cell, DKD_KVECTOR *kvec);
    void   compute_stereographic_projection(bool use_hexagonal);
};

#endif