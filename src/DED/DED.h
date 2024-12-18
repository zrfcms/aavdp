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
#include "./MC.h"

#define DED_KMAGNITUDE_LIMIT 1.0e-3
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
};

struct DKD_KNODE{
    double k[3];
    double kn;
    double intensity;
    int i, j;
    int hemisphere;//Northern = 1, Southern = -1
    DKD_KNODE *next=nullptr;
};

class DKD_KVECTOR
{
public:
    DKD_KVECTOR(CELL *cell, int npx, int npy);//Rosca-Lambert
    ~DKD_KVECTOR();
    int    numk=0;
    DKD_KNODE *khead=nullptr;
    DKD_KNODE *ktail=nullptr;
private:
    void   add_k_node(CELL *cell, double xy[2], double kn, int i=0, int j=0, bool southern_flag=false);
};

class DKD
{
public:
    int    numpx=0, numpy=0;
    double ***screenK0=nullptr;
    double **screenI=nullptr;
    double intensity_max=0.0, intensity_min=1.0e8;

    double **screenNI=nullptr, **screenSI=nullptr;
    double intensity_maxN=0.0, intensity_minN=1.0e8;
    double intensity_maxS=0.0, intensity_minS=1.0e8;
    DKD(CELL *cell, MC *mc, BETHE *bethe, double xaxis[3], double yaxis[3], double zaxis[3], double ratiox, double ratioy, double Kmag_max, int npx, int npy, char *projection);
    DKD(CELL *cell, BETHE *bethe, double xaxis[3], double yaxis[3], double zaxis[3], double ratiox, double ratioy, double voltage, double fthick, double Kmag_max, int npx, int npy, char *projection);
    DKD(CELL *cell, MC *mc, BETHE *bethe, double Kmag_max, char *projection);
    DKD(CELL *cell, BETHE *bethe, double voltage, double fthick, double Kmag_max, int nump, char *projection);
    ~DKD();
    void   dkd(char* dkd_path, char background);
    void   dkd(char* dkd_path, double vmax, double vmin, char background);
private:
    void   compute_Kikuchi_sphere_projection(CELL *cell, double xaxis[3], double yaxis[3], double zaxis[3], double ratiox, double ratioy, double kn, char *projection);
    void   compute_Kikuchi_intensity_projection(double &intens, complex<double> ***Sgh, complex<double> **Lgh, int napos, int nstrong, int npos);
    void   compute_dynamic_matrix(complex<double> **dynmat, CELL *cell, DED_GVECTOR *gvec);
    void   compute_Sgh_matrices(complex<double> ***Sgh, CELL *cell, DED_GVECTOR *gvec);
    void   compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double *EWF, int IZMAX, double Z, double DZ, double KN, int NS);
    void   compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double Z, double KN, int NS);
    void   compute_kvector_projection(CELL *cell, DKD_KVECTOR *kvec, char *projection, bool use_hexagonal);
};

#endif