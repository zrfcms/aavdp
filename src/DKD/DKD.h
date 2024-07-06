#ifndef _AAVDP_DKD_H_
#define _AAVDP_DKD_H_
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
#include "DKD_MC.h"

using namespace std;

struct BETHE{
    double c1=4.0, c2=8.0, c3=50.0;
    double c_sg=1.0;
};

struct DKD_KNODE{
    double k[3];
    double kn;
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
    double **karray;
    double *knarray;
    int    **kijarray;
private:
    DKD_KNODE *khead=nullptr;
    DKD_KNODE *ktail=nullptr;
    FILE   *fp=nullptr;
    void   add_k_vector(CELL *cell, double xy[2], double kn, int i=0, int j=0, bool southern_flag=false);
    void   free_k_node(DKD_KNODE *khead);
};

struct DKD_GNODE{
    double hkl[3];//Miller indices
    complex<double> Ug, qg;
    double sg;//excitation error
    bool   is_double_diffrac;
    bool   is_weak, is_strong;
    DKD_GNODE *next=nullptr;
    DKD_GNODE *nextw=nullptr;//weak
    DKD_GNODE *nexts=nullptr;//strong
};

class DKD_GVECTOR
{
public:
    DKD_GVECTOR(CELL *cell, BETHE bethe, double k[3], double fn[3], double cutoff);
    ~DKD_GVECTOR();
    int    numg=0;
    int    nstrong=0, nweak=0;
    DKD_GNODE *ghead=nullptr;
    DKD_GNODE *headw=nullptr;
private:
    DKD_GNODE *gtail=nullptr;
    void   add_g_vector(double hkl[3], complex<double> Ug, complex<double> qg, double sg, bool is_double_diffrac=false);
};

class DKD
{
public:
    double dmin=0.1;//smallest interplanar spacing to consider
    BETHE  bethe;//Bethe parameters

    int    numEbin=0, nump=0;
    double ***mLPNH=nullptr, ***mLPSH=nullptr; //The modified Lambert Projection Northern/Southern Hemisphere
    double ***mSPNH=nullptr, ***mSPSH=nullptr; //The master Stereographic Projection Northern/Southern Hemisphere
    DKD(DKD_MC *mc, CELL *cell, double dmin, double c1, double c2, double c3, double c_sg);
    ~DKD();
    void   img(char *img_path, double dimension=6, int resolution=512);
private:
    FILE   *fp=nullptr;
    double **lambdaE=nullptr;
    void   compute_scattering_probability(CELL *cell, DKD_MC *mc);
    void   compute_dynamic_matrix(complex<double> **dynmat, CELL *cell, DKD_GVECTOR *gvec);
    void   compute_Sgh_matrices(complex<double> ***Sgh, CELL *cell, DKD_GVECTOR *gvec);
    void   compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double *EWF, int IZMAX, double Z, double DZ, double KN, int NS);
    void   compute_Lambert_projection(int iE, bool is_hexagonal);
    void   compute_stereographic_projection(int iE, bool is_hexagonal);
};

#endif