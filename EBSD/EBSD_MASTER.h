#ifndef __EBSD_MASTER_H__
#define __EBSD_MASTER_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <complex>
#include <complex.h>
#include "../include/lapacke/lapacke.h"
#include "EBSD_MC.h"
#include "EBSD_CELL.h"
#include "EBSD_HDF5.h"
#include "EBSD_MATH.h"

using namespace std;
extern complex<double> to_complex(lapack_complex_double c);
extern lapack_complex_double to_lapack_complex(complex<double> c);

struct BETHE{
    double c1=8.0, c2=50.0, c3=100.0;
    double c_sg=1.00;
};

struct KNODE{
    double k[3];
    double kn;
    int i, j;
    int hemisphere;//Northern = 1, Southern = -1
    KNODE *next=nullptr;
};

class EBSD_KVECTOR
{
public:
    EBSD_KVECTOR(EBSD_CELL *cell, int nump);
    ~EBSD_KVECTOR();
    char   mode[20]="Rosca-Lambert";
    int    numk=0;
    double **karray;
    double *knarray;
    int    **kijarray;
private:
    KNODE  *khead=nullptr;
    KNODE  *ktail=nullptr;
    void   add_k_vector(EBSD_CELL *cell, double xy[2], double kn, int i=0, int j=0, bool southern_flag=false);
    void   free_k_node(KNODE *khead);
};

struct GNODE{
    int    number;//reflection number
    double hkl[3];//Miller indices
    complex<double> Ug, qg;
    double sg;//excitation error
    int    family_number;
    bool   is_double_diffrac;
    bool   is_weak, is_strong;
    GNODE  *next=nullptr;
    GNODE  *nextw=nullptr;//weak
    GNODE  *nexts=nullptr;//strong
};

class EBSD_GVECTOR
{
public:
    EBSD_GVECTOR(EBSD_CELL *cell, BETHE bethe, double k[3], double fn[3], double cutoff);
    ~EBSD_GVECTOR();
    int    numg=0;
    int    nstrong=0, nweak=0;
    GNODE  *ghead=nullptr;
    GNODE  *headw=nullptr;
private:
    GNODE  *gtail=nullptr;
    void   add_g_vector(double hkl[3], complex<double> Ug, complex<double> qg, double sg, bool is_double_diffrac=false);
};

class EBSD_MASTER
{
public:
    double dmin=0.05;//smallest d-spacing to consider
    int    nump=500;//number of master pattern pixels
    BETHE  bethe;//Bethe parameters

    double ****mLPNH, ****mLPSH; //The modified Lambert Projection Northern/Southern Hemisphere
    double ***mSPNH, ***mSPSH;
    EBSD_MASTER(const char *file_path);
    ~EBSD_MASTER();
    void   compute_master_pattern(const char* file_path);
    bool   write_parameters_into_hdf5(const char* file_path, int offset=0);
private:
    int    numEbin;
    int    izmax;
    int    npos;
    int    nset;
    int    width, height;
    double *Ebins=nullptr;
    double *depths=nullptr;//thickness integration for each energy
    double **lambdaE=nullptr;
    bool   read_parameters_from_nml(const char *file_path);
    void   precompute_master_pattern(EBSD_MC *mc, EBSD_CELL *cell);
    void   compute_dynamic_matrix(complex<double> **dynmat, EBSD_CELL *cell, EBSD_GVECTOR *gvec);
    void   compute_Sgh_matrices(complex<double> ***Sgh, EBSD_CELL *cell, EBSD_GVECTOR *gvec);
    void   compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double *EWF, int ZMAX, double THICK, double DEPTHSTEP, double KN, int NS);
    double get_Lambert_interpolation(double xyz[3], double ****mat, bool hexagonal_flag=false);
};

#endif