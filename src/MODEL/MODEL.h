#ifndef __AAVDP_MODEL_H__
#define __AAVDP_MODEL_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include "CONST_XRD.h"
#include "CONST_NED.h"
#include "CONST_KED.h"
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
#include "../QSPG/spglib/pointgroup.h"
#include "../QSPG/spglib/spg_database.h"
#include "../QB/QB.h"
#include "../QSPG/QSPG.h"
using namespace std;

#define PRE_CONST_RI1 0.5772157
#define DURCH_NUMBER 21
const double DURCH_TABLE[DURCH_NUMBER]={
1.000000,1.005051,1.010206,1.015472,1.020852,
1.026355,1.031985,1.037751,1.043662,1.049726,
1.055956,1.062364,1.068965,1.075780,1.082830,
1.090140,1.097737,1.105647,1.113894,1.122497,
1.131470};
double WK_get_Debye_Waller_factor(double G, double U);
double WK_get_electron_scattering_factor(double G, const double A[4], const double B[4]);
double WK_get_core_excitation_factor(double G, int Z, double V);
double WK_EI(double X);
double WK_IH2(double X);
double WK_IH1(double X1, double X2, double X3);
double WK_I1(double BI, double BJ, double G);
double WK_I2(double BI, double BJ, double G, double U);
double WK_get_absorptive_form_factor(double G, double U, const double A[4], const double B[4]);
complex<double> WK_get_scattering_amplitude(double G, double U, int Z, double V);

#define SYMPREC 0.001
#define SYMPREC2 0.000001

#define G_TO_WK 0.6283185307179586 //0.1*TWO_PI
#define DW_TO_WK 1.2665147955292222 //100.0/(8.0*pow(PI, 2))
#define PRE_CONST_V 0.04787801
#define PRE_CONST_U 0.664840340614319 //2.0*m0*e/h**2*1.0E-18

#define AVOGADRO_CONSTANT 6.02214076e23
#define PLANK_CONSTANT 6.62607015e-34 //J·s
#define ELECTRON_CHARGE 1.602176634e-19 //Coulomb
#define ELECTRON_REST_MASS 9.1093837090e-31 //kg
#define LIGHT_VELOCITY 299792458.0 //m/s

struct FOURIER
{
    //computation method (WK=Weickenmeier-Kohl), including absorptive form factors
    //Weickenmeier-Kohl method, working with reciprocal space in Angstrom scaled by a factor of 2*PI
    //Therefore, the G in [nm^-1] must be scaled by G_TO_WK=2*PI/10, and the DW in [nm^2] must be scaled by DW_TO_WK=100.0/8*PI^2 
    double Vmod, Vpmod;//modulus of Vg and Vgprime [V]
    double Umod, Upmod;//modulus of Ug and Ugprime [nm^-2]
    double Vang, Vpang;//phase factors of Vg and Vgprime [rad]
    complex<double> Vg, Ug;//potential Fourier coefficient [V] and scaled potential Fourier coefficient [nm^-2]
    double sig, sigp;//extinction distance and absorption length
    complex<double> qg;//interaction parameter for Darwin-Howie-Whelan equations [nm^-1]
};

struct FOURIER0
{
    double Vmod, Vpmod;//V0, mean inner potential
    double Umod, Upmod;//modulus of Ug and Ugprime [nm^-2]
    complex<double> Vg, Ug;//potential Fourier coefficient [V] and scaled potential Fourier coefficient [nm^-2]
    double sig, sigp;//sigp, Normal absorption length
    complex<double> qg;//interaction parameter for Darwin-Howie-Whelan equations [nm^-1]
    double gamma;//relativistic correction factor 
    double sigma;//interaction constant
    double voltage;//relativistic acceleration voltage
    double lambda;//wave length
};

class CELL
{
public:
    char   crystal_system[20];
    char   centering;//P, F, I, A, B, C, R, forming Bravais lattice combined with crystal system
    int    hall_number;
    int    space_group;//number between 1 and 230 using the International Tables for Crystallography, Volume A
    char   space_group_symbol[11];
    int    point_group;
    char   point_group_symbol[6];
    int    sampling_type;
    bool   is_trigonal=false;
    bool   is_hexagonal=false;//determine whether to use hexagonal indices
    bool   use_hexagonal=false;

    double a0, b0, c0, alpha, beta, gamma;
    double vol;
    double dmt[3][3], rmt[3][3];//direct and reciprocal metric tensor
    double dsm[3][3], rsm[3][3];//direct and reciprocal space matrix
    
    int    nsymmetry;
    int    sym_rotations[192][3][3];
    double sym_translations[192][3];
    int    npointsym;
    double point_dmats[48][3][3];//point group matrices in direct space
    double point_rmats[48][3][3];//point group matrices in reciprocal space

    int    npos=0, napos=0;//total number of asymmetric atoms
    int    *apos_type;
    int    *apos_multi=nullptr;//multiplicity for each asymmetric atomic position
    double ***apos_pos=nullptr;//fractional coordinates ranging between 0 and 1 in an asymmetric unit
    int    *apos_Z=nullptr;
    double *apos_M=nullptr;
    double *apos_DW=nullptr;
    double *apos_occupation=nullptr;

    int    ntype=0;
    char   **type_name=nullptr;

    double ave_M, ave_Z, density;
    CELL(const char *cell_path, const char types[][10], const double DWs[], double voltage);
    ~CELL();
    void   logging();
    double dot(double v1[3], double v2[3], char space);
    void   cross(double c_v[3], double v1[3], double v2[3], char space);
    double length(double v[3], char space);
    double angle(double v1[3], double v2[3], char space);
    void   normalize(double n_v[3], double v[3], char space);
    void   reciprocal_to_cartesian(double c_v[3], double v[3]);//from reciprocal
    void   cartesian_to_reciprocal(double r_v[3], double v[3]);//from cartesian
    void   direct_to_reciprocal(double r_v[3], double v[3]);
    void   direct_to_cartesian(double c_v[3], double v[3]);
    void   compute_equivalent_reciprocal_vectors(double equiv[48][3], int &nequiv, double g[3], char space='r');
    void   apply_point_group_symmetry(int equiv[48][3], int &nequiv, int px, int py, int hemisphere, int impx, int impy);

    int    HKL[3];
    FOURIER0 fouri0;
    FOURIER fouri;
    complex<double> ****LUTSgh=nullptr;
    complex<double> ***LUTUg=nullptr, ***LUTqg=nullptr;
    bool   ***is_double_diffrac=nullptr;//double diffraction reflection
    double get_interplanar_spacing(double g[3]);
    void   compute_reflection_range(double dmin);
    double get_excitation_error(double g[3], double k[3], double fn[3]);
    bool   is_centering_allowed(double g[3]);
    void   compute_Bloch_wave_coefficients(bool is_initial=true);
    void   compute_Fourier_coefficients(double voltage, bool is_initial=true);
    void   update_Fourier_coefficient(double voltage, double g[3], bool is_shift=false);
    void   update_Fourier_coefficient0(double voltage);
private:
    void   set_sampling_type(int sgnum, int pgnum);
    void   compute_lattice_matrices(double lat[3][3]);
    void   compute_symmetry_matrices(int sgnum, int (*rots)[3][3], double (*trans)[3], int noperation);
    void   compute_asymmetric_atomic_positions(double (*atom_pos)[3], int *atom_type, int natom);
    void   compute_atomic_density();
};

class MODEL
{
public:
    int    natom=0;
    double **atom_pos=nullptr;
    int    *atom_type=nullptr;
    double *atom_DW=nullptr;
    int    *atom_Z=nullptr;
    double *atom_occupation=nullptr;
    int    ntype=0;
    int    *type_index=nullptr;
    char   **type_name=nullptr;
    double vol;

    double dimension[3];
    double dimensionK[3];
    bool   is_periodic[3];
    bool   is_orthogonal=true;

    char   radiation[10];
    double lambda;
    MODEL(const char *model_path, const char types[][10], const double DWs[]);
    ~MODEL();
    void   reciprocal_to_cartesian(double c_g[3], double r_g[3]);
    void   direct_to_reciprocal(double r_g[3], double d_g[3]);
    double get_reciprocal_vector_length(double g[3]);
    void   compute_reciprocal_spacing(double spacing[3], double spacing_ratio[3]);
    virtual double get_diffraction_intensity(double theta, double g[3]);
    virtual double get_diffraction_intensity(double G, double g[3], bool is_zero);
private:
    double dsm[3][3], rsm[3][3];//direct and reciprocal space matrix
    double dmt[3][3], rmt[3][3];//direct and reciprocal metric tensor
    void   set_lattice(double mat[3][3]);
};

class XMODEL:public MODEL
{
public:
    XMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda);
    ~XMODEL();
    void   update_lorentzP_type(int lp_type);
    double get_diffraction_intensity(double theta, double g[3]);
    double get_diffraction_intensity(double G, double g[3], bool is_zero);
private:
    int    lorentzP_type=3;
    double get_Debye_Waller_factor(double S, double DW);
    double get_atomic_scattering_factor(double S, const double A[4], const double B[4], const double C);
    complex<double> get_atomic_structure_factor(double theta, double g[3]);
};

class NMODEL:public MODEL
{
public:
    NMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda);
    ~NMODEL();
    void   update_lorentzP_type(int lp_type);
    double get_diffraction_intensity(double theta, double g[3]);
    double get_diffraction_intensity(double G, double g[3], bool is_zero);
private:
    int    lorentzP_type=2;
    double get_Debye_Waller_factor(double S, double DW);
    complex<double> get_atomic_structure_factor(double theta, double g[3]);
};

class EMODEL:public MODEL
{
public:
    EMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda);
    ~EMODEL();
    double get_diffraction_intensity(double theta, double g[3]);
    double get_diffraction_intensity(double G, double g[3], bool is_zero);
private:
    double get_atomic_scattering_factor(double S, const double A[5], const double B[5]);
    complex<double> get_atomic_structure_factor(double theta, double g[3]);
};

#define VG_TO_WK 6.283185307179586 //TWO_PI
#define VDW_TO_WK 0.012665147955292222 //1.0/(8.0*pow(PI, 2))
#define VPRE_CONST_V 47.87801

class VMODEL:public MODEL
{
public:
    VMODEL(const char *model_path, const char types[][10], const double DWs[], double mvoltage);
    ~VMODEL();
    double get_diffraction_intensity(double theta, double g[3]);
    double get_diffraction_intensity(double G, double g[3], bool is_zero);
private:
    double voltage;//in kV
    void   set_wavelength();
};

#endif