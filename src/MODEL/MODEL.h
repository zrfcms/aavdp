#ifndef __AAVDP_MODEL_H__
#define __AAVDP_MODEL_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include "CONST.h"
#include "../QSPG/spglib/pointgroup.h"
#include "../QSPG/spglib/spg_database.h"
#include "../QB/QB.h"
#include "../QSPG/QSPG.h"
#include "../MATH/MATH.h"
#include "../HDF5/HDF5.h"

#define G_TO_WK 0.6283185307179586 //0.1*TWO_PI
#define DW_TO_WK 1.2665147955292222 //100.0/(8.0*pow(PI, 2))
#define PRE_CONST_V 0.04787801
#define PRE_CONST_U 0.664840340614319 //2.0*m0*e/h**2*1.0E-18
#define PRE_CONST_RI1 0.5772157
#define ELECTRON_CHARGE 1.602176634e-19 //Coulomb
#define ELECTRON_REST_MASS 9.1093837090e-31 //kg
#define ELECTRON_REST_ENERGY 0.511e6 //eV
#define LIGHT_VELOCITY 299792458.0 //m/s
#define AVOGADRO_CONSTANT 6.02214076e23
#define PLANK_CONSTANT 6.62607015e-34 //J·s

#define XLAMBDA 1.54180
#define NLAMBDA 1.54180
#define ELAMBDA 0.02510

using namespace std;

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
    int    point_group;
    char   point_group_symbol[6];
    int    space_group;//number between 1 and 230 using the International Tables for Crystallography, Volume A
    char   space_group_symbol[11];
    int    hall_number;
    char   centering;//P, F, I, A, B, C, R, forming Bravais lattice combined with crystal system
    int    sampling_type;
    bool   is_second_setting=false;
    bool   is_trigonal=false;
    bool   is_hexagonal=false;//determine whether to use hexagonal indices

    double lattice[3][3];
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

    double ave_M, ave_Z, density;
    CELL(const char *cell_path, const char types[][10], const double DWs[]);
    CELL(const char *hdf5_path);
    ~CELL();
    void   logging();
    void   hdf5(const char* hdf5_path);
    double dot(double v1[3], double v2[3], char space='r');
    double length(double v[3], char space='r');
    double angle(double v1[3], double v2[3], char space='r');
    void   normalize(double n_v[3], double v[3], char space='r');
    void   cartesian(double c_v[3], double v[3]);//from reciprocal
    void   reciprocal(double r_v[3], double v[3]);//from cartesian
    void   compute_equivalent_reciprocal_vectors(double equiv[48][3], int &nequiv, double g[3], char space='r');
    void   apply_point_group_symmetry(int equiv[48][3], int &nequiv, int px, int py, int pz, int nump);

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
    void   compute_asymmetric_atomic_positions(double (*atom_pos)[3], int *atom_type, int natom);
    void   compute_atomic_density();
    void   set_sampling_type();
    void   compute_lattice_matrices();
    void   compute_point_symmetry_matrices();

    complex<double> get_scattering_amplitude(double G, double U, int Z, double V);
    double get_Debye_Waller_factor(double G, double U);
    double get_electron_scattering_factor(double G, const double A[4], const double B[4]);
    double get_core_excitation_factor(double G, int Z, double V);//Image formation by inelastically scattered electrons in electron microscopy, 1976, H. Rose, p139-158
    double get_absorptive_form_factor(double G, double U, const double A[4], const double B[4]);
    double EI(double X);
    double IH1(double X1, double X2, double X3);
    double IH2(double X);
    double I1(double BI, double BJ, double G);
    double I2(double BI, double BJ, double G, double U);
};

class MODEL
{
public:
    double dimension[3];
    double dimensionK[3];
    bool   is_periodic[3];
    int    natom=0;
    double **atom_pos=nullptr;
    int    *atom_type=nullptr;
    double *atom_DW=nullptr;
    int    ntype=0;
    int    *type_index=nullptr;
    double lambda;
    bool   is_orthogonal=true;
    MODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda);
    ~MODEL();
    double get_reciprocal_vector_length(double g[3]);
    void   reciprocal_to_cartesian(double c_g[3], double r_g[3]);
    void   compute_reciprocal_spacing(double spacing[3], double spacing_ratio[3]);
private:
    double dsm[3][3], rsm[3][3];//direct and reciprocal space matrix
    double dmt[3][3], rmt[3][3];//direct and reciprocal metric tensor
    void   set_lattice(double mat[3][3]);
};

class XMODEL:public MODEL
{
public:
    XMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda=XLAMBDA);
    ~XMODEL();
    double get_diffraction_intensity(double theta, double g[3], bool is_lorentz);
private:
    double get_Debye_Waller_factor(double S, double DW);
    //complex<double> get_atomic_scattering_factor(double S, int type);
    double get_atomic_scattering_factor(double S, const double A[4], const double B[4], const double C);
    complex<double> get_atomic_structure_factor(double theta, double g[3]);
};

class NMODEL:public MODEL
{
public:
    NMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda=NLAMBDA);
    ~NMODEL();
    double get_diffraction_intensity(double theta, double g[3], bool is_lorentz);
private:
    double get_Debye_Waller_factor(double S, double DW);
    complex<double> get_atomic_structure_factor(double theta, double g[3]);
};

class EMODEL:public MODEL
{
public:
    EMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda=ELAMBDA);
    ~EMODEL();
    double get_diffraction_intensity(double theta, double g[3]);
private:
    double get_atomic_scattering_factor(double S, const double A[5], const double B[5]);
    complex<double> get_atomic_structure_factor(double theta, double g[3]);
};

#endif