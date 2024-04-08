#ifndef __EBSD_CELL_H__
#define __EBSD_CELL_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include "EBSD_HDF5.h"
#include "EBSD_MATH.h"
#include "EBSD_CONST.h"

using namespace std;
extern void matrix_copy(double outmat[4][4], double inmat[4][4]);

struct FOURIER
{
    char   method[5]="WK";//computation method (WK=Weickenmeier-Kohl), including absorptive form factors
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

class EBSD_CELL
{
public:
    //parameters required
    int    crystal_system;//cubic, tetragonal, orthorhombic, hexagonal, rhombohedral/trigonal, monoclinic, triclinic ranging from 1 to 7
    double a, b, c, alpha, beta, gamma;//six lattice parameters of the unit cell, respectively in [nm] and [degree]
    int    napos;//total number of asymmetric atoms
    double **apos=nullptr;//fractional coordinates ranging between 0 and 1 in an asymmetric unit
    int    *apos_Z=nullptr;//atomic number for each asymmetric atom
    double *apos_occupation=nullptr; //site occupation for each asymmetric atom
    double *apos_DW=nullptr; //Debye-Waller factor for each asymmetric atom depending on the temperature, in [nm^2]
    int    space_group;//number between 1 and 230 using the International Tables for Crystallography, Volume A
    bool   second_setting_flag=false;   
    bool   trigonal_flag=false;
    bool   hexagonal_flag=false;//determine whether to use hexagonal indices

    int    point_group;
    char   generator_string[40];
    char   space_group_symbol[11];
    char   centering_symbol;//P, F, I, A, B, C, R, forming Bravais lattice combined with crystal system
    bool   centrosymmetric_flag=false;//to be continued
    bool   nonsymmorphic_flag=false;//to be continued
    int    sampling_type;

    int    *apos_multi=nullptr;//multiplicity for each asymmetric atomic position
    double ***pos=nullptr;
    double ave_Z; //average atomic number
    double ave_M; //average atomic weight in g/mole
    double density;//density in g/cm^3
    EBSD_CELL(const char* file_path, bool logging_flag=true);
    ~EBSD_CELL();
    bool   write_parameters_into_hdf5(const char* file_path);
    void   compute_equivalent_atomic_positions(double equiv[48][3], int &nequiv, const double atom_apos[3], bool reduce_symmetry_flag=true);
    void   compute_equivalent_reciprocal_vectors(double equiv[48][3], int &nequiv, double g[3], char space='r');
    void   apply_point_group_symmetry(int equiv[48][3], int &nequiv, int px, int py, int pz, int nump);

    double dot(double hkl1[3], double hkl2[3], char space='r');
    void   cross(double c_hkl[3], double hkl1[3], double hkl2[3], char inspace='r', char outspace='r');
    double length(double hkl[3], char space='r');
    double angle(double hkl1[3], double hkl2[3], char space='r');
    void   normalize(double n_v[3], double v[3], char space='c');
    void   reciprocal(double r_v[3], double v[3]);
    void   cartesian(double c_v[3], double v[3]);

    double lambda;
    int    HKL[3];
    FOURIER0 fouri0;
    complex<double> ***LUTUg=nullptr, ***LUTqg=nullptr;
    bool   ***is_double_diffrac=nullptr;//double diffraction reflection
    complex<double> ****LUTSgh=nullptr;
    double get_interplanar_spacing(double g[3]);
    void   compute_reflection_range(double dmin);
    bool   is_centering_allowed(double g[3]);
    double get_excitation_error(double g[3], double k[3], double fn[3]);
    void   compute_Fourier_coefficient(FOURIER *fouri, double voltage, double g[3], bool shiftflag=false);
    void   compute_relativistic_wavelength(double voltage);
    void   compute_Fourier_coefficients(double voltage, bool initflag=true);
    void   compute_Bloch_wave_coefficients(int initflag=true);
private:
    double vol;//cell volume [nm^3] 
    double dmt[3][3], rmt[3][3];//direct and reciprocal metric tensor
    double dsm[3][3], rsm[3][3];//direct and reciprocal space matrix
    double tridsm[3][3];//second direct space matrix with x-[100] and z-[111] for rhombohedral/trigonal case, to be continued
    int    ngenerator;//number of generator matrices
    int    nsymmetry;//number of non-zero symmetry matrices
    double symmetry_mats[192][4][4];//all symmetry matrices for a given spacegroup
    int    npointsym;
    double point_dmats[48][3][3];//point group matrices in direct space
    double point_rmats[48][3][3];//point group matrices in reciprocal space
    bool   read_parameters_from_nml(const char* file_path);
    bool   read_parameters_from_hdf5(const char* file_path);
    void   set_sampling_type();
    void   compute_space_matrices();
    void   compute_trigonal_space_matrix(double angle);
    void   compute_symmetry_matrices();
    void   compute_symmetry_matrix_by_one_generator(double symmat[4][4], const char generator[4], const int multiplier=1);
    void   compute_symmetry_matrix_by_symmetry_matrix(double symmat[4][4], const double symmati[4][4], const double symmatj[4][4]);
    void   reduce_symmetry_matrix(double symmat[4][4], const double eps=0.0005);
    bool   is_symmetry_matrix_new(double symmat[4][4], const double eps=0.0005);
    void   compute_point_symmetry_matrices();
    void   compute_atomic_positions();
    void   compute_atomic_density();

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

#endif