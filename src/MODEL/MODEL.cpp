#include "MODEL.h"

XMODEL::XMODEL(const char *model_path, const char types[][10], double mlambda)
{
    lambda=mlambda;
    QB_tools QB;
    QB_init(&QB);
    QB_read_lmp(&QB, model_path);
    dimension[0]=QB.boundx; dimension[1]=QB.boundy; dimension[2]=QB.boundz;
    is_periodic[0]=QB.px; is_periodic[1]=QB.py; is_periodic[2]=QB.pz;
    ntype=QB.TypeNumber;
    callocate(&type_index, ntype, -1);
    for(int i=0;i<ntype;i++){
        for(int j=0;j<X_TYPE_NUMBER;j++){
            if(0==strcmp(types[i], X_TYPE[j])){
                type_index[i]=j;
                break;
            }
        }
    }
    natom=QB.TotalNumber;
    mallocate(&atom_type, natom);
    mallocate_2d(&atom_pos, natom, 3);
    for(int i=0;i<natom;i++){
        atom_pos[i][0]=QB.atom[i].x; 
        atom_pos[i][1]=QB.atom[i].y; 
        atom_pos[i][2]=QB.atom[i].z;
        atom_type[i]=QB.atom[i].type;
    }
    QB_free_atom(&QB);
}

XMODEL::~XMODEL()
{
    if(ntype>0){
        deallocate(type_index);
        ntype=0;
    }
    if(natom>0){
        deallocate_2d(atom_pos, natom);
        deallocate(atom_type);
        natom=0;
    }
}

void XMODEL::update_incident_wavelength(double mlambda)
{
    lambda=mlambda;
}

double XMODEL::get_atomic_scattering_factor(double S, const double A[4], const double B[4], const double C)
{
    double res=0.0;
    for(int i=0;i<4;i++){
        res+=A[i]*exp(-B[i]*S*S);//S, scattering vector
    }
    res+=C;
    return res;
}

complex<double> XMODEL::get_atomic_structure_factor(double theta, double k[3])
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    for(int i=0;i<natom;i++){
        int    t=type_index[atom_type[i]-1];
        const double *A=X_A[t], *B=X_B[t], C=X_C[t];
        double q=2*PI*(k[0]*atom_pos[i][0]+k[1]*atom_pos[i][1]+k[2]*atom_pos[i][2]);
        res+=get_atomic_scattering_factor(S, A, B, C)*complex<double>(cos(q), sin(q));
    }
    return res;
}

double XMODEL::get_diffraction_intensity(double theta, double k[3])
{
    complex<double> F=get_atomic_structure_factor(theta, k);
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    return res;
}

EMODEL::EMODEL(const char *model_path, const char types[][10], double mlambda)
{
    lambda=mlambda;
    QB_tools QB;
    QB_init(&QB);
    QB_read_lmp(&QB, model_path);
    dimension[0]=QB.boundx; dimension[1]=QB.boundy; dimension[2]=QB.boundz;
    is_periodic[0]=QB.px; is_periodic[1]=QB.py; is_periodic[2]=QB.pz;
    ntype=QB.TypeNumber;
    callocate(&type_index, ntype, -1);
    for(int i=0;i<ntype;i++){
        for(int j=0;j<E_TYPE_NUMBER;j++){
            if(0==strcmp(types[i], E_TYPE[j])){
                type_index[i]=j;
                break;
            }
        }
    }
    natom=QB.TotalNumber;
    mallocate(&atom_type, natom);
    mallocate_2d(&atom_pos, natom, 3);
    for(int i=0;i<natom;i++){
        atom_pos[i][0]=QB.atom[i].x; 
        atom_pos[i][1]=QB.atom[i].y; 
        atom_pos[i][2]=QB.atom[i].z;
        atom_type[i]=QB.atom[i].type;
    }
    QB_free_atom(&QB);
}

EMODEL::~EMODEL()
{
    if(ntype>0){
        deallocate(type_index);
        ntype=0;
    }
    if(natom>0){
        deallocate_2d(atom_pos, natom);
        deallocate(atom_type);
        natom=0;
    }
}

void EMODEL::update_incident_wavelength(double mlambda)
{
    lambda=mlambda;
}

double EMODEL::get_atomic_scattering_factor(double S, const double A[5], const double B[5])
{
    double res=0.0;
    for(int i=0;i<5;i++){
        res+=A[i]*exp(-B[i]*S*S);//S, scattering vector
    }
    return res;
}

complex<double> EMODEL::get_atomic_structure_factor(double theta, double k[3])
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    if(0.0<=S<=2.0){
        for(int i=0;i<natom;i++){
            int    t=type_index[atom_type[i]-1];
            const double *A=E_A1[t], *B=E_B1[t];
            double q=2*PI*(k[0]*atom_pos[i][0]+k[1]*atom_pos[i][1]+k[2]*atom_pos[i][2]);
            res+=get_atomic_scattering_factor(S, A, B)*complex<double>(cos(q), sin(q));
        }
    }else if(2.0<S<=6.0){
        for(int i=0;i<natom;i++){
            int    t=type_index[atom_type[i]-1];
            const double *A=E_A1[t], *B=E_B1[t];
            double q=2*PI*(k[0]*atom_pos[i][0]+k[1]*atom_pos[i][1]+k[2]*atom_pos[i][2]);
            res+=get_atomic_scattering_factor(S, A, B)*complex<double>(cos(q), sin(q));
        }
    }else{
        res+=complex<double>(0.0, 0.0);
    }
    return res;
}

double EMODEL::get_diffraction_intensity(double theta, double k[3])
{
    complex<double> F=get_atomic_structure_factor(theta, k);
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    return res;
}

void EMODEL::compute_reciprocal_spacing(double spacing[3], double spacing_ratio[3])
{
    double dimensionK[3];
    if((!is_periodic[0])&&(!is_periodic[1])&&(!is_periodic[2])){
        for(int i=0;i<3;i++){
          dimensionK[i]=1.0;
        }
    }else{
        double averageK=0.0;
        for(int i=0;i<3;i++){
            if(is_periodic[i]){
                dimensionK[i]=1.0/dimension[i];
                averageK+=dimensionK[i];
            }
        }
        averageK=averageK/double(is_periodic[0]+is_periodic[1]+is_periodic[2]);
        for(int i=0;i<3;i++){
            if(!is_periodic[i]){
                dimensionK[i]=averageK;
            }
        }
    }
    for(int i=0;i<3;i++){
        spacing[i]=dimensionK[i]*spacing_ratio[i];
    }
}