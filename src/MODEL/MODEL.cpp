#include "MODEL.h"

CELL::CELL(const char *cell_path, const char types[][10], double mlambda)
{
    lambda=mlambda;
    QB_tools QB;
    QB_init(&QB);
    QB_read_lmp(&QB, cell_path);

    double lx=QB.mat[0][0], ly=QB.mat[1][1], lz=QB.mat[2][2];
    double xy=QB.mat[1][0], xz=QB.mat[2][0], yz=QB.mat[2][1];
    a=lx; b=sqrt(pow(ly,2)+pow(xy,2)); c=sqrt(pow(lz,2)+pow(xz,2)+pow(yz,2));
    alpha=acos((xy*xz+ly*yz)/(b*c))/DEG_TO_RAD; beta=acos(xz/c)/DEG_TO_RAD; gamma=acos(xy/b)/DEG_TO_RAD;
    double calpha=cos(DEG_TO_RAD*alpha), cbeta=cos(DEG_TO_RAD*beta), cgamma=cos(DEG_TO_RAD*gamma);
    double salpha=sin(DEG_TO_RAD*alpha), sbeta=sin(DEG_TO_RAD*beta), sgamma=sin(DEG_TO_RAD*gamma);
    double tgamma=tan(DEG_TO_RAD*gamma);
    double det=a*b*c*a*b*c*(1.0-calpha*calpha-cbeta*cbeta-cgamma*cgamma+2.0*calpha*cbeta*cgamma);
    vol=sqrt(det);
    if(vol<1e-6){
        printf("[ERROR] Suspiciously small volume with six lattice parameters, %.5f, %.5f, %.5f, %.5f, %.5f, and %.5f.", a, b, c, alpha, beta, gamma);
        exit(EXIT_FAILURE);
    }
    dmt[0][0]=a*a; dmt[0][1]=a*b*cgamma; dmt[0][2]=a*c*cbeta;
    dmt[1][0]=dmt[0][1]; dmt[1][1]=b*b; dmt[1][2]=b*c*calpha;
    dmt[2][0]=dmt[0][2]; dmt[2][1]=dmt[1][2]; dmt[2][2]=c*c;
    rmt[0][0]=b*c*salpha*b*c*salpha/det; rmt[0][1]=a*b*c*c*(calpha*cbeta-cgamma)/det; rmt[0][2]=a*b*b*c*(cgamma*calpha-cbeta)/det;
    rmt[1][0]=rmt[0][1]; rmt[1][1]=a*c*sbeta*a*c*sbeta/det; rmt[1][2]=a*a*b*c*(cbeta*cgamma-calpha)/det;
    rmt[2][0]=rmt[0][2]; rmt[2][1]=rmt[1][2]; rmt[2][2]=a*b*sgamma*a*b*sgamma/det;
    dsm[0][0]=a; dsm[0][1]=b*cgamma; dsm[0][2]=c*cbeta;
    dsm[1][0]=0.0; dsm[1][1]=b*sgamma; dsm[1][2]=-c*(cbeta*cgamma-calpha)/sgamma;
    dsm[2][0]=0.0; dsm[2][1]=0.0; dsm[2][2]=vol/(a*b*sgamma);
    rsm[0][0]=1.0/a; rsm[0][1]=0.0; rsm[0][2]=0.0;
    rsm[1][0]=-1.0/(a*tgamma); rsm[1][1]=1.0/(b*sgamma); rsm[1][2]=0.0;
    rsm[2][0]=b*c*(cgamma*calpha-cbeta)/(vol*sgamma); rsm[2][1]=a*c*(cbeta*cgamma-calpha)/(vol*sgamma); rsm[2][2]=(a*b*sgamma)/vol;

    ntype=QB.TypeNumber;
    natom=QB.TotalNumber;
    mallocate(&atom_type, natom);
    mallocate_2d(&atom_pos, natom, 3);
    for(int i=0;i<natom;i++){
        atom_pos[i][0]=QB.atom[i].x; 
        atom_pos[i][1]=QB.atom[i].y; 
        atom_pos[i][2]=QB.atom[i].z;
        atom_type[i]=QB.atom[i].type;
    }
    for(int i=0;i<natom;i++){
        cartesian_to_direct(atom_pos[i], atom_pos[i]);
    }
    QB_free_atom(&QB);

    lambda=mlambda;
    callocate(&type_index, ntype, -1);
    for(int i=0;i<ntype;i++){
        for(int j=0;j<E_TYPE_NUMBER;j++){
            if(0==strcmp(types[i], E_TYPE[j])){
                type_index[i]=j;
                break;
            }
        }
    }
}

CELL::~CELL()
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

void CELL::cartesian_to_direct(double d_v[3], double c_v[3])
{
    vector_rotate(d_v, dsm, c_v);
}

double CELL::get_reciprocal_vector_length(double g[3])
{
    double r_g[3];
    vector_rotate(r_g, rmt, g);
    double len=vector_dot(g, r_g);
    len=sqrt(len);
    return len;
}

void CELL::compute_reciprocal_spacing(double spacing[3], double spacing_ratio[3])
{
    double dimension[3]={a, b, c};
    for(int i=0;i<3;i++){
        spacing[i]=1.0/dimension[i]*spacing_ratio[i];
    }
}

double CELL::get_atomic_scattering_factor(double S, const double A[4], const double B[4], const double C)//xrd
{
    double res=0.0;
    for(int i=0;i<4;i++){
        res+=A[i]*exp(-B[i]*S*S);//S, scattering vector
    }
    res+=C;
    return res;
}

complex<double> CELL::get_atomic_structure_factor(double theta, double g[3], bool is_lorentz_flag)//xrd
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    for(int i=0;i<natom;i++){
        int    t=type_index[atom_type[i]-1];
        const double *A=X_A[t], *B=X_B[t], C=X_C[t];
        double q=2*PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
        res+=get_atomic_scattering_factor(S, A, B, C)*complex<double>(cos(q), sin(q));
    }
    double Lp=sqrt((1+cos(2*theta)*cos(2*theta))/(sin(theta)*sin(theta)*cos(theta)));
    res*=complex<double>(Lp, 0.0);
    return res;
}

double CELL::get_diffraction_intensity(double theta, double g[3], bool is_lorentz_flag)
{
    complex<double> F=get_atomic_structure_factor(theta, g, is_lorentz_flag);
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    return res;
}

double CELL::get_diffraction_intensity(complex<double> F)
{
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    return res;
}

MODEL::MODEL(const char *model_path, const char types[][10], double mlambda, char mmode)
{
    lambda=mlambda;
    mode=mmode;
    const char (*TYPE)[10]=nullptr;
    int TYPE_NUM=0;
    switch(mode){
    case 'x':
        TYPE=X_TYPE;
        TYPE_NUM=X_TYPE_NUMBER;
        break;
    case 'e':
        TYPE=E_TYPE;
        TYPE_NUM=E_TYPE_NUMBER;
    default:
        printf("[ERROR] Unrecognized diffraction mode %c.", mode);
        exit(EXIT_FAILURE);
    }

    QB_tools QB;
    QB_init(&QB);
    QB_read_lmp(&QB, model_path);
    dimension[0]=QB.boundx; dimension[1]=QB.boundy; dimension[2]=QB.boundz;
    is_periodic[0]=QB.px; is_periodic[1]=QB.py; is_periodic[2]=QB.pz;
    matrix_lattice(QB.mat);
    ntype=QB.TypeNumber;
    callocate(&type_index, ntype, -1);
    index_type(types, TYPE, TYPE_NUM);
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

void MODEL::index_type(const char types[][10], const char TYPE[][10], int TYPE_NUM)
{
    for(int i=0;i<ntype;i++){
        for(int j=0;j<TYPE_NUM;j++){
            if(0==strcmp(types[i], TYPE[j])){
                type_index[i]=j;
                break;
            }
        }
    }
}

void MODEL::matrix_lattice(double mat[3][3])
{
    double a[3], b[3], c[3];
    double bc[3], ca[3], ab[3];
    a[0]=mat[0][0]; a[1]=mat[1][0]; a[2]=mat[2][0];
    b[0]=mat[0][1]; b[1]=mat[1][1]; b[2]=mat[2][1];
    c[0]=mat[0][2]; c[1]=mat[1][2]; c[2]=mat[2][2];
    vector_cross(bc, b, c);
    vector_cross(ca, c, a);
    vector_cross(ab, a, b);
    double vol=vector_dot(a, bc);

    matrix_copy(dsm, mat);
    rsm[0][0]=bc[0]; rsm[1][0]=ca[0]; rsm[2][0]=ab[0];
    rsm[0][1]=bc[1]; rsm[1][1]=ca[1]; rsm[2][1]=ab[1];
    rsm[0][2]=bc[2]; rsm[1][2]=ca[2]; rsm[2][2]=ab[2];
    matrix_constant(rsm, 1.0/vol, rsm);
    double t_mat[3][3];
    matrix_transpose(t_mat, dsm);
    matrix_multiply(dmt, t_mat, dsm);
    matrix_transpose(t_mat, rsm);
    matrix_multiply(rmt, t_mat, rsm);
}

MODEL::~MODEL()
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

double MODEL::get_reciprocal_vector_length(double g[3])
{
    double r_g[3];
    vector_rotate(r_g, rmt, g);
    double len=vector_dot(g, r_g);
    len=sqrt(len);
    return len;
}

void MODEL::compute_reciprocal_spacing(double spacing[3], double spacing_ratio[3])
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

double MODEL::get_atomic_scattering_factor(double S, const double A[4], const double B[4], const double C)//xrd
{
    double res=0.0;
    for(int i=0;i<4;i++){
        res+=A[i]*exp(-B[i]*S*S);//S, scattering vector
    }
    res+=C;
    return res;
}

double MODEL::get_atomic_scattering_factor(double S, const double A[5], const double B[5])
{
    double res=0.0;
    for(int i=0;i<5;i++){
        res+=A[i]*exp(-B[i]*S*S);//S, scattering vector
    }
    return res;
}

complex<double> MODEL::get_atomic_structure_factor(double theta, double g[3], bool is_lorentz_flag)//xrd
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    for(int i=0;i<natom;i++){
        int    t=type_index[atom_type[i]-1];
        const double *A=X_A[t], *B=X_B[t], C=X_C[t];
        double q=2*PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
        res+=get_atomic_scattering_factor(S, A, B, C)*complex<double>(cos(q), sin(q));
    }
    double Lp=sqrt((1+cos(2*theta)*cos(2*theta))/(sin(theta)*sin(theta)*cos(theta)));
    res*=complex<double>(Lp, 0.0);
    return res;
}

complex<double> MODEL::get_atomic_structure_factor(double theta, double g[3])
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    if(0.0<=S<=2.0){
        for(int i=0;i<natom;i++){
            int    t=type_index[atom_type[i]-1];
            const double *A=E_A1[t], *B=E_B1[t];
            double q=2*PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
            res+=get_atomic_scattering_factor(S, A, B)*complex<double>(cos(q), sin(q));
        }
    }else if(2.0<S<=6.0){
        for(int i=0;i<natom;i++){
            int    t=type_index[atom_type[i]-1];
            const double *A=E_A1[t], *B=E_B1[t];
            double q=2*PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
            res+=get_atomic_scattering_factor(S, A, B)*complex<double>(cos(q), sin(q));
        }
    }else{
        res+=complex<double>(0.0, 0.0);
    }
    return res;
}

double MODEL::get_diffraction_intensity(double theta, double g[3], bool is_lorentz_flag)
{
    complex<double> F=get_atomic_structure_factor(theta, g, is_lorentz_flag);
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    return res;
}

double MODEL::get_diffraction_intensity(double theta, double g[3])
{
    complex<double> F=get_atomic_structure_factor(theta, g);
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    return res;
}

double MODEL::get_diffraction_intensity(complex<double> F)
{
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    return res;
}

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

complex<double> XMODEL::get_atomic_structure_factor(double theta, double g[3])
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    for(int i=0;i<natom;i++){
        int    t=type_index[atom_type[i]-1];
        const double *A=X_A[t], *B=X_B[t], C=X_C[t];
        double q=2*PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
        res+=get_atomic_scattering_factor(S, A, B, C)*complex<double>(cos(q), sin(q));
    }
    return res;
}

double XMODEL::get_diffraction_intensity(double theta, double g[3])
{
    complex<double> F=get_atomic_structure_factor(theta, g);
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

complex<double> EMODEL::get_atomic_structure_factor(double theta, double g[3])
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    if(0.0<=S<=2.0){
        for(int i=0;i<natom;i++){
            int    t=type_index[atom_type[i]-1];
            const double *A=E_A1[t], *B=E_B1[t];
            double q=2*PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
            res+=get_atomic_scattering_factor(S, A, B)*complex<double>(cos(q), sin(q));
        }
    }else if(2.0<S<=6.0){
        for(int i=0;i<natom;i++){
            int    t=type_index[atom_type[i]-1];
            const double *A=E_A1[t], *B=E_B1[t];
            double q=2*PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
            res+=get_atomic_scattering_factor(S, A, B)*complex<double>(cos(q), sin(q));
        }
    }else{
        res+=complex<double>(0.0, 0.0);
    }
    return res;
}

double EMODEL::get_diffraction_intensity(double theta, double g[3])
{
    complex<double> F=get_atomic_structure_factor(theta, g);
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