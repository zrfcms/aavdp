#include "MODEL.h"

double WK_get_Debye_Waller_factor(double G, double U)
{
    return exp(-0.5*U*U*G*G);
}

double WK_get_electron_scattering_factor(double G, const double A[4], const double B[4])
{
    double F=0.0;
    double S2=G*G/FOUR_PI/FOUR_PI;//S: scattering vector
    for(int i=0;i<4;i++){
        double TEMP=B[i]*S2;
        if(TEMP<0.1){
            F+=A[i]*B[i]*(1.0-0.5*TEMP);
        }else if(TEMP>20.0){
            F+=A[i]/S2;
        }else{
            F+=A[i]*(1.0-exp(-TEMP))/S2;
        }
    }
    return F;
}

double WK_get_core_excitation_factor(double G, int Z, double V)
{
    double F, A0=0.5289;
    double K0=0.5068*sqrt(1022.0*V+V*V);
    double DE=6*Z*1.0e-3;//Energy loss
    double THETAB=G/(2.0*K0); //Bragg angle
    double THETAE=DE/(2.0*V)*(2.0*V+1022.0)/(V+1022.0);//Angle corresponding to energy loss
    double THETAN=1.0/(K0*0.885*A0/pow(Z, 0.333333));

    double OMEGA=2.0*THETAB/THETAN;
    double KAPPA=THETAE/THETAN;
    double O2=OMEGA*OMEGA;
    double K2=KAPPA*KAPPA;
    double X1=OMEGA/((1.0+O2)*sqrt(O2+4.0*K2))*log((OMEGA+sqrt(O2+4.0*K2))/(2.0*KAPPA));
    double X2=1.0/sqrt((1.0+O2)*(1.0+O2)+4.0*K2*O2)*log((1.0+2.0*K2+O2+sqrt((1.0+O2)*(1.0+O2)+4.0*K2*O2))/(2.0*KAPPA*sqrt(1.0+K2)));
    double X3;
    if(OMEGA>1e-2){
        X3=1.0/(OMEGA*sqrt(O2+4.0*(1.0+K2)))*log((OMEGA+sqrt(O2+4.0*(1.0+K2)))/(2.0*sqrt(1.0+K2)));
    }else{
        X3=1.0/(4.0*(1.0+K2));
    }
    F=4.0/(A0*A0)*TWO_PI/(K0*K0)*2*Z/(THETAN*THETAN)*(-X1+X2-X3);
    return F;
}

double WK_EI(double X)
{
    double A1=8.57332, A2=18.05901, A3=8.63476, A4=0.26777, B1=9.57332, B2=25.63295, B3=21.09965, B4=3.95849;
    if(X>60.0){
        printf("[ERROR] The input X of the function EI exceeds 60.0");
        exit(1);
    }
    if(X<-60.0) return 0.0;
    double ABSX=fabs(X);
    if(X<-1.0){
        return -(A4+ABSX*(A3+ABSX*(A2+ABSX*(A1+ABSX))))/(B4+ABSX*(B3+ABSX*(B2+ABSX*(B1+ABSX))))*exp(-ABSX)/ABSX;
    }else{
        double REI=0.577216+log(ABSX);
        int I=1;
        double SI=X;
        double SUMS=SI;
        while(fabs(SI/X)>1.0e-6){
            SI=SI*X*I/((I+1)*(I+1));
            SUMS+=SI;
            I++;
        };
        REI+=SUMS;
        return REI;
    }
}

double WK_IH2(double X)
{
    double INVX=1.0/X;
    int    I=int(200.0*INVX);
    if(I<0||I>DURCH_NUMBER-1){
        printf("[ERROR] Unrecognized index %d in searching durch table.\n", I);
        exit(1);
    }
    double D1=DURCH_TABLE[I], D2=DURCH_TABLE[I+1];
    double RIH2=D1+200.0*(D2-D1)*(INVX-0.5e-3*I);
    return RIH2;
}

double WK_IH1(double X1, double X2, double X3)
{
    double RIH1;
    if(X2<=20.0&&X3<=20.0){
        RIH1=exp(-X1)*(WK_EI(X2)-WK_EI(X3));
        return RIH1;
    }
    if(X2>20.0){
        RIH1=exp(X2-X1)*WK_IH2(X2)/X2;
    }else{
        RIH1=exp(-X1)*WK_EI(X2);
    }
    if(X3>20.0){
        RIH1-=exp(X3-X1)*WK_IH2(X3)/X3;
    }else{
        RIH1-=exp(-X1)*WK_EI(X3);
    }
    return RIH1;
}

double WK_I1(double BI, double BJ, double G)
{
    double RI1;
    double G2=G*G;
    double EPS=fmax(BI, BJ)*G2;
    if(EPS<=0.1){
        double RI1=PI*(BI*log((BI+BJ)/BI)+BJ*log((BI+BJ)/BJ));
        if(fabs(G)<1.0e-12){
            return RI1;
        }
        double BI2=BI*BI, BJ2=BJ*BJ;
        double TEMP=0.5*BI2*log(BI/(BI+BJ))+0.5*BJ2*log(BJ/(BI+BJ));
        TEMP+=0.75*(BI2+BJ2)-0.25*(BI+BJ)*(BI+BJ);
        TEMP+=-0.5*(BI-BJ)*(BI-BJ);
        RI1+=PI*G2*TEMP;
        return RI1;
    }
    double BIG2=BI*G2, BJG2=BJ*G2;
    RI1=2.0*PRE_CONST_RI1+log(BIG2)+log(BJG2)-2.0*WK_EI(-BI*BJ*G2/(BI+BJ));
    RI1+=WK_IH1(BIG2, BIG2*BI/(BI+BJ), BIG2);
    RI1+=WK_IH1(BJG2, BJG2*BJ/(BI+BJ), BJG2);
    RI1=RI1*PI/G2;
    return RI1;
}

double WK_I2(double BI, double BJ, double G, double U)
{
    double RI2;
    double U2=U*U, G2=G*G;
    double HALFU2=0.5*U2, QUARTERU2=0.25*U2;
    double BIUH=BI+HALFU2, BJUH=BJ+HALFU2;
    double BIU=BI+U2, BJU=BJ+U2;
    double EPS=fmax(fmax(BI, BJ), U2)*G2;

    if(EPS<=0.1){
        double TEMP1=(BI+U2)*log((BI+BJ+U2)/(BI+U2));
        TEMP1+=BJ*log((BI+BJ+U2)/(BJ+U2));
        TEMP1+=U2*log(U2/(BJ+U2));
        if(G==0.0){
            RI2=PI*TEMP1;
        }else{
            double TEMP2=0.5*HALFU2*HALFU2*log( BIU*BJU/(U2*U2));
            TEMP2+=0.5*BIUH*BIUH*log(BIU/(BIUH+BJUH));
            TEMP2+=0.5*BJUH*BJUH*log(BJU/(BIUH+BJUH));
            TEMP2+=0.25*BIU*BIU+0.5*BI*BI;
            TEMP2+=0.25*BJU*BJU+0.5*BJ*BJ;
            TEMP2+=-0.25*(BIUH+BJUH)*(BIUH+BJUH);
            TEMP2+=-0.5*pow((BI*BIU-BJ*BJU)/(BIUH+BJUH), 2);
            TEMP2+=-HALFU2*HALFU2;
            RI2=PI*(TEMP1+G2*TEMP2);
        }
    }else{
        RI2=WK_EI(-HALFU2*G2*BIUH/BIU)+WK_EI(-HALFU2*G2*BJUH/BJU);
        RI2=2.0*(RI2-WK_EI(-BIUH*BJUH*G2/(BIUH+BJUH))-WK_EI(-QUARTERU2*G2));
        RI2+=WK_IH1(HALFU2*G2, QUARTERU2*G2, QUARTERU2*U2*G2/BIU);
        RI2+=WK_IH1(HALFU2*G2, QUARTERU2*G2, QUARTERU2*U2*G2/BJU);
        RI2+=WK_IH1(BIUH*G2, BIUH*BIUH*G2/(BIUH+BJUH), BIUH*BIUH*G2/BIU);
        RI2+=WK_IH1(BJUH*G2, BJUH*BJUH*G2/(BIUH+BJUH), BJUH*BJUH*G2/BJU);
        RI2=RI2*PI/G2;
    }
    return RI2;
}

double WK_get_absorptive_form_factor(double G, double U, const double A[4], const double B[4])
{
    double F=0.0;
    double AF[4], BF[4];
    double FP2=FOUR_PI*FOUR_PI;
    double DWF=WK_get_Debye_Waller_factor(G, U);
    for(int i=0;i<4;i++){
        AF[i]=A[i]*FP2; BF[i]=B[i]/FP2;
    }
    for(int j=0;j<4;j++){
        F+=AF[j]*AF[j]*(DWF*WK_I1(BF[j], BF[j], G)-WK_I2(BF[j], BF[j], G, U));
        for(int i=0;i<j;i++){
            F+=2.0*AF[j]*AF[i]*(DWF*WK_I1(BF[i], BF[j], G)-WK_I2(BF[i], BF[j], G, U));
        }
    }
    return F;
}

complex<double> WK_get_scattering_amplitude(double G, double U, int Z, double V)
{
    complex<double> F(0.0, 0.0);
    double FR, FI;
    double GAMMA=(V+511.0)/511.0;
    double K0=0.5068*sqrt(1022.0*V+V*V);
    double DWF=WK_get_Debye_Waller_factor(G, U);//Debye-Waller amplitude is included
    const double *A=E_WK_AA[Z-1], *B=E_WK_BB[Z-1];
    FR=FOUR_PI*DWF*WK_get_electron_scattering_factor(G, A, B);//Elastic part
    FI=DWF*WK_get_core_excitation_factor(G, Z, V)+WK_get_absorptive_form_factor(G, U, A, B);//Inelastic part
    FR=FR*GAMMA; FI=FI*GAMMA*GAMMA/K0;//acceleration voltage (relativistic scattering) is included
    F.real(FR); F.imag(FI);
    return F;
}

CELL::CELL(const char *cell_path, const char types[][10], const double DWs[], double voltage)
{
    QB_tools QB;
    QB_init(&QB);
    int state=QB_read_file(&QB, cell_path);
    if(QB_F_UNKNOWN==state){
        printf("[ERROR] Unrecognized file %s\n", cell_path);
        exit(1);
    }
    QSPG_refined(&QB, &QB, SYMPREC);
    ntype=QB.TypeNumber;
    callocate_2d(&type_name, ntype, 10, '\0');
    if(types[0]!=nullptr&&types[0][0]!='\0'){
        for(int i=0;i<ntype;i++){
            strcpy(type_name[i], types[i]);
        }
    }else if(QB.ele.n==ntype){
        for(int i=0;i<ntype;i++){
            strcpy(type_name[i], QB.ele.t[i].name);
        }
    }
    QSPGDataset *QD;
    SpacegroupType SG;
    Pointgroup PG;
    QSPG_get_symmetry(&QB, &QD, SYMPREC);
    SG=spgdb_get_spacegroup_type(QD->hall_number);
    PG=ptg_get_pointgroup(SG.pointgroup_number);

    const char CRYSTAL_SYS[][20]={"None", "Triclinic", "Monoclinic", "Orthorhombic", "Tetragonal", "Trigonal", "Hexagonal", "Cubic"};
    strcpy(crystal_system, CRYSTAL_SYS[PG.holohedry]);
    const char CENTERING[]={'0', 'P', 'I', 'F', 'A', 'B', 'C', '0', 'R'};
    int centeringi=SG.centering; centering=CENTERING[centeringi];
    hall_number=QD->hall_number;
    space_group=QD->spacegroup_number;
    strcpy(space_group_symbol, QD->international_symbol);
    point_group=SG.pointgroup_number;
    strcpy(point_group_symbol, QD->pointgroup_symbol);
    set_sampling_type(space_group, point_group);
    if(PG.holohedry==TRIGO) is_trigonal=true;
    if(PG.holohedry==HEXA) is_hexagonal=true;
    if(is_trigonal||is_hexagonal) use_hexagonal=true;

    compute_lattice_matrices(QD->std_lattice);
    compute_symmetry_matrices(space_group, QD->rotations, QD->translations, QD->n_operations);
    compute_asymmetric_atomic_positions(QD->std_positions, QD->std_types, QD->n_std_atoms);
    callocate(&apos_Z, napos, -1);
    callocate(&apos_M, napos, 0.0);
    callocate(&apos_DW, napos, 0.0);
    callocate(&apos_occupation, napos, 1.0);
    for(int i=0;i<napos;i++){
        int typei=apos_type[i]-1;
        npos+=apos_multi[i];
        apos_DW[i]=DWs[typei]*0.01;
        for(int j=0;j<E_TYPE_NUMBER;j++){
            if(0==strcmp(type_name[typei], E_TYPE[j])){
                apos_Z[i]=j+1;
                apos_M[i]=E_MASS[j];
                break;
            }
        }
    }
    compute_atomic_density();
    QSPGDataset_free_dataset(QD);
    QB_free_atom(&QB);
    update_Fourier_coefficient0(voltage);
    logging();
}

void CELL::set_sampling_type(int sgnum, int pgnum)
{
    sampling_type=PG_SAMPLING_TYPE[pgnum-1];
    if(-1==sampling_type){
        if(sgnum>=143&&sgnum<=167){
            if(sgnum>=143&&sgnum<=146) sampling_type=10;
            else if(sgnum>=147&&sgnum<=155) sampling_type=12;
            else if(sgnum>=156&&sgnum<=161) sampling_type=14;
            else if(sgnum>=162&&sgnum<=163) sampling_type=17;
            else if(sgnum>=164&&sgnum<=167) sampling_type=16;
        }else{
            if(14==pgnum){
                if(sgnum>=115&&sgnum<=120){
                    sampling_type=6;
                }else{
                    sampling_type=8;
                }
            }else if(26==pgnum){
                if(187==sgnum||188==sgnum){
                    sampling_type=16;
                }else{
                    sampling_type=17;
                }
            }
        }
    } 
}

void CELL::compute_lattice_matrices(double lat[3][3])
{
    double lattice[3][3];
    matrix_constant(lattice, 0.1, lat);
    double lx=lattice[0][0], ly=lattice[1][1], lz=lattice[2][2];
    double xy=lattice[0][1], xz=lattice[0][2], yz=lattice[1][2];
    a0=lx; b0=sqrt(ly*ly+xy*xy); c0=sqrt(lz*lz+xz*xz+yz*yz);
    alpha=acos((xy*xz+ly*yz)/(b0*c0))*RAD_TO_DEG; beta=acos(xz/c0)*RAD_TO_DEG; gamma=acos(xy/b0)*RAD_TO_DEG;
    double a[3], b[3], c[3];
    double bc[3], ca[3], ab[3];
    matrix_copy(dsm, lattice);//*point
    matrix_transpose(dsm, dsm);
    vector_copy(a, dsm[0]);
    vector_copy(b, dsm[1]);
    vector_copy(c, dsm[2]);
    vector_cross(bc, b, c);
    vector_cross(ca, c, a);
    vector_cross(ab, a, b);
    vol=vector_dot(a, bc);
    vector_constant(bc, 1.0/vol, bc);
    vector_constant(ca, 1.0/vol, ca);
    vector_constant(ab, 1.0/vol, ab);
    vector_copy(rsm[0], bc);
    vector_copy(rsm[1], ca);
    vector_copy(rsm[2], ab);
    matrix_transpose(rsm, rsm);
    double t_mat[3][3];
    matrix_transpose(t_mat, dsm);
    matrix_multiply(dmt, dsm, t_mat);
    matrix_transpose(t_mat, rsm);
    matrix_multiply(rmt, t_mat, rsm);
}

void CELL::compute_symmetry_matrices(int sgnum, int (*rots)[3][3], double (*trans)[3], int noperation)
{
    nsymmetry=noperation;
    for(int i=0;i<noperation;i++){
        matrix_copy(sym_rotations[i], rots[i]);
        vector_copy(sym_translations[i], trans[i]);
    }
    int    num=0;
    int    rotations[192][3][3]; 
    double translations[192][3];
    if(sgnum!=SG_SYMMORPHIC_NUMBER[sgnum-1]){
        sgnum=SG_SYMMORPHIC_NUMBER[space_group-1];
        Symmetry *sym=spgdb_get_spacegroup_operations(SG_HALL_NUMBER_START[sgnum-1]);
        for(int i=0;i<sym->size;i++){
            matrix_copy(rotations[i], sym->rot[i]);
            vector_copy(translations[i], sym->trans[i]);
        }
        num=sym->size;
    }else{
        for(int i=0;i<nsymmetry;i++){
            matrix_copy(rotations[i], sym_rotations[i]);
            vector_copy(translations[i], sym_translations[i]);
        }
        num=nsymmetry;
    }
    int count_point=0;
    for(int i=0;i<num;i++){
        double temp1=vector_dot(translations[i], translations[i]);
        if(temp1<0.1){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    point_dmats[count_point][j][k]=rotations[i][j][k];
                }
            }
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    double temp2=0.0;
                    for(int m=0;m<3;m++){
                        for(int n=0;n<3;n++){
                            temp2+=dmt[j][m]*rotations[i][m][n]*rmt[n][k];
                        }
                    }
                    point_rmats[count_point][j][k]=temp2;
                }
            }
            count_point++;
        }
    }
    npointsym=count_point;
}

void CELL::compute_asymmetric_atomic_positions(double (*atom_pos)[3], int *atom_type, int natom)
{
    int *atom_mapping_index, *apos_index;
    callocate(&atom_mapping_index, natom, -1);
    callocate(&apos_index, natom, -1);
    int count=0, count_apos=0;
    while(count<natom){
        if(atom_mapping_index[count]!=-1){
            count++;
            continue;
        }
        apos_index[count_apos]=count;
        for(int i=0;i<nsymmetry;i++){
            double pos[3];
            bool is_found=false;
            vector_rotate(pos, sym_rotations[i], atom_pos[count]);
            vector_plus(pos, pos, sym_translations[i]);
            for(int j=0;j<natom;j++){
                if(atom_type[j]!=atom_type[count]){
                    continue;
                }
                double vec[3];
                vector_difference(vec, pos, atom_pos[j]);
                while(vec[0]<-0.5) vec[0]+=1;
                while(vec[1]<-0.5) vec[1]+=1;
                while(vec[2]<-0.5) vec[2]+=1;
                while(vec[0]>0.5) vec[0]-=1;
                while(vec[1]>0.5) vec[1]-=1;
                while(vec[2]>0.5) vec[2]-=1;
                if((vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])<SYMPREC2){
                    atom_mapping_index[j]=count;
                    is_found=true;
                    break;
                }
            }
            if(!is_found){
                printf("[WARNING] Position (%.8f, %.8f, %.8f) arising from asymmetric position (%.8f, %.8f, %.8f) (type %d) is beyond the unit cell", 
                       pos[0], pos[1], pos[2], atom_pos[count][0], atom_pos[count][1], atom_pos[count][2], atom_type[count]);
            }
        }
        count++; count_apos++;
    }
    for(int i=0;i<natom;i++){
        if(atom_mapping_index[i]==-1){
            printf("[WARNING] Position (%.8f, %.8f, %.8f) in the unit cell is not mapped to any asymmetric position", 
                    atom_pos[count][0], atom_pos[count][1], atom_pos[count][2]);
        }
    }

    napos=count_apos;
    callocate(&apos_type, napos, 0);
    callocate(&apos_multi, napos, 0);
    callocate_3d(&apos_pos, napos, nsymmetry, 3, 0.0);
    for(int i=0;i<napos;i++){
        apos_type[i]=atom_type[apos_index[i]];
        for(int j=0;j<natom;j++){
            if(apos_index[i]==atom_mapping_index[j]){
                vector_copy(apos_pos[i][apos_multi[i]], atom_pos[j]);
                apos_multi[i]++;
            }
        }
    }
    deallocate(atom_mapping_index);
    deallocate(apos_index);
}

void CELL::compute_atomic_density()
{
    double num_M=0.0, num_Z=0.0;
    ave_M=0.0; ave_Z=0.0;
    for(int i=0;i<napos;i++){
        int Z=apos_Z[i];
        double M=apos_M[i];
        ave_M+=M*apos_multi[i]*apos_occupation[i]; ave_Z+=Z*apos_multi[i];
        num_M+=apos_multi[i]*apos_occupation[i]; num_Z+=apos_multi[i];
    }
    density=ave_M/(vol*1.0e-21*AVOGADRO_CONSTANT);
    ave_M/=num_M; ave_Z/=num_Z;
}

void CELL::logging()
{
    printf("[INFO] LATTICE INFORMATION\n");
    printf("[INFO] Crystal system     : %s\n", crystal_system);
    printf("[INFO] Centering vector   : %c\n", centering);
    printf("[INFO] a [nm]             : %.5f\n", a0);
    printf("[INFO] b [nm]             : %.5f\n", b0);
    printf("[INFO] c [nm]             : %.5f\n", c0);
    printf("[INFO] alpha [deg]        : %.5f\n", alpha);
    printf("[INFO] beta  [deg]        : %.5f\n", beta);
    printf("[INFO] gamma [deg]        : %.5f\n", gamma);
    printf("[INFO] volume[nm^3]       : %.5f\n", vol);
    printf("[INFO] SYMMETRY INFORMATION\n");
    printf("[INFO] Hall number        : %d\n", hall_number);
    printf("[INFO] Space group        : %d\n", space_group);
    printf("[INFO] Space group symbol : %s\n", space_group_symbol);
    printf("[INFO] Point group        : %d\n", point_group);
    printf("[INFO] Point group symbol : %s\n", point_group_symbol);
    printf("[INFO] Sampling type      : %d\n", sampling_type);
    printf("[INFO] Number of symmetry operations      : %d\n", nsymmetry);
    printf("[INFO] Number of point symmetry operations: %d\n", npointsym);
    printf("[INFO] ATOMIC INFORMATION\n");
    printf("[INFO] Number of asymmetric atomic positions: %d\n", napos);
    for(int i=0;i<napos;i++){
        printf("[INFO] General atomic position, atomic number, atomic mass, multiplicity, site occupation, Debye-Waller factor: %d, %d, %.5f, %d, %.5f, %.5f\n", i+1, apos_Z[i], apos_M[i], apos_multi[i], apos_occupation[i], apos_DW[i]);
        printf("[INFO] Equivalent atomic positions (x, y, z):\n");
        for(int j=0;j<apos_multi[i];j++){
            printf("[INFO] > %.5f, %.5f, %.5f\n", apos_pos[i][j][0], apos_pos[i][j][1], apos_pos[i][j][2]);
        }
    }
    printf("[INFO] Number of atomic positions: %d\n", npos);
    printf("[INFO] Density [in g/cm^3], atomic number averaged, atomic mass averaged [g/mol] = %.5f, %.5f, %.5f\n", density, ave_Z, ave_M);
    printf("[INFO] FOURIER INFORMATION\n");
    printf("[INFO] Mean inner potential [V]:                %.5f\n", fouri0.Vmod);
    printf("[INFO] Relativistic correction factor [V]:      %.5f\n", fouri0.gamma);
    printf("[INFO] Relativistic Accelerating Potential [V]: %.5f\n", fouri0.voltage);
    printf("[INFO] Electron Wavelength [nm]:                %.5f\n", fouri0.lambda);
    printf("[INFO] Interaction constant [V nm^-1]:          %.5f\n", fouri0.sigma);
    printf("[INFO] Normal absorption length [nm]:           %.5f\n", fouri0.sigp);
}

CELL::~CELL()
{
    if(napos!=0){
        deallocate(apos_type);
        deallocate(apos_multi);
        deallocate_3d(apos_pos, napos, nsymmetry);
        deallocate(apos_Z);
        deallocate(apos_M);
        deallocate(apos_DW);
        deallocate(apos_occupation);
    }
    if(ntype!=0) deallocate_2d(type_name, ntype);
}

double CELL::dot(double v1[3], double v2[3], char space)
{
    double temp[3], res;
    switch(space)
    {
    case 'r':
        vector_rotate(temp, rmt, v2);
        res=vector_dot(v1, temp);
        break;
    case 'd':
        vector_rotate(temp, dmt, v2);
        res=vector_dot(v1, temp);
    default:
        printf("[ERROR] Unrecognized space %c in dot computation.", space);
        exit(1);
    }
    return res;
}

void CELL::cross(double c_v[3], double v1[3], double v2[3], char space)
{
    double temp[3], res;
    switch(space)
    {
    case 'r':
        vector_cross(c_v, v1, v2);
        vector_transform(c_v, c_v, dmt);
        break;
    case 'd':
        vector_cross(c_v, v1, v2);
        vector_transform(c_v, c_v, rmt);
        break;
    default:
        printf("[ERROR] Unrecognized space %c in cross computation", space);
        exit(1);
    }
}

double CELL::length(double v[3], char space)
{
    return sqrt(dot(v, v, space));
}

double CELL::angle(double v1[3], double v2[3], char space)
{
    double dot12=dot(v1, v2, space), len1=length(v1, space), len2=length(v2, space);
    if(0.0==len1||0.0==len2){
        printf("[ERROR] Zero length for vector [%d %d %d] or [%d %d %d] in angle computation.", 
               int(v1[0]), int(v1[1]), int(v1[2]), int(v2[0]), int(v2[1]), int(v2[2]));
        exit(1);
    }
    double res=dot12/(len1*len2);
    if(res>=1.0) return 0.0;
    else if(res<=-1.0) return PI;
    else return acos(res);
}

void CELL::normalize(double n_v[3], double v[3], char space)
{
    double len_v=length(v, space);
    if(0.0!=len_v){
        vector_constant(n_v, 1.0/len_v, v);
    }else{
        n_v[0]=0.0; n_v[1]=0.0; n_v[2]=0.0; 
    }
}

void CELL::reciprocal_to_cartesian(double c_v[3], double v[3])
{
    vector_rotate(c_v, rsm, v);
}

void CELL::cartesian_to_reciprocal(double r_v[3], double v[3])
{
    vector_rotate(r_v, dsm, v);
}

void CELL::direct_to_reciprocal(double r_v[3], double v[3])
{
    vector_transform(r_v, v, dmt);
}

void CELL::direct_to_cartesian(double c_v[3], double v[3])
{
    vector_rotate(c_v, dsm, v);
}

void CELL::compute_equivalent_reciprocal_vectors(double equiv[48][3], int &nequiv, double g[3], char space)
{
    double s[3]={0.0};
    vector_copy(equiv[0], g);
    nequiv=1;
    for(int i=1;i<npointsym;i++){
        for(int j=0;j<3;j++){
            if('r'==space){
                s[j]=vector_dot(point_rmats[i][j], g);
            }else if('d'==space){
                s[j]=vector_dot(point_dmats[i][j], g);
            }else{
                printf("[ERROR] Unrecognized space %c in star computation.", space);
                exit(1);
            }
        }
        bool is_new=true;
        for(int j=0;j<nequiv;j++){
            double diff=0.0;
            for(int k=0;k<3;k++){
                diff+=fabs(equiv[j][k]-s[k]);
            }
            if(diff<=1e-4){
                is_new=false;
            }
        }
        if(is_new){
            vector_copy(equiv[nequiv], s);
            nequiv++;
        }
    }
}

void CELL::apply_point_group_symmetry(int equiv[48][3], int &nequiv, int px, int py, int hemisphere, int impx, int impy)
{
    double xy[2]={double(px)/double(impx), double(py)/double(impy)};
    double xyz[3]; int ierr;
    if(use_hexagonal){
        compute_sphere_from_hexagonal_Lambert(xyz, ierr, xy);
    }else{
        compute_sphere_from_square_Lambert(xyz, ierr, xy);
    }
    if(hemisphere<0) xyz[2]=-xyz[2];
    vector_normalize(xyz, xyz);
    double kstar[3];
    cartesian_to_reciprocal(kstar, xyz);
    int iequiv; double vtmp[48][3];
    switch(sampling_type)
    {
    case 3:
        iequiv=2;
        vtmp[0][0]= kstar[0]; vtmp[0][1]= kstar[1]; vtmp[0][2]= kstar[2];
        vtmp[1][0]=-kstar[0]; vtmp[1][1]= kstar[1]; vtmp[1][2]=-kstar[2];
        break;
    case 6:
        iequiv=8;
        vtmp[0][0]= kstar[0]; vtmp[0][1]= kstar[1]; vtmp[0][2]= kstar[2];
        vtmp[1][0]=-kstar[0]; vtmp[1][1]=-kstar[1]; vtmp[1][2]= kstar[2];
        vtmp[2][0]=-kstar[0]; vtmp[2][1]= kstar[1]; vtmp[2][2]=-kstar[2];
        vtmp[3][0]=-kstar[0]; vtmp[3][1]=-kstar[1]; vtmp[3][2]=-kstar[2];
        vtmp[4][0]= kstar[0]; vtmp[4][1]=-kstar[1]; vtmp[4][2]=-kstar[2];
        vtmp[5][0]= kstar[0]; vtmp[5][1]= kstar[1]; vtmp[5][2]=-kstar[2];
        vtmp[6][0]= kstar[0]; vtmp[6][1]=-kstar[1]; vtmp[6][2]= kstar[2];
        vtmp[7][0]=-kstar[0]; vtmp[7][1]= kstar[1]; vtmp[7][2]= kstar[2];
        break;
    case 8:
        iequiv=8;
        vtmp[0][0]= kstar[0]; vtmp[0][1]= kstar[1]; vtmp[0][2]= kstar[2];
        vtmp[1][0]=-kstar[0]; vtmp[1][1]=-kstar[1]; vtmp[1][2]= kstar[2];
        vtmp[2][0]=-kstar[0]; vtmp[2][1]= kstar[1]; vtmp[2][2]=-kstar[2];
        vtmp[3][0]= kstar[0]; vtmp[3][1]=-kstar[1]; vtmp[3][2]=-kstar[2];
        vtmp[4][0]= kstar[1]; vtmp[4][1]=-kstar[0]; vtmp[4][2]=-kstar[2];
        vtmp[5][0]=-kstar[1]; vtmp[5][1]= kstar[0]; vtmp[5][2]=-kstar[2];
        vtmp[6][0]=-kstar[1]; vtmp[6][1]=-kstar[0]; vtmp[6][2]= kstar[2];
        vtmp[7][0]= kstar[1]; vtmp[7][1]= kstar[0]; vtmp[7][2]= kstar[2];
        break;
    case 9:
        iequiv=16;
        vtmp[0][0]= kstar[0];  vtmp[0][1]= kstar[1];  vtmp[0][2]= kstar[2];
        vtmp[1][0]=-kstar[0];  vtmp[1][1]=-kstar[1];  vtmp[1][2]= kstar[2];
        vtmp[2][0]=-kstar[0];  vtmp[2][1]= kstar[1];  vtmp[2][2]=-kstar[2];
        vtmp[3][0]=-kstar[0];  vtmp[3][1]=-kstar[1];  vtmp[3][2]=-kstar[2];
        vtmp[4][0]= kstar[0];  vtmp[4][1]=-kstar[1];  vtmp[4][2]=-kstar[2];
        vtmp[5][0]= kstar[0];  vtmp[5][1]= kstar[1];  vtmp[5][2]=-kstar[2];
        vtmp[6][0]= kstar[0];  vtmp[6][1]=-kstar[1];  vtmp[6][2]= kstar[2];
        vtmp[7][0]=-kstar[0];  vtmp[7][1]= kstar[1];  vtmp[7][2]= kstar[2];
        
        vtmp[8][0]=-kstar[1];  vtmp[8][1]= kstar[0];  vtmp[8][2]= kstar[2];
        vtmp[9][0]= kstar[1];  vtmp[9][1]=-kstar[0];  vtmp[9][2]= kstar[2];
        vtmp[10][0]= kstar[1]; vtmp[10][1]= kstar[0]; vtmp[10][2]=-kstar[2];
        vtmp[11][0]= kstar[1]; vtmp[11][1]=-kstar[0]; vtmp[11][2]=-kstar[2];
        vtmp[12][0]=-kstar[1]; vtmp[12][1]=-kstar[0]; vtmp[12][2]=-kstar[2];
        vtmp[13][0]=-kstar[1]; vtmp[13][1]= kstar[0]; vtmp[13][2]=-kstar[2];
        vtmp[14][0]= kstar[1]; vtmp[14][1]= kstar[0]; vtmp[14][2]= kstar[2];
        vtmp[15][0]=-kstar[1]; vtmp[15][1]=-kstar[0]; vtmp[15][2]= kstar[2];
        break;
    default:
        compute_equivalent_reciprocal_vectors(vtmp, iequiv, kstar);
        break;
    }
    double eps=-1.0e-4;
    for(int i=0;i<iequiv;i++){
        reciprocal_to_cartesian(xyz, vtmp[i]);
        vector_normalize(xyz, xyz);
        if(xyz[2]<eps){
            equiv[i][2]=-1;
        }else{
            equiv[i][2]=1;
        }
        if(use_hexagonal){
            compute_hexagonal_Lambert(xy, ierr, xyz);
        }else{
            compute_square_Lambert(xy, ierr, xyz);
        }
        xy[0]*=double(impx); xy[1]*=double(impy);
        equiv[i][0]=int(round(xy[0])); equiv[i][1]=int(round(xy[1]));
    }
    nequiv=iequiv;
}

double CELL::get_interplanar_spacing(double g[3])
{
    double dotgg=dot(g, g, 'r');
    if(dotgg<=0.0){
        printf("[ERROR] Zero reciprocal vector in interplanar spacing compution.");
        exit(1);
    }
    return 1.0/sqrt(dotgg);
}

void CELL::compute_reflection_range(double dmin)
{
    for(int i=0;i<3;i++){
        double dhkl;
        double hkl[3]={0.0, 0.0, 0.0};
        HKL[i]=0;
        do{
            HKL[i]++;
            hkl[i]=double(HKL[i]);
            dhkl=get_interplanar_spacing(hkl);
        }while(dhkl>=dmin);
    }
}

double CELL::get_excitation_error(double g[3], double k[3], double fn[3])
{
    double kpg[3]={k[0]+g[0], k[1]+g[1], k[2]+g[2]};
    double tkpg[3]={2.0*k[0]+g[0], 2.0*k[1]+g[1], 2.0*k[2]+g[2]};
    double q1=length(kpg, 'r'), q2=angle(kpg, fn, 'r');
    double res=-1.0*dot(g, tkpg, 'r')/(2.0*q1*cos(q2));
    return res;
}

bool CELL::is_centering_allowed(double g[3])
{
//Determine whether or not the reciprocal vector is allowed by the lattice centering
    bool is_allowed=true;
    int res;
    switch(centering)
    {
    case 'P':
        break;//all reflections allowed for a primitive lattice
    case 'F':
        res=int(g[0]+100)%2+int(g[1]+100)%2+int(g[2]+100)%2;
        if(1==res||2==res) is_allowed=false;
        break;
    case 'I':
        res=int(g[0]+g[1]+g[2]+100)%2;
        if(1==res) is_allowed=false;
        break;
    case 'A':
        res=int(g[1]+g[2]+100)%2;
        if(1==res) is_allowed=false;
        break;
    case 'B':
        res=int(g[0]+g[2]+100)%2;
        if(1==res) is_allowed=false;
        break;
    case 'C':
        res=int(g[0]+g[1]+100)%2;
        if(1==res) is_allowed=false;
        break;
    case 'R':
        if(is_hexagonal){
            res=int(-g[0]+g[1]+g[2]+90)%3;
            if(0!=res) is_allowed=false;
        }
        break;
    default:
        printf("[ERROR] Unrecognized centering symbol %c (not P, F, I, A, B, C, or R).", centering);
        exit(1);
    }
    return is_allowed;
}

void CELL::compute_Bloch_wave_coefficients(bool is_initial)
{
    int imh=HKL[0]*2, imk=HKL[1]*2, iml=HKL[2]*2; 
    if(is_initial){
        callocate_4d(&LUTSgh, napos, imh*2+1, imk*2+1, iml*2+1, complex<double>(0.0, 0.0));
    }
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                double g[3]={double(ih), double(ik), double(il)};
                if(is_centering_allowed(g)){
                    int ix=ih+imh, iy=ik+imk, iz=il+iml;
                    for(int i=0;i<napos;i++){
                        LUTSgh[i][ix][iy][iz]=complex<double>(0.0, 0.0);
                        double DBWF=pow(apos_Z[i], 2)*apos_occupation[i]*exp(-0.25*apos_DW[i]*dot(g, g, 'r'));
                        for(int j=0;j<apos_multi[i];j++){
                            double q=TWO_PI*(g[0]*apos_pos[i][j][0]+g[1]*apos_pos[i][j][1]+g[2]*apos_pos[i][j][2]);
                            complex<double> F(cos(q), sin(q));
                            LUTSgh[i][ix][iy][iz]+=F*complex<double>(DBWF, 0.0);
                        }
                    }
                }
            }
        }
    }
}

void CELL::compute_Fourier_coefficients(double voltage, bool is_initial)
{
    int imh=2*HKL[0], imk=2*HKL[1], iml=2*HKL[2];
    if(is_initial){
        int size1=imh*2+1, size2=imk*2+1, size3=iml*2+1;
        callocate_3d(&LUTUg, size1, size2, size3, complex<double>(0.0, 0.0));
        callocate_3d(&LUTqg, size1, size2, size3, complex<double>(0.0, 0.0));
        callocate_3d(&is_double_diffrac, size1, size2, size3, false);
    }
    LUTUg[imh][imk][iml]=fouri0.Ug; LUTqg[imh][imk][iml]=fouri0.qg;

    double g[3];
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                g[0]=double(ih); g[1]=double(ik); g[2]=double(il);
                if(is_centering_allowed(g)){
                    update_Fourier_coefficient(voltage, g, true);
                    int i=imh+ih, j=imk+ik, k=iml+il;
                    LUTUg[i][j][k]=fouri.Ug; LUTqg[i][j][k]=fouri.qg;
                    if(abs(fouri.Ug)<=1e-5){
                        is_double_diffrac[i][j][k]=true;
                    }
                }
            }
        }
    }
}

void CELL::update_Fourier_coefficient0(double voltage)
{
    double g[3]={0.0, 0.0, 0.0};
    update_Fourier_coefficient(voltage, g);
    double temp=1000.0*ELECTRON_CHARGE*voltage/ELECTRON_REST_MASS/pow(LIGHT_VELOCITY, 2);
    fouri0.Vmod=fouri.Vmod; fouri0.Vpmod=fouri.Vpmod;
    fouri0.Umod=fouri.Umod; fouri0.Upmod=fouri.Upmod;
    fouri0.Vg=fouri.Vg; fouri0.Ug=fouri.Ug;
    fouri0.gamma=1.0+temp;
    fouri0.voltage=voltage*(1.0+0.5*temp)*1000.0+fouri0.gamma*fouri0.Vmod;
    fouri0.lambda=1.0e9*PLANK_CONSTANT/sqrt(2.0*ELECTRON_REST_MASS*ELECTRON_CHARGE*fouri0.voltage);
    fouri0.sigma=1.0e-18*TWO_PI*ELECTRON_REST_MASS*ELECTRON_CHARGE*fouri0.gamma*fouri0.lambda/pow(PLANK_CONSTANT, 2);
    update_Fourier_coefficient(voltage, g, true);
    fouri0.sig=fouri.sig; fouri0.sigp=fouri.sigp; fouri0.qg=fouri.qg;
}

void CELL::update_Fourier_coefficient(double voltage, double g[3], bool is_shift)
{
    complex<double> FV(0.0, 0.0), FVp(0.0, 0.0);
    double GWK;
    if(0.0==(g[0]*g[0]+g[1]*g[1]+g[2]*g[2])){
        GWK=0.0;
    }else{
        GWK=length(g, 'r')*G_TO_WK;
    }
    for(int i=0;i<napos;i++){//loop over atoms in the asymmetric unit
        double UWK=sqrt(apos_DW[i]*DW_TO_WK)+1.0e-12;//root-mean-square value of the atomic vibration amplitude or thermal displacement in [nm]
        int    Z=apos_Z[i];
        complex<double> Fscat=WK_get_scattering_amplitude(GWK, UWK, Z, voltage)*complex<double>(apos_occupation[i], 0.0);
        complex<double> Fpos(0.0, 0.0);
        for(int j=0;j<apos_multi[i];j++){//loop over atoms in orbit
            double q=-1.0*TWO_PI*(g[0]*apos_pos[i][j][0]+g[1]*apos_pos[i][j][1]+g[2]*apos_pos[i][j][2]);
            complex<double> jFpos(cos(q), sin(q));
            Fpos+=jFpos;
        }
        FV+=Fscat.real()*Fpos; FVp+=Fscat.imag()*Fpos;

    }
    double constV=PRE_CONST_V/vol/FOUR_PI, constU=PRE_CONST_U;
    fouri.Vmod=constV*abs(FV); fouri.Vang=arg(FV);
    fouri.Vpmod=constV*abs(FVp); fouri.Vpang=arg(FVp);
    fouri.Umod=constU*fouri.Vmod; fouri.Upmod=constU*fouri.Vpmod;
    fouri.Vg.real(constV*(FV.real()-FVp.imag())); fouri.Vg.imag(constV*(FV.imag()+FVp.real()));
    fouri.Ug=constU*fouri.Vg;
    if(fabs(fouri.Umod)>0.0){
        fouri.sig=1.0/fabs(fouri.Umod)/fouri0.lambda;
    }else{
        fouri.sig=1.0e8;
    }
    if(fabs(fouri.Upmod)>0.0){
        fouri.sigp=1.0/fabs(fouri.Upmod)/fouri0.lambda;
    }else{
        fouri.sigp=1.0e8;
    }
    if(is_shift){
        fouri.qg.real(cos(fouri.Vang)/fouri.sig-sin(fouri.Vpang)/fouri.sigp); 
        fouri.qg.imag(cos(fouri.Vpang)/fouri.sigp+sin(fouri.Vang)/fouri.sig);
    }else{
        double betag=fouri.Vpang-fouri.Vang;
        fouri.qg.real(1.0/fouri.sig-sin(betag)/fouri.sigp); 
        fouri.qg.imag(cos(betag)/fouri.sigp);
    }
}

MODEL::MODEL(const char *model_path, const char types[][10], const double DWs[])
{
    QB_tools QB;
    QB_init(&QB);
    int state=QB_read_file(&QB, model_path);
    if(QB_F_UNKNOWN==state){
        printf("[ERROR] Unrecognized file %s\n", model_path);
        exit(1);
    }
    set_lattice(QB.mat);
    is_periodic[0]=QB.px; is_periodic[1]=QB.py; is_periodic[2]=QB.pz;
    ntype=QB.TypeNumber;
    natom=QB.TotalNumber;
    mallocate_2d(&atom_pos, natom, 3);
    mallocate(&atom_type, natom);
    mallocate(&atom_DW, natom);
    for(int i=0;i<natom;i++){
        atom_pos[i][0]=QB.atom[i].x; 
        atom_pos[i][1]=QB.atom[i].y; 
        atom_pos[i][2]=QB.atom[i].z;
        atom_type[i]=QB.atom[i].type;
        atom_DW[i]=DWs[atom_type[i]-1];
    }
    callocate(&atom_Z, natom, -1);
    callocate(&atom_occupation, natom, 1.0);
    callocate(&type_index, ntype, -1);
    callocate_2d(&type_name, ntype, 10, '\0');
    if(QB.ele.n==ntype){
        for(int i=0;i<ntype;i++){
            strcpy(type_name[i], QB.ele.t[i].name);
        }
    }
    QB_free_atom(&QB);
}

void MODEL::set_lattice(double mat[3][3])
{
    double lx=mat[0][0], ly=mat[1][1], lz=mat[2][2];
    double a[3]={mat[0][0], mat[0][1], mat[0][2]}, b[3]={mat[1][0], mat[1][1], mat[1][2]}, c[3]={mat[2][0], mat[2][1], mat[2][2]};
    double la=vector_length(a), lb=vector_length(b), lc=vector_length(c);
    if(fabs(lx-la)<1.0e-6&&fabs(ly-lb)<1.0e-6&&fabs(lz-lc)<1.0e-6){
        is_orthogonal=true;
        vol=lx*ly*lz;
        double E[3][3]={{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        matrix_copy(dsm, E); matrix_copy(rsm, E);
        matrix_copy(dmt, E); matrix_copy(rmt, E);
        dimension[0]=lx; dimension[1]=ly; dimension[2]=lz;
        dimensionK[0]=1.0/lx; dimensionK[1]=1.0/ly; dimensionK[2]=1.0/lz;
    }else{
        is_orthogonal=false;
        double bc[3], ca[3], ab[3];
        vector_cross(bc, b, c);
        vector_cross(ca, c, a);
        vector_cross(ab, a, b);
        vol=vector_dot(a, bc);
        vector_constant(bc, 1.0/vol, bc);
        vector_constant(ca, 1.0/vol, ca);
        vector_constant(ab, 1.0/vol, ab);
        vector_copy(dsm[0], a);
        vector_copy(dsm[1], b);
        vector_copy(dsm[2], c);
        vector_copy(rsm[0], bc);
        vector_copy(rsm[1], ca);
        vector_copy(rsm[2], ab);
        for(int i=0;i<3;i++){
            vector_normalize(dsm[i], dsm[i]);
            vector_normalize(rsm[i], rsm[i]);
        }
        matrix_transpose(dsm, dsm);
        matrix_transpose(rsm, rsm);
        double t_mat[3][3];
        matrix_transpose(t_mat, dsm);
        matrix_multiply(dmt, dsm, t_mat);
        matrix_transpose(t_mat, rsm);
        matrix_multiply(rmt, t_mat, rsm);
        dimension[0]=la; dimension[1]=lb; dimension[2]=lc;
        dimensionK[0]=vector_length(bc); dimensionK[1]=vector_length(ca); dimensionK[2]=vector_length(ab);
    }
}

MODEL::~MODEL()
{
    if(natom>0){
        deallocate_2d(atom_pos, natom);
        deallocate(atom_type);
        deallocate(atom_DW);
        deallocate(atom_Z);
        deallocate(atom_occupation);
    }
    if(ntype>0){
        deallocate(type_index);
        deallocate_2d(type_name, ntype);
    }
}

void MODEL::reciprocal_to_cartesian(double c_g[3], double r_g[3])
{
    vector_rotate(c_g, rsm, r_g);
}

void MODEL::direct_to_reciprocal(double r_g[3], double d_g[3])
{
    vector_transform(r_g, d_g, dmt);
}

double MODEL::get_reciprocal_vector_length(double g[3])
{
    double r_g[3];
    vector_rotate(r_g, rmt, g);
    double dot=vector_dot(g, r_g);
    return sqrt(dot);
}

void MODEL::compute_reciprocal_spacing(double spacing[3], double spacing_ratio[3])
{
    double dimK[3];
    if((!is_periodic[0])&&(!is_periodic[1])&&(!is_periodic[2])){
        for(int i=0;i<3;i++){
          dimK[i]=1.0;
        }
    }else{
        double aveK=0.0;
        for(int i=0;i<3;i++){
            if(is_periodic[i]){
                dimK[i]=dimensionK[i];
                aveK+=dimK[i];
            }
        }
        aveK=aveK/double(is_periodic[0]+is_periodic[1]+is_periodic[2]);
        for(int i=0;i<3;i++){
            if(!is_periodic[i]){
                dimK[i]=aveK;
            }        
        }
    }
    for(int i=0;i<3;i++){
        spacing[i]=dimK[i]*spacing_ratio[i];
    }
}

double MODEL::get_diffraction_intensity(double theta, double g[3])
{
    return -1.0;
}

double MODEL::get_diffraction_intensity(double G, double g[3], bool is_zero)
{
    return -1.0;
}

XMODEL::XMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda):MODEL(model_path, types, DWs)
{
    lambda=mlambda;
    strcpy(radiation, "X-ray");
    if(types[0]!=nullptr&&types[0][0]!='\0'){
        for(int i=0;i<ntype;i++){
            strcpy(type_name[i], types[i]);
        }
    }
    for(int i=0;i<ntype;i++){
        int j;
        for(j=0;j<X_TYPE_NUMBER;j++){
            if(0==strcmp(type_name[i], X_TYPE[j])){
                type_index[i]=j;
                break;
            }
        }
        if(j==X_TYPE_NUMBER){
            printf("[ERROR] Unrecognized type %s for x-ray diffraction", type_name[i]);
            exit(1);
        }
    }
}

XMODEL::~XMODEL()
{

}

double XMODEL::get_Debye_Waller_factor(double S, double DW)
{
    return exp(-DW*S*S);//S, scattering vector
}

double XMODEL::get_atomic_scattering_factor(double S, const double A[4], const double B[4], const double C)//xrd
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
    for(int i=0;i<natom;++i){
        int    t=type_index[atom_type[i]-1];
        const double *A=X_A[t], *B=X_B[t], C=X_C[t];
        double q=TWO_PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
        res+=get_Debye_Waller_factor(S, atom_DW[i])*get_atomic_scattering_factor(S, A, B, C)*complex<double>(cos(q), sin(q));
    }
    return res;
}

void XMODEL::update_lorentzP_type(int lp_type)
{
    lorentzP_type=lp_type;
}

double XMODEL::get_diffraction_intensity(double theta, double g[3])
{
    complex<double> F=get_atomic_structure_factor(theta, g);
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    switch(lorentzP_type)
    {
    case 0:
        break;
    case 1:
        res*=((1+cos(2*theta)*cos(2*theta))/2);
        break;
    case 2:
        res*=((1+cos(2*theta)*cos(2*theta))/2/sin(2*theta));
        break;
    case 3:
        res*=((1+cos(2*theta)*cos(2*theta))/(sin(theta)*sin(theta)*cos(theta)));
        break;
    default:
        printf("[ERROR] Unrecognized type %d of lorentz polarization factor", lorentzP_type);
        exit(1);
    }
    return res;
}

double XMODEL::get_diffraction_intensity(double G, double g[3], bool is_zero)
{
    return -1.0;
}

NMODEL::NMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda):MODEL(model_path, types, DWs)
{
    lambda=mlambda;
    strcpy(radiation, "neutron");
    if(types[0]!=nullptr&&types[0][0]!='\0'){
        for(int i=0;i<ntype;i++){
            strcpy(type_name[i], types[i]);
        }
    }
    for(int i=0;i<ntype;i++){
        int j;
        for(j=0;j<N_TYPE_NUMBER;j++){
            if(0==strcmp(type_name[i], N_TYPE[j])){
                type_index[i]=j;
                break;
            }
        }
        if(j==N_TYPE_NUMBER){
            printf("[ERROR] Unrecognized type %s for neutron diffraction", type_name[i]);
            exit(1);
        }
    }
}

NMODEL::~NMODEL()
{

}

double NMODEL::get_Debye_Waller_factor(double S, double DW)
{
    return exp(-DW*S*S);//S, scattering vector
}

complex<double> NMODEL::get_atomic_structure_factor(double theta, double g[3])//xrd
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    for(int i=0;i<natom;i++){
        int    t=type_index[atom_type[i]-1];
        double q=TWO_PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
        res+=get_Debye_Waller_factor(S, atom_DW[i])*N_N[t]*complex<double>(cos(q), sin(q));
    }
    return res;
}

void NMODEL::update_lorentzP_type(int lp_type)
{
    lorentzP_type=lp_type;
}

double NMODEL::get_diffraction_intensity(double theta, double g[3])
{
    complex<double> F=get_atomic_structure_factor(theta, g);
    double res=(F.real()*F.real()+F.imag()*F.imag())/natom;
    switch(lorentzP_type)
    {
    case 0:
        break;
    case 1:
        res*=(1.0/sin(2*theta));
        break;
    case 2:
        res*=(1.0/(sin(theta)*sin(theta)*cos(theta)));
        break;
    default:
        printf("[ERROR] Unrecognized type %d of lorentz polarization factor", lorentzP_type);
        exit(1);
    }
    return res;
}

double NMODEL::get_diffraction_intensity(double G, double g[3], bool is_zero)
{
    return -1.0;
}

EMODEL::EMODEL(const char *model_path, const char types[][10], const double DWs[], double mlambda):MODEL(model_path, types, DWs)
{
    lambda=mlambda;
    strcpy(radiation, "electron");
    if(types[0]!=nullptr&&types[0][0]!='\0'){
        for(int i=0;i<ntype;i++){
            strcpy(type_name[i], types[i]);
        }
    }
    for(int i=0;i<ntype;i++){
        int j;
        for(j=0;j<E_TYPE_NUMBER;j++){
            if(0==strcmp(type_name[i], E_TYPE[j])){
                type_index[i]=j;
                break;
            }
        }
        if(j==E_TYPE_NUMBER){
            printf("[ERROR] Unrecognized type %s for electron diffraction", type_name[i]);
            exit(1);
        }
    }
    for(int i=0;i<natom;i++){
        int index=type_index[atom_type[i]-1];
        atom_Z[i]=index+1;
    }
}

EMODEL::~EMODEL()
{

}

double EMODEL::get_atomic_scattering_factor(double S, const double A[5], const double B[5])//sed
{
    double res=0.0;
    for(int i=0;i<5;i++){
        res+=A[i]*exp(-B[i]*S*S);//S, scattering vector
    }
    return res;
}

complex<double> EMODEL::get_atomic_structure_factor(double theta, double g[3])//xrd
{
    complex<double> res(0.0, 0.0);
    double S=sin(theta)/lambda;
    if(0.0<=S<=2.0){
        for(int i=0;i<natom;i++){
            int    t=type_index[atom_type[i]-1];
            const double *A=E_GP_A1[t], *B=E_GP_B1[t];
            double q=TWO_PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
            res+=get_atomic_scattering_factor(S, A, B)*complex<double>(cos(q), sin(q));
        }
    }else if(2.0<S<=6.0){
        for(int i=0;i<natom;i++){
            int    t=type_index[atom_type[i]-1];
            const double *A=E_GP_A2[t], *B=E_GP_B2[t];
            double q=TWO_PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
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

double EMODEL::get_diffraction_intensity(double G, double g[3], bool is_zero)
{
    double theta=asin(0.5*lambda*G);
    return get_diffraction_intensity(theta, g);
}

VMODEL::VMODEL(const char *model_path, const char types[][10], const double DWs[], double mvoltage):MODEL(model_path, types, DWs)
{
    voltage=mvoltage;
    strcpy(radiation, "electron");
    if(types[0]!=nullptr&&types[0][0]!='\0'){
        for(int i=0;i<ntype;i++){
            strcpy(type_name[i], types[i]);
        }
    }
    for(int i=0;i<ntype;i++){
        int j;
        for(j=0;j<E_TYPE_NUMBER;j++){
            if(0==strcmp(type_name[i], E_TYPE[j])){
                type_index[i]=j;
                break;
            }
        }
        if(j==E_TYPE_NUMBER){
            printf("[ERROR] Unrecognized type %s for electron diffraction", type_name[i]);
            exit(1);
        }
    }
    for(int i=0;i<natom;i++){
        int index=type_index[atom_type[i]-1];
        atom_Z[i]=index+1;
    }
    set_wavelength();
}

VMODEL::~VMODEL()
{

}

void VMODEL::set_wavelength()
{
    double g[3]={0.0};
    get_diffraction_intensity(0.0, g, true);
}

double VMODEL::get_diffraction_intensity(double theta, double g[3])
{
    return -1.0;
}

double VMODEL::get_diffraction_intensity(double G, double g[3], bool is_zero)
{
    complex<double> FV(0.0, 0.0), FVp(0.0, 0.0);
    double GWK;
    if(is_zero){
        GWK=0.0;
        g[0]=g[1]=g[2]=0.0;
    }else{
        GWK=G*VG_TO_WK;
    }
    for(int i=0;i<natom;i++){
        double UWK=sqrt(atom_DW[i]*VDW_TO_WK)+1.0e-12;//avoiding UWK=0.0
        int    Z=atom_Z[i];
        complex<double> Fscat=WK_get_scattering_amplitude(GWK, UWK, Z, voltage)*complex<double>(atom_occupation[i], 0.0);
        double q=-1.0*TWO_PI*(g[0]*atom_pos[i][0]+g[1]*atom_pos[i][1]+g[2]*atom_pos[i][2]);
        complex<double> Fpos(cos(q), sin(q));
        FV+=Fscat.real()*Fpos; FVp+=Fscat.imag()*Fpos;
    }
    double constV=VPRE_CONST_V/vol/FOUR_PI;
    if(is_zero){
        double temp=1000.0*ELECTRON_CHARGE*voltage/ELECTRON_REST_MASS/pow(LIGHT_VELOCITY, 2.0);
        double mvoltage=1000.0*voltage*(1.0+0.5*temp)+constV*abs(FV)*(1.0+temp);
        lambda=1.0e10*PLANK_CONSTANT/sqrt(2.0*ELECTRON_REST_MASS*ELECTRON_CHARGE*mvoltage);
    }
    complex<double> Vg(constV*(FV.real()-FVp.imag()), constV*(FV.imag()+FVp.real()));
    return abs(Vg)*abs(Vg);
}
