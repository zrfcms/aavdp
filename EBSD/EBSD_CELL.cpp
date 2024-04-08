#include "EBSD_CELL.h"

void matrix_copy(double outmat[4][4], double inmat[4][4])
{
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            outmat[i][j]=inmat[i][j];
        }
    }
}

EBSD_CELL::EBSD_CELL(const char* file_path, bool logging_flag)
{
    size_t len=strlen(file_path);
    bool flag;
    if(len>=4&&strcmp(file_path+len-4, ".nml")==0){
        flag=read_parameters_from_nml(file_path);
    }else if(len>=3&&strcmp(file_path+len-3, ".h5")==0){
        flag=read_parameters_from_hdf5(file_path);
    }else if(len>=5&&strcmp(file_path+len-5, ".hdf5")==0){
        flag=read_parameters_from_hdf5(file_path);
    }else{
        printf("Error! Unrecognized file %s.", file_path);
        exit(EXIT_FAILURE);  
    }
    if(!flag){
        printf("Error! Unrecognized parameters in file %s.", file_path);
        exit(EXIT_FAILURE);
    }

    if(crystal_system==4) hexagonal_flag=true;
    if(crystal_system==5&&second_setting_flag) hexagonal_flag=true;
    if(crystal_system==5) trigonal_flag=true;

    int point_group=1;
    for(int i=0;i<POINT_GROUP_NUMBER;i++){
        if(PG_FIRST_SPACE_GROUP[i]<=point_group) point_group=i+1;
    }
    strcpy(generator_string, GENERATOR_STRING[space_group-1]);
    strcpy(space_group_symbol, SPACE_GROUP_SYMBOL[space_group-1]);
    centering_symbol=space_group_symbol[0];
    if(generator_string[0]=='1') centrosymmetric_flag=true;
    for(int i=0;i<SYMMORPHIC_SPACE_GROUP_NUMBER;i++){
        if(point_group==SYMMORPHIC_SPACE_GROUP[i]) nonsymmorphic_flag=true;
    }

    set_sampling_type();
    compute_space_matrices();
    compute_symmetry_matrices();
    compute_point_symmetry_matrices();
    compute_atomic_positions();
    compute_atomic_density();

    if(logging_flag){
        printf("-->Crystal Structure Information<--\n");
        printf("  a [nm]             : %.5f\n", a);
        printf("  b [nm]             : %.5f\n", b);
        printf("  b [nm]             : %.5f\n", c);
        printf("  alpha [deg]        : %.5f\n", alpha);
        printf("  beta  [deg]        : %.5f\n", beta);
        printf("  gamma [deg]        : %.5f\n", gamma);
        printf("  Volume [nm^3]      : %.8f\n", vol);
        printf("  Space group        : %d\n", space_group);
        printf("  Space group symbol : %s\n", space_group_symbol);
        printf("  Generator String   : %s\n", generator_string);
        if(second_setting_flag){
            if(crystal_system!=5){
                printf("  Using second origin setting\n");
            }else{
                printf("  Using rhombohedral/trigonal parameters\n");
            }
        }
        if(centrosymmetric_flag){
            printf("  Structure is centrosymmetric\n");
        }else{
            printf("  Structure is non-centrosymmetric\n");
        }
        printf("  Number of generators                : %d\n", ngenerator);
        printf("  Number of symmetry matrices         : %d\n", nsymmetry);
        printf("  Number of point symmetry matrices   : %d\n", npointsym);
        printf("  Number of asymmetric atomic positions : %d\n", napos);
        for(int i=0;i<napos;i++){
            printf("  General atomic position, atomic number, multiplicity                        : %d, %d, %d\n", i+1, apos_Z[i], apos_multi[i]);
            printf("  Equivalent atomic positions  (x, y, z, site occupation, Debye-Waller factor): \n");
            for(int j=0;j<apos_multi[i];j++){
                printf("         > %.5f, %.5f, %.5f, %.5f, %.5f\n", pos[i][j][0], pos[i][j][1], pos[i][j][2], apos_occupation[i], apos_DW[i]);
            }
        }
        printf("  Density, atomic number averaged, atomic mass averaged = %.5f, %.5f, %.5f\n", density, ave_Z, ave_M);
    }
}

EBSD_CELL::~EBSD_CELL()
{
    deallocate_3d(pos, napos, nsymmetry);
    deallocate_2d(apos, napos);
    deallocate(apos_Z);
    deallocate(apos_occupation);
    deallocate(apos_DW);
    deallocate(apos_multi);
    int size1=4*HKL[0]+1, size2=4*HKL[1]+1;
    if(LUTUg!=nullptr){
        deallocate_3d(LUTUg, size1, size2);
    }
    if(LUTqg!=nullptr){
        deallocate_3d(LUTqg, size1, size2);
    }
    if(is_double_diffrac!=nullptr){
        deallocate_3d(is_double_diffrac, size1, size2);
    }
    if(LUTSgh!=nullptr){
        deallocate_4d(LUTSgh, napos, size1, size2);
    }
}

bool EBSD_CELL::read_parameters_from_nml(const char* file_path)
{
    FILE *fp=fopen(file_path, "r");
    if(fp==NULL){
        printf("Error! Unable to open file %s.\n", file_path);
    }
    fseek(fp, 0, SEEK_SET);

    char keys[][MAX_LENTH_IN_NAME]={"lattice", "space_group", "atoms"}; int keyi=0;
    char readbuff[MAX_LENTH_IN_NAME];
    for(int i=0;(keyi<3)&&(i<MAX_WORD_NUMBER);i++){
        fscanf(fp, "%s", readbuff);
        if('#'==readbuff[0]){
            while(('\n'!=fgetc(fp))&&(!feof(fp)));
            continue;
        }
        if(0==strcmp(readbuff, keys[0])){
            fscanf(fp, "%d", &crystal_system);
            alpha=beta=gamma=90.0;
            switch(crystal_system){
                case 1: fscanf(fp, "%lf", &a); b=c=a; break;
                case 2: fscanf(fp, "%lf", &a); fscanf(fp, "%lf", &c); b=a; break;
                case 3: fscanf(fp, "%lf", &a); fscanf(fp, "%lf", &b); fscanf(fp, "%lf", &c); break;
                case 4: fscanf(fp, "%lf", &a); fscanf(fp, "%lf", &c); b=a; gamma=120.0; break;
                case 5: fscanf(fp, "%lf", &a); fscanf(fp, "%lf", &alpha); b=c=a; beta=gamma=alpha; break;
                case 6: fscanf(fp, "%lf", &a); fscanf(fp, "%lf", &b); fscanf(fp, "%lf", &c); fscanf(fp, "%lf", &beta); break;
                case 7: fscanf(fp, "%lf", &a); fscanf(fp, "%lf", &b); fscanf(fp, "%lf", &c); fscanf(fp, "%lf", &alpha); fscanf(fp, "%lf", &beta); fscanf(fp, "%lf", &gamma); break;
                default: printf("Error! Unrecognized crystal system %d.", crystal_system); exit(EXIT_FAILURE);
            }
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, keys[1])){
            int space_group_setting;
            fscanf(fp, "%d", &space_group);
            fscanf(fp, "%d", &space_group_setting);
            if(1==space_group_setting) second_setting_flag=false;
            else if(2==space_group_setting) second_setting_flag=true;
            else{
                printf("Error! Unrecognized space group setting %d.", space_group_setting);
                exit(EXIT_FAILURE);
            }
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, keys[2])){
            fscanf(fp, "%d", &napos);
            mallocate_2d(&apos, napos, 3);
            mallocate(&apos_Z, napos); mallocate(&apos_occupation, napos); mallocate(&apos_DW, napos); 
            for(int j=0;j<napos;j++){
                fscanf(fp, "%d", &apos_Z[j]);
                if((apos_Z[j]<1)||(apos_Z[j]>TYPE_NUMBER)){
                    printf("Error! Unrecognized atom number %d (not between 1 and 98).", apos_Z[j]);
                    exit(EXIT_FAILURE); 
                }
                for(int k=0;k<3;k++){
                    fscanf(fp, "%lf", &apos[j][k]);
                }
                fscanf(fp, "%lf", &apos_occupation[j]);
                if((apos_occupation[j]<0.0)||(apos_occupation[j]>1.0)){
                    printf("Error! Unrecognized site occupation %.5f (not between 0.0 and 1.0).", apos_occupation[j]);
                    exit(EXIT_FAILURE); 
                }
                fscanf(fp, "%lf", &apos_DW[j]);
            }
            keyi++;
            continue;
        }
    }
    if(3==keyi) return true;
    else return false;
}

bool EBSD_CELL::read_parameters_from_hdf5(const char* file_path)
{
    size_t size_T;
    double *lattice;
    size_t size_A[2];
    double **atomdata; 
    int space_group_setting;
    EBSD_HDF5 hdf;
    hdf.open(file_path);
    hdf.read("/CrystalData/CrystalSystem", crystal_system);
    hdf.read_array("/CrystalData/LatticeParameters", &lattice, size_T);
    hdf.read("/CrystalData/Natomtypes", napos);
    hdf.read_array_2d("/CrystalData/AtomData", &atomdata, size_A[0], size_A[1]);
    hdf.read_array("/CrystalData/Atomtypes", &apos_Z, size_T);
    hdf.read("/CrystalData/SpaceGroupNumber", space_group);
    hdf.read("/CrystalData/SpaceGroupSetting", space_group_setting);
    a=lattice[0]; b=lattice[1]; c=lattice[2]; alpha=lattice[3]; beta=lattice[4]; gamma=lattice[5];
    mallocate_2d(&apos, napos, 3); 
    mallocate(&apos_occupation, napos);
    mallocate(&apos_DW, napos);
    for(int i=0;i<napos;i++){
        for(int j=0;j<3;j++){
            apos[i][j]=atomdata[j][i];
        }
        apos_occupation[i]=atomdata[3][i];
        apos_DW[i]=atomdata[4][i];
    }
    second_setting_flag=1==space_group_setting?false:true;
    hdf.close();
    return true; 
}

void EBSD_CELL::set_sampling_type()
{
    sampling_type=PG_SAMPLING_TYPE[point_group-1];
    if(-1==sampling_type||14==point_group||26==point_group){
        if(trigonal_flag){
            if(space_group>=143&&space_group<=146) sampling_type=10;
            if(147==space_group||148==space_group) sampling_type=12;
            if(149==space_group||150==space_group||151==space_group) sampling_type=12;
            if(152==space_group||153==space_group||154==space_group||155==space_group) sampling_type=12;
            if(156==space_group||157==space_group||158==space_group) sampling_type=14;
            if(159==space_group||160==space_group||161==space_group) sampling_type=14;
            if(162==space_group||163==space_group) sampling_type=17;
            if(164==space_group||165==space_group||166==space_group||167==space_group) sampling_type=12;
        }else{
            if(14==point_group){
                if(space_group>=115&&space_group<=120){
                    sampling_type=6;
                }else{
                    sampling_type=8;
                }
            }else if(26==point_group){
                if(187==space_group||188==space_group){
                    sampling_type=16;
                }else{
                    sampling_type=17;
                }
            }
        }
    }
}

bool EBSD_CELL::write_parameters_into_hdf5(const char* file_path)
{
    EBSD_HDF5 hdf;
    hdf.open(file_path);
    hdf.write_group("/CrystalData");
    double **atomdata;
    mallocate_2d(&atomdata, napos, 5);
    for(int i=0;i<napos;i++){
        for(int j=0;j<3;j++) atomdata[i][j]=apos[i][j];
        atomdata[i][3]=apos_occupation[i];
        atomdata[i][4]=apos_DW[i];
    }
    hdf.write_array_2d("/CrystalData/AtomData", atomdata, napos, 5);
    hdf.write_array("/CrystalData/Atomtypes", apos_Z, napos);
    char str[]="NULL";
    hdf.write("/CrystalData/CreationDate", str);
    hdf.write("/CrystalData/CreationTime", str);
    hdf.write("/CrystalData/Creator", str);
    hdf.write("/CrystalData/CrystalSystem", crystal_system);
    double lattice[6]={a, b, c, alpha, beta, gamma};
    hdf.write_array("/CrystalData/LatticeParameters", lattice, 6);
    hdf.write("/CrystalData/Natomtypes", napos);
    hdf.write("/CrystalData/ProgramName", str);
    hdf.write("/CrystalData/Source", str);
    hdf.write("/CrystalData/SpaceGroupNumber", space_group);
    int space_group_setting=second_setting_flag?2:1;
    hdf.write("/CrystalData/SpaceGroupSetting", space_group_setting);
    hdf.close();
    return true;
}

void EBSD_CELL::compute_space_matrices()
{
    double calpha=cos(DEG_TO_RAD*alpha);
    double cbeta=cos(DEG_TO_RAD*beta);
    double cgamma=cos(DEG_TO_RAD*gamma);
    double salpha=sin(DEG_TO_RAD*alpha);
    double sbeta=sin(DEG_TO_RAD*beta);
    double sgamma=sin(DEG_TO_RAD*gamma);
    double tgamma=tan(DEG_TO_RAD*gamma);
    double det=a*b*c*a*b*c*(1.0-calpha*calpha-cbeta*cbeta-cgamma*cgamma+2.0*calpha*cbeta*cgamma);
    vol=sqrt(det);
    if(vol<1e-6){
        printf("Error! Suspiciously small volume with six lattice parameters, %.5f, %.5f, %.5f, %.5f, %.5f, and %.5f.", a, b, c, alpha, beta, gamma);
        exit(EXIT_FAILURE);
    }
    dmt[0][0]=a*a; 
    dmt[0][1]=a*b*cgamma; 
    dmt[0][2]=a*c*cbeta;
    dmt[1][0]=dmt[0][1]; 
    dmt[1][1]=b*b;  
    dmt[1][2]=b*c*calpha;
    dmt[2][0]=dmt[0][2]; 
    dmt[2][1]=dmt[1][2];  
    dmt[2][2]=c*c;
    rmt[0][0]=b*c*salpha*b*c*salpha/det; 
    rmt[0][1]=a*b*c*c*(calpha*cbeta-cgamma)/det; 
    rmt[0][2]=a*b*b*c*(cgamma*calpha-cbeta)/det;
    rmt[1][0]=rmt[0][1];              
    rmt[1][1]=a*c*sbeta*a*c*sbeta/det;                   
    rmt[1][2]=a*a*b*c*(cbeta*cgamma-calpha)/det;
    rmt[2][0]=rmt[0][2];              
    rmt[2][1]=rmt[1][2];
    rmt[2][2]=a*b*sgamma*a*b*sgamma/det;
    dsm[0][0]=a;   
    dsm[0][1]=b*cgamma; 
    dsm[0][2]=c*cbeta;
    dsm[1][0]=0.0; 
    dsm[1][1]=b*sgamma; 
    dsm[1][2]=-c*(cbeta*cgamma-calpha)/sgamma;
    dsm[2][0]=0.0; 
    dsm[2][1]=0.0;      
    dsm[2][2]=vol/(a*b*sgamma);
    rsm[0][0]=1.0/a;                                  
    rsm[0][1]=0.0;                                    
    rsm[0][2]=0.0;
    rsm[1][0]=-1.0/(a*tgamma);                        
    rsm[1][1]=1.0/(b*sgamma);                         
    rsm[1][2]=0.0;
    rsm[2][0]=b*c*(cgamma*calpha-cbeta)/(vol*sgamma); 
    rsm[2][1]=a*c*(cbeta*cgamma-calpha)/(vol*sgamma); 
    rsm[2][2]=(a*b*sgamma)/vol;
    if(trigonal_flag){
        compute_trigonal_space_matrix(alpha);
    }
}

void EBSD_CELL::compute_trigonal_space_matrix(double angle)
{
    double x1=0.5/cos(DEG_TO_RAD*0.5*angle);
    double x2=sqrt(1.0-x1*x1);
    double y1=2.0*sin(DEG_TO_RAD*0.5*angle)/sqrt(3.0);
    double y2=sqrt(1.0-y1*y1);
    double matx[3][3]={{1.0, 0.0, 0.0}, {0.0, x1, -x2}, {0.0, x2, x1}};
    double maty[3][3]={{y1, 0.0, -y2}, {0.0, 1.0, 0.0}, {y2, 0.0, y1}};
    double matxt[3][3], matyt[3][3];
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            matxt[i][j]=matx[j][i]; matyt[i][j]=maty[j][i];
        }
    }
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            tridsm[i][j]=0.0;
            for(int k=0;k<3;k++){
                tridsm[i][j]+=matyt[i][k]*matxt[k][j];
            } 
        }
    }
}

void EBSD_CELL::compute_symmetry_matrices()
{
    double symmat[4][4]; char operators[4];
    //symmetry matrices from generator string
    ngenerator=1; 
    compute_symmetry_matrix_by_one_generator(symmat, "aOOO");//the first symmetry matrix
    matrix_copy(symmetry_mats[0], symmat);
    ngenerator+=(int(generator_string[1])-POINT_SYMMETRY_NUMBER);//the subsequent symmetry matrices described by four characters, int get corresponding ASCII value
    int geni;
    for(int i=1;i<ngenerator;i++){
        for(int j=0;j<4;j++){
            geni=2+4*(i-1)+j;
            operators[j]=generator_string[geni];
        }
        compute_symmetry_matrix_by_one_generator(symmat, operators);
        matrix_copy(symmetry_mats[i], symmat);
    }
    if(generator_string[0]=='1'){
        compute_symmetry_matrix_by_one_generator(symmat, "hOOO");//the symmetry matrix corresponding to the inversion operator
        matrix_copy(symmetry_mats[ngenerator++], symmat);
    }
    if(second_setting_flag&&generator_string[++geni]!='0'){//the translation corresponding to the second setting
        operators[0]='a';
        for(int i=1;i<4;i++){
            operators[i]=generator_string[++geni];
        }
        for(int i=1;i<ngenerator;i++){
            double mat1[4][4], mat2[4][4];
            matrix_copy(mat1, symmetry_mats[i]);
            compute_symmetry_matrix_by_one_generator(symmat, operators, -1);//translate to first setting
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    mat2[j][k]=0.0;
                    for(int n=0;n<4;n++){
                        mat2[j][k]+=mat1[j][n]*symmat[n][k];
                    } 
                }
            }
            compute_symmetry_matrix_by_one_generator(symmat, operators);//translate back to second setting
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    mat1[j][k]=0.0;
                    for(int n=0;n<4;n++){
                        mat1[j][k]+=symmat[j][n]*mat2[n][k];
                    } 
                }
            }
            reduce_symmetry_matrix(mat1, 1e-12);
            matrix_copy(symmetry_mats[i], mat1);
        }
    }
    //symmetry matrices from symmetry matrices
    nsymmetry=ngenerator;
    for(int i=0;i<ngenerator;i++){
        compute_symmetry_matrix_by_symmetry_matrix(symmat, symmetry_mats[i], symmetry_mats[i]);
        if(is_symmetry_matrix_new(symmat)){
            matrix_copy(symmetry_mats[nsymmetry++], symmat);
        }
    }
    int i=0;
    while(i<nsymmetry){
        int j=i+1;
        while(j<nsymmetry){
            compute_symmetry_matrix_by_symmetry_matrix(symmat, symmetry_mats[j], symmetry_mats[i]);
            if(is_symmetry_matrix_new(symmat)){
                matrix_copy(symmetry_mats[nsymmetry++], symmat);
                if(nsymmetry>=MAX_MULTIPLICITY) i=j=nsymmetry; 
            }
            j++;
        }
        i++;
    }
    for(int i=0;i<nsymmetry;i++){
        for(int j=0;j<3;j++){
            symmetry_mats[i][j][3]=fmod(symmetry_mats[i][j][3], 1.0);
        }
    }
}

void EBSD_CELL::compute_symmetry_matrix_by_one_generator(double symmat[4][4], const char generator[4], const int multiplier)
{
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            symmat[i][j]=0.0;
        }
    }
    symmat[3][3]=1.0;

    switch(generator[0])//point symmetry operation corresponding to 3x3 submatrix
    {
        case 'a': symmat[0][0]= 1.0; symmat[1][1]= 1.0; symmat[2][2]= 1.0; break; 
        case 'b': symmat[0][0]=-1.0; symmat[1][1]=-1.0; symmat[2][2]= 1.0; break;
        case 'c': symmat[0][0]=-1.0; symmat[1][1]= 1.0; symmat[2][2]=-1.0; break;
        case 'd': symmat[0][2]= 1.0; symmat[1][0]= 1.0; symmat[2][1]= 1.0; break;
        case 'e': symmat[0][1]= 1.0; symmat[1][0]= 1.0; symmat[2][2]=-1.0; break;
        case 'f': symmat[0][1]=-1.0; symmat[1][0]=-1.0; symmat[2][2]=-1.0; break;
        case 'g': symmat[0][1]=-1.0; symmat[1][0]= 1.0; symmat[2][2]= 1.0; break;
        case 'h': symmat[0][0]=-1.0; symmat[1][1]=-1.0; symmat[2][2]=-1.0; break;
        case 'i': symmat[0][0]= 1.0; symmat[1][1]= 1.0; symmat[2][2]=-1.0; break;
        case 'j': symmat[0][0]= 1.0; symmat[1][1]=-1.0; symmat[2][2]= 1.0; break;
        case 'k': symmat[0][1]=-1.0; symmat[1][0]=-1.0; symmat[2][2]= 1.0; break;
        case 'l': symmat[0][1]= 1.0; symmat[1][0]= 1.0; symmat[2][2]= 1.0; break;
        case 'm': symmat[0][1]= 1.0; symmat[1][0]=-1.0; symmat[2][2]=-1.0; break;
        case 'n': symmat[0][1]=-1.0; symmat[1][0]= 1.0; symmat[1][1]=-1.0; symmat[2][2]=1.0; break;
        default: printf("Error! Unrecognized label %c (not between a and n).", generator[0]); exit(EXIT_FAILURE);
    }

    for(int i=1;i<4;i++){//translational component (forward or reverse)
        switch(generator[i])
        {
            case 'A': symmat[i-1][3]=multiplier/6.0; break; 
            case 'B': symmat[i-1][3]=multiplier/4.0; break;
            case 'C': symmat[i-1][3]=multiplier/3.0; break;
            case 'D': symmat[i-1][3]=multiplier/2.0; break;
            case 'E': symmat[i-1][3]=multiplier*2.0/3.0; break;
            case 'F': symmat[i-1][3]=multiplier*3.0/4.0; break;
            case 'G': symmat[i-1][3]=multiplier*5.0/6.0; break;
            case 'O': symmat[i-1][3]=0.0; break;
            case 'X': symmat[i-1][3]=-multiplier*3.0/8.0; break;
            case 'Y': symmat[i-1][3]=-multiplier/4.0; break;
            case 'Z': symmat[i-1][3]=-multiplier/8.0; break;
            default: printf("Error! Unrecognized label %c (not in the 11 uppercase letters).", generator[i]); exit(EXIT_FAILURE);
        }
    }
}

void EBSD_CELL::compute_symmetry_matrix_by_symmetry_matrix(double symmat[4][4], const double symmati[4][4], const double symmatj[4][4])
{
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            symmat[i][j]=0.0;
            for(int k=0;k<4;k++){
                symmat[i][j]+=symmati[i][k]*symmatj[k][j];
            }
        }
    }
    reduce_symmetry_matrix(symmat);
}

void EBSD_CELL::reduce_symmetry_matrix(double symmat[4][4], const double eps)
{
    for(int i=0;i<3;i++){//reduce the translations to the fundamental unit cell
        if(fabs(symmat[i][3])<eps) symmat[i][3]=0.0;
        if(symmat[i][3]<0.0){
            while(symmat[i][3]<0.0){
                symmat[i][3]+=1.0;
            }
        }
        if(symmat[i][3]>=1.0){
            while(symmat[i][3]>=1.0){
                symmat[i][3]-=1.0;
            }
        }
        if(fabs(symmat[i][3]-1.0)<eps) symmat[i][3]=0.0;
    }
}

bool EBSD_CELL::is_symmetry_matrix_new(double symmat[4][4], const double eps)
{
    int k=0, n=0;
    while((k<nsymmetry)&&(n<12)){
        n=0; 
        for(int i=0;i<3;i++){
            for(int j=0;j<4;j++){
                if(fabs(symmat[i][j]-symmetry_mats[k][i][j])<eps) n++;
            }
        }
        k++;
    }
    if(n!=12){
        return true;
    }else{
        return false;
    }
}

void EBSD_CELL::compute_point_symmetry_matrices()
{
    npointsym=0;
    for(int i=0;i<nsymmetry;i++){
        double temp1=symmetry_mats[i][0][3]*symmetry_mats[i][0][3]+symmetry_mats[i][1][3]*symmetry_mats[i][1][3]+symmetry_mats[i][2][3]*symmetry_mats[i][2][3];
        if(temp1<0.1){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    point_dmats[npointsym][j][k]=symmetry_mats[i][j][k];
                }
            }
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    double temp2=0.0;
                    for(int m=0;m<3;m++){
                        for(int n=0;n<3;n++){
                            temp2+=dmt[j][m]*symmetry_mats[i][m][n]*rmt[n][k];
                        }
                    }
                    point_rmats[npointsym][j][k]=temp2;
                }
            }
            npointsym++;
        }
    }
}

void EBSD_CELL::compute_atomic_positions()
{
    mallocate(&apos_multi, napos);
    mallocate_3d(&pos, napos, nsymmetry, 3);
    double equiv[48][3];
    for(int i=0;i<napos;i++){
        int nequiv;
        compute_equivalent_atomic_positions(equiv, nequiv, apos[i]);
        apos_multi[i]=nequiv;
        for(int j=0;j<nequiv;j++){
            for(int k=0;k<3;k++){
                pos[i][j][k]=equiv[j][k];
            }
        }
    }
}

void EBSD_CELL::compute_atomic_density()
{
    double num_M=0.0, num_Z=0.0;
    ave_M=0.0; ave_Z=0.0;
    for(int i=0;i<napos;i++){
        int Z=apos_Z[i];
        double M=ATOM_MASS[Z-1];
        ave_M+=M*apos_multi[i]*apos_occupation[i]; ave_Z+=Z*apos_multi[i];
        num_M+=apos_multi[i]*apos_occupation[i]; num_Z+=apos_multi[i];
    }
    density=ave_M/(vol*1.0e-21*AVOGADRO_CONSTANT);
    ave_M/=num_M; ave_Z/=num_Z;
}

void EBSD_CELL::compute_equivalent_atomic_positions(double equiv[48][3], int &nequiv, const double atom_apos[3], bool reduce_symmetry_flag)
{
    nequiv=0;
    for(int i=0;i<3;i++){
        equiv[nequiv][i]=atom_apos[i];
    }
    nequiv++;
    double atom_pos[3];
    for(int i=1;i<nsymmetry;i++){
        for(int j=0;j<3;j++){
            atom_pos[j]=symmetry_mats[i][j][3];
            for(int k=0;k<3;k++){
                atom_pos[j]+=symmetry_mats[i][j][k]*atom_apos[k];
            }
        }
        for(int j=0;j<3;j++){
            if(fabs(atom_pos[j])<1e-6) atom_pos[j]=0.0;
        }
        if(reduce_symmetry_flag){
            for(int j=0;j<3;j++){
                atom_pos[j]=fmod(atom_pos[j]+100.0, 1.0);
                if(fabs(atom_pos[j])<1e-6) atom_pos[j]=0.0;
            }
        }
        bool is_new=true;
        for(int i=0;i<nequiv;i++){
            double diff=0.0;
            for(int j=0;j<3;j++){
                diff+=fabs(equiv[i][j]-atom_pos[j]);
            }
            if(diff<1e-4) is_new=false;
        }
        if(is_new){
            for(int i=0;i<3;i++){
                equiv[nequiv][i]=atom_pos[i];
            }
            nequiv++;
        }
    }
}

void EBSD_CELL::compute_equivalent_reciprocal_vectors(double equiv[48][3], int &nequiv, double g[3], char space)
{
    double s[3];
    equiv[0][0]=g[0]; equiv[0][1]=g[1]; equiv[0][2]=g[2];
    nequiv=1;
    for(int i=1;i<npointsym;i++){
        for(int j=0;j<3;j++){
            s[j]=0.0;
            for(int k=0;k<3;k++){
                if('r'==space){
                    s[j]+=point_rmats[i][j][k]*g[k];
                }else if('d'==space){
                    s[j]+=point_dmats[i][j][k]*g[k];
                }else{
                    printf("Error! Unrecognized space %c in star computation.", space);
                    exit(EXIT_FAILURE);
                }
            }
        }
        bool is_new=true;
        for(int j=0;j<nequiv;i++){
            double diff=0.0;
            for(int k=0;k<3;k++){
                diff+=fabs(equiv[j][k]-s[k]);
            }
            if(diff<=1e-4){
                is_new=false;
            }
        }
        if(is_new){
            equiv[nequiv][0]=s[0]; equiv[nequiv][1]=s[1]; equiv[nequiv][2]=s[2];
            nequiv++;
        }
    }
}

void EBSD_CELL::apply_point_group_symmetry(int equiv[48][3], int &nequiv, int px, int py, int pz, int nump)
{
    double xy[2]={double(px)/double(nump), double(py)/double(nump)};
    double xyz[3]; int ierr;
    if(hexagonal_flag){
        compute_sphere_from_hexagonal_Lambert(xyz, ierr, xy);
    }else{
        compute_sphere_from_square_Lambert(xyz, ierr, xy);
    }
    if(pz<0) xyz[2]=-xyz[2];
    normalize(xyz, xyz);
    double kstar[3];
    reciprocal(kstar, xyz);
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
        cartesian(xyz, vtmp[i]);
        normalize(xyz, xyz);
        if(xyz[2]<eps){
            equiv[i][2]=-1;
        }else{
            equiv[i][2]=1;
        }
        if(hexagonal_flag){
            compute_hexagonal_Lambert(xy, ierr, xyz);
        }else{
            compute_square_Lambert(xy, ierr, xyz);
        }
        xy[0]*=nump; xy[1]*=nump;
        equiv[i][0]=int(round(xy[0])); equiv[i][1]=int(round(xy[1]));
    }
    nequiv=iequiv;
}

double EBSD_CELL::dot(double hkl1[3], double hkl2[3], char space)
{
    double temp[3], res;
    switch(space)
    {
        case 'r':
            temp[0]=rmt[0][0]*hkl2[0]+rmt[0][1]*hkl2[1]+rmt[0][2]*hkl2[2];
            temp[1]=rmt[1][0]*hkl2[0]+rmt[1][1]*hkl2[1]+rmt[1][2]*hkl2[2];
            temp[2]=rmt[2][0]*hkl2[0]+rmt[2][1]*hkl2[1]+rmt[2][2]*hkl2[2];
            res=hkl1[0]*temp[0]+hkl1[1]*temp[1]+hkl1[2]*temp[2];
            break;
        case 'c':
            res=hkl1[0]*hkl2[0]+hkl1[1]*hkl2[1]+hkl1[2]*hkl2[2];
            break;
        default:
            printf("Error! Unrecognized space %c in dot computation.", space);
            exit(EXIT_FAILURE);
    }
    return res;
}

void EBSD_CELL::cross(double c_hkl[3], double hkl1[3], double hkl2[3], char inspace, char outspace)
{
    switch(inspace)
    {
        case 'c':
            c_hkl[0]=hkl1[1]*hkl2[2]-hkl1[2]*hkl2[1];
            c_hkl[1]=hkl1[2]*hkl2[0]-hkl1[0]*hkl2[2];
            c_hkl[2]=hkl1[0]*hkl2[1]-hkl1[1]*hkl2[0];
            break;
        default:
            printf("Error! Unrecognized space %c in dot computation.", inspace);
            exit(EXIT_FAILURE);
    }

}

double EBSD_CELL::length(double hkl[3], char space)
{
    return sqrt(dot(hkl, hkl, space));
}

double EBSD_CELL::angle(double hkl1[3], double hkl2[3], char space)
{
    double dot12=dot(hkl1, hkl2, space), len1=length(hkl1, space), len2=length(hkl2, space);
    if(0.0==len1||0.0==len2){
        printf("Error! Zero length for vector [%d %d %d] or [%d %d %d] in angle computation.", 
               hkl1[0], hkl1[1], hkl1[2], hkl2[0], hkl2[1], hkl2[2]);
        exit(EXIT_FAILURE);
    }
    double res=dot12/(len1*len2);
    if(res>=1.0) return 0.0;
    else if(res<=-1.0) return PI;
    else return acos(res);
}

void EBSD_CELL::normalize(double n_v[3], double v[3], char space)
{
    double len_v=length(v, space);
    if(0.0!=len_v){
        n_v[0]=v[0]/len_v; n_v[1]=v[1]/len_v; n_v[2]=v[2]/len_v;
    }else{
        n_v[0]=0.0; n_v[1]=0.0; n_v[2]=0.0; 
    }
}

void EBSD_CELL::reciprocal(double r_v[3], double v[3])//from cartesian
{
    double vc[3]={v[0], v[1], v[2]};
    r_v[0]=dsm[0][0]*vc[0]+dsm[1][0]*vc[1]+dsm[2][0]*vc[2];
    r_v[1]=dsm[0][1]*vc[0]+dsm[1][1]*vc[1]+dsm[2][1]*vc[2];
    r_v[2]=dsm[0][2]*vc[0]+dsm[1][2]*vc[1]+dsm[2][2]*vc[2];
}

void EBSD_CELL::cartesian(double c_v[3], double v[3])//from reciprocal
{
    double vc[3]={v[0], v[1], v[2]};
    c_v[0]=rsm[0][0]*vc[0]+rsm[0][1]*vc[1]+rsm[0][2]*vc[2];
    c_v[1]=rsm[1][0]*vc[0]+rsm[1][1]*vc[1]+rsm[1][2]*vc[2];
    c_v[2]=rsm[2][0]*vc[0]+rsm[2][1]*vc[1]+rsm[2][2]*vc[2];
}

double EBSD_CELL::get_interplanar_spacing(double g[3])
{
    double dotgg=dot(g, g);
    if(dotgg<=0.0){
        printf("Error! Zero reciprocal vector in interplanar spacing compution.");
        exit(EXIT_FAILURE);
    }
    return 1.0/sqrt(dotgg);
}

void EBSD_CELL::compute_reflection_range(double dmin)
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

bool EBSD_CELL::is_centering_allowed(double g[3])
{
//Determine whether or not the reciprocal vector is allowed by the lattice centering
    bool is_allowed=true;
    int res;
    switch(centering_symbol)
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
            if(hexagonal_flag){
                res=int(-g[0]+g[1]+g[2]+90)%3;
                if(0!=res) is_allowed=false;
            }
        default:
            printf("Error! Unrecognized centering symbol %c (not P, F, I, A, B, C, or R).", centering_symbol);
            exit(EXIT_FAILURE);
    }
    return is_allowed;
}

double EBSD_CELL::get_excitation_error(double g[3], double k[3], double fn[3])
{
    double kpg[3]={k[0]+g[0], k[1]+g[1], k[2]+g[2]};
    double tkpg[3]={2.0*k[0]+g[0], 2.0*k[1]+g[1], 2.0*k[2]+g[2]};
    double q1=length(kpg), q2=angle(kpg, fn);
    double res=-1.0*dot(g, tkpg)/(2.0*q1*cos(q2));
    return res;
}

void EBSD_CELL::compute_Fourier_coefficient(FOURIER *fouri, double voltage, double g[3], bool shiftflag)
{
    complex<double> FV(0.0, 0.0), FVp(0.0, 0.0);
    double GWK;
    if(0.0==(g[0]*g[0]+g[1]*g[1]+g[2]*g[2])){
        GWK=0.0;
    }else{
        GWK=length(g)*G_TO_WK;
    }
    for(int i=0;i<napos;i++){//loop over atoms in the asymmetric unit
        double UWK=sqrt(apos_DW[i]*DW_TO_WK);//root-mean-square value of the atomic vibration amplitude or thermal displacement in [nm]
        int    Z=apos_Z[i];
        complex<double> Fscat=get_scattering_amplitude(GWK, UWK, Z, voltage)*complex<double>(apos_occupation[i], 0.0);
        complex<double> Fpos(0.0, 0.0);
        for(int j=0;j<apos_multi[i];j++){//loop over atoms in orbit
            double q=-1.0*TWO_PI*(g[0]*pos[i][j][0]+g[1]*pos[i][j][1]+g[2]*pos[i][j][2]);
            complex<double> jFpos(cos(q), sin(q));
            Fpos+=jFpos;
        }
        FV+=Fscat.real()*Fpos; FVp+=Fscat.imag()*Fpos;
    }
    double constV=PRE_CONST_V/vol/FOUR_PI, constU=PRE_CONST_U;
    fouri->Vmod=constV*abs(FV); fouri->Vang=arg(FV);
    fouri->Vpmod=constV*abs(FVp); fouri->Vpang=arg(FVp);
    fouri->Umod=constU*fouri->Vmod; fouri->Upmod=constU*fouri->Vpmod;
    fouri->Vg.real(constV*(FV.real()-FVp.imag())); fouri->Vg.imag(constV*(FV.imag()+FVp.real()));
    fouri->Ug=constU*fouri->Vg;
    if(fabs(fouri->Umod)>0.0){
        fouri->sig=1.0/fabs(fouri->Umod)/lambda;
    }else{
        fouri->sig=1.0e8;
    }
    if(fabs(fouri->Upmod)>0.0){
        fouri->sigp=1.0/fabs(fouri->Upmod)/lambda;
    }else{
        fouri->sigp=1.0e8;
    }
    if(shiftflag){
        fouri->qg.real(cos(fouri->Vang)/fouri->sig-sin(fouri->Vpang)/fouri->sigp); 
        fouri->qg.imag(cos(fouri->Vpang)/fouri->sigp+sin(fouri->Vang)/fouri->sig);
    }else{
        double betag=fouri->Vpang-fouri->Vang;
        fouri->qg.real(1.0/fouri->sig-sin(betag)/fouri->sigp); 
        fouri->qg.imag(cos(betag)/fouri->sigp);
    }
}

void EBSD_CELL::compute_relativistic_wavelength(double voltage)
{
    FOURIER fouri;
    double g[3]={0.0, 0.0, 0.0};
    compute_Fourier_coefficient(&fouri, voltage, g);
    double temp=1000.0*ELECTRON_CHARGE*voltage/ELECTRON_REST_MASS/pow(LIGHT_VELOCITY, 2);
    fouri0.Vmod=fouri.Vmod; fouri0.Vpmod=fouri.Vpmod;
    fouri0.Umod=fouri.Umod; fouri0.Upmod=fouri.Upmod;
    fouri0.Vg=fouri.Vg; fouri0.Ug=fouri.Ug;
    fouri0.gamma=1.0+temp;
    fouri0.voltage=voltage*(1.0+0.5*temp)*1000.0+fouri0.gamma*fouri0.Vmod;
    fouri0.lambda=lambda=1.0e9*PLANK_CONSTANT/sqrt(2.0*ELECTRON_REST_MASS*ELECTRON_CHARGE*fouri0.voltage);
    fouri0.sigma=1.0e-18*TWO_PI*ELECTRON_REST_MASS*ELECTRON_CHARGE*fouri0.gamma*fouri0.lambda/pow(PLANK_CONSTANT, 2);
    compute_Fourier_coefficient(&fouri, voltage, g, true);
    fouri0.sig=fouri.sig; fouri0.sigp=fouri.sigp; fouri0.qg=fouri.qg;
}

void EBSD_CELL::compute_Fourier_coefficients(double voltage, bool initflag)
{
    printf("-->Fourier Coefficients Information<--\n");
    compute_relativistic_wavelength(voltage);
    printf("  Mean inner potential in [V] %.5f\n", fouri0.Vmod);
    printf("  Wavelength corrected for refraction\n");
    printf("  Relativistic correction factor [gamma] in [V] %.5f\n", fouri0.gamma);
    printf("  Relativistic Accelerating Potential [V] %.5f\n", fouri0.voltage);
    printf("  Electron Wavelength [nm] %.5f\n", fouri0.lambda);
    printf("  Interaction constant [V nm^-1] %.5f\n", fouri0.sigma);
    printf("  Normal absorption length [nm] = %.5f\n", fouri0.sigp);

    int imh=2*HKL[0], imk=2*HKL[1], iml=2*HKL[2];
    if(initflag){
        int size1=imh*2+1, size2=imk*2+1, size3=iml*2+1;
        callocate_3d(&LUTUg, size1, size2, size3, 0.0+0.0i);
        callocate_3d(&LUTqg, size1, size2, size3, 0.0+0.0i);
        callocate_3d(&is_double_diffrac, size1, size2, size3, false);
    }
    LUTUg[imh][imk][iml]=fouri0.Ug; LUTqg[imh][imk][iml]=fouri0.qg;

    FOURIER fouri; double g[3];
    printf("Generating Fourier coefficient lookup table ... ");
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                g[0]=double(ih); g[1]=double(ik); g[2]=double(il);
                if(is_centering_allowed(g)){
                    compute_Fourier_coefficient(&fouri, voltage, g, true);
                    int i=imh+ih, j=imk+ik, k=iml+il;
                    LUTUg[i][j][k]=fouri.Ug; LUTqg[i][j][k]=fouri.qg;
                    if(abs(fouri.Ug)<=1e-5){
                        is_double_diffrac[i][j][k]=true;
                    }
                }
            }
        }
    }
    printf("Done!\n");
}

void EBSD_CELL::compute_Bloch_wave_coefficients(int initflag)
{
    int imh=HKL[0]*2, imk=HKL[1]*2, iml=HKL[2]*2; 
    if(initflag){
        callocate_4d(&LUTSgh, napos, imh*2+1, imk*2+1, iml*2+1, 0.0+0.0i);
    }
    printf("Generating Sgh coefficient lookup table ... ");
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                double g[3]={double(ih), double(ik), double(il)};
                if(is_centering_allowed(g)){
                    int ix=ih+imh, iy=ik+imk, iz=il+iml;
                    for(int i=0;i<napos;i++){
                        LUTSgh[i][ix][iy][iz]=complex<double>(0.0, 0.0);
                        double DBWF=pow(apos_Z[i], 2)*apos_occupation[i]*exp(-0.25*apos_DW[i]*dot(g, g));
                        for(int j=0;j<apos_multi[i];j++){
                            double q=TWO_PI*(g[0]*pos[i][j][0]+g[1]*pos[i][j][1]+g[2]*pos[i][j][2]);
                            complex<double> F(cos(q), sin(q));
                            LUTSgh[i][ix][iy][iz]+=F*complex<double>(DBWF, 0.0);
                        }
                    }
                }
            }
        }
    }
    printf("Done!\n");
}

complex<double> EBSD_CELL::get_scattering_amplitude(double G, double U, int Z, double V)
{
    complex<double> F(0.0, 0.0);
    double FR, FI;
    double GAMMA=(V+511.0)/511.0;
    double K0=0.5068*sqrt(1022.0*V+V*V);
    double DWF=get_Debye_Waller_factor(G, U);//Debye-Waller amplitude is included
    const double *A=AA[Z-1], *B=BB[Z-1];
    FR=FOUR_PI*DWF*get_electron_scattering_factor(G, A, B);//Elastic part
    FI=DWF*get_core_excitation_factor(G, Z, V)+get_absorptive_form_factor(G, U, A, B);//Inelastic part
    FR=FR*GAMMA; FI=FI*GAMMA*GAMMA/K0;//acceleration voltage (relativistic scattering) is included
    F.real(FR); F.imag(FI);
    return F;
}

double EBSD_CELL::get_Debye_Waller_factor(double G, double U)
{
    return exp(-0.5*U*U*G*G);
}

double EBSD_CELL::get_electron_scattering_factor(double G, const double A[4], const double B[4])
{
    double F=0.0;
    double S2=G*G/FOUR_PI/FOUR_PI;//S: scattering vector
    for(int i=0;i<4;i++){
        double TEMP=B[i]*S2;
        if(TEMP<0.1){
            F+=A[i]*B[i]*(1.0-0.5*TEMP);
        }else if(TEMP>20){
            F+=A[i]/S2;
        }else{
            F+=A[i]*(1.0-exp(-TEMP))/S2;
        }
    }
    return F;
}

double EBSD_CELL::get_core_excitation_factor(double G, int Z, double V)
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

double EBSD_CELL::get_absorptive_form_factor(double G, double U, const double A[4], const double B[4])
{
    double F=0.0;
    double AF[4], BF[4];
    double FP2=FOUR_PI*FOUR_PI;
    double DWF=exp(-0.5*U*U*G*G);
    for(int i=0;i<4;i++){
        AF[i]=A[i]*FP2; BF[i]=B[i]/FP2;
    }
    for(int j=0;j<4;j++){
        F+=AF[j]*AF[j]*(DWF*I1(BF[j], BF[j], G)-I2(BF[j], BF[j], G, U));
        for(int i=0;i<j;i++){
            F+=2.0*AF[j]*AF[i]*(DWF*I1(BF[i], BF[j], G)-I2(BF[i], BF[j], G, U));
        }
    }
    return F;
}

double EBSD_CELL::EI(double X)
{
    double A1=8.57332, A2=18.05901, A3=8.63476, A4=0.26777, B1=9.57332, B2=25.63295, B3=21.09965, B4=3.95849;
    if(X>60.0){
        printf("Error! the input of the function EI X > 60.");
        exit(1);
    }
    if(X<-60.0) return 0.0;
    double ABSX=abs(X);
    if(X<-1.0){
        return -(A4+ABSX*(A3+ABSX*(A2+ABSX*(A1+ABSX))))/(B4+ABSX*(B3+ABSX*(B2+ABSX*(B1+ABSX))))*exp(-ABSX)/ABSX;
    }else{
        double REI=0.577216+log(ABSX);
        int I=1;
        double SI=X;
        double SUMS=SI;
        while(abs(SI/X)>1.0e-6){
            SI=SI*X*I/((I+1)*(I+1));
            SUMS+=SI;
            I++;
        };
        REI+=SUMS;
        return REI;
    }
}

double EBSD_CELL::IH1(double X1, double X2, double X3)
{
    double RIH1;
    if(X2<=20.0&&X3<=20.0){
        RIH1=exp(-X1)*(EI(X2)-EI(X3));
        return RIH1;
    }
    if(X2>20.0){
        RIH1=exp(X2-X1)*IH2(X2)/X2;
    }else{
        RIH1=exp(-X1)*EI(X2);
    }
    if(X3>20.0){
        RIH1-=exp(X3-X1)*IH2(X3)/X3;
    }else{
        RIH1-=exp(-X1)*EI(X3);
    }
    return RIH1;
}

double EBSD_CELL::IH2(double X)
{
    double INVX=1.0/X;
    int    I=int(200.0*INVX);
    if(I<0||I>DURCH_NUMBER-1){
        printf("Error! Unrecognized index %d in searching durch table.\n", I);
        exit(EXIT_FAILURE);
    }
    double D1=DURCH_TABLE[I], D2=DURCH_TABLE[I+1];
    double RIH2=D1+200.0*(D2-D1)*(INVX-0.5e-3*I);
    return RIH2;
}

double EBSD_CELL::I1(double BI, double BJ, double G)
{
    double RI1;
    double G2=G*G;
    double EPS=fmax(BI, BJ)*G2;
    if(EPS<=0.1){
        double RI1=PI*(BI*log((BI+BJ)/BI)+BJ*log((BI+BJ)/BJ));
        if(fabs(G)<ZERO){
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
    RI1=2.0*PRE_CONST_RI1+log(BIG2)+log(BJG2)-2.0*EI(-BI*BJ*G2/(BI+BJ));
    RI1+=IH1(BIG2, BIG2*BI/(BI+BJ), BIG2);
    RI1+=IH1(BJG2, BJG2*BJ/(BI+BJ), BJG2);
    RI1=RI1*PI/G2;
    return RI1;
}

double EBSD_CELL::I2(double BI, double BJ, double G, double U)
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
        RI2=EI(-HALFU2*G2*BIUH/BIU)+EI(-HALFU2*G2*BJUH/BJU);
        RI2=2.0*(RI2-EI(-BIUH*BJUH*G2/(BIUH+BJUH))-EI(-QUARTERU2*G2));
        RI2+=IH1(HALFU2*G2, QUARTERU2*G2, QUARTERU2*U2*G2/BIU);
        RI2+=IH1(HALFU2*G2, QUARTERU2*G2, QUARTERU2*U2*G2/BJU);
        RI2+=IH1(BIUH*G2, BIUH*BIUH*G2/(BIUH+BJUH), BIUH*BIUH*G2/BIU);
        RI2+=IH1(BJUH*G2, BJUH*BJUH*G2/(BIUH+BJUH), BJUH*BJUH*G2/BJU);
        RI2=RI2*PI/G2;
    }
    return RI2;
}