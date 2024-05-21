#include "MODEL.h"

CELL::CELL(const char *cell_path, const char types[][10], const double DWs[])
{
    QB_tools QB;
    QB_init(&QB);
    QB_read_lmp(&QB, cell_path);
    QSPG_refined(&QB, &QB, SYMPREC);
    QSPGDataset *QD;
    SpacegroupType SG;
    Pointgroup PG;
    QSPG_get_symmetry(&QB, &QD, SYMPREC);
    SG=spgdb_get_spacegroup_type(QD->hall_number);
    PG=ptg_get_pointgroup(SG.pointgroup_number);

    point_group=SG.pointgroup_number;
    strcpy(point_group_symbol, QD->pointgroup_symbol);
    space_group=QD->spacegroup_number;
    strcpy(space_group_symbol, QD->international_symbol);
    hall_number=QD->hall_number;
    matrix_copy(lattice, QD->std_lattice);
    const char CENTERING[]={'0', 'P', 'I', 'F', 'A', 'B', 'C', '0', 'R'};
    int centeringi=SG.centering; centering=CENTERING[centeringi];
    if(QD->choice[0]=='2') is_second_setting=true;
    if(PG.holohedry==TRIGO) is_trigonal=true;
    if(PG.holohedry==HEXA) is_hexagonal=true;
    if(is_second_setting&&PG.holohedry==TRIGO) is_hexagonal=true;//second setting

    nsymmetry=QD->n_operations;
    for(int i=0;i<nsymmetry;i++){
        matrix_copy(sym_rotations[i], QD->rotations[i]);
        vector_copy(sym_translations[i], QD->translations[i]);
    }
    compute_asymmetric_atomic_positions(QD->std_positions, QD->std_types, QD->n_std_atoms);
    callocate(&apos_Z, napos, -1);
    callocate(&apos_M, napos, 0.0);
    callocate(&apos_DW, napos, 0.0);
    callocate(&apos_occupation, napos, 1.0);
    for(int i=0;i<napos;i++){
        int typei=apos_type[i]-1;
        npos+=apos_multi[i];
        apos_DW[i]=DWs[typei];
        for(int j=0;j<TYPE_NUMBER;j++){
            if(0==strcmp(types[typei], ATOM_SYMBOL[j])){
                apos_Z[i]=j+1;
                apos_M[i]=ATOM_MASS[j];
                break;
            }
        }
    }
    compute_atomic_density();
    QSPGDataset_free_dataset(QD);
    QB_free_atom(&QB);
    logging();
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
    double t_mat[3][3], c_v[3];
    matrix_transpose(t_mat, lattice);
    vector_cross(c_v, t_mat[1], t_mat[2]);
    double vol=vector_dot(t_mat[0], c_v);
    double num_M=0.0, num_Z=0.0;
    ave_M=0.0; ave_Z=0.0;
    for(int i=0;i<napos;i++){
        int Z=apos_Z[i];
        double M=apos_M[i];
        ave_M+=M*apos_multi[i]*apos_occupation[i]; ave_Z+=Z*apos_multi[i];
        num_M+=apos_multi[i]*apos_occupation[i]; num_Z+=apos_multi[i];
    }
    density=ave_M/(vol*1.0e-21*AVOGADRO_CONSTANT)*1.0e3;
    ave_M/=num_M; ave_Z/=num_Z;
}

CELL::CELL(const char *hdf5_path)
{
    double **lat, **tra, ***pos; int ***rot;
    char cen[2];
    int setting, trigonal, hexagonal;
    size_t size1, size2, size3;
    HDF5 hdf;
    hdf.open(hdf5_path);
    hdf.read("/CrystalStructure/PointGroupNumber", point_group);
    hdf.read("/CrystalStructure/PointGroupSymbol", point_group_symbol);
    hdf.read("/CrystalStructure/SpaceGroupNumber", space_group);
    hdf.read("/CrystalStructure/SpaceGroupSymbol", space_group_symbol);
    hdf.read("/CrystalStructure/HallNumber", hall_number);
    hdf.read("/CrystalStructure/CenteringVector", cen);
    hdf.read("/CrystalStructure/IsSecondSetting", setting);
    hdf.read("/CrystalStructure/IsTrigonal", trigonal);
    hdf.read("/CrystalStructure/IsHexagonal", hexagonal);
    hdf.read_array_2d("/CrystalStructure/Lattice", &lat, size1, size2);
    hdf.read("/CrystalStructure/SymmetryNumber", nsymmetry);
    hdf.read_array_3d("/CrystalStructure/SymmetryRotations", &rot, size1, size2, size3);
    hdf.read_array_2d("/CrystalStructure/SymmetryTranslations", &tra, size1, size2);
    hdf.read("/CrystalStructure/AsymmetricNumber", napos);
    hdf.read("/CrystalStructure/PositionNumber", npos);
    hdf.read_array("/CrystalStructure/AsymmetricType", &apos_type, size1);
    hdf.read_array("/CrystalStructure/AsymmetricMultiplicity", &apos_multi, size1);
    hdf.read_array_3d("/CrystalStructure/AsymmetricEquivalentPositions", &pos, size1, size2, size3);
    hdf.read_array("/CrystalStructure/AsymmetricAtomicNumber", &apos_Z, size1);
    hdf.read_array("/CrystalStructure/AsymmetricAtomicMass", &apos_M, size1);
    hdf.read_array("/CrystalStructure/AsymmetricDebyeWallerFactor", &apos_DW, size1);
    hdf.read_array("/CrystalStructure/AsymmetricSiteOccupation", &apos_occupation, size1);
    hdf.read("/CrystalStructure/AveragedAtomicMass", ave_M);
    hdf.read("/CrystalStructure/AveragedAtomicNumber", ave_Z);
    hdf.read("/CrystalStructure/Density", density);
    hdf.close();
    mallocate_3d(&apos_pos, napos, nsymmetry, 3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            lattice[i][j]=lat[i][j];
        }
    }
    for(int i=0;i<192;i++){
        for(int j=0;j<3;j++){
            sym_translations[i][j]=tra[i][j];
            for(int k=0;k<3;k++){
                sym_rotations[i][j][k]=rot[j][k][i];
            }
        }
    }
    for(int i=0;i<napos;i++){
        for(int j=0;j<nsymmetry;j++){
            for(int k=0;k<3;k++){
                apos_pos[i][j][k]=pos[j][k][i];
            }
        }
    }
    centering=cen[0];
    is_second_setting=setting==1?true:false;
    is_trigonal=trigonal==1?true:false;
    is_hexagonal=hexagonal==1?true:false;
    logging();
}

void CELL::logging()
{
    printf("-->Crystal Structure Information<--\n");
    printf("Point group        : %d\n", point_group);
    printf("Point group symbol : %s\n", point_group_symbol);
    printf("Space group        : %d\n", space_group);
    printf("Space group symbol : %s\n", space_group_symbol);
    printf("Hall number        : %d\n", hall_number);
    printf("Centering vector   : %c\n", centering);
    printf("Second setting? ");
    if(is_second_setting) printf("Yes\n");
    else printf("No\n");
    printf("Trigonal? ");
    if(is_trigonal) printf("Yes\n");
    else printf("No\n");
    printf("Hexagonal? ");
    if(is_hexagonal) printf("Yes\n");
    else printf("No\n");

    printf("Lattice:\n");
    for(int i=0;i<3;i++){
        printf("> %.5f, %.5f, %.5f\n", lattice[i][0], lattice[i][1], lattice[i][2]);
    }
    printf("Number of symmetry operations: %d\n", nsymmetry);
    printf("Number of asymmetric atomic positions: %d\n", napos);
    for(int i=0;i<napos;i++){
        printf("General atomic position, atomic number, atomic mass, multiplicity, site occupation, Debye-Waller factor: %d, %d, %.5f, %d, %.5f, %.5f\n", i+1, apos_Z[i], apos_M[i], apos_multi[i], apos_occupation[i], apos_DW[i]);
        printf("Equivalent atomic positions (x, y, z):\n");
        for(int j=0;j<apos_multi[i];j++){
            printf("> %.5f, %.5f, %.5f\n", apos_pos[i][j][0], apos_pos[i][j][1], apos_pos[i][j][2]);
        }
    }
    printf("Number of atomic positions: %d\n", npos);
    printf("Density [in g/cm^3], atomic number averaged, atomic mass averaged [g/mol] = %.5f, %.5f, %.5f\n", density, ave_Z, ave_M);
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
}

void CELL::hdf5(const char* hdf5_path)
{
    double **lat; mallocate_2d(&lat, 3, 3);
    int ***rot; mallocate_3d(&rot, 3, 3, 192);
    double **tra; mallocate_2d(&tra, 192, 3);
    double ***pos; mallocate_3d(&pos, nsymmetry, 3, napos);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            lat[i][j]=lattice[i][j];
        }
    }
    for(int i=0;i<192;i++){
        for(int j=0;j<3;j++){
            tra[i][j]=sym_translations[i][j];
            for(int k=0;k<3;k++){
                rot[j][k][i]=sym_rotations[i][j][k];
            }
        }
    }
    for(int i=0;i<napos;i++){
        for(int j=0;j<nsymmetry;j++){
            for(int k=0;k<3;k++){
                pos[j][k][i]=apos_pos[i][j][k];
            }
        }
    }
    char cen[]={centering, '\0'};
    HDF5 hdf;
    hdf.open(hdf5_path);
    hdf.write_group("/CrystalStructure");
    hdf.write("/CrystalStructure/PointGroupNumber", point_group);
    hdf.write("/CrystalStructure/PointGroupSymbol", point_group_symbol);
    hdf.write("/CrystalStructure/SpaceGroupNumber", space_group);
    hdf.write("/CrystalStructure/SpaceGroupSymbol", space_group_symbol);
    hdf.write("/CrystalStructure/HallNumber", hall_number);
    hdf.write("/CrystalStructure/CenteringVector", cen);
    hdf.write("/CrystalStructure/IsSecondSetting", is_second_setting);
    hdf.write("/CrystalStructure/IsTrigonal", is_trigonal);
    hdf.write("/CrystalStructure/IsHexagonal", is_hexagonal);
    hdf.write_array_2d("/CrystalStructure/Lattice", lat, 3, 3);
    hdf.write("/CrystalStructure/SymmetryNumber", nsymmetry);
    hdf.write_array_3d("/CrystalStructure/SymmetryRotations", rot, 3, 3, 192);
    hdf.write_array_2d("/CrystalStructure/SymmetryTranslations", tra, 192, 3);
    hdf.write("/CrystalStructure/AsymmetricNumber", napos);
    hdf.write("/CrystalStructure/PositionNumber", npos);
    hdf.write_array("/CrystalStructure/AsymmetricType", apos_type, napos);
    hdf.write_array("/CrystalStructure/AsymmetricMultiplicity", apos_multi, napos);
    hdf.write_array_3d("/CrystalStructure/AsymmetricEquivalentPositions", pos, nsymmetry, 3, napos);
    hdf.write_array("/CrystalStructure/AsymmetricAtomicNumber", apos_Z, napos);
    hdf.write_array("/CrystalStructure/AsymmetricAtomicMass", apos_M, napos);
    hdf.write_array("/CrystalStructure/AsymmetricDebyeWallerFactor", apos_DW, napos);
    hdf.write_array("/CrystalStructure/AsymmetricSiteOccupation", apos_occupation, napos);
    hdf.write("/CrystalStructure/AveragedAtomicMass", ave_M);
    hdf.write("/CrystalStructure/AveragedAtomicNumber", ave_Z);
    hdf.write("/CrystalStructure/Density", density);
    hdf.close();
}

void CELL::set_sampling_type()
{
    sampling_type=PG_SAMPLING_TYPE[point_group-1];
    if(-1==sampling_type||14==point_group||26==point_group){
        if(is_trigonal){
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

void CELL::compute_lattice_matrices()
{
    double a[3], b[3], c[3];
    double bc[3], ca[3], ab[3];
    matrix_constant(dsm, 0.1, lattice);//*point
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
    matrix_multiply(dmt, t_mat, dsm);
    matrix_transpose(t_mat, rsm);
    matrix_multiply(rmt, t_mat, rsm);
}

void CELL::compute_point_symmetry_matrices()
{
    int count_point=0;
    for(int i=0;i<nsymmetry;i++){
        double temp1=vector_dot(sym_translations[i], sym_translations[i]);
        if(temp1<0.1){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    point_dmats[count_point][j][k]=sym_rotations[i][j][k];
                }
            }
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    double temp2=0.0;
                    for(int m=0;m<3;m++){
                        for(int n=0;n<3;n++){
                            temp2+=dmt[j][m]*sym_rotations[i][m][n]*rmt[n][k];
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

double CELL::dot(double v1[3], double v2[3], char space)
{
    double temp[3], res;
    switch(space)
    {
    case 'r':
        temp[0]=rmt[0][0]*v2[0]+rmt[0][1]*v2[1]+rmt[0][2]*v2[2];
        temp[1]=rmt[1][0]*v2[0]+rmt[1][1]*v2[1]+rmt[1][2]*v2[2];
        temp[2]=rmt[2][0]*v2[0]+rmt[2][1]*v2[1]+rmt[2][2]*v2[2];
        res=v1[0]*temp[0]+v1[1]*temp[1]+v1[2]*temp[2];
        break;
    case 'd':
        temp[0]=dmt[0][0]*v2[0]+dmt[0][1]*v2[1]+dmt[0][2]*v2[2];
        temp[1]=dmt[1][0]*v2[0]+dmt[1][1]*v2[1]+dmt[1][2]*v2[2];
        temp[2]=dmt[2][0]*v2[0]+dmt[2][1]*v2[1]+dmt[2][2]*v2[2];
        res=v1[0]*temp[0]+v1[1]*temp[1]+v1[2]*temp[2];
    default:
        printf("[ERROR] Unrecognized space %c in dot computation.", space);
        exit(EXIT_FAILURE);
    }
    return res;
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
               v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
        exit(EXIT_FAILURE);
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
        n_v[0]=v[0]/len_v; n_v[1]=v[1]/len_v; n_v[2]=v[2]/len_v;
    }else{
        n_v[0]=0.0; n_v[1]=0.0; n_v[2]=0.0; 
    }
}

void CELL::cartesian(double c_v[3], double v[3])
{
    double vc[3]={v[0], v[1], v[2]};
    c_v[0]=rsm[0][0]*vc[0]+rsm[0][1]*vc[1]+rsm[0][2]*vc[2];
    c_v[1]=rsm[1][0]*vc[0]+rsm[1][1]*vc[1]+rsm[1][2]*vc[2];
    c_v[2]=rsm[2][0]*vc[0]+rsm[2][1]*vc[1]+rsm[2][2]*vc[2];
}

void CELL::reciprocal(double r_v[3], double v[3])
{
    double vc[3]={v[0], v[1], v[2]};
    r_v[0]=dsm[0][0]*vc[0]+dsm[1][0]*vc[1]+dsm[2][0]*vc[2];
    r_v[1]=dsm[0][1]*vc[0]+dsm[1][1]*vc[1]+dsm[2][1]*vc[2];
    r_v[2]=dsm[0][2]*vc[0]+dsm[1][2]*vc[1]+dsm[2][2]*vc[2];
}

void CELL::compute_equivalent_reciprocal_vectors(double equiv[48][3], int &nequiv, double g[3], char space)
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
                    printf("[ERROR] Unrecognized space %c in star computation.", space);
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

void CELL::apply_point_group_symmetry(int equiv[48][3], int &nequiv, int px, int py, int pz, int nump)
{
    double xy[2]={double(px)/double(nump), double(py)/double(nump)};
    double xyz[3]; int ierr;
    if(is_hexagonal){
        compute_sphere_from_hexagonal_Lambert(xyz, ierr, xy);
    }else{
        compute_sphere_from_square_Lambert(xyz, ierr, xy);
    }
    if(pz<0) xyz[2]=-xyz[2];
    vector_normalize(xyz, xyz);
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
        vector_normalize(xyz, xyz);
        if(xyz[2]<eps){
            equiv[i][2]=-1;
        }else{
            equiv[i][2]=1;
        }
        if(is_hexagonal){
            compute_hexagonal_Lambert(xy, ierr, xyz);
        }else{
            compute_square_Lambert(xy, ierr, xyz);
        }
        xy[0]*=nump; xy[1]*=nump;
        equiv[i][0]=int(round(xy[0])); equiv[i][1]=int(round(xy[1]));
    }
    nequiv=iequiv;
}

double CELL::get_interplanar_spacing(double g[3])
{
    double dotgg=dot(g, g);
    if(dotgg<=0.0){
        printf("[ERROR] Zero reciprocal vector in interplanar spacing compution.");
        exit(EXIT_FAILURE);
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
    double q1=length(kpg), q2=angle(kpg, fn);
    double res=-1.0*dot(g, tkpg)/(2.0*q1*cos(q2));
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
    default:
        printf("[ERROR] Unrecognized centering symbol %c (not P, F, I, A, B, C, or R).", centering);
        exit(EXIT_FAILURE);
    }
    return is_allowed;
}

void CELL::compute_Bloch_wave_coefficients(bool is_initial)
{
    int imh=HKL[0]*2, imk=HKL[1]*2, iml=HKL[2]*2; 
    if(is_initial){
        callocate_4d(&LUTSgh, napos, imh*2+1, imk*2+1, iml*2+1, complex<double>(0.0, 0.0));
    }
    printf("Generating Bloch wave coefficient lookup table ... ");
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
                            double q=TWO_PI*(g[0]*apos_pos[i][j][0]+g[1]*apos_pos[i][j][1]+g[2]*apos_pos[i][j][2]);
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

void CELL::compute_Fourier_coefficients(double voltage, bool is_initial)
{
    printf("-->Fourier Coefficients Information<--\n");
    update_Fourier_coefficient0(voltage);
    printf("Mean inner potential in [V] %.5f\n", fouri0.Vmod);
    printf("Wavelength corrected for refraction\n");
    printf("Relativistic correction factor [gamma] in [V] %.5f\n", fouri0.gamma);
    printf("Relativistic Accelerating Potential [V] %.5f\n", fouri0.voltage);
    printf("Electron Wavelength [nm] %.5f\n", fouri0.lambda);
    printf("Interaction constant [V nm^-1] %.5f\n", fouri0.sigma);
    printf("Normal absorption length [nm] = %.5f\n", fouri0.sigp);

    int imh=2*HKL[0], imk=2*HKL[1], iml=2*HKL[2];
    if(is_initial){
        int size1=imh*2+1, size2=imk*2+1, size3=iml*2+1;
        callocate_3d(&LUTUg, size1, size2, size3, complex<double>(0.0, 0.0));
        callocate_3d(&LUTqg, size1, size2, size3, complex<double>(0.0, 0.0));
        callocate_3d(&is_double_diffrac, size1, size2, size3, false);
    }
    LUTUg[imh][imk][iml]=fouri0.Ug; LUTqg[imh][imk][iml]=fouri0.qg;

    double g[3];
    printf("Generating Fourier coefficient lookup table ... ");
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                g[0]=double(ih); g[1]=double(ik); g[2]=double(il);
                if(is_centering_allowed(g)){
                    update_Fourier_coefficient(voltage, g, true);
                    int i=imh+ih, j=imk+ik, k=iml+il;
                    LUTUg[i][j][k]=fouri.Ug; LUTqg[i][j][k]=fouri.qg;
                    if(fabs(fouri.Ug)<=1e-5){
                        is_double_diffrac[i][j][k]=true;
                    }
                }
            }
        }
    }
    printf("Done!\n");
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
        GWK=length(g)*G_TO_WK;
    }
    for(int i=0;i<napos;i++){//loop over atoms in the asymmetric unit
        double UWK=sqrt(apos_DW[i]*DW_TO_WK);//root-mean-square value of the atomic vibration amplitude or thermal displacement in [nm]
        int    Z=apos_Z[i];
        complex<double> Fscat=get_scattering_amplitude(GWK, UWK, Z, voltage)*complex<double>(apos_occupation[i], 0.0);
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

complex<double> CELL::get_scattering_amplitude(double G, double U, int Z, double V)
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

double CELL::get_Debye_Waller_factor(double G, double U)
{
    return exp(-0.5*U*U*G*G);
}

double CELL::get_electron_scattering_factor(double G, const double A[4], const double B[4])
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

double CELL::get_core_excitation_factor(double G, int Z, double V)
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

double CELL::get_absorptive_form_factor(double G, double U, const double A[4], const double B[4])
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

double CELL::EI(double X)
{
    double A1=8.57332, A2=18.05901, A3=8.63476, A4=0.26777, B1=9.57332, B2=25.63295, B3=21.09965, B4=3.95849;
    if(X>60.0){
        printf("[ERROR] the input of the function EI X > 60.");
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

double CELL::IH1(double X1, double X2, double X3)
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

double CELL::IH2(double X)
{
    double INVX=1.0/X;
    int    I=int(200.0*INVX);
    if(I<0||I>DURCH_NUMBER-1){
        printf("[ERROR] Unrecognized index %d in searching durch table.\n", I);
        exit(EXIT_FAILURE);
    }
    double D1=DURCH_TABLE[I], D2=DURCH_TABLE[I+1];
    double RIH2=D1+200.0*(D2-D1)*(INVX-0.5e-3*I);
    return RIH2;
}

double CELL::I1(double BI, double BJ, double G)
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

double CELL::I2(double BI, double BJ, double G, double U)
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
    is_periodic[0]=QB.px; is_periodic[1]=QB.py; is_periodic[2]=QB.pz;
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
    QSPG_refined(&QB, &QB, SYMPREC);
    quatify_lattice(QB.mat);
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

void MODEL::quatify_lattice(double mat[3][3])
{
    double lx=mat[0][0], ly=mat[1][1], lz=mat[2][2];
    double xy=mat[1][0], xz=mat[2][0], yz=mat[2][1];
    dimension[0]=lx; dimension[1]=sqrt(xy*xy+ly*ly); dimension[2]=sqrt(xz*xz+yz*yz+lz*lz);

    double a[3]={lx, 0.0, 0.0}, b[3]={xy, ly, 0.0}, c[3]={xz, yz, lz};
    double bc[3], ca[3], ab[3];
    vector_cross(bc, b, c);
    vector_cross(ca, c, a);
    vector_cross(ab, a, b);
    double vol=vector_dot(a, bc);
    vector_constant(bc, 1.0/vol, bc);
    vector_constant(ca, 1.0/vol, ca);
    vector_constant(ab, 1.0/vol, ab);
    dimensionK[0]=vector_length(bc); dimensionK[1]=vector_length(ca); dimensionK[2]=vector_length(ab);
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
    matrix_multiply(dmt, t_mat, dsm);
    matrix_transpose(t_mat, rsm);
    matrix_multiply(rmt, t_mat, rsm);
}

MODEL::~MODEL()
{
    if(ntype>0){
        deallocate(type_index);
    }
    if(natom>0){
        deallocate_2d(atom_pos, natom);
        deallocate(atom_type);
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

void MODEL::reciprocal_to_cartesian(double c_g[3], double r_g[3])
{
    vector_rotate(c_g, rsm, r_g);
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
        double F=get_atomic_scattering_factor(S, A, B, C);
        res+=F*complex<double>(cos(q), sin(q));
    }
    if(is_lorentz_flag){
        double Lp=sqrt((1+cos(2*theta)*cos(2*theta))/(sin(theta)*sin(theta)*cos(theta)));
        res*=Lp;
    }
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