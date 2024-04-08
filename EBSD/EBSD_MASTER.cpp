#include "EBSD_MASTER.h"

complex<double> to_complex(lapack_complex_double c){
    complex<double> res;
    res.real(creal(c));
    res.imag(cimag(c));
    return res;
}

lapack_complex_double to_lapack_complex(complex<double> c){
    lapack_complex_double res=c.real()+c.imag()*I;
    return res;
}

EBSD_KVECTOR::EBSD_KVECTOR(EBSD_CELL *cell, int nump)
{
    double delta=1.0/double(nump);
    double kn=1.0/cell->lambda;
    double xy[2];
    add_k_vector(cell, xy, kn);
    switch(cell->sampling_type)
    {
    case 7:
        for(int i=0;i<=nump;i++){
            for(int j=0;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 8:
        for(int i=0;i<=nump;i++){
            for(int j=-i;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j);
            }
        }
        break;
    case 9:
        for(int i=0;i<=nump;i++){
            for(int j=0;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j);
            }
        }
        break;
    default:
        printf("Error! Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
        exit(EXIT_FAILURE);
    }
    mallocate_2d(&karray, numk, 3); mallocate_2d(&kijarray, numk, 3);
    mallocate(&knarray, numk);
    KNODE *temp=khead;
    for(int ik=0;ik<numk;ik++){
        karray[ik][0]=temp->k[0]; karray[ik][1]=temp->k[1]; karray[ik][2]=temp->k[2];
        knarray[ik]=temp->kn;
        kijarray[ik][0]=temp->i; kijarray[ik][1]=temp->j; kijarray[ik][2]=temp->hemisphere;
        temp=temp->next;
    }
    free_k_node(khead);
}

EBSD_KVECTOR::~EBSD_KVECTOR()
{
    deallocate_2d(karray, numk);
    deallocate(knarray);
    deallocate_2d(kijarray, numk);
}

void EBSD_KVECTOR::add_k_vector(EBSD_CELL *cell, double xy[2], double kn, int i, int j, bool southern_flag)
{
    if(ktail==nullptr){
        khead=ktail=new KNODE;
        double kstar[3]={0.0, 0.0, kn};//c* as the center of the Rosca-Lambert projection
        cell->reciprocal(kstar, kstar);
        khead->k[0]=kstar[0]; khead->k[1]=kstar[1]; khead->k[2]=kstar[2];
        khead->kn=kn;
        khead->i=i; khead->j=j;
        khead->hemisphere=1;//Northern hemisphere
        numk++;
    }else{
        double kstar[3], r_kstar[3];
        int ierr;
        ktail->next=new KNODE;
        ktail=ktail->next;
        if(cell->hexagonal_flag){
            compute_sphere_from_hexagonal_Lambert(kstar, ierr, xy);
        }else{
            compute_sphere_from_square_Lambert(kstar, ierr, xy);
        }
        normalize_vector(kstar, kstar);
        kstar[0]*=kn; kstar[1]*=kn; kstar[2]*=kn;
        cell->reciprocal(r_kstar, kstar);
        ktail->k[0]=r_kstar[0]; ktail->k[1]=r_kstar[1]; ktail->k[2]=r_kstar[2];
        ktail->kn=kn;
        ktail->i=i; ktail->j=j;
        ktail->hemisphere=1;
        numk++;
        if(southern_flag){
            kstar[0]=-kstar[0]; kstar[1]=-kstar[1]; kstar[2]=-kstar[2];
            cell->reciprocal(r_kstar, kstar);
            ktail->next=new KNODE;
            ktail=ktail->next;
            ktail->k[0]=r_kstar[0]; ktail->k[1]=r_kstar[1]; ktail->k[2]=r_kstar[2];
            ktail->kn=kn;
            ktail->i=i; ktail->j=j;
            ktail->hemisphere=-1;
            numk++;
        }
    }
}

void EBSD_KVECTOR::free_k_node(KNODE *khead)
{
    KNODE *cur=khead;
    while(cur!=nullptr){
        KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

EBSD_GVECTOR::EBSD_GVECTOR(EBSD_CELL *cell, BETHE bethe, double k[3], double fn[3], double cutoff)
{
    double g[3]={0.0, 0.0, 0.0};
    int imh=cell->HKL[0], imk=cell->HKL[1], iml=cell->HKL[2];
    int cx=imh*2, cy=imk*2, cz=iml*2;
    complex<double> Ug=cell->LUTUg[cx][cy][cz];
    complex<double> qg=cell->LUTqg[cx][cy][cz];
    double sg=0.0;
    bool   is_double_diffrac=cell->is_double_diffrac[cx][cy][cz];
    add_g_vector(g, Ug, qg, sg, is_double_diffrac);
    double kn=1.0/cell->lambda;
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                if(abs(ih)+abs(ik)+abs(il)!=0){//avoid double counting the origin
                    g[0]=double(ih); g[1]=double(ik); g[2]=double(il);
                    double dhkl=cell->get_interplanar_spacing(g);
                    if(cell->is_centering_allowed(g)&&(dhkl>cutoff)){
                        int ix=ih+cx, iy=ik+cy, iz=il+cz;
                        Ug=cell->LUTUg[ix][iy][iz];
                        qg=cell->LUTqg[ix][iy][iz];
                        sg=cell->get_excitation_error(g, k, fn);
                        is_double_diffrac=cell->is_double_diffrac[ix][iy][iz];
                        if(is_double_diffrac){
                            if(fabs(sg)<=bethe.c_sg){
                                add_g_vector(g, Ug, qg, sg, is_double_diffrac);
                            }
                        }else{
                            double rg=kn*fabs(sg)/fabs(Ug);
                            if(rg<=bethe.c3){
                                add_g_vector(g, Ug, qg, sg, is_double_diffrac);
                            }
                        }
                    }
                }
            }
        }
    }

    GNODE *gtemp=ghead->next;
    int **glist, num=0;
    mallocate_2d(&glist, numg, 3);
    while(gtemp!=nullptr){
        for(int j=0;j<3;j++){
            glist[num][j]=int(gtemp->hkl[j]);
        }
        num++;
        gtemp=gtemp->next;
    }

    GNODE *tailw=nullptr, *tails=nullptr;
    gtemp=ghead->next;
    gtemp->is_strong=true;
    gtemp->is_weak=false;
    tails=gtemp;
    nstrong++;
    for(int i=1;i<num;i++){
        gtemp=gtemp->next;
        double min=MAX_MIN;
        double sgp=kn*fabs(gtemp->sg);
        for(int j=0;j<num;j++){
            double temp;
            int ix=glist[i][0]-glist[j][0]+cx, iy=glist[i][1]-glist[j][1]+cy, iz=glist[i][2]-glist[j][2]+cz;
            if(cell->is_double_diffrac[ix][iy][iz]){
                temp=1.0e4;
            }else{
                temp=sgp/abs(cell->LUTUg[ix][iy][iz]);
            }
            if(temp<min) min=temp;
        }
        if(min>bethe.c2){//ignore the reflection
            gtemp->is_weak=false;
            gtemp->is_strong=false;
            continue;
        }
        if((min>bethe.c1)&&(min<=bethe.c2)){//weak
            if(0==nweak){
                headw=tailw=gtemp;
            }else{
                tailw->nextw=gtemp;
                tailw=gtemp;
            }
            gtemp->is_weak=true;
            gtemp->is_strong=false;
            nweak++;
            continue;
        }
        if(min<=bethe.c1){//strong
            tails->nexts=gtemp;
            tails=gtemp;
            gtemp->is_weak=false;
            gtemp->is_strong=true;
            nstrong++;
        }
    }
    deallocate_2d(glist, numg);
    //printf("%d %d %d %d\n", num, numg, nstrong, nweak);
}

EBSD_GVECTOR::~EBSD_GVECTOR()
{
    GNODE *cur=ghead;
    while(cur!=nullptr){
        GNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void EBSD_GVECTOR::add_g_vector(double hkl[3], complex<double> Ug, complex<double> qg, double sg, bool is_double_diffrac)
{
    if(ghead==nullptr&&gtail==nullptr){
        ghead=gtail=new GNODE;
    }
    gtail->next=new GNODE;
    gtail=gtail->next;
    gtail->number=numg;
    gtail->hkl[0]=hkl[0]; gtail->hkl[1]=hkl[1]; gtail->hkl[2]=hkl[2];
    gtail->Ug=Ug;
    gtail->qg=qg;
    gtail->sg=sg;
    gtail->family_number=0;
    gtail->is_double_diffrac=is_double_diffrac;
    numg++;
}

EBSD_MASTER::EBSD_MASTER(const char *file_path)
{
    size_t len=strlen(file_path);
    bool flag;
    if(len>=4&&strcmp(file_path+len-4, ".nml")==0){
        flag=read_parameters_from_nml(file_path);
    }else{
        printf("Error! Unrecognized file %s.", file_path);
        exit(EXIT_FAILURE);
    }
    if(!flag){
        printf("Error! Unrecognized parameters in file %s.", file_path);
        exit(EXIT_FAILURE);
    }
}

EBSD_MASTER::~EBSD_MASTER()
{
    
}

bool EBSD_MASTER::read_parameters_from_nml(const char *file_path)
{
    FILE *fp=fopen(file_path, "r");
    if(fp==NULL){
        printf("Error! Unable to open file %s.\n", file_path);
    }
    fseek(fp, 0, SEEK_SET);

    char keys[][MAX_LENTH_IN_NAME]={"smallest_spacing", "pixel_number", "Bethe_parameters"}; int keyi=0;
    char readbuff[MAX_LENTH_IN_NAME];
    for(int i=0;(keyi<3)&&(i<MAX_WORD_NUMBER);i++){
        fscanf(fp, "%s", readbuff);
        if('#'==readbuff[0]){
            while(('\n'!=fgetc(fp))&&(!feof(fp)));
            continue;
        }
        if(0==strcmp(readbuff, keys[0])){
            fscanf(fp, "%lf", &dmin);
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, keys[1])){
            fscanf(fp, "%d", &nump);
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, keys[2])){
            fscanf(fp, "%lf", &bethe.c1);
            fscanf(fp, "%lf", &bethe.c2);
            fscanf(fp, "%lf", &bethe.c3);
            keyi++;
            continue;
        }
    }
    if(3==keyi) return true;
    else return false;
}

bool EBSD_MASTER::write_parameters_into_hdf5(const char* file_path, int offset)
{
    double BetheParameters[4]={bethe.c1, bethe.c2, bethe.c3, bethe.c_sg};
    double *EkeVs=Ebins;
    int lastEnergy=1;
    int numEbins=numEbin;
    char xtalname[]="NULL";
    size_t offsets_LP[4]={0, offset, 0, 0};
    size_t offsets_SP[4]={offset, 0, 0};

    EBSD_HDF5 hdf;
    hdf.open(file_path);
    hdf.write_group("/EMData");
    hdf.write_group("/EMData/EBSDmaster");
    hdf.write_array("/EMData/EBSDmaster/BetheParameters", BetheParameters, 4);
    hdf.write_array("/EMData/EBSDmaster/EkeVs", EkeVs, numEbins);
    hdf.write("/EMData/EBSDmaster/lastEnergy", lastEnergy);
    size_t size_LP[4]={nset, 1, width, height}, size_SP[3]={1, width, height};
    hdf.write_hyper_array_4d("/EMData/EBSDmaster/mLPNH", mLPNH, size_LP, offsets_LP, nset, numEbin, width, height);
    hdf.write_hyper_array_4d("/EMData/EBSDmaster/mLPSH", mLPSH, size_LP, offsets_LP, nset, numEbin, width, height);
    hdf.write_hyper_array_3d("/EMData/EBSDmaster/mSPNH", mSPNH, size_SP, offsets_SP, numEbin, width, height);
    hdf.write_hyper_array_3d("/EMData/EBSDmaster/mSPSH", mSPSH, size_SP, offsets_SP, numEbin, width, height);
    hdf.write("/EMData/EBSDmaster/numEbins", numEbins);
    hdf.write("/EMData/EBSDmaster/xtalname", xtalname);
    hdf.close();
    return true;
}

void EBSD_MASTER::compute_master_pattern(const char* file_path)
{
    EBSD_MC mc(file_path);
    EBSD_CELL cell(file_path);
    precompute_master_pattern(&mc, &cell);

    width=height=2*nump+1;
    nset=cell.napos;
    callocate_4d(&mLPNH, cell.napos, 1, width, height, 0.0);
    callocate_4d(&mLPSH, cell.napos, 1, width, height, 0.0);
    callocate_3d(&mSPNH, 1, width, height, 0.0);
    callocate_3d(&mSPSH, 1, width, height, 0.0);
    
    cell.compute_reflection_range(dmin);
    printf("Range of reflections along a*, b*, and c* = %d, %d, and %d.\n", cell.HKL[0], cell.HKL[1], cell.HKL[2]);
    cell.compute_Bloch_wave_coefficients();
    clock_t start, finish, time;
    start=clock();
    int iEstart=numEbin-1;
    for(int iE=iEstart;iE>=0;iE--){
        clock_t istart, ifinish, itime;
        istart=clock();
        printf("Starting computation for energy bin (in reverse order) %d of %d; ", iE+1, numEbin);
        printf("energy [keV] = %.2f.\n", Ebins[iE]);

        //set the accelerating voltage and compute the Fourier coefficients
        double voltage=Ebins[iE];
        if(iE==iEstart){
            cell.compute_Fourier_coefficients(voltage);
        }else{
            cell.compute_Fourier_coefficients(voltage, false);
        }

        //create the incident beam direction list
        EBSD_KVECTOR kvec(&cell, nump);
        printf("Independent beam directions to be considered = %d\n", kvec.numk);

        //create the master reflection list
        int sum_strong=0, sum_weak=0;
        for(int ik=0;ik<kvec.numk;ik++){
            double *kk=kvec.karray[ik];
            double *fn=kk;
            EBSD_GVECTOR gvec(&cell, bethe, kk, fn, dmin);
            sum_strong+=gvec.nstrong; sum_weak+=gvec.nweak;

            complex<double> **dynmat;
            complex<double> **Lgh, ***Sgh;
            callocate_2d(&dynmat, gvec.nstrong, gvec.nstrong, 0.0+0.0i);
            callocate_2d(&Lgh, gvec.nstrong, gvec.nstrong, 0.0+0.0i);
            callocate_3d(&Sgh, cell.napos, gvec.nstrong, gvec.nstrong, 0.0+0.0i);
            //generate the dynamical matrix and solve the dynamical eigenvalue equation
            compute_dynamic_matrix(dynmat, &cell, &gvec);
            compute_Sgh_matrices(Sgh, &cell, &gvec);
            compute_Lgh_matrix(Lgh, dynmat, lambdaE[iE], izmax, depths[iE], mc.depthstep, kvec.knarray[ik], gvec.nstrong);
            double *sval;
            callocate(&sval, cell.napos, 0.0);
            for(int i=0;i<cell.napos;i++){
                complex<double> iSgh[gvec.nstrong][gvec.nstrong];
                for(int j=0;j<gvec.nstrong;j++){
                    for(int k=0;k<gvec.nstrong;k++){
                        iSgh[j][k]=Sgh[i][j][k];
                    }
                }
                for(int j=0;j<gvec.nstrong;j++){
                    for(int k=0;k<gvec.nstrong;k++){
                        complex<double> mtmp=Lgh[j][k]*iSgh[j][k];
                        sval[i]+=mtmp.real();
                    }
                }
                sval[i]/=double(npos);
            }

            int ipx=kvec.kijarray[ik][0], ipy=kvec.kijarray[ik][1], ipz=kvec.kijarray[ik][2];
            int iequiv[48][3], nequiv;
            cell.apply_point_group_symmetry(iequiv, nequiv, ipx, ipy, ipz, nump);
            for(int i=0;i<nequiv;i++){
                int ix=iequiv[i][0]+nump, iy=iequiv[i][1]+nump;
                if(-1==iequiv[i][2]){
                    for(int j=0;j<cell.napos;j++){
                        mLPSH[j][0][ix][iy]=sval[j];
                    }
                }else if(1==iequiv[i][2]){
                    for(int j=0;j<cell.napos;j++){
                        mLPNH[j][0][ix][iy]=sval[j];
                    }
                }
            }
        
            if(0==ik%1000){
                printf("Completed beam direction %d of %d.\n", ik+1, kvec.numk);
            }
            deallocate_2d(dynmat, gvec.nstrong);
            deallocate_2d(Lgh, gvec.nstrong);
            deallocate_3d(Sgh, cell.napos, gvec.nstrong);
            deallocate(sval);
        }

        //convert the hexagonally sampled array to a square Lambert projection
        if(cell.hexagonal_flag){
            double ****auxNH, ****auxSH;
            mallocate_4d(&auxNH, cell.napos, 1, width, height);
            mallocate_4d(&auxSH, cell.napos, 1, width, height);
            for(int i=0;i<cell.napos;i++){
                for(int j=0;j<1;j++){
                    for(int k=0;k<width;k++){
                        for(int n=0;n<height;n++){
                            auxNH[i][j][k][n]=mLPNH[i][j][k][n];
                            auxSH[i][j][k][n]=mLPSH[i][j][k][n];
                        }
                    }
                }
            }
            for(int i=-nump;i<=nump;i++){
                for(int j=-nump;j<=nump;j++){
                    double xy[2]={double(i)/double(nump), double(j)/double(nump)}; 
                    double xyz[3]; int ierr;
                    compute_sphere_from_square_Lambert(xyz, ierr, xy);
                    compute_hexagonal_Lambert(xy, ierr, xyz);
                    if(0!=ierr){
                        printf("Error! Unable to compute hexagonal Lambert interpolation using (%.2f, %.2f, %.2f).\n", xyz[0], xyz[1], xyz[2]);
                        exit(EXIT_FAILURE);
                    }
                    int ix=int(xy[0]), iy=int(xy[1]);
                    int ixp=ix+1, iyp=iy+1;
                    if(ixp>nump) ixp=ix;
                    if(iyp>nump) iyp=iy;
                    double dx=xy[0]-ix, dy=xy[1]-iy;
                    double dxm=1.0-dx, dym=1.0-dy;
                    ixp+=nump; iyp+=nump; ix+=nump; iy+=nump;
                    for(int k=0;k<nset;k++){
                        mLPNH[k][0][i][j]=auxNH[k][0][ix][iy]*dxm*dym+auxNH[k][0][ixp][iy]*dx*dym+
                                          auxNH[k][0][ix][iyp]*dxm*dy+auxNH[k][0][ixp][iyp]*dx*dy;
                        mLPSH[k][0][i][j]=auxSH[k][0][ix][iy]*dxm*dym+auxSH[k][0][ixp][iy]*dx*dym+
                                          auxSH[k][0][ix][iyp]*dxm*dy+auxSH[k][0][ixp][iyp]*dx*dy;
                    }
                }
            }
            deallocate_4d(auxNH, cell.napos, 1, width);
            deallocate_4d(auxSH, cell.napos, 1, width);
        }
        //make sure that the outer pixel rim of the mLPSH patterns is identical to that of the mLPNH array
        int num=nump*2;
        for(int i=0;i<=num;i++){
            for(int j=0;j<cell.napos;j++){
                mLPSH[j][0][i][0]=mLPNH[j][0][i][0];
                mLPSH[j][0][i][num]=mLPNH[j][0][i][num];
                mLPSH[j][0][0][i]=mLPNH[j][0][0][i];
                mLPSH[j][0][num][i]=mLPNH[j][0][num][i];
            }
        }
        for(int i=-nump;i<=nump;i++){
            for(int j=-nump;j<=nump;j++){
                double xy[2]={double(i)/double(nump), double(i)/double(nump)};
                double xyz[3]; int ierr;
                compute_sphere_from_stereographic_projection(xyz, ierr, xy);
                normalize_vector(xyz, xyz);
                int ix=i+nump, iy=j+nump;
                if(0!=ierr){
                    mSPNH[0][ix][iy]=0.0; mSPSH[0][ix][iy]=0.0;
                }else{
                    mSPNH[0][ix][iy]=get_Lambert_interpolation(xyz, mLPNH); 
                    mSPSH[0][ix][iy]=get_Lambert_interpolation(xyz, mLPSH);
                }
            }
        }
        double max_LPNH=0.0, max_LPSH=0.0, min_LPNH=1.0e5, min_LPSH=1.0e5;
        double max_SPNH=0.0, max_SPSH=0.0, min_SPNH=1.0e5, min_SPSH=1.0e5;
        for(int i=-nump;i<nump;i++){
            for(int j=-nump;j<nump;j++){
                int ix=i+nump, iy=j+nump;
                if(max_LPNH<mLPNH[0][0][ix][iy]) max_LPNH=mLPNH[0][0][ix][iy];
                if(min_LPNH>mLPNH[0][0][ix][iy]) min_LPNH=mLPNH[0][0][ix][iy];
                if(max_LPSH<mLPSH[0][0][ix][iy]) max_LPSH=mLPSH[0][0][ix][iy];
                if(min_LPSH>mLPSH[0][0][ix][iy]) min_LPSH=mLPSH[0][0][ix][iy];
            }
        }
        printf("mLPNH max: %.5f min: %.5f, mLPSH max: %.5f min: %.5f\n", max_LPNH, min_LPNH, max_LPSH, min_LPSH);
        printf("-> Average number of strong reflections = %d.\n", int(round(double(sum_strong)/double(kvec.numk))));
        printf("-> Average number of weak reflections = %d.\n", int(round(double(sum_weak)/double(kvec.numk))));
        ifinish=clock(); itime=ifinish-istart;
        printf("Execution time [s]: %.2f.\n", itime);
        write_parameters_into_hdf5(file_path, iE);
        printf("Intermediate data stored in file %s.\n", file_path);
    }
    finish=clock(); time=finish-start;
    printf("Execution time [s]: %.2f.\n", time);
    printf("Intermediate data stored in file %s.\n", file_path);
}

void EBSD_MASTER::precompute_master_pattern(EBSD_MC *mc, EBSD_CELL *cell)
{
    numEbin=mc->numEbin;
    mallocate(&Ebins, mc->numEbin);
    for(int i=0;i<mc->numEbin;i++){
        Ebins[i]=mc->Ehistmin+i*mc->Ebin;
    }
    mallocate(&depths, mc->numEbin);
    izmax=0;
    for(int iE=0;iE<mc->numEbin;iE++){
        for(int ix=0;ix<mc->num_z;ix++){
            for(int iy=0;iy<mc->num_z;iy++){
                int sum_z=0;
                for(int i=0;i<mc->numzbin;i++){
                    sum_z+=mc->accum_z[iE][i][ix][iy];
                }
                int iz=1, isum_z=mc->accum_z[iE][0][ix][iy];
                while(isum_z<0.99*double(sum_z)){
                    isum_z+=mc->accum_z[iE][iz][ix][iy];
                    iz++;
                }
                if(iz>izmax) izmax=iz;
            }
        }
        depths[iE]=double(izmax)*mc->depthstep;
    }
    mallocate_2d(&lambdaE, mc->numEbin, izmax);
    for(int iE=0;iE<mc->numEbin;iE++){
        cell->compute_relativistic_wavelength(Ebins[iE]);
        for(int iz=0;iz<izmax;iz++){
            int sum_z=0;
            for(int ix=0;ix<mc->num_z;ix++){
                for(int iy=0;iy<mc->num_z;iy++){
                    sum_z+=mc->accum_z[iE][iz][ix][iy];
                }
            }
            lambdaE[iE][iz]=double(sum_z)/double(mc->num_e)*exp(TWO_PI*double(iz)*mc->depthstep/cell->fouri0.sigp);
        }
    }
    npos=0;
    for(int i=0;i<cell->napos;i++){
        npos+=cell->apos_multi[i];
    }
}

void EBSD_MASTER::compute_dynamic_matrix(complex<double> **dynmat, EBSD_CELL *cell, EBSD_GVECTOR *gvec)
{
    GNODE *rtemp=gvec->ghead->next;
    int imh=2*cell->HKL[0], imk=2*cell->HKL[1], iml=2*cell->HKL[2];
    int ir=0;
    while(rtemp!=nullptr){
        GNODE *ctemp=gvec->ghead->next;
        int ic=0;
        while(ctemp!=nullptr){
            int ih, ik, il;
            if(ic!=ir){
                if(0!=gvec->nweak){
                    complex<double> wsum(0.0, 0.0);
                    GNODE *wtemp=gvec->headw;
                    while(wtemp!=nullptr){
                        complex<double> Ughp, Uhph;
                        complex<double> temp(1.0/wtemp->sg, 0.0);
                        ih=rtemp->hkl[0]-wtemp->hkl[0]+imh;
                        ik=rtemp->hkl[1]-wtemp->hkl[1]+imk;
                        il=rtemp->hkl[2]-wtemp->hkl[2]+iml;
                        Ughp=cell->LUTUg[ih][ik][il];
                        ih=wtemp->hkl[0]-ctemp->hkl[0]+imh;
                        ik=wtemp->hkl[1]-ctemp->hkl[1]+imk;
                        il=wtemp->hkl[2]-ctemp->hkl[2]+iml;
                        Uhph=cell->LUTUg[ih][ik][il];
                        wsum+=Ughp*Uhph*temp;
                        wtemp=wtemp->nextw;
                    }
                    ih=rtemp->hkl[0]-ctemp->hkl[0]+imh;
                    ik=rtemp->hkl[1]-ctemp->hkl[1]+imk;
                    il=rtemp->hkl[2]-ctemp->hkl[2]+iml;
                    complex<double> temp(0.5*cell->lambda, 0.0);
                    dynmat[ir][ic]=cell->LUTUg[ih][ik][il]-temp*wsum;
                }else{
                    ih=rtemp->hkl[0]-ctemp->hkl[0]+imh;
                    ik=rtemp->hkl[1]-ctemp->hkl[1]+imk;
                    il=rtemp->hkl[2]-ctemp->hkl[2]+iml;
                    dynmat[ir][ic]=cell->LUTUg[ih][ik][il];
                }
            }else{
                if(0!=gvec->nweak){
                    double wsgsum=0.0;
                    GNODE *wtemp=gvec->headw;
                    while(wtemp!=nullptr){
                        complex<double> Ughp;
                        ih=rtemp->hkl[0]-wtemp->hkl[0]+imh;
                        ik=rtemp->hkl[1]-wtemp->hkl[1]+imk;
                        il=rtemp->hkl[2]-wtemp->hkl[2]+iml;
                        Ughp=cell->LUTUg[ih][ik][il];
                        wsgsum+=fabs(Ughp)*fabs(Ughp)/wtemp->sg;
                        wtemp=wtemp->nextw;
                    }
                    wsgsum*=cell->fouri0.lambda/2.0;
                    dynmat[ir][ir].real(2.0*rtemp->sg/cell->fouri0.lambda-wsgsum); dynmat[ir][ir].imag(cell->fouri0.Upmod);
                }else{
                    dynmat[ir][ir].real(2.0*rtemp->sg/cell->fouri0.lambda); dynmat[ir][ir].imag(cell->fouri0.Upmod);
                }
            }
            ctemp=ctemp->nexts;
            ic++;
        }
        rtemp=rtemp->nexts;
        ir++;
    }
}

void EBSD_MASTER::compute_Sgh_matrices(complex<double>*** Sgh, EBSD_CELL *cell, EBSD_GVECTOR* gvec)
{
    GNODE *rtemp=gvec->ghead->next;
    int ir=0;
    int imh=2*cell->HKL[0], imk=2*cell->HKL[1], iml=2*cell->HKL[2];
    while(rtemp!=nullptr){
        GNODE *ctemp=gvec->ghead->next;
        int ic=0;
        while(ctemp!=nullptr){
            int ih, ik, il;
            ih=ctemp->hkl[0]-rtemp->hkl[0]+imh;
            ik=ctemp->hkl[1]-rtemp->hkl[1]+imk;
            il=ctemp->hkl[2]-rtemp->hkl[2]+iml;
            for(int i=0;i<nset;i++){
                Sgh[i][ir][ic]=cell->LUTSgh[i][ih][ik][il];
            }
            ic++;
            ctemp=ctemp->nexts;
        }
        ir++;
        rtemp=rtemp->nexts;
    }
}

void EBSD_MASTER::compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double *EWF, int ZMAX, double THICK, double DEPTHSTEP, double KN, int NS)
{
    int INFO;
    char JOBVL='N', JOBVR='V';
    int LDA=NS, LWORK=-1, LWMAX=5000; 
    lapack_complex_double *W; mallocate(&W, NS); 
    lapack_complex_double *VL, *CG; mallocate(&VL, NS*NS); mallocate(&CG, NS*NS);
    lapack_complex_double *WORK; mallocate(&WORK, LWMAX); 
    double *RWORK; mallocate(&RWORK, 2*NS);
    lapack_complex_double *A; mallocate(&A, LDA*NS);
    for(int i=0;i<LDA;i++){
        for(int j=0;j<NS;j++){
            A[i*LDA+j]=to_lapack_complex(DMAT[i][j]);
        }
    }
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//determine the optimal LWORK, i.e., WORK[0]
    LWORK=min(LWMAX, int(creal(WORK[0]))); LDA=NS;
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//call the eigenvalue solver
    if(0!=INFO){
        printf("Error! Unable to return INFO as zero when calling the eigenvalue solver.");
        exit(EXIT_FAILURE);
    }
    deallocate(VL);
    deallocate(WORK); deallocate(RWORK);

    lapack_complex_double *CGINV; mallocate(&CGINV, NS*NS);
    for(int i=0;i<NS*NS;i++) CGINV[i]=CG[i];
    int *IPIV; mallocate(&IPIV, NS); 
    LDA=NS; 
    INFO=LAPACKE_zgetrf(LAPACK_ROW_MAJOR, LDA, NS, CGINV, LDA, IPIV);
    int MILWORK=-1;
    lapack_complex_double *GETMILWORK; mallocate(&GETMILWORK, LWMAX); 
    LDA=NS; 
    INFO=LAPACKE_zgetri_work(LAPACK_ROW_MAJOR, NS, CGINV, LDA, IPIV, GETMILWORK, MILWORK);
    MILWORK=int(creal(GETMILWORK[0]));
    lapack_complex_double *MIWORK; callocate(&MIWORK, MILWORK, 0.0+0.0*I);
    LDA=NS;
    INFO=LAPACKE_zgetri_work(LAPACK_ROW_MAJOR, NS, CGINV, LDA, IPIV, MIWORK, MILWORK);
    deallocate(IPIV); deallocate(MIWORK);

    complex<double> *NW; callocate(&NW, NS, 0.0+0.0i);
    complex<double> F(2.0*KN, 0.0);
    for(int i=0;i<NS;i++){
        NW[i]=to_complex(W[i]);
        NW[i]=NW[i]/F;
    }
    complex<double> **NCGINV; callocate_2d(&NCGINV, NS, NS, 0.0+0.0i);
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCGINV[i][j]=to_complex(CGINV[i*NS+j]);
        }
    }
    complex<double> **NCG; callocate_2d(&NCG, NS, NS, 0.0+0.0i);
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCG[i][j]=to_complex(CG[i*NS+j]);
        }
    }
    complex<double> **IJK; callocate_2d(&IJK, NS, NS, 0.0+0.0i);
    double TPI=TWO_PI*DEPTHSTEP, DZT=DEPTHSTEP/THICK;
    for(int j=0;j<NS;j++){
        for(int k=0;k<NS;k++){
            complex<double> sumq(0.0, 0.0);
            complex<double> q(TPI*(NW[j].imag()+NW[k].imag()), TPI*(NW[j].real()-NW[k].real()));
            if(q.real()<0.0) q=-q;
            for(int iz=0;iz<ZMAX;iz++){
                sumq+=EWF[iz]*exp(-double(iz)*q);
            }
            IJK[j][k]=conj(NCGINV[j][0])*sumq*NCGINV[k][0];
            IJK[j][k]*=DZT;
        }
    }
    complex<double> **NCGCONJ; callocate_2d(&NCGCONJ, NS, NS, 0.0+0.0i);
    complex<double> **NCGT; callocate_2d(&NCGT, NS, NS, 0.0+0.0i);
    complex<double> **TEMP; callocate_2d(&TEMP, NS, NS, 0.0+0.0i);
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCGCONJ[i][j]=conj(NCG[i][j]);
        }
    }
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCGT[i][j]=NCG[j][i];
        }
    }
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            for(int k=0;k<NS;k++){
                TEMP[i][j]+=NCGCONJ[i][k]*IJK[k][j];
            } 
        }
    }
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            Lgh[i][j]=complex<double>(0.0, 0.0);
            for(int k=0;k<NS;k++){
                Lgh[i][j]+=TEMP[i][k]*NCGT[k][j];
            } 
        }
    }
    deallocate_2d(NCGCONJ, NS); deallocate_2d(NCGT, NS); deallocate_2d(TEMP, NS);
}

double EBSD_MASTER::get_Lambert_interpolation(double xyz[3], double ****mat, bool hexagonal_flag)
{
    int ix, iy, ixp, iyp;
    double dx, dy, dxm, dym;
    compute_Lambert_interpolation(xyz, nump, hexagonal_flag, ix, iy, ixp, iyp, dx, dy, dxm, dym);
    double res=0.0;
    for(int i=0;i<nset;i++){
        res+=mat[i][0][ix][iy]*dxm*dym+mat[i][0][ixp][iy]*dx*dym+
             mat[i][0][ix][iyp]*dxm*dy+mat[i][0][ixp][iyp]*dx*dy;
    }
    return res;
}

