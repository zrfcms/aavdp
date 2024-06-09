#include "DKD.h"

bool is_in_hexagonal_grid(double xy[2])
{
    double x=fabs(xy[0]-0.5*xy[1]);
    double y=fabs(SQRT_HALF_3*xy[1]);
    if(x>1.0||y>SQRT_HALF_3) return false;
    if(x+y*SQRT_3_INVERSE>1.0) return false;
    return true;
}

DKD_KVECTOR::DKD_KVECTOR(CELL *cell, int npx, int npy)
{
    double delta=1.0/double(npx);
    double kn=1.0/cell->fouri0.lambda;
    double xy[2];
    add_k_vector(cell, xy, kn);
    switch(cell->sampling_type)
    {
    case 1://triclinic 1
        for(int j=-npy;j<=npy;j++){
            for(int i=-npx;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 2://triclinic -1
        for(int j=-npy;j<=npy;j++){
            for(int i=-npx;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j);
            }
        }
        break;
    case 3://monoclinic 2
        for(int j=-npy;j<=npy;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 4://monoclinic m
        for(int j=0;j<=npy;j++){
            for(int i=-npx;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 5://monoclinic 2/m, orthorhombic 222, mm2, tetragonal 4, -4
        for(int j=0;j<=npy;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 6://orthorhombic mmm, tetragonal 4/m, 422, -4m2, cubic m-3, 432
        for(int j=0;j<=npy;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j);
            }
        }
        break;
    case 7://tetragonal 4mm
        for(int i=0;i<=npx;i++){
            for(int j=0;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 8://tetragonal -42m, cubic -43m
        for(int i=0;i<=npx;i++){
            for(int j=-i;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j);
            }
        }
        break;
    case 9://tetragonal 4/mmm, cubic m-3m
        for(int i=0;i<=npx;i++){
            for(int j=0;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_vector(cell, xy, kn, i, j);
            }
        }
        break;
    case 10://hexagonal 3
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j, true);
                }
            }
        }
        break;
    case 11://rhombohedral 3
        printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
        exit(EXIT_FAILURE);
    case 12://hexagonal -3, 321, -6 [not implemented: rhombohedral 32]
        if(cell->is_trigonal&&cell->is_second_setting){
            printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
            exit(EXIT_FAILURE);
        }
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j);
                }
            }
        }
    case 13://hexagonal 312 [not implemented: rhombohedral -3]
        if(cell->is_trigonal&&cell->is_second_setting){
            printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
            exit(EXIT_FAILURE);
        }
        for(int j=0;j>=-npx;j--){
            for(int i=j/2;i<npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j);
                }
            }
        }
        for(int i=0;i<=npx;i++){
            for(int j=0;j<i/2;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j);
                }
            }
        }
        break;
    case 14://hexagonal 3m1 [not implemented: rhombohedral 3m]
        if(cell->is_trigonal&&cell->is_second_setting){
            printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
            exit(EXIT_FAILURE);
        }
        for(int j=1;j<=npx;j++){
            for(int i=j;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j, true);
                }
            }
        }
        break;
    case 15://hexagonal 31m, 6
        for(int j=1;j<=npx;j++){
            for(int i=j;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j, true);
                }
            }
        }
        break;
    case 16://hexagonal -3m1, 622, -6m2 [not implemented: rhombohedral -3m]
        if(cell->is_trigonal&&cell->is_second_setting){
            printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
            exit(EXIT_FAILURE);
        }
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)<(PI/6.0-1.0e-4)))){
                    continue;
                }else if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j);
                }
            }
        }
        break;
    case 17://hexagonal -31m, 6/m, -62m
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)>(PI/3.0+1.0e-4)))){
                    continue;
                }else if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j);
                }
            }
        }
        break;
    case 18://hexagonal 6mm
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)>(PI/6.0+1.0e-4)))){
                    continue;
                }else if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j, true);
                }
            }
        }
        break;
    case 19:
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)>(PI/6.0+1.0e-4)))){
                    continue;
                }else if(is_in_hexagonal_grid(xy)){
                    add_k_vector(cell, xy, kn, i, j);
                }
            }
        }
        break;
    default:
        printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
        exit(EXIT_FAILURE);
    }
    mallocate_2d(&karray, numk, 3); mallocate_2d(&kijarray, numk, 3);
    mallocate(&knarray, numk);
    DKD_KNODE *temp=khead;
    for(int ik=0;ik<numk;ik++){
        karray[ik][0]=temp->k[0]; karray[ik][1]=temp->k[1]; karray[ik][2]=temp->k[2];
        knarray[ik]=temp->kn;
        kijarray[ik][0]=temp->i; kijarray[ik][1]=temp->j; kijarray[ik][2]=temp->hemisphere;
        temp=temp->next;
    }
    free_k_node(khead);
}

DKD_KVECTOR::~DKD_KVECTOR()
{
    deallocate_2d(karray, numk);
    deallocate(knarray);
    deallocate_2d(kijarray, numk);
}

void DKD_KVECTOR::add_k_vector(CELL *cell, double xy[2], double kn, int i, int j, bool southern_flag)
{
    if(ktail==nullptr){
        khead=ktail=new DKD_KNODE;
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
        ktail->next=new DKD_KNODE;
        ktail=ktail->next;
        if(cell->is_hexagonal){
            compute_sphere_from_hexagonal_Lambert(kstar, ierr, xy);
        }else{
            compute_sphere_from_square_Lambert(kstar, ierr, xy);
        }
        vector_normalize(kstar, kstar);
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
            ktail->next=new DKD_KNODE;
            ktail=ktail->next;
            ktail->k[0]=r_kstar[0]; ktail->k[1]=r_kstar[1]; ktail->k[2]=r_kstar[2];
            ktail->kn=kn;
            ktail->i=i; ktail->j=j;
            ktail->hemisphere=-1;
            numk++;
        }
    }
}

void DKD_KVECTOR::free_k_node(DKD_KNODE *khead)
{
    DKD_KNODE *cur=khead;
    while(cur!=nullptr){
        DKD_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

DKD_GVECTOR::DKD_GVECTOR(CELL *cell, BETHE bethe, double k[3], double fn[3], double cutoff)
{
    int imh=cell->HKL[0], imk=cell->HKL[1], iml=cell->HKL[2];
    int cx=imh*2, cy=imk*2, cz=iml*2;
    double g[3]={0.0, 0.0, 0.0};
    complex<double> Ug=cell->LUTUg[cx][cy][cz];
    complex<double> qg=cell->LUTqg[cx][cy][cz];
    double sg=0.0;
    bool   is_double_diffrac=cell->is_double_diffrac[cx][cy][cz];
    add_g_vector(g, Ug, qg, sg, is_double_diffrac);
    double kn=1.0/cell->fouri0.lambda;
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

    DKD_GNODE *gtemp=ghead->next;
    int **glist, num=0;
    mallocate_2d(&glist, numg, 3);
    while(gtemp!=nullptr){
        for(int j=0;j<3;j++){
            glist[num][j]=int(gtemp->hkl[j]);
        }
        num++;
        gtemp=gtemp->next;
    }

    DKD_GNODE *tailw=nullptr, *tails=nullptr;
    gtemp=ghead->next;
    gtemp->is_strong=true;
    gtemp->is_weak=false;
    tails=gtemp;
    nstrong++;
    for(int i=1;i<num;i++){
        gtemp=gtemp->next;
        double imin=1.0e8;
        double sgp=kn*fabs(gtemp->sg);
        for(int j=0;j<num;j++){
            double temp;
            int ix=glist[i][0]-glist[j][0]+cx, iy=glist[i][1]-glist[j][1]+cy, iz=glist[i][2]-glist[j][2]+cz;
            if(cell->is_double_diffrac[ix][iy][iz]){
                temp=1.0e4;
            }else{
                temp=sgp/fabs(cell->LUTUg[ix][iy][iz]);
            }
            if(temp<imin) imin=temp;
        }
        if(imin>bethe.c2){//ignore the reflection
            gtemp->is_weak=false;
            gtemp->is_strong=false;
            continue;
        }
        if((imin>bethe.c1)&&(imin<=bethe.c2)){//weak
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
        if(imin<=bethe.c1){//strong
            tails->nexts=gtemp;
            tails=gtemp;
            gtemp->is_weak=false;
            gtemp->is_strong=true;
            nstrong++;
        }
    }
    deallocate_2d(glist, numg);
}

DKD_GVECTOR::~DKD_GVECTOR()
{
    DKD_GNODE *cur=ghead;
    while(cur!=nullptr){
        DKD_GNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void DKD_GVECTOR::add_g_vector(double hkl[3], complex<double> Ug, complex<double> qg, double sg, bool is_double_diffrac)
{
    if(ghead==nullptr&&gtail==nullptr){
        ghead=gtail=new DKD_GNODE;
    }
    gtail->next=new DKD_GNODE;
    gtail=gtail->next;
    gtail->hkl[0]=hkl[0]; gtail->hkl[1]=hkl[1]; gtail->hkl[2]=hkl[2];
    gtail->Ug=Ug;
    gtail->qg=qg;
    gtail->sg=sg;
    gtail->is_double_diffrac=is_double_diffrac;
    numg++;
}

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

DKD::DKD(const char *hdf5_path, double dmin, double c1, double c2, double c3, double c_sg)
{
    DKD_MC mc(hdf5_path);
    nump=mc.nump; numEbin=mc.numEbin;

    CELL cell(hdf5_path);
    int npos=cell.npos; napos=cell.napos;
    cell.compute_reflection_range(dmin);
    cell.compute_Bloch_wave_coefficients();
    printf("Range of reflections along a*, b*, and c* = %d, %d, and %d.\n", cell.HKL[0], cell.HKL[1], cell.HKL[2]);

    clock_t start, finish;
    start=clock();
    int iEstart=numEbin-1, imp=nump/2;
    callocate_4d(&mLPNH, cell.napos, 1, nump, nump, 0.0);
    callocate_4d(&mLPSH, cell.napos, 1, nump, nump, 0.0);
    callocate_3d(&mSPNH, 1, nump, nump, 0.0);
    callocate_3d(&mSPSH, 1, nump, nump, 0.0);
    for(int iE=iEstart;iE>=0;iE--){
        clock_t istart, ifinish;
        istart=clock();
        printf("Starting computation for energy bin (in reverse order) %d of %d; ", iE+1, numEbin);
        //set the accelerating voltage and compute the Fourier coefficients
        double voltage=mc.Ebins[iE];
        printf("energy [keV] = %.2f.\n", voltage);
        if(iE==iEstart){
            cell.compute_Fourier_coefficients(voltage);
        }else{
            cell.compute_Fourier_coefficients(voltage, false);
        }

        double *lambdaE; mallocate(&lambdaE, mc.izmax);
        for(int iz=0;iz<mc.izmax;iz++){
            int sum_z=0;
            for(int ix=0;ix<mc.numpz;ix++){
                for(int iy=0;iy<mc.numpz;iy++){
                    sum_z+=mc.accum_z[iE][iz][ix][iy];
                }
            }
            lambdaE[iz]=double(sum_z)/double(mc.num_e)*exp(TWO_PI*double(iz)*mc.depthstep/cell.fouri0.sigp);
        }

        //create the incident beam direction list
        DKD_KVECTOR kvec(&cell, imp, imp);
        printf("Independent beam directions to be considered = %d\n", kvec.numk);

        //create the master reflection list
        int sum_strong=0, sum_weak=0;
        for(int ik=0;ik<kvec.numk;ik++){
            double *kk=kvec.karray[ik];
            double *fn=kk;
            DKD_GVECTOR gvec(&cell, bethe, kk, fn, dmin);
            sum_strong+=gvec.nstrong; sum_weak+=gvec.nweak;

            complex<double> **dynmat;
            complex<double> **Lgh, ***Sgh;
            callocate_2d(&dynmat, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
            callocate_2d(&Lgh, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
            callocate_3d(&Sgh, cell.napos, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
            //generate the dynamical matrix and solve the dynamical eigenvalue equation
            compute_dynamic_matrix(dynmat, &cell, &gvec);
            compute_Sgh_matrices(Sgh, &cell, &gvec);
            compute_Lgh_matrix(Lgh, dynmat, lambdaE, mc.izmax, mc.depths[iE], mc.depthstep, kvec.knarray[ik], gvec.nstrong);
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
            cell.apply_point_group_symmetry(iequiv, nequiv, ipx, ipy, ipz, imp);
            for(int i=0;i<nequiv;i++){
                int ix=iequiv[i][0]+imp, iy=iequiv[i][1]+imp;
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
        if(cell.is_hexagonal){
            double ****auxNH, ****auxSH;
            mallocate_4d(&auxNH, cell.napos, 1, nump, nump);
            mallocate_4d(&auxSH, cell.napos, 1, nump, nump);
            for(int i=0;i<cell.napos;i++){
                for(int j=0;j<1;j++){
                    for(int k=0;k<nump;k++){
                        for(int n=0;n<nump;n++){
                            auxNH[i][j][k][n]=mLPNH[i][j][k][n];
                            auxSH[i][j][k][n]=mLPSH[i][j][k][n];
                        }
                    }
                }
            }
            for(int i=-imp;i<=imp;i++){
                for(int j=-imp;j<=imp;j++){
                    double xy[2]={double(i)/double(imp), double(j)/double(imp)}; 
                    double xyz[3]; int ierr;
                    compute_sphere_from_square_Lambert(xyz, ierr, xy);
                    compute_hexagonal_Lambert(xy, ierr, xyz);
                    if(0!=ierr){
                        printf("[ERROR] Unable to compute hexagonal Lambert interpolation using (%.2f, %.2f, %.2f).\n", xyz[0], xyz[1], xyz[2]);
                        exit(EXIT_FAILURE);
                    }
                    xy[0]*=double(imp); xy[1]*=double(imp);
                    int ix=int(xy[0]), iy=int(xy[1]);
                    int ixp=ix+1, iyp=iy+1;
                    if(ixp>imp) ixp=ix;
                    if(iyp>imp) iyp=iy;
                    double dx=xy[0]-ix, dy=xy[1]-iy;
                    double dxm=1.0-dx, dym=1.0-dy;
                    ixp+=imp; iyp+=imp; ix+=imp; iy+=imp;
                    int idx=i+imp, idy=j+imp;
                    for(int k=0;k<napos;k++){
                        mLPNH[k][0][idx][idy]=auxNH[k][0][ix][iy]*dxm*dym+auxNH[k][0][ixp][iy]*dx*dym+
                                          auxNH[k][0][ix][iyp]*dxm*dy+auxNH[k][0][ixp][iyp]*dx*dy;
                        mLPSH[k][0][idx][idy]=auxSH[k][0][ix][iy]*dxm*dym+auxSH[k][0][ixp][iy]*dx*dym+
                                          auxSH[k][0][ix][iyp]*dxm*dy+auxSH[k][0][ixp][iyp]*dx*dy;
                    }
                }
            }
            deallocate_4d(auxNH, cell.napos, 1, nump);
            deallocate_4d(auxSH, cell.napos, 1, nump);
        }
        //make sure that the outer pixel rim of the mLPSH patterns is identical to that of the mLPNH array
        int num=imp*2;
        for(int i=0;i<=num;i++){
            for(int j=0;j<cell.napos;j++){
                mLPSH[j][0][i][0]=mLPNH[j][0][i][0];
                mLPSH[j][0][i][num]=mLPNH[j][0][i][num];
                mLPSH[j][0][0][i]=mLPNH[j][0][0][i];
                mLPSH[j][0][num][i]=mLPNH[j][0][num][i];
            }
        }
        for(int i=-imp;i<=imp;i++){
            for(int j=-imp;j<=imp;j++){
                double xy[2]={double(i)/double(imp), double(i)/double(imp)};
                double xyz[3]; int ierr;
                compute_sphere_from_stereographic_projection(xyz, ierr, xy);
                vector_normalize(xyz, xyz);
                int ix=i+imp, iy=j+imp;
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
        for(int i=-imp;i<imp;i++){
            for(int j=-imp;j<imp;j++){
                int ix=i+imp, iy=j+imp;
                if(max_LPNH<mLPNH[0][0][ix][iy]) max_LPNH=mLPNH[0][0][ix][iy];
                if(min_LPNH>mLPNH[0][0][ix][iy]) min_LPNH=mLPNH[0][0][ix][iy];
                if(max_LPSH<mLPSH[0][0][ix][iy]) max_LPSH=mLPSH[0][0][ix][iy];
                if(min_LPSH>mLPSH[0][0][ix][iy]) min_LPSH=mLPSH[0][0][ix][iy];
            }
        }
        printf("mLPNH max: %.5f min: %.5f, mLPSH max: %.5f min: %.5f\n", max_LPNH, min_LPNH, max_LPSH, min_LPSH);
        printf("-> Average number of strong reflections = %d.\n", int(round(double(sum_strong)/double(kvec.numk))));
        printf("-> Average number of weak reflections = %d.\n", int(round(double(sum_weak)/double(kvec.numk))));
        ifinish=clock();
        printf("Execution time [s]: %.2f.\n", double(ifinish-istart)/CLOCKS_PER_SEC);
        hdf5(hdf5_path, iE);
        printf("Intermediate data stored in file %s.\n", hdf5_path);
    }
    finish=clock();
    printf("Execution time [s]: %.2f.\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("Intermediate data stored in file %s.\n", hdf5_path);
}

void DKD::compute_dynamic_matrix(complex<double> **dynmat, CELL *cell, DKD_GVECTOR *gvec)
{
    DKD_GNODE *rtemp=gvec->ghead->next;
    int imh=2*cell->HKL[0], imk=2*cell->HKL[1], iml=2*cell->HKL[2];
    double k0_2=2.0/cell->fouri0.lambda, k0_2i=0.5*cell->fouri0.lambda;
    for(int ir=0;ir<gvec->nstrong;ir++){
        DKD_GNODE *ctemp=gvec->ghead->next;
        for(int ic=0;ic<gvec->nstrong;ic++){
            int ih, ik, il;
            if(ic!=ir){
                complex<double> wsum(0.0, 0.0);
                DKD_GNODE *wtemp=gvec->headw;
                for(int iw=0;iw<gvec->nweak;iw++){
                    complex<double> Ughp, Uhph;
                    ih=rtemp->hkl[0]-wtemp->hkl[0]+imh;
                    ik=rtemp->hkl[1]-wtemp->hkl[1]+imk;
                    il=rtemp->hkl[2]-wtemp->hkl[2]+iml;
                    Ughp=cell->LUTUg[ih][ik][il];
                    ih=wtemp->hkl[0]-ctemp->hkl[0]+imh;
                    ik=wtemp->hkl[1]-ctemp->hkl[1]+imk;
                    il=wtemp->hkl[2]-ctemp->hkl[2]+iml;
                    Uhph=cell->LUTUg[ih][ik][il];
                    wsum+=Ughp*Uhph*complex<double>(1.0/wtemp->sg, 0.0);
                    wtemp=wtemp->nextw;
                }
                wsum*=complex<double>(k0_2i, 0.0);
                ih=rtemp->hkl[0]-ctemp->hkl[0]+imh;
                ik=rtemp->hkl[1]-ctemp->hkl[1]+imk;
                il=rtemp->hkl[2]-ctemp->hkl[2]+iml;
                dynmat[ir][ic]=cell->LUTUg[ih][ik][il]-wsum;
            }else{
                double wsum=0.0;
                DKD_GNODE *wtemp=gvec->headw;
                for(int iw=0;iw<gvec->nweak;iw++){
                    complex<double> Ughp;
                    ih=rtemp->hkl[0]-wtemp->hkl[0]+imh;
                    ik=rtemp->hkl[1]-wtemp->hkl[1]+imk;
                    il=rtemp->hkl[2]-wtemp->hkl[2]+iml;
                    Ughp=cell->LUTUg[ih][ik][il];
                    wsum+=fabs(Ughp)*fabs(Ughp)/wtemp->sg;
                    wtemp=wtemp->nextw;
                }
                wsum*=k0_2i;
                dynmat[ir][ir]=complex<double>(k0_2*rtemp->sg-wsum, cell->fouri0.Upmod);
            }
            ctemp=ctemp->nexts;
        }
        rtemp=rtemp->nexts;
    }
}

void DKD::compute_Sgh_matrices(complex<double>*** Sgh, CELL *cell, DKD_GVECTOR* gvec)
{
    int imh=2*cell->HKL[0], imk=2*cell->HKL[1], iml=2*cell->HKL[2];
    DKD_GNODE *rtemp=gvec->ghead->next;
    for(int ir=0;ir<gvec->nstrong;ir++){
        DKD_GNODE *ctemp=gvec->ghead->next;
        for(int ic=0;ic<gvec->nstrong;ic++){
            int ih, ik, il;
            ih=ctemp->hkl[0]-rtemp->hkl[0]+imh;
            ik=ctemp->hkl[1]-rtemp->hkl[1]+imk;
            il=ctemp->hkl[2]-rtemp->hkl[2]+iml;
            for(int i=0;i<napos;i++){
                Sgh[i][ir][ic]=cell->LUTSgh[i][ih][ik][il];
            }
            ctemp=ctemp->nexts;
        }
        rtemp=rtemp->nexts;
    }
}

void DKD::compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double *EWF, int IZMAX, double Z, double DZ, double KN, int NS)
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
    LWORK=std::min(LWMAX, int(creal(WORK[0]))); LDA=NS;
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//call the eigenvalue solver
    if(0!=INFO){
        printf("[ERROR] Unable to return INFO as zero when calling the eigenvalue solver.");
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
    double TPI=TWO_PI*DZ, DZT=DZ/Z;
    for(int j=0;j<NS;j++){
        for(int k=0;k<NS;k++){
            complex<double> sumq(0.0, 0.0);
            complex<double> q(TPI*(NW[j].imag()+NW[k].imag()), TPI*(NW[j].real()-NW[k].real()));
            if(q.real()<0.0) q=-q;
            for(int iz=0;iz<IZMAX;iz++){
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

double DKD::get_Lambert_interpolation(double xyz[3], double ****mat, bool is_hexagonal)
{
    int ix, iy, ixp, iyp;
    double dx, dy, dxm, dym;
    compute_Lambert_interpolation(xyz, nump/2, is_hexagonal, ix, iy, ixp, iyp, dx, dy, dxm, dym);
    double res=0.0;
    for(int i=0;i<napos;i++){
        res+=mat[i][0][ix][iy]*dxm*dym+mat[i][0][ixp][iy]*dx*dym+
             mat[i][0][ix][iyp]*dxm*dy+mat[i][0][ixp][iyp]*dx*dy;
    }
    return res;
}

void DKD::hdf5(const char *hdf5_path, int offset)
{
    size_t offsets_LP[4]={0, offset, 0, 0}, size_LP[4]={napos, 1, nump, nump};
    size_t offsets_SP[4]={offset, 0, 0}, size_SP[3]={1, nump, nump};
    HDF5 hdf;
    hdf.open(hdf5_path);
    hdf.write_group("/DynamicKikuchiDiffraction");
    hdf.write("/DynamicKikuchiDiffraction/SmallestInterplanarSpacing", dmin);
    hdf.write("/DynamicKikuchiDiffraction/PatternPixelNumber", nump);
    hdf.write_group("/DynamicKikuchiDiffraction/BetheParameters");
    hdf.write("/DynamicKikuchiDiffraction/BetheParameters/c1", bethe.c1);
    hdf.write("/DynamicKikuchiDiffraction/BetheParameters/c2", bethe.c2);
    hdf.write("/DynamicKikuchiDiffraction/BetheParameters/c3", bethe.c3);
    hdf.write("/DynamicKikuchiDiffraction/BetheParameters/c_sg", bethe.c_sg);
    hdf.write("/DynamicKikuchiDiffraction/AsymmetricNumber", napos);
    hdf.write("/DynamicKikuchiDiffraction/EnergyBinNumber", numEbin);
    hdf.write_hyper_array_4d("/DynamicKikuchiDiffraction/modifiedLPNH", mLPNH, size_LP, offsets_LP, napos, numEbin, nump, nump);
    hdf.write_hyper_array_4d("/DynamicKikuchiDiffraction/modifiedLPSH", mLPSH, size_LP, offsets_LP, napos, numEbin, nump, nump);
    hdf.write_hyper_array_3d("/DynamicKikuchiDiffraction/masterSPNH", mSPNH, size_SP, offsets_SP, numEbin, nump, nump);
    hdf.write_hyper_array_3d("/DynamicKikuchiDiffraction/masterSPSH", mSPSH, size_SP, offsets_SP, numEbin, nump, nump);
    hdf.close();
}

DKD::DKD(const char *hdf5_path)
{
    size_t size1, size2, size3, size4;
    HDF5 hdf;
    hdf.open(hdf5_path);
    hdf.read("/DynamicKikuchiDiffraction/SmallestInterplanarSpacing", dmin);
    hdf.read("/DynamicKikuchiDiffraction/PatternPixelNumber", nump);
    hdf.read("/DynamicKikuchiDiffraction/BetheParameters/c1", bethe.c1);
    hdf.read("/DynamicKikuchiDiffraction/BetheParameters/c2", bethe.c2);
    hdf.read("/DynamicKikuchiDiffraction/BetheParameters/c3", bethe.c3);
    hdf.read("/DynamicKikuchiDiffraction/BetheParameters/c_sg", bethe.c_sg);
    hdf.read("/DynamicKikuchiDiffraction/AsymmetricNumber", napos);
    hdf.read("/DynamicKikuchiDiffraction/EnergyBinNumber", numEbin);
    hdf.read_array_4d("/DynamicKikuchiDiffraction/modifiedLPNH", &mLPNH, size1, size2, size3, size4);
    hdf.read_array_4d("/DynamicKikuchiDiffraction/modifiedLPSH", &mLPSH, size1, size2, size3, size4);
    hdf.read_array_3d("/DynamicKikuchiDiffraction/masterSPNH", &mSPNH, size1, size2, size3);
    hdf.read_array_3d("/DynamicKikuchiDiffraction/masterSPSH", &mSPSH, size1, size2, size3);
    hdf.close();
}

DKD::~DKD()
{

}

void DKD::img(const char *img_path, double dimension, int resolution)
{
    char name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
    split_path(name, ext, img_path);
    if(0!=strcmp(ext, ".png")){
        printf("[ERROR] Unrecognized extension %s.", ext);
        exit(EXIT_FAILURE);
    }
    char png_path[PATH_CHAR_NUMBER]; 
    char exts[4][EXT_CHAR_NUMBER]; strcpy(exts[3], ext);
    for(int i=0;i<napos;i++){
        for(int j=0;j<numEbin;j++){
            int_to_str(exts[1], i); strcat(exts[1], "-"); int_to_str(exts[2], j);
            strcpy(exts[0], ".DKD_LPNH."); 
            strcpy(png_path, name); merge_path(png_path, exts, 4);
            image_array(png_path, mLPNH[i][j], nump, nump, dimension, dimension, resolution);
            printf("Image data for the modified lambert projection of northern hemisphere at apos %d, energy %d stored in %s.\n", i, j, png_path);
            strcpy(exts[0], ".DKD_LPSH."); 
            strcpy(png_path, name); merge_path(png_path, exts, 4);
            image_array(png_path, mLPSH[i][j], nump, nump, dimension, dimension, resolution);
            printf("Image data for the modified lambert projection of sorthern hemisphere at apos %d, energy %d stored in %s.\n", i, j, png_path);
        }
    }
    strcpy(exts[2], "");
    for(int i=0;i<numEbin;i++){
        strcpy(exts[0], ".DKD_SPNH."); 
        strcpy(png_path, name); merge_path(png_path, exts, 4);
        image_array(png_path, mSPNH[i], nump, nump, dimension, dimension, resolution);
        printf("Image data for the master stereographic projection of northern hemisphere at energy %d stored in %s.\n", i, png_path);
        strcpy(exts[0], ".DKD_SPSH."); 
        strcpy(png_path, name); merge_path(png_path, exts, 4);
        image_array(png_path, mSPSH[i], nump, nump, dimension, dimension, resolution);
        printf("Image data for the master stereographic projection of northern hemisphere at energy %d stored in %s.\n", i, png_path);
    }
}