#include "DED.h"

DED_GVECTOR::DED_GVECTOR(CELL *cell, BETHE *bethe, double kk[3], double fn[3], double kn, double dmin)
{
    int imh=cell->HKL[0], imk=cell->HKL[1], iml=cell->HKL[2];
    int cx=imh*2, cy=imk*2, cz=iml*2;
    double g[3]={0.0, 0.0, 0.0};
    complex<double> Ug=cell->LUTUg[cx][cy][cz];
    complex<double> qg=cell->LUTqg[cx][cy][cz];
    double sg=0.0;
    bool   is_double_diffrac=cell->is_double_diffrac[cx][cy][cz];
    add_g_node(g, Ug, qg, sg, is_double_diffrac);
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                if(0==abs(ih)+abs(ik)+abs(il)) continue;
                g[0]=double(ih); g[1]=double(ik); g[2]=double(il);
                double dhkl=cell->get_interplanar_spacing(g);
                if(cell->is_centering_allowed(g)&&(dhkl>dmin)){
                    int ix=ih+cx, iy=ik+cy, iz=il+cz;
                    Ug=cell->LUTUg[ix][iy][iz];
                    qg=cell->LUTqg[ix][iy][iz];
                    sg=cell->get_excitation_error(g, kk, fn);
                    is_double_diffrac=cell->is_double_diffrac[ix][iy][iz];                    
                    if(is_double_diffrac){
                        if(fabs(sg)<=bethe->c_sg){
                            add_g_node(g, Ug, qg, sg, is_double_diffrac);
                        }
                    }else{
                        double rg=kn*fabs(sg)/abs(Ug);
                        if(rg<=bethe->c3){
                            add_g_node(g, Ug, qg, sg, is_double_diffrac); 
                        }
                    }
                }
            }
        }
    }

    DED_GNODE *gtemp=ghead;
    int **garray; callocate_2d(&garray, numg, 3, 0);
    for(int i=0;i<numg&&gtemp!=nullptr;i++){
        for(int j=0;j<3;j++){
            garray[i][j]=gtemp->hkl[j];
        }
        gtemp=gtemp->next;
    }

    DED_GNODE *tails=nullptr, *tailw=nullptr;
    gtemp=heads=tails=ghead;
    gtemp->is_strong=true;
    gtemp->is_weak=false;
    nstrong++;
    for(int i=1;i<numg;i++){
        gtemp=gtemp->next;
        double imin=1.0e8;
        double sgp=kn*fabs(gtemp->sg);
        for(int j=0;j<numg;j++){
            double temp;
            int ix=garray[i][0]-garray[j][0]+cx, iy=garray[i][1]-garray[j][1]+cy, iz=garray[i][2]-garray[j][2]+cz;
            if(cell->is_double_diffrac[ix][iy][iz]){
                temp=1.0e4;
            }else{
                temp=sgp/abs(cell->LUTUg[ix][iy][iz]);
            }
            if(temp<imin) imin=temp;
        }
        if(imin>bethe->c2){//ignore the reflection
            gtemp->is_weak=false;
            gtemp->is_strong=false;
            continue;
        }
        if((imin>bethe->c1)&&(imin<=bethe->c2)){//weak
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
        if(imin<=bethe->c1){//strong
            tails->nexts=gtemp;
            tails=gtemp;
            gtemp->is_weak=false;
            gtemp->is_strong=true;
            nstrong++;
        }
    }
    deallocate_2d(garray, numg);
}

DED_GVECTOR::~DED_GVECTOR()
{
    DED_GNODE *cur=ghead;
    while(cur!=nullptr){
        DED_GNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void DED_GVECTOR::add_g_node(double hkl[3], complex<double> Ug, complex<double> qg, double sg, bool is_double_diffrac)
{
    if(ghead==nullptr&&gtail==nullptr){
        ghead=gtail=new DED_GNODE;
        vector_copy(gtail->hkl, hkl);
        gtail->Ug=Ug;
        gtail->qg=qg;
        gtail->sg=sg;
        gtail->is_double_diffrac=is_double_diffrac;
    }else{
        gtail->next=new DED_GNODE;
        gtail=gtail->next;
        vector_copy(gtail->hkl, hkl);
        gtail->Ug=Ug;
        gtail->qg=qg;
        gtail->sg=sg;
        gtail->is_double_diffrac=is_double_diffrac;
    }
    numg++;
}

void copy_knode_data(DED_KNODE *knode1, DED_KNODE *knode2)
{
    vector_copy(knode1->hkl, knode2->hkl);
    vector_copy(knode1->K, knode2->K);
    knode1->Kmagnitude=knode2->Kmagnitude;
    knode1->intensity=knode2->intensity;
}

void swap_knode_data(DED_KNODE *knode1, DED_KNODE *knode2)
{
    DED_KNODE *ktemp=new DED_KNODE;
    copy_knode_data(ktemp, knode1);
    copy_knode_data(knode1, knode2);
    copy_knode_data(knode2, ktemp);
    delete ktemp;
}

void quick_sort(DED_KNODE *kstart, DED_KNODE *kend)
{
    if(kstart==nullptr||kend==nullptr||kstart==kend) return;
    DED_KNODE *knode1=kstart;
    DED_KNODE *knode2=kstart->next;
    double Kmagnitude=kstart->Kmagnitude;
    while(knode2!=kend->next&&knode2!=nullptr){
        if(knode2->Kmagnitude<Kmagnitude){
            knode1=knode1->next;
            if(knode1!=knode2){
                swap_knode_data(knode1, knode2);
            }
        }
        knode2=knode2->next;
    }
    swap_knode_data(knode1, kstart);
    quick_sort(kstart, knode1);
    quick_sort(knode1->next, kend);
}

DED::DED(CELL *cell, BETHE *bethe, int zone[3], int fnorm[3], double fthick, double voltage, double Kmag_max, double threshold)
{
    Kmagnitude_max=Kmag_max*10.0;
    double dmin=1.0/Kmag_max*0.1;
    double ft=fthick*0.1;
    printf("[INFO] Starting computation of dynamical electron diffraction...\n");
    cell->compute_reflection_range(dmin);
    cell->compute_Fourier_coefficients(voltage);
    printf("[INFO] Range of Miller indices along a*, b*, or c* in reciprocal space: %d %d %d\n", cell->HKL[0], cell->HKL[1], cell->HKL[2]);
    double kn=1.0/cell->fouri0.lambda;
    double kk[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    double fn[3]={double(fnorm[0]), double(fnorm[1]), double(fnorm[2])};
    cell->direct_to_reciprocal(kk, kk);
    cell->normalize(kk, kk, 'r');
    vector_constant(kk, kn, kk);
    cell->direct_to_reciprocal(fn, fn);
    cell->normalize(fn, fn, 'r');
    vector_constant(fn, kn, fn);
    DED_GVECTOR gvec(cell, bethe, kk, fn, kn, dmin);
    printf("[INFO] Total number of reflections, strong reflections, and weak reflections: %d %d %d\n", gvec.numg, gvec.nstrong, gvec.nweak);
    complex<double> **dmat;
    callocate_2d(&dmat, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
    compute_dynamic_matrix(dmat, cell, &gvec);
    double *intens;
    callocate(&intens, gvec.nstrong, 0.0);
    compute_diffraction_intensity(intens, dmat, ft, kn, gvec.nstrong);
    DED_GNODE *gtemp=gvec.heads;
    for(int i=0;i<gvec.nstrong&&gtemp!=nullptr;i++){
        if(intens[i]>DED_INTENSITY_LIMIT){
            add_k_node(cell, gtemp->hkl, intens[i]);
        }
        gtemp=gtemp->nexts;
    }
    printf("[INFO] Intensity at the transmission spot: %.8f\n", khead->intensity);
    printf("[INFO] Number of diffraction intensity (including intensity at the transmission spot): %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    filter_diffraction_intensity(threshold);
    printf("[INFO] Number of filtered diffraction intensity (including intensity at the transmission spot): %d\n", numk);
    printf("[INFO] Range of filtered diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    find_first_and_second_knearests();
    if(knearest_2==nullptr){
        printf("[INFO] The first nearest diffraction vectors R1: [%.8f %.8f %.8f]\n", knearest_1->K[0]*0.1, knearest_1->K[1]*0.1, knearest_1->K[2]*0.1);
        printf("[WARN] Unable to find the second nearest diffraction vector\n");
    }else{
        printf("[INFO] The first and second nearest diffraction vectors R1, R2: [%.8f %.8f %.8f], [%.8f %.8f %.8f] (R2/R1 %.8f and angle %.8f)\n", 
        knearest_1->K[0]*0.1, knearest_1->K[1]*0.1, knearest_1->K[2]*0.1, knearest_2->K[0]*0.1, knearest_2->K[1]*0.1, knearest_2->K[2]*0.1, knearest_2->Kmagnitude/knearest_1->Kmagnitude, vector_angle(knearest_2->K, knearest_1->K)*RAD_TO_DEG);
    }
    rotate_by_first_knearest(cell, zone);
    printf("[INFO] Ending computation of kinematic electron diffraction\n");
}

DED::~DED()
{
    DED_KNODE *cur=khead;
    while(cur!=nullptr){
        DED_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void DED::add_k_node(CELL *cell, double hkl[3], double intensity)
{
    if(khead==nullptr&&ktail==nullptr){
        khead=ktail=new DED_KNODE;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, hkl);
        cell->reciprocal_to_cartesian(ktail->K, ktail->K);
        ktail->Kmagnitude=vector_length(ktail->K);
        ktail->intensity=intensity;
    }else{
        ktail->next=new DED_KNODE;
        ktail=ktail->next;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, hkl);
        cell->reciprocal_to_cartesian(ktail->K, ktail->K);
        ktail->Kmagnitude=vector_length(ktail->K);
        ktail->intensity=intensity;
        if(intensity<intensity_min) intensity_min=intensity;
        if(intensity>intensity_max) intensity_max=intensity;
    }
    numk++;
}

void DED::compute_dynamic_matrix(complex<double> **dmat, CELL *cell, DED_GVECTOR *gvec)
{
    DED_GNODE *rtemp=gvec->ghead;
    int imh=2*cell->HKL[0], imk=2*cell->HKL[1], iml=2*cell->HKL[2];
    double k0_2=2.0/cell->fouri0.lambda, k0_2i=0.5*cell->fouri0.lambda;
    for(int ir=0;ir<gvec->nstrong;ir++){
        DED_GNODE *ctemp=gvec->ghead;
        for(int ic=0;ic<gvec->nstrong;ic++){
            int ih, ik, il;
            if(ic!=ir){
                complex<double> wsum(0.0, 0.0);
                DED_GNODE *wtemp=gvec->headw;
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
                dmat[ir][ic]=cell->LUTUg[ih][ik][il]-wsum;
            }else{
                double wsum=0.0;
                DED_GNODE *wtemp=gvec->headw;
                for(int iw=0;iw<gvec->nweak;iw++){
                    complex<double> Ughp;
                    ih=rtemp->hkl[0]-wtemp->hkl[0]+imh;
                    ik=rtemp->hkl[1]-wtemp->hkl[1]+imk;
                    il=rtemp->hkl[2]-wtemp->hkl[2]+iml;
                    Ughp=cell->LUTUg[ih][ik][il];
                    wsum+=abs(Ughp)*abs(Ughp)/wtemp->sg;
                    wtemp=wtemp->nextw;
                }
                wsum*=k0_2i;
                dmat[ir][ir]=complex<double>(k0_2*rtemp->sg-wsum, cell->fouri0.Upmod);
            }
            ctemp=ctemp->nexts;
        }
        rtemp=rtemp->nexts;
    }
}

void DED::compute_diffraction_intensity(double *INTENS, complex<double> **DMAT, double Z, double KN, int NS)
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
            A[i*LDA+j]=DMAT[i][j];
        }
    }
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//determine the optimal LWORK, i.e., WORK[0]
    LWORK=m_min(LWMAX, int(WORK[0].real())); LDA=NS;
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//call the eigenvalue solver
    if(0!=INFO){
        printf("[ERROR] Unable to return INFO as zero when calling the eigenvalue solver.");
        exit(1);
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
    MILWORK=int(GETMILWORK[0].real());
    lapack_complex_double *MIWORK; callocate(&MIWORK, MILWORK, lapack_complex_double(0.0, 0.0));
    LDA=NS;
    INFO=LAPACKE_zgetri_work(LAPACK_ROW_MAJOR, NS, CGINV, LDA, IPIV, MIWORK, MILWORK);
    deallocate(IPIV); deallocate(MIWORK);

    complex<double> *NW; callocate(&NW, NS, complex<double>(0.0, 0.0));
    complex<double> F(2.0*KN, 0.0);
    for(int i=0;i<NS;i++){
        NW[i]=W[i];
        NW[i]=NW[i]/F;
    }
    complex<double> **NCGINV; callocate_2d(&NCGINV, NS, NS, complex<double>(0.0, 0.0));
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCGINV[i][j]=CGINV[i*NS+j];
        }
    }
    complex<double> **NCG; callocate_2d(&NCG, NS, NS, complex<double>(0.0, 0.0));
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCG[i][j]=CG[i*NS+j];
        }
    }
    complex<double> TPI(0.0, TWO_PI*Z);
    for(int j=0;j<NS;j++){
        complex<double> PW(0.0, 0.0);
        for(int k=0;k<NS;k++){
            PW+=(NCGINV[k][0]*NCG[j][k]*exp(TPI*NW[k]));
        }
        INTENS[j]=abs(PW)*abs(PW);
    }
}

void DED::filter_diffraction_intensity(double threshold)
{
    double intensity_threshold=intensity_max*threshold;
    intensity_min=1.0e8;
    while(khead!=nullptr&&khead->intensity<=intensity_threshold){
        DED_KNODE *ktemp=khead;
        khead=khead->next;
        delete ktemp;
        numk--;
    }
    DED_KNODE *cur=khead;
    while(cur!=nullptr&&cur->next!=nullptr){
        if(cur->next->intensity<=intensity_threshold){
            DED_KNODE *ktemp=cur->next;
            cur->next=cur->next->next;
            delete ktemp;
            numk--;
        }else{
            if(intensity_min>cur->intensity) intensity_min=cur->intensity;
            cur=cur->next;
        }
    }
}

void DED::find_first_and_second_knearests()
{
    quick_sort(khead, ktail);
    knearest_1=khead->next;
    if(knearest_1==nullptr){
        printf("[ERROR] Unable to find the first nearest diffraction vector\n");
        exit(1);
    }
    knearest_2=knearest_1->next;
    while(knearest_2!=nullptr){
        if(knearest_2->Kmagnitude-knearest_1->Kmagnitude>DED_KMAGNITUDE_LIMIT&&vector_angle(knearest_2->K, knearest_1->K)<=HALF_PI&&vector_angle(knearest_2->K, knearest_1->K)>1.0e-6) break;
        knearest_2=knearest_2->next;
    }
}

void DED::rotate_by_first_knearest(CELL *cell, int zone[3])
{
    axes[2][0]=double(zone[0]); axes[2][1]=double(zone[1]); axes[2][2]=double(zone[2]);
    vector_normalize(axes[2], axes[2]);
    vector_copy(axes[0], knearest_1->K);
    vector_normalize(axes[0], axes[0]);
    vector_cross(axes[1], axes[0], axes[2]);
    vector_normalize(axes[1], axes[1]);
    vector_cross(axes[0], axes[1], axes[2]);
    vector_normalize(axes[0], axes[0]);
}

void DED::rotate(double x[3], double y[3])
{
    vector_copy(axes[0], x);
    vector_normalize(axes[0], axes[0]);
    vector_copy(axes[1], y);
    vector_normalize(axes[1], axes[1]);
    if(fabs(vector_dot(axes[0], axes[1]))>1.0e-6||fabs(vector_dot(axes[1], axes[2]))>1.0e-6||fabs(vector_dot(axes[0], axes[2]))>1.0e-6){
        printf("[ERROR] The orthogonality condition is not satisfied with x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f]", x[0], x[1], x[2], y[0], y[1], y[2]);
        exit(1);
    }
}

void DED::ded(char *ded_path)
{
    FILE *fp=nullptr;
    fp=fopen(ded_path,"w");
    fprintf(fp, "# h\tk\tl\tK_1\tK_2\tK_3\tx\ty\tz\tintensity\tintensity_norm (%d points, x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f], and z-[%.8f %.8f %.8f])\n", numk-1, 
            axes[0][0], axes[0][1], axes[0][2], axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
    double *pos_x=nullptr, *pos_y=nullptr, *intensity=nullptr;
    callocate(&pos_x, numk-1, 0.0); 
    callocate(&pos_y, numk-1, 0.0); 
    callocate(&intensity, numk-1, 0.0);
    DED_KNODE *ktemp=khead->next;
    double constn=100.0/intensity_max;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        double xyz[3]; vector_rotate(xyz, axes, ktemp->K);
        double intensity_norm=constn*ktemp->intensity;
        fprintf(fp, "%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", int(ktemp->hkl[0]), int(ktemp->hkl[1]), int(ktemp->hkl[2]), 
                ktemp->K[0]*0.1, ktemp->K[1]*0.1, ktemp->K[2]*0.1, xyz[0]*0.1, xyz[1]*0.1, xyz[2]*0.1, ktemp->intensity, intensity_norm);
        fflush(fp);
        pos_x[i-1]=xyz[0]*0.1; pos_y[i-1]=xyz[1]*0.1; intensity[i-1]=intensity_norm;
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s.\n", ded_path);

    char png_path[strlen(ded_path)+5];
    strcpy(png_path, ded_path); strcat(png_path, ".png");
    img(png_path, pos_x, pos_y, intensity, numk-1, Kmagnitude_max*0.1);
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void DED::ded(char *ded_path, double sigma, double dx)
{
    int    nbin=2*round(Kmagnitude_max*0.1/dx)+1;
    int    nbin_half=nbin/2;
    double *pos_x=nullptr, **intensity=nullptr;
    callocate(&pos_x, nbin, 0.0);
    callocate_2d(&intensity, nbin, nbin, 0.0);
    for(int i=0;i<=nbin_half;i++){
        pos_x[nbin_half+i]=double(i)*dx; 
        pos_x[nbin_half-i]=-double(i)*dx;
    }

    double **intensity_c=nullptr; 
    callocate_2d(&intensity_c, nbin, nbin, 0.0);
    DED_KNODE *ktemp=khead->next;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        double xyz[3]; vector_rotate(xyz, axes, ktemp->K);
        int m=round(xyz[0]*0.1/dx)+nbin_half;
        int n=round(xyz[1]*0.1/dx)+nbin_half;
        gaussian(intensity_c, pos_x, pos_x, nbin, ktemp->intensity, pos_x[m], pos_x[n], sigma);
        for(int j=0;j<nbin;j++){
            for(int k=0;k<nbin;k++){
                intensity[j][k]+=intensity_c[j][k];
                intensity_c[j][k]=0.0;
            }
        }
        ktemp=ktemp->next;
    }
    deallocate_2d(intensity_c, nbin);

    double imax=0.0, imin=1.0e8;
    for(int i=0;i<nbin;i++){
        for(int j=0;j<nbin;j++){
            if(imax<intensity[i][j]) imax=intensity[i][j];
            if(imin>intensity[i][j]) imin=intensity[i][j];
        }
    }
    printf("[INFO] Number of profiled diffraction intensity: %d\n", nbin*nbin);
    printf("[INFO] Range of profiled diffraction intensity: %.8f %.8f\n", imin, imax);

    FILE *fp=nullptr;
    fp=fopen(ded_path,"w");
    fprintf(fp, "# x\ty\tintensity\tintensity_norm (%d points, rotated by x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f], and z-[%.8f %.8f %.8f])\n", nbin*nbin,
            axes[0][0], axes[0][1], axes[0][2], axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
    double constn=100.0/imax;
    for(int i=0;i<nbin;i++){
        for(int j=0;j<nbin;j++){
            fprintf(fp, "%.8f\t%.8f\t%.8f\t", pos_x[j], pos_x[i], intensity[i][j]);
            intensity[i][j]*=constn;
            fprintf(fp, "%.8f\n", intensity[i][j]);
            fflush(fp);
        }
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", ded_path);

    char png_path[strlen(ded_path)+5];
    strcpy(png_path, ded_path); strcat(png_path, ".png");
    img(png_path, intensity, nbin, nbin, imax, imin);
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void DED::img(char *png_path, double *x, double *y, double *value, int num, double limit)
{
    double height=6.0, width=6.0;
    int tick_max=int(limit);
    int n_major_tick=2*tick_max+1;
    double *major_ticks; mallocate(&major_ticks, n_major_tick);
    for(int i=0;i<n_major_tick;i++){
        major_ticks[i]=double(-tick_max+i);
    }
    GRAPH graph(width, height, 300);
    graph.set_xlim(-limit, limit);
    graph.set_ylim(-limit, limit);
    graph.set_xticks(major_ticks, n_major_tick);
    graph.set_yticks(major_ticks, n_major_tick);
    graph.set_tick_in(false);
    graph.scatter(x, y, value, num);
    graph.draw(png_path);
}

void DED::img(char* png_path, double **value, int numpx, int numpy, double vmax, double vmin, char background)
{
    double *wdata;
    unreshape_2d(&wdata, value, numpy, numpx);
    unsigned char *pixels;
    int num=numpx*numpy;
    mallocate(&pixels, 3*num);

    double vdiff=vmax-vmin;
    unsigned char rgb_min[3]={0}, rgb_max[3]={255};
    switch(background)
    {
    case 'b':
        rgb_min[0]=rgb_min[1]=rgb_min[2]=0;
        rgb_max[0]=rgb_max[1]=rgb_max[2]=255;
        break;
    case 'w':
        rgb_min[0]=rgb_min[1]=rgb_min[2]=255;
        rgb_max[0]=rgb_max[1]=rgb_max[2]=0;
        break;
    default:
        printf("[ERROR] Unrecognized background %s\n", background);
        exit(1);
    }
    unsigned char rgb_diff[3]; vector_difference(rgb_diff, rgb_max, rgb_min);
    unsigned char rgb[3];
    for(int i=0;i<num;i++){
        if(wdata[i]>=vmax){
            vector_copy(rgb, rgb_max);
        }else if(wdata[i]<=vmin){
            vector_copy(rgb, rgb_min);
        }else{
            rgb[0]=rgb_min[0]+int((wdata[i]-vmin)/vdiff*rgb_diff[0]);
            rgb[1]=rgb_min[1]+int((wdata[i]-vmin)/vdiff*rgb_diff[1]);
            rgb[2]=rgb_min[2]+int((wdata[i]-vmin)/vdiff*rgb_diff[2]);
        }
        pixels[i*3]=rgb[0]; pixels[i*3+1]=rgb[1]; pixels[i*3+2]=rgb[2];
    }
    image_pixels(png_path, pixels, numpx, numpy);
}

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
    double xy[2]={0.0, 0.0};
    add_k_node(cell, xy, kn);
    switch(cell->sampling_type)
    {
    case 1://triclinic 1
        for(int j=-npy;j<=npy;j++){
            for(int i=-npx;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 2://triclinic -1
        for(int j=-npy;j<=npy;j++){
            for(int i=-npx;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j);
            }
        }
        break;
    case 3://monoclinic 2
        for(int j=-npy;j<=npy;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 4://monoclinic m
        for(int j=0;j<=npy;j++){
            for(int i=-npx;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 5://monoclinic 2/m, orthorhombic 222, mm2, tetragonal 4, -4
        for(int j=0;j<=npy;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 6://orthorhombic mmm, tetragonal 4/m, 422, -4m2, cubic m-3, 432
        for(int j=0;j<=npy;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j);
            }
        }
        break;
    case 7://tetragonal 4mm
        for(int i=0;i<=npx;i++){
            for(int j=0;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j, true);
            }
        }
        break;
    case 8://tetragonal -42m, cubic -43m
        for(int i=0;i<=npx;i++){
            for(int j=-i;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j);
            }
        }
        break;
    case 9://tetragonal 4/mmm, cubic m-3m
        for(int i=0;i<=npx;i++){
            for(int j=0;j<=i;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                add_k_node(cell, xy, kn, i, j);
            }
        }
        break;
    case 10://hexagonal 3
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_node(cell, xy, kn, i, j, true);
                }
            }
        }
        break;
    case 11://rhombohedral 3
        printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
        exit(1);
    case 12://hexagonal -3, 321, -6 [not implemented: rhombohedral 32]
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_node(cell, xy, kn, i, j);
                }
            }
        }
        break;
    case 13://hexagonal 312 [not implemented: rhombohedral -3]
        for(int j=0;j>=-npx;j--){
            for(int i=j/2;i<npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_node(cell, xy, kn, i, j);
                }
            }
        }
        for(int i=0;i<=npx;i++){
            for(int j=0;j<i/2;j++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_node(cell, xy, kn, i, j);
                }
            }
        }
        break;
    case 14://hexagonal 3m1 [not implemented: rhombohedral 3m]
        for(int j=1;j<=npx;j++){
            for(int i=j;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_node(cell, xy, kn, i, j, true);
                }
            }
        }
        break;
    case 15://hexagonal 31m, 6
        for(int j=1;j<=npx;j++){
            for(int i=j;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                if(is_in_hexagonal_grid(xy)){
                    add_k_node(cell, xy, kn, i, j, true);
                }
            }
        }
        break;
    case 16://hexagonal -3m1, 622, -6m2 [not implemented: rhombohedral -3m]
        for(int j=0;j<=npx;j++){
            for(int i=0;i<=npx;i++){
                xy[0]=i*delta; xy[1]=j*delta;
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)<(PI/6.0-1.0e-4)))){
                    continue;
                }else if(is_in_hexagonal_grid(xy)){
                    add_k_node(cell, xy, kn, i, j);
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
                    add_k_node(cell, xy, kn, i, j);
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
                    add_k_node(cell, xy, kn, i, j, true);
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
                    add_k_node(cell, xy, kn, i, j);
                }
            }
        }
        break;
    default:
        printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
        exit(1);
    }
}

DKD_KVECTOR::~DKD_KVECTOR()
{
    DKD_KNODE *cur=khead;
    while(cur!=nullptr){
        DKD_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void DKD_KVECTOR::add_k_node(CELL *cell, double xy[2], double kn, int i, int j, bool southern_flag)
{
    if(ktail==nullptr){
        khead=ktail=new DKD_KNODE;
        double kstar[3]={0.0, 0.0, kn}, r_kstar[3];//c* as the center of the Rosca-Lambert projection
        cell->cartesian_to_reciprocal(r_kstar, kstar);
        khead->k[0]=r_kstar[0]; khead->k[1]=r_kstar[1]; khead->k[2]=r_kstar[2];
        khead->kn=kn;
        khead->i=i; khead->j=j;
        khead->hemisphere=1;//Northern hemisphere
        numk++;
    }else{
        double kstar[3], r_kstar[3];
        int ierr;
        ktail->next=new DKD_KNODE;
        ktail=ktail->next;
        if(cell->use_hexagonal){
            compute_sphere_from_hexagonal_Lambert(kstar, ierr, xy);
        }else{
            compute_sphere_from_square_Lambert(kstar, ierr, xy);
        }
        vector_normalize(kstar, kstar);
        kstar[0]*=kn; kstar[1]*=kn; kstar[2]*=kn;
        cell->cartesian_to_reciprocal(r_kstar, kstar);
        ktail->k[0]=r_kstar[0]; ktail->k[1]=r_kstar[1]; ktail->k[2]=r_kstar[2];
        ktail->kn=kn;
        ktail->i=i; ktail->j=j;
        ktail->hemisphere=1;
        numk++;
        if(southern_flag){
            kstar[0]=-kstar[0]; kstar[1]=-kstar[1]; kstar[2]=-kstar[2];
            cell->cartesian_to_reciprocal(r_kstar, kstar);
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

DKD::DKD(CELL *cell, MC *mc, BETHE *bethe, double Kmag_max, char *projection)
{
    numpx=numpy=mc->numpE;
    double dmin=1.0/Kmag_max*0.1;
    cell->compute_reflection_range(dmin);
    cell->compute_Bloch_wave_coefficients();
    printf("[INFO] Range of Miller indices along a*, b*, or c* in reciprocal space: %d %d %d\n", cell->HKL[0], cell->HKL[1], cell->HKL[2]);

    clock_t start, finish; start=clock();
    double voltage=mc->EkeV;
    cell->compute_Fourier_coefficients(voltage);

    DKD_KVECTOR kvec(cell, numpx/2, numpy/2);
    printf("[INFO] Independent beam directions to be considered = %d\n", kvec.numk);

    int sum_strong=0, sum_weak=0;
    DKD_KNODE *ktemp=kvec.khead;
    for(int ik=0;ik<kvec.numk;ik++){
        DED_GVECTOR gvec(cell, bethe, ktemp->k, ktemp->k, ktemp->kn, dmin);
        sum_strong+=gvec.nstrong; sum_weak+=gvec.nweak;

        complex<double> ***Sgh;
        callocate_3d(&Sgh, cell->napos, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
        compute_Sgh_matrices(Sgh, cell, &gvec);

        complex<double> **dmat, **Lgh;
        callocate_2d(&dmat, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
        callocate_2d(&Lgh, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
        compute_dynamic_matrix(dmat, cell, &gvec);
        compute_Lgh_matrix(Lgh, dmat, mc->weight, mc->izmax, mc->z, mc->dz, ktemp->kn, gvec.nstrong);
        
        ktemp->intensity=0.0;
        for(int i=0;i<cell->napos;i++){
            for(int j=0;j<gvec.nstrong;j++){
                for(int k=0;k<gvec.nstrong;k++){
                    complex<double> temp=Lgh[j][k]*Sgh[i][j][k];
                    ktemp->intensity+=temp.real();
                }
            }
        }
        ktemp->intensity/=double(cell->npos);
        ktemp=ktemp->next;

        if(0==(ik+1)%1000){
            printf("[INFO] Completed beam direction %d of %d.\n", ik+1, kvec.numk);
        }
        deallocate_2d(Lgh, gvec.nstrong);
        deallocate_2d(dmat, gvec.nstrong);
        deallocate_3d(Sgh, cell->napos, gvec.nstrong);
    }
    finish=clock();
    printf("[INFO] Execution time [s]: %.2f.\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("[INFO] Average number of strong reflections = %d.\n", int(round(double(sum_strong)/double(kvec.numk))));
    printf("[INFO] Average number of weak reflections = %d.\n", int(round(double(sum_weak)/double(kvec.numk))));
    compute_kvector_projection(cell, &kvec, projection, cell->use_hexagonal);
    printf("[INFO] Range of intensity on the projection of northern hemisphere: %.8f, %.8f\n", intensity_minN, intensity_maxN);
    printf("[INFO] Range of intensity on the projection of southern hemisphere: %.8f, %.8f\n", intensity_minS, intensity_maxS);

}

DKD::DKD(CELL *cell, BETHE *bethe, double voltage, double fthick, double Kmag_max, int nump, char *projection)
{
    numpx=numpy=nump;
    double dmin=1.0/Kmag_max*0.1;
    double ft=fthick*0.1;
    cell->compute_reflection_range(dmin);
    cell->compute_Bloch_wave_coefficients();
    printf("[INFO] Range of Miller indices along a*, b*, or c* in reciprocal space: %d %d %d\n", cell->HKL[0], cell->HKL[1], cell->HKL[2]);

    clock_t start, finish; start=clock();
    cell->compute_Fourier_coefficients(voltage);

    DKD_KVECTOR kvec(cell, numpx/2, numpy/2);
    printf("[INFO] Independent beam directions to be considered = %d\n", kvec.numk);

    int sum_strong=0, sum_weak=0;
    DKD_KNODE *ktemp=kvec.khead;
    for(int ik=0;ik<kvec.numk;ik++){
        DED_GVECTOR gvec(cell, bethe, ktemp->k, ktemp->k, ktemp->kn, dmin);
        sum_strong+=gvec.nstrong; sum_weak+=gvec.nweak;

        complex<double> ***Sgh;
        callocate_3d(&Sgh, cell->napos, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
        compute_Sgh_matrices(Sgh, cell, &gvec);

        complex<double> **dmat, **Lgh;
        callocate_2d(&dmat, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
        callocate_2d(&Lgh, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
        compute_dynamic_matrix(dmat, cell, &gvec);
        compute_Lgh_matrix(Lgh, dmat, ft, ktemp->kn, gvec.nstrong);
        
        ktemp->intensity=0.0;
        for(int i=0;i<cell->napos;i++){
            for(int j=0;j<gvec.nstrong;j++){
                for(int k=0;k<gvec.nstrong;k++){
                    complex<double> temp=Lgh[j][k]*Sgh[i][j][k];
                    ktemp->intensity+=temp.real();
                }
            }
        }
        ktemp->intensity/=double(cell->npos);
        ktemp=ktemp->next;

        if(0==(ik+1)%1000){
            printf("[INFO] Completed beam direction %d of %d.\n", ik+1, kvec.numk);
        }
        deallocate_2d(Lgh, gvec.nstrong);
        deallocate_2d(dmat, gvec.nstrong);
        deallocate_3d(Sgh, cell->napos, gvec.nstrong);
    }
    finish=clock();
    printf("[INFO] Execution time [s]: %.2f.\n", double(finish-start)/CLOCKS_PER_SEC);
    printf("[INFO] Average number of strong reflections = %d.\n", int(round(double(sum_strong)/double(kvec.numk))));
    printf("[INFO] Average number of weak reflections = %d.\n", int(round(double(sum_weak)/double(kvec.numk))));
    compute_kvector_projection(cell, &kvec, projection, cell->use_hexagonal);
    printf("[INFO] Range of intensity on the projection of northern hemisphere: %.8f, %.8f\n", intensity_minN, intensity_maxN);
    printf("[INFO] Range of intensity on the projection of southern hemisphere: %.8f, %.8f\n", intensity_minS, intensity_maxS);

}

void DKD::compute_dynamic_matrix(complex<double> **dmat, CELL *cell, DED_GVECTOR *gvec)
{
    DED_GNODE *rtemp=gvec->ghead;
    int imh=2*cell->HKL[0], imk=2*cell->HKL[1], iml=2*cell->HKL[2];
    double k0_2=2.0/cell->fouri0.lambda, k0_2i=0.5*cell->fouri0.lambda;
    for(int ir=0;ir<gvec->nstrong;ir++){
        DED_GNODE *ctemp=gvec->ghead;
        for(int ic=0;ic<gvec->nstrong;ic++){
            int ih, ik, il;
            if(ic!=ir){
                complex<double> wsum(0.0, 0.0);
                DED_GNODE *wtemp=gvec->headw;
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
                dmat[ir][ic]=cell->LUTUg[ih][ik][il]-wsum;
            }else{
                double wsum=0.0;
                DED_GNODE *wtemp=gvec->headw;
                for(int iw=0;iw<gvec->nweak;iw++){
                    complex<double> Ughp;
                    ih=rtemp->hkl[0]-wtemp->hkl[0]+imh;
                    ik=rtemp->hkl[1]-wtemp->hkl[1]+imk;
                    il=rtemp->hkl[2]-wtemp->hkl[2]+iml;
                    Ughp=cell->LUTUg[ih][ik][il];
                    wsum+=abs(Ughp)*abs(Ughp)/wtemp->sg;
                    wtemp=wtemp->nextw;
                }
                wsum*=k0_2i;
                dmat[ir][ir]=complex<double>(k0_2*rtemp->sg-wsum, cell->fouri0.Upmod);
            }
            ctemp=ctemp->nexts;
        }
        rtemp=rtemp->nexts;
    }
}

void DKD::compute_Sgh_matrices(complex<double>*** Sgh, CELL *cell, DED_GVECTOR* gvec)
{
    int imh=2*cell->HKL[0], imk=2*cell->HKL[1], iml=2*cell->HKL[2];
    DED_GNODE *rtemp=gvec->ghead;
    for(int ir=0;ir<gvec->nstrong;ir++){
        DED_GNODE *ctemp=gvec->ghead;
        for(int ic=0;ic<gvec->nstrong;ic++){
            int ih, ik, il;
            ih=ctemp->hkl[0]-rtemp->hkl[0]+imh;
            ik=ctemp->hkl[1]-rtemp->hkl[1]+imk;
            il=ctemp->hkl[2]-rtemp->hkl[2]+iml;
            for(int i=0;i<cell->napos;i++){
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
            A[i*LDA+j]=DMAT[i][j];
        }
    }
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//determine the optimal LWORK, i.e., WORK[0]
    LWORK=m_min(LWMAX, int(WORK[0].real())); LDA=NS;
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//call the eigenvalue solver
    if(0!=INFO){
        printf("[ERROR] Unable to return INFO as zero when calling the eigenvalue solver.");
        exit(1);
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
    MILWORK=int(GETMILWORK[0].real());
    lapack_complex_double *MIWORK; callocate(&MIWORK, MILWORK, lapack_complex_double(0.0, 0.0));
    LDA=NS;
    INFO=LAPACKE_zgetri_work(LAPACK_ROW_MAJOR, NS, CGINV, LDA, IPIV, MIWORK, MILWORK);
    deallocate(IPIV); deallocate(MIWORK);

    complex<double> *NW; callocate(&NW, NS, complex<double>(0.0, 0.0));
    complex<double> F(2.0*KN, 0.0);
    for(int i=0;i<NS;i++){
        NW[i]=W[i];
        NW[i]=NW[i]/F;
    }
    complex<double> **NCGINV; callocate_2d(&NCGINV, NS, NS, complex<double>(0.0, 0.0));
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCGINV[i][j]=CGINV[i*NS+j];
        }
    }
    complex<double> **NCG; callocate_2d(&NCG, NS, NS, complex<double>(0.0, 0.0));
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCG[i][j]=CG[i*NS+j];
        }
    }
    complex<double> **IJK; callocate_2d(&IJK, NS, NS, complex<double>(0.0, 0.0));
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
    complex<double> **NCGCONJ; callocate_2d(&NCGCONJ, NS, NS, complex<double>(0.0, 0.0));
    complex<double> **NCGT; callocate_2d(&NCGT, NS, NS, complex<double>(0.0, 0.0));
    complex<double> **TEMP; callocate_2d(&TEMP, NS, NS, complex<double>(0.0, 0.0));
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

void DKD::compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double Z, double KN, int NS)
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
            A[i*LDA+j]=DMAT[i][j];
        }
    }
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//determine the optimal LWORK, i.e., WORK[0]
    LWORK=m_min(LWMAX, int(WORK[0].real())); LDA=NS;
    INFO=LAPACKE_zgeev_work(LAPACK_ROW_MAJOR, JOBVL, JOBVR, NS, A, LDA, W, VL, NS, CG, NS, WORK, LWORK, RWORK);//call the eigenvalue solver
    if(0!=INFO){
        printf("[ERROR] Unable to return INFO as zero when calling the eigenvalue solver.");
        exit(1);
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
    MILWORK=int(GETMILWORK[0].real());
    lapack_complex_double *MIWORK; callocate(&MIWORK, MILWORK, lapack_complex_double(0.0, 0.0));
    LDA=NS;
    INFO=LAPACKE_zgetri_work(LAPACK_ROW_MAJOR, NS, CGINV, LDA, IPIV, MIWORK, MILWORK);
    deallocate(IPIV); deallocate(MIWORK);

    complex<double> *NW; callocate(&NW, NS, complex<double>(0.0, 0.0));
    complex<double> F(2.0*KN, 0.0);
    for(int i=0;i<NS;i++){
        NW[i]=W[i];
        NW[i]=NW[i]/F;
    }
    complex<double> **NCGINV; callocate_2d(&NCGINV, NS, NS, complex<double>(0.0, 0.0));
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCGINV[i][j]=CGINV[i*NS+j];
        }
    }
    complex<double> **NCG; callocate_2d(&NCG, NS, NS, complex<double>(0.0, 0.0));
    for(int i=0;i<NS;i++){
        for(int j=0;j<NS;j++){
            NCG[i][j]=CG[i*NS+j];
        }
    }
    complex<double> **IJK; callocate_2d(&IJK, NS, NS, complex<double>(0.0, 0.0));
    double TPI=TWO_PI*Z;
    for(int j=0;j<NS;j++){
        for(int k=0;k<NS;k++){
            complex<double> sumq(0.0, 0.0);
            complex<double> q(TPI*(NW[j].imag()+NW[k].imag()), TPI*(NW[j].real()-NW[k].real()));
            if(q.real()<0.0) q=-q;
            IJK[j][k]=conj(NCGINV[j][0])*exp(-q)*NCGINV[k][0];
        }
    }

    complex<double> **NCGCONJ; callocate_2d(&NCGCONJ, NS, NS, complex<double>(0.0, 0.0));
    complex<double> **NCGT; callocate_2d(&NCGT, NS, NS, complex<double>(0.0, 0.0));
    complex<double> **TEMP; callocate_2d(&TEMP, NS, NS, complex<double>(0.0, 0.0));
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

void DKD::compute_kvector_projection(CELL *cell, DKD_KVECTOR *kvec, char *projection, bool use_hexagonal)
{
    callocate_2d(&screenNI, numpy, numpx, 0.0);
    callocate_2d(&screenSI, numpy, numpx, 0.0);
    DKD_KNODE *ktemp=kvec->khead;
    int imp=numpx/2;
    for(int i=0;i<kvec->numk;i++){
        int iequiv[48][3], nequiv;
        cell->apply_point_group_symmetry(iequiv, nequiv, ktemp->i, ktemp->j, ktemp->hemisphere, imp, imp);
        for(int i=0;i<nequiv;i++){
            int ix=iequiv[i][0]+imp, iy=iequiv[i][1]+imp;
            if(-1==iequiv[i][2]){
                screenSI[ix][iy]=ktemp->intensity;
            }else if(1==iequiv[i][2]){
                screenNI[ix][iy]=ktemp->intensity;
            }
        }
        ktemp=ktemp->next;
    }
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            if(intensity_maxN<screenNI[i][j]) intensity_maxN=screenNI[i][j];
            if(intensity_minN>screenNI[i][j]) intensity_minN=screenNI[i][j];
        }
    }
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            if(intensity_maxS<screenSI[i][j]) intensity_maxS=screenSI[i][j];
            if(intensity_minS>screenSI[i][j]) intensity_minS=screenSI[i][j];
        }
    }
    img("./test.N.png", screenNI, intensity_maxN, intensity_minN, 'b');
    img("./test.S.png", screenSI, intensity_maxS, intensity_minS, 'b');

    double **matN, **matS;
    callocate_2d(&matN, numpy, numpx, 0.0);
    callocate_2d(&matS, numpy, numpx, 0.0);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            matN[i][j]=screenNI[i][j];
            matS[i][j]=screenSI[i][j];
        }
    }
    for(int i=-imp;i<=imp;i++){
        for(int j=-imp;j<=imp;j++){
            double xy[2]={double(i)/double(imp), double(j)/double(imp)}; 
            double xyz[3]; int ierr;
            compute_sphere_from_square_Lambert(xyz, ierr, xy);
            if(use_hexagonal){
                compute_hexagonal_Lambert(xy, ierr, xyz);
            }else{
                compute_square_Lambert(xy, ierr, xyz);
            }
            if(0!=ierr){
                printf("[ERROR] Unable to compute Lambert interpolation using (%.2f, %.2f, %.2f).\n", xyz[0], xyz[1], xyz[2]);
                exit(1);
            }
            xy[0]*=double(imp); xy[1]*=double(imp);
            int ix=floor(xy[0]), iy=floor(xy[1]);
            int ixp=ix+1, iyp=iy+1;
            if(ixp>imp) ixp=ix;
            if(iyp>imp) iyp=iy;
            double dx=xy[0]-ix, dy=xy[1]-iy;
            double dxm=1.0-dx, dym=1.0-dy;
            ixp+=imp; iyp+=imp; ix+=imp; iy+=imp;
            int idx=i+imp, idy=j+imp;
            screenNI[idx][idy]=matN[ix][iy]*dxm*dym+matN[ixp][iy]*dx*dym+
                                matN[ix][iyp]*dxm*dy+matN[ixp][iyp]*dx*dy;
            screenSI[idx][idy]=matS[ix][iy]*dxm*dym+matS[ixp][iy]*dx*dym+
                                matS[ix][iyp]*dxm*dy+matS[ixp][iyp]*dx*dy;
        }
    }
    int num=numpx-1;
    for(int i=0;i<numpx;i++){
        screenSI[i][0]=screenNI[i][0];
        screenSI[i][num]=screenNI[i][num];
        screenSI[0][i]=screenNI[0][i];
        screenSI[num][i]=screenNI[num][i];
    }
    deallocate_2d(matN, numpy);
    deallocate_2d(matS, numpy);
    intensity_maxN=0.0; intensity_minN=1.0e8;
    intensity_maxS=0.0; intensity_minS=1.0e8;
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            if(intensity_maxN<screenNI[i][j]) intensity_maxN=screenNI[i][j];
            if(intensity_minN>screenNI[i][j]) intensity_minN=screenNI[i][j];
        }
    }
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            if(intensity_maxS<screenSI[i][j]) intensity_maxS=screenSI[i][j];
            if(intensity_minS>screenSI[i][j]) intensity_minS=screenSI[i][j];
        }
    }

    if(0==strcmp(projection, "stereo")){
        callocate_2d(&matN, numpy, numpx, 0.0);
        callocate_2d(&matS, numpy, numpx, 0.0);
        for(int i=0;i<numpy;i++){
            for(int j=0;j<numpx;j++){
                matN[i][j]=screenNI[i][j];
                matS[i][j]=screenSI[i][j];
            }
        }
        intensity_maxN=0.0; intensity_minN=1.0e8;
        intensity_maxS=0.0; intensity_minS=1.0e8;
        for(int i=-imp;i<=imp;i++){
            for(int j=-imp;j<=imp;j++){
                double xy[2]={double(i)/double(imp), double(j)/double(imp)};
                double xyz[3]; int ierr;
                compute_sphere_from_stereographic_projection(xyz, ierr, xy);
                vector_normalize(xyz, xyz);
                int ii=i+imp, jj=j+imp;
                if(0!=ierr){
                    screenNI[ii][jj]=0.0; screenSI[ii][jj]=0.0;
                }else{
                    compute_square_Lambert(xy, ierr, xyz);
                    if(0!=ierr){
                        printf("[ERROR] Unable to compute Lambert interpolation using (%.2f, %.2f, %.2f).\n", xyz[0], xyz[1], xyz[2]);
                        exit(1);
                    }
                    xy[0]*=imp; xy[1]*=imp;
                    int ix=int(imp+xy[0])-imp, iy=int(imp+xy[1])-imp;
                    int ixp=ix+1, iyp=iy+1;
                    if(ixp>imp) ixp=ix;
                    if(iyp>imp) iyp=iy;
                    if(ix<-imp) ix=ixp;
                    if(iy<-imp) iy=iyp;
                    double dx=xy[0]-ix, dy=xy[1]-iy; 
                    double dxm=1.0-dx, dym=1.0-dy;
                    ixp+=imp; iyp+=imp; ix+=imp; iy+=imp;
                    screenNI[ii][jj]=matN[ix][iy]*dxm*dym+matN[ixp][iy]*dx*dym+
                                matN[ix][iyp]*dxm*dy+matN[ixp][iyp]*dx*dy;
                    screenSI[ii][jj]=matS[ix][iy]*dxm*dym+matS[ixp][iy]*dx*dym+
                                matS[ix][iyp]*dxm*dy+matS[ixp][iyp]*dx*dy;
                    if(intensity_maxN<screenNI[ii][jj]) intensity_maxN=screenNI[ii][jj];
                    if(intensity_minN>screenNI[ii][jj]) intensity_minN=screenNI[ii][jj];
                    if(intensity_maxS<screenSI[ii][jj]) intensity_maxS=screenSI[ii][jj];
                    if(intensity_minS>screenSI[ii][jj]) intensity_minS=screenSI[ii][jj];
                }
            }
        }
    }
}

DKD::~DKD()
{
    if(numpy!=0){
        deallocate_2d(screenNI, numpy);
        deallocate_2d(screenSI, numpy);
    }
}

void DKD::dkd(char* dkd_path, char background)
{
    FILE *fp=nullptr;
    fp=fopen(dkd_path,"w");
    fprintf(fp, "KIKUCHI_IMAGE_SIZE\n");
    fprintf(fp, "%d %d\n", numpx, numpy);
    fprintf(fp, "KIKUCHI_IMAGE_NVALUE\n");
    fflush(fp);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            fprintf(fp, "%.8f\n", screenNI[i][j]);
            fflush(fp);
        }
    }
    fprintf(fp, "KIKUCHI_IMAGE_SVALUE\n");
    fflush(fp);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            fprintf(fp, "%.8f\n", screenSI[i][j]);
            fflush(fp);
        }
    }
    printf("[INFO] Information for Kikuchi pattern stored in %s\n", dkd_path);
    char png_path[strlen(dkd_path)+5];
    strcpy(png_path, dkd_path); strcat(png_path, ".N.png");
    img(png_path, screenNI, intensity_maxN, intensity_minN, background);
    printf("[INFO] Image for Kikuchi pattern stored in %s\n", png_path);
    strcpy(png_path, dkd_path); strcat(png_path, ".S.png");
    img(png_path, screenSI, intensity_maxS, intensity_minS, background);
    printf("[INFO] Image for Kikuchi pattern stored in %s\n", png_path);
}

void DKD::dkd(char* dkd_path, double vmax, double vmin, char background)
{
    FILE *fp=nullptr;
    fp=fopen(dkd_path,"w");
    fprintf(fp, "KIKUCHI_IMAGE_SIZE\n");
    fprintf(fp, "%d %d\n", numpx, numpy);
    fprintf(fp, "KIKUCHI_IMAGE_NVALUE\n");
    fflush(fp);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            fprintf(fp, "%.8f\n", screenNI[i][j]);
            fflush(fp);
        }
    }
    fprintf(fp, "KIKUCHI_IMAGE_SVALUE\n");
    fflush(fp);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            fprintf(fp, "%.8f\n", screenSI[i][j]);
            fflush(fp);
        }
    }
    printf("[INFO] Information for Kikuchi pattern stored in %s\n", dkd_path);
    char png_path[strlen(dkd_path)+5];
    strcpy(png_path, dkd_path); strcat(png_path, ".N.png");
    img(png_path, screenNI, vmax, vmin, background);
    printf("[INFO] Image for Kikuchi pattern stored in %s\n", png_path);
    strcpy(png_path, dkd_path); strcat(png_path, ".S.png");
    img(png_path, screenSI, vmax, vmin, background);
    printf("[INFO] Image for Kikuchi pattern stored in %s\n", png_path);
}

void DKD::img(char* png_path, double **value, double vmax, double vmin, char background)
{
    double *wdata;
    unreshape_2d(&wdata, value, numpy, numpx);
    unsigned char *pixels;
    int num=numpx*numpy;
    mallocate(&pixels, 3*num);

    double vdiff=vmax-vmin;
    unsigned char rgb_min[3]={0}, rgb_max[3]={255};
    switch(background)
    {
    case 'b':
        rgb_min[0]=rgb_min[1]=rgb_min[2]=0;
        rgb_max[0]=rgb_max[1]=rgb_max[2]=255;
        break;
    case 'w':
        rgb_min[0]=rgb_min[1]=rgb_min[2]=255;
        rgb_max[0]=rgb_max[1]=rgb_max[2]=0;
        break;
    default:
        printf("[ERROR] Unrecognized background %s\n", background);
        exit(1);
    }
    unsigned char rgb_diff[3]; vector_difference(rgb_diff, rgb_max, rgb_min);
    unsigned char rgb[3];
    for(int i=0;i<num;i++){
        if(wdata[i]>=vmax){
            vector_copy(rgb, rgb_max);
        }else if(wdata[i]<=vmin){
            vector_copy(rgb, rgb_min);
        }else{
            rgb[0]=rgb_min[0]+int((wdata[i]-vmin)/vdiff*rgb_diff[0]);
            rgb[1]=rgb_min[1]+int((wdata[i]-vmin)/vdiff*rgb_diff[1]);
            rgb[2]=rgb_min[2]+int((wdata[i]-vmin)/vdiff*rgb_diff[2]);
        }
        pixels[i*3]=rgb[0]; pixels[i*3+1]=rgb[1]; pixels[i*3+2]=rgb[2];
    }
    image_pixels(png_path, pixels, numpx, numpy);
}