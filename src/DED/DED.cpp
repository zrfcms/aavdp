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

DED::DED(CELL *cell, BETHE *bethe, int zone[3], int fnorm[3], double fthick, double voltage, double dmin, double threshold)
{
    printf("[INFO] Starting computation of dynamical electron diffraction...\n");
    cell->compute_reflection_range(dmin);
    cell->compute_Fourier_coefficients(voltage);
    printf("[INFO] Range of Miller index h along a* in reciprocal space: %d %d\n", 0, cell->HKL[0]);
    printf("[INFO] Range of Miller index k along b* in reciprocal space: %d %d\n", 0, cell->HKL[1]);
    printf("[INFO] Range of Miller index l along c* in reciprocal space: %d %d\n", 0, cell->HKL[2]);

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
    compute_diffraction_intensity(intens, dmat, fthick, kn, gvec.nstrong);
    // DED_GNODE *gtemp=gvec.heads;
    // for(int i=0;i<gvec.nstrong&&gtemp!=nullptr;i++){
    //     if(intens[i]>DED_INTENSITY_LIMIT){
    //         add_k_node(cell, gtemp->hkl, intens[i]);
    //     }
    //     gtemp=gtemp->nexts;
    // }
    DED_GNODE *gtemp=gvec.ghead;
    for(int i=0;i<gvec.numg&&gtemp!=nullptr;i++){
        cell->update_Fourier_coefficient(voltage, gtemp->hkl);
        double intensity=abs(cell->fouri.Vg*cell->fouri.Vg);
        if(intensity>DED_INTENSITY_LIMIT){
            add_k_node(cell, gtemp->hkl, intensity);
        }
        gtemp=gtemp->next;
    }
    double intensity_threshold=intensity_max*threshold;
    filter(khead, intensity_threshold);
    quick_sort(khead, ktail);
    // if(numk<3){
    //     printf("[ERROR] Not enough diffraction vectors to construct the diffraction pattern\n");
    //     exit(1);
    // }
    knearest_1=khead->next;
    knearest_2=knearest_1->next;
    while(knearest_2!=nullptr){
        if(knearest_2->Kmagnitude-knearest_1->Kmagnitude>KMAGNITUDE_LIMIT&&vector_angle(knearest_2->K, knearest_1->K)<HALF_PI&&vector_angle(knearest_2->K, knearest_1->K)>ZERO_LIMIT) break;
        knearest_2=knearest_2->next;
    }
    printf("[INFO] Intensity at the transmission spot: %.8f\n", khead->intensity);
    printf("[INFO] Number of diffraction vectors: %d\n", numk);
    printf("[INFO] Range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
    if(knearest_2==nullptr){
        printf("[WARNING] Unable to find the second nearest diffraction vector\n");
    }else{
        printf("[INFO] The first and second nearest diffraction vectors R1, R2: [%d %d %d], [%d %d %d] (R2/R1 %.8f and angle %.8f)\n", 
        int(knearest_1->hkl[0]), int(knearest_1->hkl[1]), int(knearest_1->hkl[2]), int(knearest_2->hkl[0]), int(knearest_2->hkl[1]), int(knearest_2->hkl[2]), knearest_2->Kmagnitude/knearest_1->Kmagnitude, vector_angle(knearest_2->K, knearest_1->K)*RAD_TO_DEG);
        uvw[2][0]=double(zone[0]); uvw[2][1]=double(zone[1]); uvw[2][2]=double(zone[2]);
        vector_copy(uvw[0], knearest_1->hkl);
        cell->cross(uvw[1], uvw[0], uvw[2], 'r');
        cell->reciprocal_to_cartesian(axes[0], uvw[0]);
        vector_normalize(axes[0], axes[0]);
        cell->reciprocal_to_cartesian(axes[1], uvw[1]);
        vector_normalize(axes[1], axes[1]);
        cell->reciprocal_to_cartesian(axes[2], uvw[2]);
        vector_normalize(axes[2], axes[2]);
        vector_cross(axes[1], axes[0], axes[2]);
        for(int i=0;i<3;i++){
            printf("%.5f %.5f %.5f\n", uvw[i][0], uvw[i][1], uvw[i][2]);
            printf("%.5f %.5f %.5f\n", axes[i][0], axes[i][1], axes[i][2]);
        }
    }
    printf("[INFO] Ending computation of dynamical electron diffraction\n");
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

void DED::filter(DED_KNODE *kstart, double threshold)
{
    intensity_min=1.0e8;
    while(kstart!=nullptr&&kstart->intensity<=threshold){
        DED_KNODE *ktemp=kstart;
        kstart=kstart->next;
        delete ktemp;
        numk--;
    }
    DED_KNODE *cur=kstart;
    while(cur!=nullptr&&cur->next!=nullptr){
        if(cur->next->intensity<=threshold){
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

void DED::ded(char *ded_path)
{
    FILE *fp=nullptr;
    fp=fopen(ded_path,"w");
    fprintf(fp, "# h\tk\tl\tx\ty\tintensity\tintensity_norm (%d points, x-[%d %d %d], y-[%d %d %d], and z-[%d %d %d])\n", 
            numk, int(uvw[0][0]), int(uvw[0][1]), int(uvw[0][2]), int(uvw[1][0]), int(uvw[1][1]), int(uvw[2][2]), int(uvw[2][0]), int(uvw[2][1]), int(uvw[2][2]));
    double *pos_x=nullptr, *pos_y=nullptr, *intensity=nullptr;
    callocate(&pos_x, numk, 0.0); 
    callocate(&pos_y, numk, 0.0); 
    callocate(&intensity, numk, 0.0);
    DED_KNODE *ktemp=khead;
    fprintf(fp, "%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\n", 0, 0, 0, 0.0, 0.0, ktemp->intensity, 100.0);
    fflush(fp);
    pos_x[0]=pos_y[0]=0.0; intensity[0]=100.0;
    ktemp=ktemp->next;

    double constn=100.0/intensity_max;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        double xyz[3]; vector_rotate(xyz, axes, ktemp->K);
        double intensity_norm=constn*ktemp->intensity;
        fprintf(fp, "%.8f\t%.8f\t%.8f\t%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", ktemp->K[0], ktemp->K[1], ktemp->K[2], int(ktemp->hkl[0]), int(ktemp->hkl[1]), int(ktemp->hkl[2]), xyz[0], xyz[1], xyz[2], ktemp->intensity, intensity_norm);
        fflush(fp);
        pos_x[i]=xyz[0]; pos_y[i]=xyz[1]; intensity[i]=intensity_norm;
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s.\n", ded_path);

    char png_path[PATH_CHAR_NUMBER];
    strcpy(png_path, ded_path); strcat(png_path, ".png");
    img(png_path, pos_x, pos_y, intensity, numk);
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void DED::rotate(CELL *cell, int x[3], int y[3])
{
    uvw[0][0]=double(x[0]); uvw[0][1]=double(x[1]); uvw[0][2]=double(x[2]);
    uvw[1][0]=double(y[0]); uvw[1][1]=double(y[1]); uvw[1][2]=double(y[2]);
    cell->reciprocal_to_cartesian(axes[0], uvw[0]);
    vector_normalize(axes[0], axes[0]);
    cell->reciprocal_to_cartesian(axes[1], uvw[1]);
    vector_normalize(axes[1], axes[1]);
    if(fabs(vector_dot(axes[0], axes[1]))>1.0e-6||fabs(vector_dot(axes[1], axes[2]))>1.0e-6||fabs(vector_dot(axes[0], axes[2]))>1.0e-6){
        printf("[ERROR] The orthogonality condition is not satisfied with x-[%d %d %d], y-[%d %d %d]", x[0], x[1], x[2], y[0], y[1], y[2]);
        exit(1);
    }
}

void DED::img(char *png_path, double *x, double *y, double *value, int num)
{
    double height=6.0, width=6.0;
    int tick_max=int(ceil(ktail->Kmagnitude*10.0));
    int n_major_tick=2*tick_max+1;
    double Kmagnitude_max=double(tick_max)*0.1;
    double *major_ticks; mallocate(&major_ticks, n_major_tick);
    for(int i=0;i<n_major_tick;i++){
        major_ticks[i]=double(-tick_max+i);
    }
    GRAPH graph(width, height, 300);
    graph.set_xlim(-Kmagnitude_max, Kmagnitude_max);
    graph.set_ylim(-Kmagnitude_max, Kmagnitude_max);
    graph.set_xticks(major_ticks, n_major_tick);
    graph.set_yticks(major_ticks, n_major_tick);
    graph.set_tick_in(false);
    graph.scatter(x, y, value, num);
    graph.draw(png_path);
}

void DED::copy_knode_data(DED_KNODE *knode1, DED_KNODE *knode2)
{
    vector_copy(knode1->hkl, knode2->hkl);
    vector_copy(knode1->K, knode2->K);
    knode1->Kmagnitude=knode2->Kmagnitude;
    knode1->intensity=knode2->intensity;
}

void DED::swap_knode_data(DED_KNODE *knode1, DED_KNODE *knode2)
{
    DED_KNODE *ktemp=new DED_KNODE;
    copy_knode_data(ktemp, knode1);
    copy_knode_data(knode1, knode2);
    copy_knode_data(knode2, ktemp);
    delete ktemp;
}

void DED::quick_sort(DED_KNODE *kstart, DED_KNODE *kend)
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

DKD_KVECTOR::DKD_KVECTOR(CELL *cell, int impx, int impy)
{
    deltax=1.0/double(impx); deltay=1.0/double(impy);
    double kn=1.0/cell->fouri0.lambda;
    add_k_node(cell, 0, 0, kn);
    switch(cell->sampling_type)
    {
    case 1://triclinic 1
        for(int j=-impy;j<=impy;j++){
            for(int i=-impx;i<=impx;i++){
                add_k_node(cell, i, j, kn, true);
            }
        }
        break;
    case 2://triclinic -1
        for(int j=-impy;j<=impy;j++){
            for(int i=-impx;i<=impx;i++){
                add_k_node(cell, i, j, kn);
            }
        }
        break;
    case 3://monoclinic 2
        for(int j=-impy;j<=impy;j++){
            for(int i=0;i<=impx;i++){
                add_k_node(cell, i, j, kn, true);
            }
        }
        break;
    case 4://monoclinic m
        for(int j=0;j<=impy;j++){
            for(int i=-impx;i<=impx;i++){
                add_k_node(cell, i, j, kn, true);
            }
        }
        break;
    case 5://monoclinic 2/m, orthorhombic 222, mm2, tetragonal 4, -4
        for(int j=0;j<=impy;j++){
            for(int i=0;i<=impx;i++){
                add_k_node(cell, i, j, kn, true);
            }
        }
        break;
    case 6://orthorhombic mmm, tetragonal 4/m, 422, -4m2, cubic m-3, 432
        for(int j=0;j<=impy;j++){
            for(int i=0;i<=impx;i++){
                add_k_node(cell, i, j, kn);
            }
        }
        break;
    case 7://tetragonal 4mm
        for(int i=0;i<=impx;i++){
            for(int j=0;j<=i;j++){
                add_k_node(cell, i, j, kn, true);
            }
        }
        break;
    case 8://tetragonal -42m, cubic -43m
        for(int i=0;i<=impx;i++){
            for(int j=-i;j<=i;j++){
                add_k_node(cell, i, j, kn);
            }
        }
        break;
    case 9://tetragonal 4/mmm, cubic m-3m
        for(int i=0;i<=impx;i++){
            for(int j=0;j<=i;j++){
                add_k_node(cell, i, j, kn);
            }
        }
        break;
    case 10://hexagonal 3
        for(int j=0;j<=impx;j++){
            for(int i=0;i<=impx;i++){
                add_k_node(cell, i, j, kn, true);
            }
        }
        break;
    case 11://rhombohedral 3
        printf("[ERROR] Unrecognized sampling type %d in k-vector computation.", cell->sampling_type);
        exit(1);
    case 12://hexagonal -3, 321, -6 [not implemented: rhombohedral 32]
        for(int j=0;j<=impx;j++){
            for(int i=0;i<=impx;i++){
                add_k_node(cell, i, j, kn);
            }
        }
        break;
    case 13://hexagonal 312 [not implemented: rhombohedral -3]
        for(int j=0;j>=-impx;j--){
            for(int i=j/2;i<impx;i++){
                add_k_node(cell, i, j, kn);
            }
        }
        for(int i=0;i<=impx;i++){
            for(int j=0;j<i/2;j++){
                add_k_node(cell, i, j, kn);
            }
        }
        break;
    case 14://hexagonal 3m1 [not implemented: rhombohedral 3m]
        for(int j=1;j<=impx;j++){
            for(int i=j;i<=impx;i++){
                add_k_node(cell, i, j, kn, true);
            }
        }
        break;
    case 15://hexagonal 31m, 6
        for(int j=1;j<=impx;j++){
            for(int i=j;i<=impx;i++){
                add_k_node(cell, i, j, kn, true);
            }
        }
        break;
    case 16://hexagonal -3m1, 622, -6m2 [not implemented: rhombohedral -3m]
        for(int j=0;j<=impx;j++){
            for(int i=0;i<=impx;i++){
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)<(PI/6.0-1.0e-4)))){
                    continue;
                }else{
                    add_k_node(cell, i, j, kn);
                }
            }
        }
        break;
    case 17://hexagonal -31m, 6/m, -62m
        for(int j=0;j<=impx;j++){
            for(int i=0;i<=impx;i++){
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)>(PI/3.0+1.0e-4)))){
                    continue;
                }else{
                    add_k_node(cell, i, j, kn);
                }
            }
        }
        break;
    case 18://hexagonal 6mm
        for(int j=0;j<=impx;j++){
            for(int i=0;i<=impx;i++){
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)>(PI/6.0+1.0e-4)))){
                    continue;
                }else{
                    add_k_node(cell, i, j, kn, true);
                }
            }
        }
        break;
    case 19:
        for(int j=0;j<=impx;j++){
            for(int i=0;i<=impx;i++){
                double x=double(i)-double(j)/2.0;
                double y=double(j)*SQRT_HALF_3;
                if((x<0.0)||((x>=0.0)&&(atan2(y, x)>(PI/6.0+1.0e-4)))){
                    continue;
                }else{
                    add_k_node(cell, i, j, kn);
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

void DKD_KVECTOR::add_k_node(CELL *cell, int i, int j, double kn, bool southern_flag)
{
    if(khead==nullptr&&ktail==nullptr){
        khead=ktail=new DKD_KNODE;
        double kstar[3]={0.0, 0.0, kn}, r_kstar[3];//c* as the center of the Rosca-Lambert projection
        cell->cartesian_to_reciprocal(r_kstar, kstar);
        vector_copy(khead->k, r_kstar);
        khead->kn=kn;
        khead->i=i; khead->j=j;
        khead->hemisphere=1;//Northern hemisphere
        numk++;
    }else{
        double xy[2]={double(i)*deltax, double(j)*deltay};
        if(cell->sampling_type>=10){//case 10-19 or use_hexagonal ?
            double x=fabs(xy[0]-0.5*xy[1]), y=fabs(SQRT_HALF_3*xy[1]);
            if(x>1.0||y>SQRT_HALF_3) return;
            if(x+y*SQRT_3_INVERSE>1.0) return;
        }
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
        vector_constant(kstar, kn, kstar);
        cell->cartesian_to_reciprocal(r_kstar, kstar);
        vector_copy(ktail->k, r_kstar);
        ktail->kn=kn;
        ktail->i=i; ktail->j=j;
        ktail->hemisphere=1;
        numk++;
        if(southern_flag){
            vector_constant(kstar, -1.0, kstar);
            cell->cartesian_to_reciprocal(r_kstar, kstar);
            ktail->next=new DKD_KNODE;
            ktail=ktail->next;
            vector_copy(ktail->k, r_kstar);
            ktail->kn=kn;
            ktail->i=i; ktail->j=j;
            ktail->hemisphere=-1;
            numk++;
        }
    }
}

void DKD_KVECTOR::add_k_intensity(DKD_KNODE *knode, complex<double> **Lgh, complex<double> ***Sgh, int nstrong, int napos, int npos)
{
    knode->intensity=0.0;
    for(int i=0;i<napos;i++){
        for(int j=0;j<nstrong;j++){
            for(int k=0;k<nstrong;k++){
                complex<double> temp=Lgh[j][k]*Sgh[i][j][k];
                knode->intensity+=temp.real();
            }
        }
    }
    knode->intensity/=double(npos);
    if(intensity_max<knode->intensity) intensity_max=knode->intensity;
    if(intensity_min>knode->intensity) intensity_min=knode->intensity;
}

DKD::DKD(CELL *cell, DED_MC *mc, BETHE *bethe, double voltage, double dmin, int nump)
{
    printf("[INFO] Starting computation of dynamical electron diffraction...\n");
    cell->compute_reflection_range(dmin);
    printf("[INFO] Range of Miller index h along a* in reciprocal space: %d %d\n", 0, cell->HKL[0]);
    printf("[INFO] Range of Miller index k along b* in reciprocal space: %d %d\n", 0, cell->HKL[1]);
    printf("[INFO] Range of Miller index l along c* in reciprocal space: %d %d\n", 0, cell->HKL[2]);
    printf("[INFO] Generating Fourier coefficient lookup table ...");
    cell->compute_Fourier_coefficients(voltage);
    printf("[INFO] Done!\n");
    printf("[INFO] Generating Bloch wave coefficient lookup table ...");
    cell->compute_Bloch_wave_coefficients();
    printf("[INFO] Done!\n");
    
    numpx=numpy=nump;
    impx=impy=nump/2;
    DKD_KVECTOR kvec(cell, impx, impy);
    printf("[INFO] Independent beam directions to be considered = %d\n", kvec.numk);
    
    int sum_strong=0, sum_weak=0;
    clock_t start, finish; start=clock();
    DKD_KNODE *ktemp=kvec.khead;
    for(int ik=0;ik<kvec.numk&&ktemp!=nullptr;ik++){
        double *kk=ktemp->k, *fn=kk;
        double kn=ktemp->kn;
        DED_GVECTOR gvec(cell, bethe, kk, fn, kn, dmin);
        sum_strong+=gvec.nstrong; sum_weak+=gvec.nweak;

        int nstrong=gvec.nstrong;
        complex<double> ***Sgh;
        callocate_3d(&Sgh, cell->napos, nstrong, nstrong, complex<double>(0.0, 0.0));
        compute_Sgh_matrices(Sgh, cell, &gvec);
        complex<double> **dmat, **Lgh;
        callocate_2d(&dmat, nstrong, nstrong, complex<double>(0.0, 0.0));
        callocate_2d(&Lgh, nstrong, nstrong, complex<double>(0.0, 0.0));
        compute_dynamic_matrix(dmat, cell, &gvec);
        compute_Lgh_matrix(Lgh, dmat, mc->weight, mc->izmax, mc->z, mc->dz, kn, nstrong);
        kvec.add_k_intensity(ktemp, Lgh, Sgh, nstrong, cell->napos, cell->npos);
        deallocate_2d(Lgh, nstrong);
        deallocate_2d(dmat, nstrong);
        deallocate_3d(Sgh, cell->napos, nstrong);
        ktemp=ktemp->next;
        if(0==(ik+1)%1000) printf("[INFO] Completed beam direction %d of %d.\n", ik+1, kvec.numk);
    }
    finish=clock();
    printf("[INFO] Average number of strong reflections = %d.\n", int(round(double(sum_strong)/double(kvec.numk))));
    printf("[INFO] Average number of weak reflections = %d.\n", int(round(double(sum_weak)/double(kvec.numk))));
    printf("[INFO] Execution time [s]: %.2f.\n", double(finish-start)/CLOCKS_PER_SEC);
    compute_Lambert_projection(cell, &kvec);
    compute_stereographic_projection(cell->use_hexagonal);
}

void DKD::compute_dynamic_matrix(complex<double> **dmat, CELL *cell, DED_GVECTOR *gvec)
{
    DED_GNODE *rtemp=gvec->ghead->next;
    int imh=2*cell->HKL[0], imk=2*cell->HKL[1], iml=2*cell->HKL[2];
    double k0_2=2.0/cell->fouri0.lambda, k0_2i=0.5*cell->fouri0.lambda;
    for(int ir=0;ir<gvec->nstrong;ir++){
        DED_GNODE *ctemp=gvec->ghead->next;
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
    DED_GNODE *rtemp=gvec->ghead->next;
    for(int ir=0;ir<gvec->nstrong;ir++){
        DED_GNODE *ctemp=gvec->ghead->next;
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

void DKD::compute_Lambert_projection(CELL *cell, DKD_KVECTOR *kvec)
{
    callocate_2d(&mLPNH, numpx, numpy, 0.0);
    callocate_2d(&mLPSH, numpx, numpy, 0.0);
    DKD_KNODE *ktemp=kvec->khead;
    for(int ik=0;ik<kvec->numk&&ktemp!=nullptr;ik++){
        int iequiv[48][3], nequiv;
        cell->apply_point_group_symmetry(iequiv, nequiv, ktemp->i, ktemp->j, ktemp->hemisphere, impx, impy);
        for(int i=0;i<nequiv;i++){
            int ix=iequiv[i][0]+impx, iy=iequiv[i][1]+impy;
            if(-1==iequiv[i][2]){
                mLPSH[ix][iy]=ktemp->intensity;
            }else if(1==iequiv[i][2]){
                mLPNH[ix][iy]=ktemp->intensity;
            }
        }
    }

    double **matN, **matS;
    callocate_2d(&matN, numpx, numpy, 0.0);
    callocate_2d(&matS, numpx, numpy, 0.0);
    for(int i=0;i<numpx;i++){
        for(int j=0;j<numpy;j++){
            matN[i][j]=mLPNH[i][j];
            matS[i][j]=mLPSH[i][j];
        }
    }
    for(int i=-impx;i<=impx;i++){
        for(int j=-impy;j<=impy;j++){
            double xy[2]={double(i)/double(impx), double(j)/double(impy)}; 
            double xyz[3]; int ierr;
            compute_sphere_from_square_Lambert(xyz, ierr, xy);
            if(cell->use_hexagonal){
                compute_hexagonal_Lambert(xy, ierr, xyz);
            }else{
                compute_square_Lambert(xy, ierr, xyz);
            }
            if(0!=ierr){
                printf("[ERROR] Unable to compute Lambert interpolation using (%.2f, %.2f, %.2f).\n", xyz[0], xyz[1], xyz[2]);
                exit(1);
            }
            xy[0]*=double(impx); xy[1]*=double(impy);
            int ix=floor(xy[0]), iy=floor(xy[1]);
            int ixp=ix+1, iyp=iy+1;
            if(ixp>impx) ixp=ix;
            if(iyp>impy) iyp=iy;
            double dx=xy[0]-ix, dy=xy[1]-iy;
            double dxm=1.0-dx, dym=1.0-dy;
            ixp+=impx; iyp+=impy; ix+=impx; iy+=impy;
            int idx=i+impx, idy=j+impy;
            mLPNH[idx][idy]=matN[ix][iy]*dxm*dym+matN[ixp][iy]*dx*dym+
                                matN[ix][iyp]*dxm*dy+matN[ixp][iyp]*dx*dy;
            mLPSH[idx][idy]=matS[ix][iy]*dxm*dym+matS[ixp][iy]*dx*dym+
                                matS[ix][iyp]*dxm*dy+matS[ixp][iyp]*dx*dy;
        }
    }
    deallocate_2d(matN, numpx);
    deallocate_2d(matS, numpx);
    for(int i=0;i<numpx;i++){
        mLPSH[i][0]=mLPNH[i][0];
        mLPSH[i][numpy-1]=mLPNH[i][numpy-1];
    }
    for(int i=0;i<numpy;i++){
        mLPSH[0][i]=mLPNH[0][i];
        mLPSH[numpx-1][i]=mLPNH[numpx-1][i];
    }
}

void DKD::compute_stereographic_projection(bool use_hexagonal)
{
    callocate_2d(&mSPNH, numpx, numpy, 0.0);
    callocate_2d(&mSPSH, numpx, numpy, 0.0);
    for(int i=-impx;i<=impx;i++){
        for(int j=-impy;j<=impy;j++){
            double xy[2]={double(i)/double(impx), double(j)/double(impy)};
            double xyz[3]; int ierr;
            compute_sphere_from_stereographic_projection(xyz, ierr, xy);
            vector_normalize(xyz, xyz);
            int ii=i+impx, jj=j+impy;
            if(0!=ierr){
                mSPNH[ii][jj]=0.0; mSPSH[ii][jj]=0.0;
            }else{
                if(use_hexagonal){
                    compute_hexagonal_Lambert(xy, ierr, xyz);
                }else{
                    compute_square_Lambert(xy, ierr, xyz);
                }
                if(0!=ierr){
                    printf("[ERROR] Unable to compute Lambert interpolation using (%.2f, %.2f, %.2f).\n", xyz[0], xyz[1], xyz[2]);
                    exit(1);
                }
                xy[0]*=impx; xy[1]*=impy;
                int ix=int(impx+xy[0])-impx, iy=int(impy+xy[1])-impy;
                int ixp=ix+1, iyp=iy+1;
                if(ixp>impx) ixp=ix;
                if(iyp>impy) iyp=iy;
                if(ix<-impx) ix=ixp;
                if(iy<-impy) iy=iyp;
                double dx=xy[0]-ix, dy=xy[1]-iy; 
                double dxm=1.0-dx, dym=1.0-dy;
                ixp+=impx; iyp+=impy; ix+=impx; iy+=impy;
                mSPNH[ii][jj]=mLPNH[ix][iy]*dxm*dym+mLPNH[ixp][iy]*dx*dym+
                                     mLPNH[ix][iyp]*dxm*dy+mLPNH[ixp][iyp]*dx*dy;
                mSPSH[ii][jj]=mLPSH[ix][iy]*dxm*dym+mLPSH[ixp][iy]*dx*dym+
                                     mLPSH[ix][iyp]*dxm*dy+mLPSH[ixp][iyp]*dx*dy;
            }
        }
    }
}

DKD::~DKD()
{
    if(numpx!=0&&numpy!=0){
        deallocate_2d(mLPNH, numpx);
        deallocate_2d(mLPSH, numpx);
        deallocate_2d(mSPNH, numpx);
        deallocate_2d(mSPSH, numpx);
    }
}

void DKD::img(char *LPNH_path, char *LPSH_path, char *SPNH_path, char *SPSH_path, double dimension, int resolution)
{
    image_array(LPNH_path, mLPNH, numpx, numpy, dimension, dimension, resolution);
    printf("[INFO] Image data for the modified lambert projection of northern hemisphere stored in %s\n", LPNH_path);
    image_array(LPSH_path, mLPSH, numpx, numpy, dimension, dimension, resolution);
    printf("[INFO] Image data for the modified lambert projection of sorthern hemisphere stored in %s\n", LPSH_path);
    image_array(SPNH_path, mSPNH, numpx, numpy, dimension, dimension, resolution);
    printf("[INFO] Image data for the master stereographic projection of northern hemisphere stored in %s\n", SPNH_path);
    image_array(SPSH_path, mSPSH, numpx, numpy, dimension, dimension, resolution);
    printf("[INFO] Image data for the master stereographic projection of northern hemisphere stored in %s\n", SPSH_path);
}