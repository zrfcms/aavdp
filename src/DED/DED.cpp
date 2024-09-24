#include "DED.h"

DED_KVECTOR::DED_KVECTOR(CELL *cell, double kz[3], double g1[3], double precangle)
{
    // double kstar[3], gx[3], gy[3];
    // vector_copy(kstar, kz);
    // vector_copy(gx, g1);
    // cell->direct_to_reciprocal(kstar, kstar);
    // printf("Reciprocal lattice vectors are [%.2f %.2f %.2f]\n", kstar[0], kstar[1], kstar[2]);
    // vector_cross(gy, gx, kstar);
    // cell->direct_to_reciprocal(gy, gy);
    // cell->normalize(kstar, kstar, 'r');
    // double gx_len=cell->length(gx, 'r');
    // cell->normalize(gx, gx, 'r');
    // cell->normalize(gy, gy, 'r');
    // printf("Reciprocal lattice vectors are [%.2f %.2f %.2f] along x-axis and [%.2f %.2f %.2f] along y-axis, with normal [%.2f %.2f %.2f]\n", gx[0], gx[1], gx[2], gy[0], gy[1], gy[2], kstar[0], kstar[1], kstar[2]);

    // double kn=1.0/cell->fouri0.lambda;
    // double rg=2.0*kn*sin(precangle);
    // double drg=2.0*kn*sin(precangle-0.25*1.0e-3)/double(ncircle)/gx_len, dth=TWO_PI/double(pcircle);
    // for(int i=-ncircle;i<=ncircle;i++){
    //     for(int j=0;j<pcircle;j++){
    //         double irg=rg+drg*i;
    //         double kt[3], ktx[3], kty[3];
    //         vector_constant(ktx, -irg*cos(j*dth), gx);
    //         vector_constant(kty, -irg*sin(j*dth), gy);
    //         vector_plus(kt, ktx, kty);
    //         double kr[3];
    //         double dt=cell->dot(kt, kt, 'r');
    //         vector_constant(kr, sqrt(kn*kn-dt*dt), kstar);
    //         vector_plus(kr, kt, kr);
    //         double kn=cell->dot(kr, kstar, 'r');
    //         add_k_node(kr, kt, kn);
    //     }
    // }
    // callocate_2d(&karray, numk, 3, 0.0); callocate(&knarray, numk, 0.0);
    // DED_KNODE *ktemp=khead;
    // FILE *fp=fopen("kvectors.txt", "w");
    // for(int ik=0;ik<numk&&ktemp!=nullptr;ik++){
    //     fprintf(fp, "%.8f\t%.8f\t%.8f\n", ktemp->k[0], ktemp->k[1], ktemp->k[2]);
    //     fflush(fp);
    //     vector_copy(karray[ik], ktemp->k);
    //     knarray[ik]=ktemp->kn;
    //     ktemp=ktemp->next;
    // }
    // free_k_node();
    // fclose(fp);
    numk=1;
    callocate_2d(&karray, 1, 3, 0.0); callocate(&knarray, 1, 0.0);
    double kn=1.0/cell->fouri0.lambda;
    double kstar[3]={0.0, 0.0, kn}, r_kstar[3];
    cell->cartesian_to_reciprocal(r_kstar, kstar);
    vector_copy(karray[0], kstar);
    knarray[0]=kn;
}

DED_KVECTOR::~DED_KVECTOR()
{
    deallocate_2d(karray, numk);
    deallocate(knarray);
}

void DED_KVECTOR::add_k_node(double k[3], double kt[3], double kn)
{
    if(khead==nullptr&&ktail==nullptr){
        khead=ktail=new DED_KNODE;
        vector_copy(ktail->k, k);
        vector_copy(ktail->kt, kt);
        ktail->kn=kn;      
    }else{
        ktail->next=new DED_KNODE;
        ktail=ktail->next;
        vector_copy(ktail->k, k);
        vector_copy(ktail->kt, kt);
        ktail->kn=kn;
    }
    numk++;
}

void DED_KVECTOR::free_k_node()
{
    DED_KNODE *cur=khead;
    while(cur!=nullptr){
        DED_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

DED_GVECTOR::DED_GVECTOR(CELL *cell, double kz[3], double fn[3], double cutoff, double precangle, double goffset)
{
    double kn=1.0/cell->fouri0.lambda;
    double kstar[3];
    cell->direct_to_reciprocal(kstar, kz);
    cell->normalize(kstar, kstar, 'r');
    vector_constant(kstar, kn, kstar);
    printf("kstar %.2f %.2f %.2f\n", kstar[0], kstar[1], kstar[2]);
    double ng=kn*cos(precangle);
    double rg=2.0*kn*sin(precangle);

    int imh=cell->HKL[0], imk=cell->HKL[1], iml=cell->HKL[2];
    int cx=imh*2, cy=imk*2, cz=iml*2;

    double g[3]={0.0, 0.0, 0.0};
    complex<double> Ug=cell->LUTUg[cx][cy][cz];
    complex<double> qg=cell->LUTqg[cx][cy][cz];
    double sg=0.0, sig=0.0;
    bool is_double_diffrac=false;
    add_g_node(g, Ug, qg, sg, sig, is_double_diffrac);
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                if(abs(ih)+abs(ik)+abs(il)!=0){//avoid double counting the origin
                    g[0]=double(ih); g[1]=double(ik); g[2]=double(il);
                    double dhkl=cell->get_interplanar_spacing(g);
                    if(cell->is_centering_allowed(g)&&(dhkl>cutoff)){
                        double dt=cell->dot(g, kstar, 'r');
                        double gp[3], gv[3];
                        vector_constant(gp, dt, kstar);
                        vector_difference(gv, g, gp);
                        double gv_len=cell->length(gv, 'r');
                        double gp_len=cell->length(gp, 'r');
                        if(dt<=0.0) gp_len=-gp_len;
                        double y=gv_len*rg, z=ng*ng-gv_len*gv_len;
                        double upper_bound=goffset+ng-sqrt(z-y);
                        double lower_bound=-goffset+ng-sqrt(z+y);
                        if(gp_len>=lower_bound&&gp_len<=upper_bound){
                            int ix=ih+cx, iy=ik+cy, iz=il+cz;
                            Ug=cell->LUTUg[ix][iy][iz];
                            qg=cell->LUTqg[ix][iy][iz];
                            sg=cell->get_excitation_error(g, kstar, kstar);
                            sig=0.0;
                            is_double_diffrac=cell->is_double_diffrac[ix][iy][iz];
                            if(is_double_diffrac){
                                add_g_node(g, Ug, qg, sg, sig, is_double_diffrac);
                            }else{
                                add_g_node(g, Ug, qg, sg, sig, is_double_diffrac);
                            }
                        }
                    }
                }
            }
        }
    }
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

void DED_GVECTOR::filter_g_node(CELL *cell, DED_BETHE *bethe, double k[3], double fn[3])
{
    double kn=1.0/cell->fouri0.lambda;
    int imh=cell->HKL[0], imk=cell->HKL[1], iml=cell->HKL[2];
    int cx=imh*2, cy=imk*2, cz=iml*2;

    DED_GNODE *gtemp=ghead->next;
    int num=0;
    int **garray; callocate_2d(&garray, numg, 3, 0);
    while(gtemp!=nullptr){
        for(int i=0;i<3;i++){
            garray[num][i]=int(gtemp->hkl[i]);
        }
        num++;
        gtemp=gtemp->next;
    }

    DED_GNODE *tails=nullptr, *tailw=nullptr;
    gtemp=ghead->next;
    gtemp->is_strong=true;
    gtemp->is_weak=false;
    heads=tails=gtemp;
    nstrong++;
    for(int i=1;i<num;i++){
        gtemp=gtemp->next;
        double imin=1.0e8;
        double sg=cell->get_excitation_error(gtemp->hkl, k, fn);
        //double sgp=kn*fabs(sg);
        double sgp=kn*fabs(gtemp->sg);
        for(int j=0;j<num;j++){
            double temp;
            int ix=garray[i][0]-garray[j][0]+cx, iy=garray[i][1]-garray[j][1]+cy, iz=garray[i][2]-garray[j][2]+cz;
            if(cell->is_double_diffrac[ix][iy][iz]){
                temp=1.0e4;
            }else{
                temp=sgp/abs(cell->LUTUg[ix][iy][iz]);
            }
            if(temp<imin) imin=temp;
        }
        //if(abs(int(gtemp->hkl[0]))==1&&abs(int(gtemp->hkl[1]))==1&&abs(int(gtemp->hkl[2]))==1) printf("%.2f\n", imin);
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

void DED_GVECTOR::add_g_node(double hkl[3], complex<double> Ug, complex<double> qg, double sg, double sig, bool is_double_diffrac)
{
    if(ghead==nullptr&&gtail==nullptr){
        ghead=gtail=new DED_GNODE;
    }
    gtail->next=new DED_GNODE;
    gtail=gtail->next;
    gtail->hkl[0]=hkl[0]; gtail->hkl[1]=hkl[1]; gtail->hkl[2]=hkl[2];
    printf("%.2f %.2f %.2f\n", hkl[0], hkl[1], hkl[2]);
    gtail->Ug=Ug;
    gtail->qg=qg;
    gtail->sg=sg;
    gtail->sig=sig;
    gtail->is_double_diffrac=is_double_diffrac;
    numg++;
}

void DED_GVECTOR::unfilter_g_node()
{
    nstrong=nweak=0;
    DED_GNODE *gtemp=ghead;
    for(int i=0;i<numg&&gtemp!=nullptr;i++){
        gtemp->is_strong=gtemp->is_weak=false;
        gtemp->nexts=gtemp->nextw=nullptr;
        gtemp=gtemp->next;
    }
    heads=headw=nullptr;
}

DED::DED(CELL *cell, int zone[3], int fnorm[3], double thickness, double dmin, double voltage, double precangle, double c1, double c2, double c3, double c_sg)
{
    cell->compute_reflection_range(dmin);
    printf("[INFO] Range of reflections along a*, b*, and c* = %d, %d, and %d.\n", cell->HKL[0], cell->HKL[1], cell->HKL[2]);
    cell->compute_Fourier_coefficients(voltage);

    double kz[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    double fn[3]={double(fnorm[0]), double(fnorm[1]), double(fnorm[2])};
    struct DED_BETHE bethe={c1, c2, c3, c_sg};
    double g1[3], g2[3];
    cell->compute_shortest_reciprocal_vectors(g1, g2, kz);
    printf("[INFO] The first shortest reciprocal lattice vector along zone-[%d %d %d]: %.2f %.2f %.2f\n", zone[0], zone[1], zone[2], g1[0], g1[1], g1[2]);
    printf("[INFO] The second shortest reciprocal lattice vector along zone-[%d %d %d]: %.2f %.2f %.2f\n", zone[0], zone[1], zone[2], g2[0], g2[1], g2[2]);
    DED_KVECTOR kvec(cell, kz, g1, precangle);
    DED_GVECTOR gvec(cell, kz, kz, dmin, precangle);
    printf("[INFO] Independent beam directions to be considered = %d\n", kvec.numk);
    printf("[INFO] Total number of reflections = %d.\n", gvec.numg);
    for(int ik=0;ik<kvec.numk;ik++){
        double *kk=kvec.karray[ik];
        double *fn=kk;
        gvec.filter_g_node(cell, &bethe, kk, fn);

        complex<double> **dmat;
        callocate_2d(&dmat, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
        compute_dynamic_matrix(dmat, cell, &gvec);

        complex<double> **Lgh;
        callocate_2d(&Lgh, gvec.nstrong, gvec.nstrong, complex<double>(0.0, 0.0));
        compute_Lgh_matrix(Lgh, dmat, thickness, kvec.knarray[ik], gvec.nstrong);
        printf("[INFO] Total number of strong reflections = %d.\n", gvec.nstrong);

        DED_GNODE *gtemp=gvec.heads;
        for(int i=0;i<gvec.nstrong&&gtemp!=nullptr;i++){
            complex<double> sumq(0.0, 0.0);
            for(int j=0;j<gvec.nstrong;j++){
                sumq+=Lgh[i][j];
            }
            gtemp->sig+=pow(fabs(sumq), 2);
            printf("%.2f %.2f %.2f %.2f\n", gtemp->hkl[0], gtemp->hkl[1], gtemp->hkl[2], gtemp->sig);
            gtemp=gtemp->nexts;
        }
        gvec.unfilter_g_node();
    }
    double DM[2][2]={{vector_dot(g2, g2), -vector_dot(g1, g2)}, {-vector_dot(g1, g2), vector_dot(g1, g1)}};
    double DD=1.0/(DM[0][0]*DM[1][1]-DM[0][1]*DM[1][0]);
    DED_GNODE *gtemp=gvec.ghead;
    for(int i=0;i<gvec.numg&&gtemp!=nullptr;i++){
        if(gtemp->sig!=0.0){
            double x1, x2;
            x1=vector_dot(gtemp->hkl, g1);
            x2=vector_dot(gtemp->hkl, g2);
            double pos[2]={DD*(DM[0][0]*x1+DM[0][1]*x2), DD*(DM[1][0]*x1+DM[1][1]*x2)};
            add_intensity_node(gtemp->hkl, pos, gtemp->sig);
        }
        gtemp=gtemp->next;
    }
}

DED::~DED()
{
    DED_INODE *cur=ihead;
    while(cur!=nullptr){
        DED_INODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
}

void DED::ded(const char *ded_path)
{
    FILE *fp=fopen(ded_path, "w");
    fprintf(fp, "# h\tk\tl\tx\ty\tintensity\tintensity_norm (%d points)\n", numi);
    double consti=100.0/intensity_max;
    DED_INODE *itemp=ihead;
    for(int i=0;i<numi&&itemp!=nullptr;i++){
        fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", itemp->hkl[0], itemp->hkl[1], itemp->hkl[2], itemp->pos[0], itemp->pos[1], itemp->intensity, consti*itemp->intensity);
        fflush(fp);
        itemp=itemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for selected-area electron pattern stored in %s.\n", ded_path);
}

void DED::add_intensity_node(double hkl[3], double pos[2], double intensity)
{
    if(ihead==nullptr&&itail==nullptr){
        ihead=itail=new DED_INODE;
        vector_copy(itail->hkl, hkl);
        itail->pos[0]=pos[0]; itail->pos[1]=pos[1];
        itail->intensity=intensity;
    }else{
        itail->next=new DED_INODE;
        itail=itail->next;
        vector_copy(itail->hkl, hkl);
        itail->pos[0]=pos[0]; itail->pos[1]=pos[1];
        itail->intensity=intensity;
    }
    if(intensity<intensity_min) intensity_min=intensity;
    if(intensity>intensity_max) intensity_max=intensity;
    numi++;
}

void DED::compute_dynamic_matrix(complex<double> **dmat, CELL *cell, DED_GVECTOR *gvec)
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
                //printf("%.2f+%.2fi\t", dmat[ir][ic].real(), dmat[ir][ic].imag());
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
                //printf("%.2f+%.2fi\t", dmat[ir][ic].real(), dmat[ir][ic].imag());
            }
            ctemp=ctemp->nexts;
        }
        //printf("\n");
        rtemp=rtemp->nexts;
    }
}

void DED::compute_Lgh_matrix(complex<double> **Lgh, complex<double> **DMAT, double Z, double KN, int NS)
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
    LWORK=std::min(LWMAX, int(WORK[0].real())); LDA=NS;
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
    if(0!=INFO){
        printf("[ERROR] Unable to return INFO as zero when calling the eigenvalue solver.");
        exit(EXIT_FAILURE);
    }
    int MILWORK=64*NS;
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
    double TPI=TWO_PI*Z;
    for(int j=0;j<NS;j++){
        for(int k=0;k<NS;k++){
            complex<double> sumq(0.0, 0.0);
            complex<double> q(TPI*(NW[j].imag()+NW[k].imag()), TPI*(NW[j].real()-NW[k].real()));
            if(q.real()<0.0) q=-q;
            Lgh[j][k]=NCG[j][k]*exp(-q)*NCGINV[k][0];
        }
    }
    // for(int j=0;j<NS;j++){
    //     complex<double> q(cos(Z*NW[j].real()), sin(Z*NW[j].real()));
    //     q*=exp(-Z*NW[j].imag());
    //     for(int k=0;k<NS;k++){
    //         Lgh[k][j]=NCG[k][j]*q*NCGINV[j][0];
    //         //printf("%.2f %.2f %.2f %.2f %.2f %.2f\n", q.real(), q.imag(), NCG[j][k].real(), NCG[j][k].imag(), NCGINV[k][0].real(), NCGINV[k][0].imag());
    //     }
    // }
}