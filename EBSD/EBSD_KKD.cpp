#include "EBSD_KKD.h"

QUATERNION quate_conjg(QUATERNION q)
{
    QUATERNION qc={q.c1, -q.c2, -q.c3, -q.c4};
    return qc;
}

QUATERNION quate_multi(QUATERNION q1, QUATERNION q2)
{
    QUATERNION qm;
    double eps=1.0;
    qm.c1=q1.c1*q2.c1-q1.c2*q2.c2-q1.c3*q2.c3-q1.c4*q2.c4;
    qm.c2=q1.c1*q2.c2+q1.c2*q2.c1+eps*(q1.c3*q2.c4-q1.c4*q2.c3);
    qm.c3=q1.c1*q2.c3+q1.c3*q2.c1+eps*(q1.c4*q2.c2-q1.c2*q2.c4);
    qm.c4=q1.c1*q2.c4+q1.c4*q2.c1+eps*(q1.c2*q2.c3-q1.c3*q2.c2); 
    return qm;
}

void quate_rotate(double r_v[3], double v[3], QUATERNION q)
{
    QUATERNION qv={0.0, v[0], v[1], v[2]};
    QUATERNION r_qv=quate_multi(q, quate_multi(qv, quate_conjg(q)));
    r_v[0]=r_qv.c2; r_v[1]=r_qv.c3; r_v[2]=r_qv.c4; 
}

QUATERNION quate_convert(double v[3], double angle)
{
    QUATERNION q;
    if(fabs(angle)<1e-12){
        q.c1=1.0; q.c2=q.c3=q.c4=0.0;
    }else{
        double c=cos(0.5*angle);
        double s=sin(0.5*angle);
        q.c1=c; q.c2=s*v[0]; q.c3=s*v[1]; q.c4=s*v[2];
    }
    return q;
}

EBSD_KKD::EBSD_KKD(const char *file_path)
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

EBSD_KKD::~EBSD_KKD()
{

}

bool EBSD_KKD::read_parameters_from_nml(const char *file_path)
{
    FILE *fp=fopen(file_path, "r");
    if(fp==NULL){
        printf("Error! Unable to open file %s.\n", file_path);
    }
    fseek(fp, 0, SEEK_SET);

    char keys[][MAX_LENTH_IN_NAME]={"smallest_spacing", "voltage", "intensity_cutoff"}; int keyi=0;
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
            fscanf(fp, "%lf", &voltage);
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, keys[2])){
            fscanf(fp, "%lf", &ci);
            keyi++;
            continue;
        }
    }
    if(3==keyi) return true;
    else return false;
}

bool EBSD_KKD::write_parameters_into_hdf5(const char* file_path)
{
    EBSD_HDF5 hdf;
    hdf.open(file_path);
    hdf.write_group("/EMData");
    hdf.write_group("/EMData/EMkinematical");
    hdf.write_array_2d("/EMData/EMkinematical/masterNH", mLPNH, npx, npy);
    hdf.write_array_2d("/EMData/EMkinematical/masterSH", mLPSH, npx, npy);
    hdf.write_array_2d("/EMData/EMkinematical/stereoNH", mSPNH, npx, npy);
    hdf.write_array_2d("/EMData/EMkinematical/stereoSH", mSPSH, npx, npy);
    hdf.close();
    return true;
}

void EBSD_KKD::compute_master_pattern(const char* file_path)
{
    EBSD_CELL cell(file_path);
    cell.compute_reflection_range(dmin);
    printf("Range of reflections along a*, b*, and c* = %d, %d, and %d.\n", cell.HKL[0], cell.HKL[1], cell.HKL[2]);

    int imh=cell.HKL[0], imk=cell.HKL[1], iml=cell.HKL[2];
    int num=(2*imh+1)*(2*imk+1)*(2*iml+1);
    double **kvec;
    double *LUTUmod, *LUTUang;
    callocate_2d(&kvec, num, 3, 0.0);
    callocate(&LUTUmod, num, 0.0);
    callocate(&LUTUang, num, 0.0);
    int count=0;
    cell.compute_relativistic_wavelength(voltage);
    for(int ih=-imh;ih<=imh;ih++){
        for(int ik=-imk;ik<=imk;ik++){
            for(int il=-iml;il<=iml;il++){
                FOURIER fouri;
                double hkl[3]={-double(ih), -double(ik), -double(il)};
                cell.compute_Fourier_coefficient(&fouri, voltage, hkl);
                if(fouri.Vmod>=ci){
                    if(0==(abs(ih)+abs(ik)+abs(il))){
                        LUTUmod[count]=0.0;
                    }else{
                        if(0.0!=fouri.Vmod){
                            LUTUmod[count]=pow(fabs(fouri.Vmod), 0.5);
                        }else{
                            LUTUmod[count]=0.0;
                        }
                    }
                    double g[3]={double(ih), double(ik), double(il)};
                    double k[3];
                    cell.cartesian(k, g);
                    cell.normalize(k, k);
                    kvec[count][0]=k[0];
                    kvec[count][1]=k[1];
                    kvec[count][2]=k[2];
                    LUTUang[count]=asin(0.5*cell.fouri0.lambda*cell.length(g));
                    count++;
                }
            }
        }
    }
    printf("Total number of entries found = %d.\n", count);

    npx=2*nump+1; npy=2*nump+1;
    callocate_2d(&mLPNH, npx, npy, 1.0);
    callocate_2d(&mLPSH, npx, npy, 1.0);
    mallocate_2d(&mSPNH, npx, npy);
    mallocate_2d(&mSPSH, npx, npy);
    double dphi=2.0*PI/double(numphi);
    double z[3]={0.0, 0.0, 1.0};
    double **kc;
    mallocate_2d(&kc, numphi, 3);
    for(int i=0;i<=count;i++){
        double ctheta=cos(LUTUang[i]);
        double stheta=sin(LUTUang[i]);
        for(int j=0;j<numphi;j++){//located at the Kossel cone
            double phi=double(j)*dphi;
            kc[j][0]=ctheta*cos(phi);
            kc[j][1]=ctheta*sin(phi);
            kc[j][2]=stheta;
        }

        double k[3]={kvec[i][0], kvec[i][1], kvec[i][2]};
        double x=cell.dot(z, k, 'c');
        if(1.0!=x){
            double ax[3];
            cell.cross(ax, k, z, 'c', 'c');
            cell.normalize(ax, ax);
            double ang=acos(x);
            if(ang<0.0){
                ax[0]=-ax[0]; ax[1]=-ax[1]; ax[2]=-ax[2];
                ang=-ang;
            }
            QUATERNION q=quate_conjg(quate_convert(ax, ang));
            for(int j=0;j<numphi;j++){
                double r_kc[3];
                quate_rotate(r_kc, kc[j], q);
                cell.normalize(r_kc, r_kc);
                kc[j][0]=r_kc[0]; kc[j][1]=r_kc[1]; kc[j][2]=r_kc[2]; 
            }
        }
        for(int j=0;j<numphi;j++){
            double xy[2]; int ierr;
            compute_square_Lambert(xy, ierr, kc[j]);
            xy[0]*=double(nump); xy[1]*=double(nump);
            if(0!=ierr){
                printf("Error! Unable to compute square Lambert interpolation using (%.2f, %.2f).\n", xy[0], xy[1]);
            }
            int ix=int(xy[0]+nump)-nump;
            int iy=int(xy[1]+nump)-nump;
            if(kc[j][2]>=0.0){
                apply_anti_aliasing(mLPNH, xy, ix, iy, LUTUmod[i]);
            }else{
                apply_anti_aliasing(mLPSH, xy, ix, iy, LUTUmod[i]);
            }
        }
    }
    deallocate_2d(kc, numphi);
    deallocate_2d(kvec, num);
    deallocate(LUTUmod);
    deallocate(LUTUang);
    double minn=MAX_MIN, mins=MAX_MIN;
    double maxn=-MAX_MIN, maxs=-MAX_MIN;
    for(int j=0;j<npx;j++){
        for(int k=0;k<npy;k++){
            if(minn>mLPNH[j][k]) minn=mLPNH[j][k];
            if(mins>mLPSH[j][k]) mins=mLPSH[j][k];
            if(maxn<mLPNH[j][k]) maxn=mLPNH[j][k];
            if(maxs<mLPSH[j][k]) maxs=mLPSH[j][k];
        }
    }
    double minm=min(minn, mins);
    for(int j=0;j<npx;j++){
        for(int k=0;k<npy;k++){
            mLPNH[j][k]=(mLPNH[j][k]-minm)/(maxn-minm);
            mLPSH[j][k]=(mLPSH[j][k]-minm)/(maxs-minm);
        }
    }
    printf("minn %.5f mins %.5f\n", minn, mins);
    printf("maxn %.5f maxs %.5f\n", maxn, maxs);
    for(int i=-nump;i<=nump;i++){
        for(int j=-nump;j<=nump;j++){
            double xy[2]={double(i)/double(nump), double(j)/double(nump)};
            double xyz[3]; int ierr;
            compute_sphere_from_stereographic_projection(xyz, ierr, xy);
            normalize_vector(xyz, xyz);
            int ix=i+nump, iy=j+nump;
            if(0!=ierr){
                mSPNH[ix][iy]=0.0;
                mSPSH[ix][iy]=0.0;
            }else{
                mSPNH[ix][iy]=get_Lambert_interpolation(xyz, mLPNH);
                mSPSH[ix][iy]=get_Lambert_interpolation(xyz, mLPSH);
            }
        }
    }
    write_parameters_into_hdf5(file_path);
    create_images("mLPNH.png", mLPNH);
    create_images("mLPSH.png", mLPSH);
    create_images("mSPNH.png", mSPNH);
    create_images("mSPSH.png", mSPSH);
}

void EBSD_KKD::apply_anti_aliasing(double **master, double xy[2], int ix, int iy, double intensity)
{
    if(abs(ix)<nump&&abs(iy)<nump){
        double dx=xy[0]-double(ix);
        double dy=xy[1]-double(iy);
        int i=ix+nump, j=iy+nump;
        if(fabs(dx)>fabs(dy)){
            double d=fabs(dx);
            if(dx<0.0){
                master[i-1][j]-=d*intensity;
                master[i][j]-=(1.0-d)*intensity;
            }else{
                master[i][j]-=(1.0-d)*intensity;
                master[i+1][j]-=d*intensity;
            }
        }else{
            double d=fabs(dy);
            if(dy<0.0){
                master[i][j-1]-=d*intensity;
                master[i][j]-=(1.0-d)*intensity;
            }else{
                master[i][j]-=(1.0-d)*intensity;
                master[i][j+1]-=d*intensity;
            }
        }
    }
}

void EBSD_KKD::test_images(const char *file_path)
{
    size_t npx=2*nump+1, npy=2*nump+1;
    EBSD_HDF5 hdf;
    hdf.open(file_path);
    hdf.read_array_2d("/EMData/EMkinematical/masterNH", &mLPNH, npx, npy);
    hdf.read_array_2d("/EMData/EMkinematical/masterSH", &mLPSH, npx, npy);
    hdf.read_array_2d("/EMData/EMkinematical/stereoNH", &mSPNH, npx, npy);
    hdf.read_array_2d("/EMData/EMkinematical/stereoSH", &mSPSH, npx, npy);
    hdf.close();
    create_images("mLPNH-ref.png", mLPNH);
    create_images("mLPSH-ref.png", mLPSH);
    create_images("mSPNH-ref.png", mSPNH);
    create_images("mSPSH-ref.png", mSPSH);
}

double EBSD_KKD::get_Lambert_interpolation(double xyz[3], double **mat, bool hexagonal_flag)
{
    int ix, iy, ixp, iyp;
    double dx, dy, dxm, dym;
    compute_Lambert_interpolation(xyz, nump, hexagonal_flag, ix, iy, ixp, iyp, dx, dy, dxm, dym);
    double res=mat[ix][iy]*dxm*dym+mat[ixp][iy]*dx*dym+mat[ix][iyp]*dxm*dy+mat[ixp][iyp]*dx*dy;
    return res;
}

void EBSD_KKD::create_images(const char* png_path, double **mat)
{
    int nrow=2*nump+1, ncol=2*nump+1;
    double minm=MAX_MIN, maxm=-MAX_MIN, diffm;
    for(int j=0;j<nrow;j++){
        for(int k=0;k<ncol;k++){
            if(minm>mat[j][k]) minm=mat[j][k];
            if(maxm<mat[j][k]) maxm=mat[j][k];
        }
    }
    printf("min %.5f max %.5f\n", minm, maxm);
    diffm=maxm-minm;
    int num=nrow*ncol;
    double *wmat;
    reshape2d(&wmat, mat, nrow, ncol);
    unsigned char *pixels;
    mallocate(&pixels, 3*num);
    for(int i=0;i<num;i++){
        pixels[i*3]=pixels[i*3+1]=pixels[i*3+2]=round((wmat[i]-minm)/diffm*255.0);
    }
    create_image(png_path, pixels, nrow, ncol);
}