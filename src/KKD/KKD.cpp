#include "KKD.h"

KKD::KKD(int zone[3], double threshold_ratio, double screenD, int screen_dpi, int screenW, int screenH)
{
    SED sed(model, spacing, Kmagnitude_max);
    radiusK=1.0/model->lambda;
    double dK_eff=4*radiusK/sqrt(double(pow(screen_dpi, 2)*pow(screenD, 2)+1));
    int screen_npx=screen_dpi*screenW, screen_npy=screen_dpi*screenH;
    if(0==zone[0]&&0==zone[1]&&0==zone[2]){
        compute_Kikuchi_sphere_projection(screen_npx);
        dK_eff*=2.0;
    }else{
        compute_Kikuchi_sphere_projection(zone, screen_npx, screen_npy, screenD, screenW, screenH);
    }
    compute_Kikuchi_intensity_projection(&sed, threshold_ratio, dK_eff);
}

KKD::~KKD()
{
    if(0!=numpx&&0!=numpy){
        radiusK=0.0;
        deallocate_3d(screenG, numpx, numpy);
        deallocate_2d(screenI, numpx);
        numpx=0; numpy=0;
    }
}

void KKD::compute_Kikuchi_sphere_projection(int screen_npx)
{
    double dp, dt;
    dp=dt=2*PI/double(screen_npx-1);
    int numpx=screen_npx, numpy=screen_npx/2;
    callocate_3d(&screenG, numpx, numpy, 3, 0.0);
    for(int i=0;i<numpx;i++){
        for(int j=0;j<numpy;j++){
            double phi=double(i)*dp;//azimuth angle, 0-2pi
            double theta=double(j)*dt;//polar angle, 0-pi
            screenG[i][j][0]=radiusK*sin(theta)*cos(phi);//spherical to cartesian
            screenG[i][j][1]=radiusK*sin(theta)*sin(phi);
            screenG[i][j][2]=radiusK*cos(theta);
        }
    }
}

void KKD::compute_Kikuchi_sphere_projection(int zone[3], int screen_npx, int screen_npy, double screenD, double screenW, double screenH)
{
    numpx=screen_npx; numpy=screen_npy;
    callocate_3d(&screenG, screen_npx, screen_npy, 3, screenD);
    double imW=screenW/2.0, imH=screenH/2.0;
    double dy=screenW/double(screen_npx-1), dz=screenH/double(screen_npy-1);
    for(int i=0;i<screen_npx;i++){
        for(int j=0;j<screen_npy;j++){
            screenG[i][j][1]=double(i)*dy-imW;
            screenG[i][j][2]=double(j)*dz-imH;
        }
    }
    double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    vector_normalize(n_zone, n_zone);
    double x=n_zone[0], y=n_zone[1], z=n_zone[2];
    double xy=sqrt(x*x+y*y);
    double Rp[3][3]={{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    if(0.0!=xy){
        Rp[0][0]=x/xy; Rp[0][1]=-y/xy;
        Rp[1][0]=y/xy; Rp[1][1]=x/xy;//rotate azimuth angle around z-axis
    }
    double Rt[3][3]={{xy, 0.0, -z}, {0.0, 1.0, 0.0}, {z, 0.0, xy}};//rotate polar angle around y-axis
    double R[3][3];
    matrix_multiply(R, Rp, Rt);
    for(int i=0;i<screen_npx;i++){
        for(int j=0;j<screen_npy;j++){
            vector_normalize(screenG[i][j], screenG[i][j]);
            screenG[i][j][0]*=radiusK;
            screenG[i][j][1]*=radiusK;
            screenG[i][j][2]*=radiusK;
            vector_rotate(screenG[i][j], screenG[i][j], R);//zone axis transformation
        }
    }
}

void KKD::compute_Kikuchi_intensity_projection(SED *sed, double threshold_ratio, double interceptK)
{
    callocate_2d(&screenI, numpx, numpy, 0.0);
    double intensity_max=fmax(fabs(sed->intensity_max), fabs(sed->intensity_min));
    double intensity_threshold=intensity_max*threshold_ratio;
    for(int i=0;i<sed->numk;i++){
        double kmag=vector_length(sed->Kvectors[i]);
        if(kmag<0.2||fabs(sed->Kintensity[i])<intensity_threshold){
            continue;
        }
        for(int j=0;j<numpx;j++){
            for(int k=0;k<numpy;k++){
                double d[3];
                vector_differance(d, screenG[j][k], sed->Kvectors[i]);
                double dmag=vector_length(d);
                double upper_bound=sqrt(radiusK*radiusK+kmag*interceptK/2.0);
                double lower_bound=sqrt(radiusK*radiusK-kmag*interceptK/2.0);
                if(dmag<=upper_bound&&dmag>=lower_bound){
                    screenI[j][k]+=sed->Kintensity[i];
                }
            }
        }
    }
}

// KKD::KKD(const char *file_path)
// {
//     size_t len=strlen(file_path);
//     bool flag;
//     if(len>=4&&strcmp(file_path+len-4, ".nml")==0){
//         flag=read_parameters_from_nml(file_path);
//     }else{
//         printf("Error! Unrecognized file %s.", file_path);
//         exit(EXIT_FAILURE);
//     }
//     if(!flag){
//         printf("Error! Unrecognized parameters in file %s.", file_path);
//         exit(EXIT_FAILURE);
//     }
// }

// bool KKD::read_parameters_from_nml(const char *file_path)
// {
//     FILE *fp=fopen(file_path, "r");
//     if(fp==NULL){
//         printf("Error! Unable to open file %s.\n", file_path);
//     }
//     fseek(fp, 0, SEEK_SET);

//     char keys[][MAX_LENTH_IN_NAME]={"smallest_spacing", "voltage", "intensity_cutoff"}; int keyi=0;
//     char readbuff[MAX_LENTH_IN_NAME];
//     for(int i=0;(keyi<3)&&(i<MAX_WORD_NUMBER);i++){
//         fscanf(fp, "%s", readbuff);
//         if('#'==readbuff[0]){
//             while(('\n'!=fgetc(fp))&&(!feof(fp)));
//             continue;
//         }
//         if(0==strcmp(readbuff, keys[0])){
//             fscanf(fp, "%lf", &dmin);
//             keyi++;
//             continue;
//         }
//         if(0==strcmp(readbuff, keys[1])){
//             fscanf(fp, "%lf", &voltage);
//             keyi++;
//             continue;
//         }
//         if(0==strcmp(readbuff, keys[2])){
//             fscanf(fp, "%lf", &threshold);
//             keyi++;
//             continue;
//         }
//     }
//     if(3==keyi) return true;
//     else return false;
// }

// bool KKD::write_parameters_into_hdf5(const char* file_path)
// {
//     EBSD_HDF5 hdf;
//     hdf.open(file_path);
//     hdf.write_group("/EMData");
//     hdf.write_group("/EMData/EMkinematical");
//     hdf.write_array_2d("/EMData/EMkinematical/masterNH", mLPNH, npx, npy);
//     hdf.write_array_2d("/EMData/EMkinematical/masterSH", mLPSH, npx, npy);
//     hdf.write_array_2d("/EMData/EMkinematical/stereoNH", mSPNH, npx, npy);
//     hdf.write_array_2d("/EMData/EMkinematical/stereoSH", mSPSH, npx, npy);
//     hdf.close();
//     return true;
// }

// void KKD::compute_master_pattern(const char* file_path)
// {
//     EBSD_CELL cell(file_path);
//     cell.compute_reflection_range(dmin);
//     printf("Range of reflections along a*, b*, and c* = %d, %d, and %d.\n", cell.HKL[0], cell.HKL[1], cell.HKL[2]);

//     int imh=cell.HKL[0], imk=cell.HKL[1], iml=cell.HKL[2];
//     int num=(2*imh+1)*(2*imk+1)*(2*iml+1);
//     double **kvec;
//     double *LUTUmod, *LUTUang;
//     callocate_2d(&kvec, num, 3, 0.0);
//     callocate(&LUTUmod, num, 0.0);
//     callocate(&LUTUang, num, 0.0);
//     int count=0;
//     cell.compute_relativistic_wavelength(voltage);
//     for(int ih=-imh;ih<=imh;ih++){
//         for(int ik=-imk;ik<=imk;ik++){
//             for(int il=-iml;il<=iml;il++){
//                 FOURIER fouri;
//                 double hkl[3]={-double(ih), -double(ik), -double(il)};
//                 cell.compute_Fourier_coefficient(&fouri, voltage, hkl);
//                 if(fouri.Vmod>=threshold){
//                     if(0==(abs(ih)+abs(ik)+abs(il))){
//                         LUTUmod[count]=0.0;
//                     }else{
//                         if(0.0!=fouri.Vmod){
//                             LUTUmod[count]=pow(fabs(fouri.Vmod), 0.5);
//                         }else{
//                             LUTUmod[count]=0.0;
//                         }
//                     }
//                     double g[3]={double(ih), double(ik), double(il)};
//                     double k[3];
//                     cell.cartesian(k, g);
//                     cell.normalize(k, k);
//                     kvec[count][0]=k[0];
//                     kvec[count][1]=k[1];
//                     kvec[count][2]=k[2];
//                     LUTUang[count]=asin(0.5*cell.fouri0.lambda*cell.length(g));
//                     count++;
//                 }
//             }
//         }
//     }
//     printf("Total number of entries found = %d.\n", count);

//     npx=2*nump+1; npy=2*nump+1;
//     callocate_2d(&mLPNH, npx, npy, 1.0);
//     callocate_2d(&mLPSH, npx, npy, 1.0);
//     mallocate_2d(&mSPNH, npx, npy);
//     mallocate_2d(&mSPSH, npx, npy);
//     double dphi=2.0*PI/double(numKc);
//     double z[3]={0.0, 0.0, 1.0};
//     double **kc;
//     mallocate_2d(&kc, numKc, 3);
//     for(int i=0;i<=count;i++){
//         double ctheta=cos(LUTUang[i]);
//         double stheta=sin(LUTUang[i]);
//         for(int j=0;j<numKc;j++){//located at the Kossel cone
//             double phi=double(j)*dphi;
//             kc[j][0]=ctheta*cos(phi);
//             kc[j][1]=ctheta*sin(phi);
//             kc[j][2]=stheta;
//         }

//         double k[3]={kvec[i][0], kvec[i][1], kvec[i][2]};
//         double x=cell.dot(z, k, 'c');
//         if(1.0!=x){
//             double axis[3];
//             cell.cross(axis, k, z, 'c', 'c');
//             cell.normalize(axis, axis);
//             double angle=acos(x);
//             if(angle<0.0){
//                 axis[0]=-axis[0]; axis[1]=-axis[1]; axis[2]=-axis[2];
//                 angle=-angle;
//             }
//             QUATERNION q=quate_conjg(quate_convert(axis, angle));
//             for(int j=0;j<numKc;j++){
//                 double r_kc[3];
//                 quate_rotate(r_kc, kc[j], q);
//                 cell.normalize(r_kc, r_kc);
//                 kc[j][0]=r_kc[0]; kc[j][1]=r_kc[1]; kc[j][2]=r_kc[2]; 
//             }
//         }
//         for(int j=0;j<numKc;j++){
//             double xy[2]; int ierr;
//             compute_square_Lambert(xy, ierr, kc[j]);
//             xy[0]*=double(nump); xy[1]*=double(nump);
//             if(0!=ierr){
//                 printf("Error! Unable to compute square Lambert interpolation using (%.2f, %.2f).\n", xy[0], xy[1]);
//             }
//             int ix=int(xy[0]+nump)-nump;
//             int iy=int(xy[1]+nump)-nump;
//             if(kc[j][2]>=0.0){
//                 apply_anti_aliasing(mLPNH, xy, ix, iy, LUTUmod[i]);
//             }else{
//                 apply_anti_aliasing(mLPSH, xy, ix, iy, LUTUmod[i]);
//             }
//         }
//     }
//     deallocate_2d(kc, numKc);
//     deallocate_2d(kvec, num);
//     deallocate(LUTUmod);
//     deallocate(LUTUang);
//     double minn=MAX_MIN, mins=MAX_MIN;
//     double maxn=-MAX_MIN, maxs=-MAX_MIN;
//     for(int j=0;j<npx;j++){
//         for(int k=0;k<npy;k++){
//             if(minn>mLPNH[j][k]) minn=mLPNH[j][k];
//             if(mins>mLPSH[j][k]) mins=mLPSH[j][k];
//             if(maxn<mLPNH[j][k]) maxn=mLPNH[j][k];
//             if(maxs<mLPSH[j][k]) maxs=mLPSH[j][k];
//         }
//     }
//     double minm=min(minn, mins);
//     for(int j=0;j<npx;j++){
//         for(int k=0;k<npy;k++){
//             mLPNH[j][k]=(mLPNH[j][k]-minm)/(maxn-minm);
//             mLPSH[j][k]=(mLPSH[j][k]-minm)/(maxs-minm);
//         }
//     }
//     printf("minn %.5f mins %.5f\n", minn, mins);
//     printf("maxn %.5f maxs %.5f\n", maxn, maxs);
//     for(int i=-nump;i<=nump;i++){
//         for(int j=-nump;j<=nump;j++){
//             double xy[2]={double(i)/double(nump), double(j)/double(nump)};
//             double xyz[3]; int ierr;
//             compute_sphere_from_stereographic_projection(xyz, ierr, xy);
//             normalize_vector(xyz, xyz);
//             int ix=i+nump, iy=j+nump;
//             if(0!=ierr){
//                 mSPNH[ix][iy]=0.0;
//                 mSPSH[ix][iy]=0.0;
//             }else{
//                 mSPNH[ix][iy]=get_Lambert_interpolation(xyz, mLPNH);
//                 mSPSH[ix][iy]=get_Lambert_interpolation(xyz, mLPSH);
//             }
//         }
//     }
//     write_parameters_into_hdf5(file_path);
//     create_images("mLPNH.png", mLPNH);
//     create_images("mLPSH.png", mLPSH);
//     create_images("mSPNH.png", mSPNH);
//     create_images("mSPSH.png", mSPSH);
// }

// void KKD::apply_anti_aliasing(double **master, double xy[2], int ix, int iy, double intensity)
// {
//     if(abs(ix)<nump&&abs(iy)<nump){
//         double dx=xy[0]-double(ix);
//         double dy=xy[1]-double(iy);
//         int i=ix+nump, j=iy+nump;
//         if(fabs(dx)>fabs(dy)){
//             double d=fabs(dx);
//             if(dx<0.0){
//                 master[i-1][j]-=d*intensity;
//                 master[i][j]-=(1.0-d)*intensity;
//             }else{
//                 master[i][j]-=(1.0-d)*intensity;
//                 master[i+1][j]-=d*intensity;
//             }
//         }else{
//             double d=fabs(dy);
//             if(dy<0.0){
//                 master[i][j-1]-=d*intensity;
//                 master[i][j]-=(1.0-d)*intensity;
//             }else{
//                 master[i][j]-=(1.0-d)*intensity;
//                 master[i][j+1]-=d*intensity;
//             }
//         }
//     }
// }

// void KKD::test_images(const char *file_path)
// {
//     size_t npx=2*nump+1, npy=2*nump+1;
//     EBSD_HDF5 hdf;
//     hdf.open(file_path);
//     hdf.read_array_2d("/EMData/EMkinematical/masterNH", &mLPNH, npx, npy);
//     hdf.read_array_2d("/EMData/EMkinematical/masterSH", &mLPSH, npx, npy);
//     hdf.read_array_2d("/EMData/EMkinematical/stereoNH", &mSPNH, npx, npy);
//     hdf.read_array_2d("/EMData/EMkinematical/stereoSH", &mSPSH, npx, npy);
//     hdf.close();
//     create_images("mLPNH-ref.png", mLPNH);
//     create_images("mLPSH-ref.png", mLPSH);
//     create_images("mSPNH-ref.png", mSPNH);
//     create_images("mSPSH-ref.png", mSPSH);
// }

// double KKD::get_Lambert_interpolation(double xyz[3], double **mat, bool hexagonal_flag)
// {
//     int ix, iy, ixp, iyp;
//     double dx, dy, dxm, dym;
//     compute_Lambert_interpolation(xyz, nump, hexagonal_flag, ix, iy, ixp, iyp, dx, dy, dxm, dym);
//     double res=mat[ix][iy]*dxm*dym+mat[ixp][iy]*dx*dym+mat[ix][iyp]*dxm*dy+mat[ixp][iyp]*dx*dy;
//     return res;
// }

// void KKD::create_images(const char* png_path, double **mat)
// {
//     int nrow=2*nump+1, ncol=2*nump+1;
//     double minm=MAX_MIN, maxm=-MAX_MIN, diffm;
//     for(int j=0;j<nrow;j++){
//         for(int k=0;k<ncol;k++){
//             if(minm>mat[j][k]) minm=mat[j][k];
//             if(maxm<mat[j][k]) maxm=mat[j][k];
//         }
//     }
//     printf("min %.5f max %.5f\n", minm, maxm);
//     diffm=maxm-minm;
//     int num=nrow*ncol;
//     double *wmat;
//     reshape2d(&wmat, mat, nrow, ncol);
//     int char *pixels;
//     mallocate(&pixels, 3*num);
//     for(int i=0;i<num;i++){
//         pixels[i*3]=pixels[i*3+1]=pixels[i*3+2]=round((wmat[i]-minm)/diffm*255.0);
//     }
//     create_image(png_path, pixels, nrow, ncol);
// }