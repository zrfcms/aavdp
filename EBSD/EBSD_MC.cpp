#include "EBSD_MC.h"

EBSD_MC::EBSD_MC(const char *file_path)
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
}

EBSD_MC::~EBSD_MC()
{
    deallocate_3d(accum_SP, size_SP[0], size_SP[1]);
    deallocate_3d(accum_E, size_E[0], size_E[1]);
    deallocate_4d(accum_z, size_z[0], size_z[1], size_z[2]);
}

bool EBSD_MC::read_parameters_from_nml(const char *file_path)
{
    FILE *fp=fopen(file_path, "r");
    if(fp==NULL){
        printf("Error! Unable to open file %s.\n", file_path);
    }
    fseek(fp, 0, SEEK_SET);

    char keys[][MAX_LENTH_IN_NAME]={"tilt_angle_omega", "tilt_angle_sigma", "energy", "depth", 
                                    "electron_number", "pixel_number"}; 
    int keyi=0;
    char readbuff[MAX_LENTH_IN_NAME];
    for(int i=0;(keyi<6)&&(i<MAX_WORD_NUMBER);i++){
        fscanf(fp, "%s", readbuff);
        if('#'==readbuff[0]){
            while(('\n'!=fgetc(fp))&&(!feof(fp)));
            continue;
        }
        if(0==strcmp(readbuff, keys[0])){
            fscanf(fp, "%lf", &omega);
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, keys[1])){
            fscanf(fp, "%d", &sigma);
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, keys[2])){
            fscanf(fp, "%lf", &Ehistmin);
            fscanf(fp, "%lf", &EkeV);
            fscanf(fp, "%lf", &Ebin);
            keyi++;
            continue;
        }
        if(0==strcmp(readbuff, keys[3])){
            fscanf(fp, "%lf", &depthmax);
            fscanf(fp, "%lf", &depthstep);
            keyi++;
            continue; 
        }
        if(0==strcmp(readbuff, keys[4])){
            fscanf(fp, "%d", &num_e);
            fscanf(fp, "%d", &num_e);
            keyi++;
            continue; 
        }
        if(0==strcmp(readbuff, keys[5])){
            fscanf(fp, "%d", &multiplier);
            keyi++;
            continue; 
        }
    }
    if(3==keyi) return true;
    else return false;
}

bool EBSD_MC::read_parameters_from_hdf5(const char *file_path)
{
    size_t size_SP5[3];
    double ***accum_SP5=nullptr;
    size_t size_E5[3];
    int    ***accum_E5=nullptr;
    size_t size_z5[4];
    int    ****accum_z5=nullptr;
    EBSD_HDF5 hdf;
    hdf.open(file_path);
    hdf.read_array_3d("/EMData/MCOpenCL/accumSP", &accum_SP5, size_SP5[0], size_SP5[1], size_SP5[2]);
    hdf.read_array_3d("/EMData/MCOpenCL/accum_e", &accum_E5, size_E5[0], size_E5[1], size_E5[2]);
    hdf.read_array_4d("/EMData/MCOpenCL/accum_z", &accum_z5, size_z5[0], size_z5[1], size_z5[2], size_z5[3]);
    hdf.read("/EMData/MCOpenCL/numEbins", numEbin);
    hdf.read("/EMData/MCOpenCL/numzbins", numzbin);
    hdf.read("/EMData/MCOpenCL/totnum_el", num_e);
    hdf.read("/NMLparameters/MCCLNameList/omega", omega);
    hdf.read("/NMLparameters/MCCLNameList/sig", sigma);
    hdf.read("/NMLparameters/MCCLNameList/Ehistmin", Ehistmin);
    hdf.read("/NMLparameters/MCCLNameList/EkeV", EkeV);
    hdf.read("/NMLparameters/MCCLNameList/Ebinsize", Ebin);
    hdf.read("/NMLparameters/MCCLNameList/depthmax", depthmax);
    hdf.read("/NMLparameters/MCCLNameList/depthstep", depthstep);
    hdf.read("/NMLparameters/MCCLNameList/multiplier", multiplier);
    hdf.read("/NMLparameters/MCCLNameList/numsx", nump);
    size_SP[0]=size_SP5[2]; size_SP[1]=size_SP5[0]; size_SP[2]=size_SP5[1];
    mallocate_3d(&accum_SP, size_SP[0], size_SP[1], size_SP[2]);
    for(int i=0;i<size_SP[0];i++){
        for(int j=0;j<size_SP[1];j++){
            for(int k=0;k<size_SP[2];k++){
                accum_SP[i][j][k]=accum_SP5[j][k][i];
            }
        }
    }
    size_E[0]=size_E5[2]; size_E[1]=size_E5[0]; size_E[2]=size_E5[1];
    mallocate_3d(&accum_E, size_E[0], size_E[1], size_E[2]);
    for(int i=0;i<size_E[0];i++){
        for(int j=0;j<size_E[1];j++){
            for(int k=0;k<size_E[2];k++){
                accum_E[i][j][k]=accum_E5[j][k][i];
            }
        }
    }
    size_z[0]=size_z5[3]; size_z[1]=size_z5[2]; size_z[2]=size_z5[0]; size_z[3]=size_z5[1];
    mallocate_4d(&accum_z, size_z[0], size_z[1], size_z[2], size_z[3]);
    for(int i=0;i<size_z[0];i++){
        for(int j=0;j<size_z[1];j++){
            for(int k=0;k<size_z[2];k++){
                for(int n=0;n<size_z[3];n++){
                    accum_z[i][j][k][n]=accum_z5[k][n][j][i];
                }
            }
        }
    }
    return true; 
}

bool EBSD_MC::write_parameters_into_hdf5(const char *file_path)
{
    EBSD_HDF5 hdf;
    hdf.open(file_path);
    hdf.write_group("/EMData");
    hdf.write_group("/EMData/MCOpenCL");
    hdf.write_array_3d("/EMData/MCOpenCL/accumSP", accum_SP, size_SP[0], size_SP[1], size_SP[2]);
    hdf.write_array_3d("/EMData/MCOpenCL/accum_e", accum_E, size_E[0], size_E[1], size_E[2]);
    hdf.write_array_4d("/EMData/MCOpenCL/accum_z", accum_z, size_z[0], size_z[1], size_z[2], size_z[3]);
    hdf.write("/EMData/MCOpenCL/multiplier", multiplier);
    hdf.write("/EMData/MCOpenCL/numEbins", numEbin);
    hdf.write("/EMData/MCOpenCL/numzbins", numzbin);
    hdf.write("/EMData/MCOpenCL/totnum_el", num_e);
    hdf.close();
    return true;
}

double EBSD_MC::compute_alphainfreepath(double E)
{
    return 3.4E-3*pow(ave_Z, 2.0/3.0)/E;
}

//Set the free path for the energy and scale it by a random number
double EBSD_MC::compute_freepath(double E, double alpha)
{
    double area_inv=CONST_IN_AREA_INV*pow(E*(E+1024.0)/ave_Z/(E+511.0), 2.0)*(alpha*(1.0+alpha));
    double lambda=ave_M/(AVOGADRO_CONSTANT*density)*area_inv;
    R=rng.uniform();
    return -lambda*log(R);
}

double EBSD_MC::compute_energyloss(double E, double step){
    double J=(9.76*ave_Z+58.5/pow(ave_Z, 0.19))*1.0E-3/1.166;
    J=1.0/J;
    double dE=78500.0*density*ave_Z/ave_M/E*log(E*J+1.0)*step;
    return dE;
}

void EBSD_MC::run(int seed) 
{
    Emin=Ehistmin-Ebin/2.0;
    numpx=nump; numpy=nump;
    numEbin=int((EkeV-Ehistmin)/Ebin)+1; 
    numzbin=int(depthmax/depthstep)+1;
    num_E=nump; num_z=(nump-1)/10+1;
    size_E[0]=numEbin; size_E[1]=size_E[2]=num_E; 
    callocate_3d(&accum_E, numEbin, num_E, num_E, 0);
    size_z[0]=numEbin; size_z[1]=numzbin; size_z[2]=size_z[3]=num_z; 
    callocate_4d(&accum_z, numEbin, numzbin, num_z, num_z, 0);
    rng.seed(seed);
    count_E=0; count_z=0;
    for(int ie=0;ie<num_e;ie++){
        if(ie%10000==0){
            printf("Completed electrons %d;\n", ie);
            printf("Back-scattered electrons hits = %d.\n", count_E);
        }

        //Set the initial energy, coordinate, and direction (cosine) for this incident electron
        double E0=EkeV;
        double xyz[3]={0.0, 0.0, 0.0};
        double cdir[3]={cos((90.0-sigma)*DEG_TO_RAD), 0.0, -sin((90.0-sigma)*DEG_TO_RAD)};

        //Set the free path and advance the coordinate
        double alpha=compute_alphainfreepath(E0);
        double step=compute_freepath(E0, alpha);
        for(int i=0;i<3;i++){
            xyz[i]+=step*CEN_TO_ANG*cdir[i];
        }

        int jt=0;
        int nx=(numpx-1)/2;
        while(jt<SCATTERING_EVENT_NUMBER){
            //Advance the energy
            double dE=compute_energyloss(E0, step);
            E0=E0-dE;
            if(E0<0.0) break; //Exit if the energy becomes low enough
              
            //Find the deflection and azimuthal angle by the scattering event
            R=rng.uniform();
            double cphi=1.0-2.0*alpha*R/(1.0+alpha-R);
            double sphi=sin(acos(cphi));
            R=rng.uniform();
            double cpsi=cos(2*PI*R);
            double spsi=sin(2*PI*R);

            //Advance the direction (cosine)
            double cdirt[3];
            if(abs(cdir[2])>ONE){
                cdirt[0]=sphi*cpsi; cdirt[1]=sphi*spsi; cdirt[2]=(cdir[2]/abs(cdir[2]))*cphi;
            }
            else{
                double dsq=sqrt(1.0-cdir[2]*cdir[2]);
                double dsqi=1.0/dsq;
                cdirt[0]=sphi*(cdir[0]*cdir[2]*cpsi-cdir[1]*spsi)*dsqi+cdir[0]*cphi;
                cdirt[1]=sphi*(cdir[1]*cdir[2]*cpsi+cdir[0]*spsi)*dsqi+cdir[1]*cphi;
                cdirt[2]=-sphi*cpsi*dsq+cdir[2]*cphi;
            }
            double lent=sqrt(cdirt[0]*cdirt[0]+cdirt[1]*cdirt[1]+cdirt[2]*cdirt[2]);
            for(int i=0;i<3;i++){
                cdirt[i]=cdirt[i]/lent;
            }

            //Advance the free path and the coordinate
            alpha=compute_alphainfreepath(E0);
            step=compute_freepath(E0, alpha);
            double xyzt[3];
            for(int i=0;i<3;i++){
                xyzt[i]=xyz[i]+step*CEN_TO_ANG*cdirt[i];
            }

            double zmax=xyzt[1]*tan(omega*DEG_TO_RAD);
            //Determine whether the electron exit the crystal
            if(xyzt[2]>zmax){
                double xyL[2];
                int ierr;
                compute_square_Lambert(xyL, ierr, cdirt);  //Coordinate in the Lambert projection
                xyL[0]=xyL[0]*nx; xyL[1]=xyL[1]*nx;
                int    pxyL[2]={int(round(xyL[1]))+nx, int(round(-xyL[0]))+nx};     //Find the nearest pixel taking the reversal of coordinate frame (x,y)->(y,-x) into account
                
                if(pxyL[0]<nump&&pxyL[1]<=nump){
                    if(E0>Emin){
                        int iE=round((E0-Ehistmin)/Ebin);
                        //first add this electron to the correct exit distance vs. energy bin (coarser than the angular plot)
                        double disz=abs(xyz[2]/cdirt[2]);   //distance from last scattering point to surface along trajectory
                        int    iz=round(disz*0.1/depthstep);
                        if((iz>0)&&(iz<numzbin)){
                            int px=round(pxyL[0]/10.0);
                            int py=round(pxyL[1]/10.0);
                            accum_z[iE][iz][px][py]+=1;
                            count_z+=1;
                        }
                        accum_E[iE][pxyL[0]][pxyL[1]]+=1;
                        count_E+=1;
                    }
                }
                break;
            }
            for(int i=0;i<3;i++){
                cdir[i]=cdirt[i]; xyz[i]=xyzt[i];
            }
            jt+=1;
        }
    }
}