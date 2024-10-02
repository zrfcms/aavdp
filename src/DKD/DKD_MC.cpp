#include "DKD_MC.h"

void compute_Lambert_Projection(double xy[2], int &ierr, double xyz[3]) 
{
    ierr=0;
    xy[0]=0.0; xy[1]=0.0;
    if(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]<1.0e-8){
        ierr=1;
    }else{
        double mag=sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
        xyz[0]=xyz[0]/mag; xyz[1]=xyz[1]/mag; xyz[2]=xyz[2]/mag;
        if(fabs(fabs(xyz[2])-1.0)>1.0e-8){
            double q;
            if((fabs(xyz[1])<=fabs(xyz[0]))&&(fabs(xyz[0])>1.0e-8)){
                q=fabs(xyz[0])/xyz[0]*sqrt(2.0*(1.0+xyz[2]));
                xy[0]=q*PI_SQRT_HALF; 
                xy[1]=q*atan(xyz[1]/xyz[0])/PI_SQRT_HALF;
            }else if(fabs(xyz[1])>1.0e-8){
                q=fabs(xyz[1])/xyz[1]*sqrt(2.0*(1.0+xyz[2]));
                xy[0]=q*atan(xyz[0]/xyz[1])/PI_SQRT_HALF; 
                xy[1]=q*PI_SQRT_HALF;
            }
        }
    }
    xy[0]=xy[0]/PI_HALF_SQRT; xy[1]=xy[1]/PI_HALF_SQRT;
}

RNG::RNG(int nseed)
{
    this->nseed=nseed;
    callocate_2d(&this->default_seeds, this->nseed, this->ns, 0);
    int *seeds;
    int seed_num=nseed*ns;
    FILE *fp;
    fp=fopen(rand_path, "rb");
    mallocate(&seeds, seed_num);
    fread(seeds, sizeof(int), seed_num, fp);
    fclose(fp);
    for(int i=0;i<this->nseed;i++){
        for(int j=0;j<this->ns;j++){
            this->default_seeds[i][j]=seeds[i*ns+j];
        }
    }
    deallocate(seeds);
}

RNG::~RNG()
{
    if(nseed>0){
        deallocate_2d(default_seeds, nseed);
    }
}

void RNG::seed(int id)
{
    if(id>=nseed){
        printf("[ERROR] The seed id %d exceeds the allowed range %d-%d", id, 0, nseed-1);
        exit(1);
    }
    for(int i=0;i<ns;i++){
        state[i]=default_seeds[id][i];
    }
}

double RNG::random(){
    int b;
    int z1=state[0], z2=state[1], z3=state[2], z4=state[3];
    b  = ((z1 << 6) ^ z1) >> 13;
    z1 = ((z1 & 4294967294U) << 18) ^ b;
    b  = ((z2 << 2) ^ z2) >> 27;
    z2 = ((z2 & 4294967288U) << 2) ^ b;
    b  = ((z3 << 13) ^ z3) >> 21;
    z3 = ((z3 & 4294967280U) << 7) ^ b;
    b  = ((z4 << 3) ^ z4) >> 12;
    z4 = ((z4 & 4294967168U) << 13) ^ b;
    state[0]=z1; state[1]=z2; state[2]=z3; state[3]=z4;
    return fabs(z1 ^ z2 ^ z3 ^ z4)/2147483647.0;
}

DKD_MC::DKD_MC(CELL *cell, double omega, double sigma, double Emax, double Emin, double Ebin, double zmax, double zstep, int num_e, int nump)
{
    this->ave_M=cell->ave_M; 
    this->ave_Z=cell->ave_Z;
    this->density=cell->density;

    this->omega=omega; this->sigma=sigma;
    this->EkeV=Emax; this->Emin=Emin; this->Ebin=Ebin;
    this->depthmax=zmax; this->depthstep=zstep;
    this->num_e=num_e; this->nump=nump;
    this->numEbin=int((EkeV-Emin)/Ebin);
    callocate(&this->Ebins, this->numEbin, 0.0);
    for(int i=0;i<this->numEbin;i++){
        this->Ebins[i]=this->EkeV-(double)i*this->Ebin;
    }
    this->numzbin=int(depthmax/depthstep)+1;
    this->numpz=(nump-1)/10+1;

    callocate_3d(&this->accum_E, this->numEbin, this->nump, this->nump, 0);
    callocate_4d(&this->accum_z, this->numEbin, this->numzbin, this->numpz, this->numpz, 0);
    rng=new RNG(nbatch);
    int ebin=int(ceil(double(num_e)/double(nbatch)));
    int count_bse=0, count_e=0;
    for(int i=0;i<this->nbatch;i++){
        if((num_e-i*ebin)<ebin) ebin=num_e-i*ebin;
        this->compute(count_bse, count_e, ebin, i);
    }
    this->compute_depth_distribution();
    for(int iE=0;iE<this->numEbin;iE++){
        int nmax=0, nmin=1e8, x, y;
        for(int ix=0;ix<this->nump;ix++){
            for(int iy=0;iy<this->nump;iy++){
                if(nmax<this->accum_E[iE][ix][iy]) nmax=this->accum_E[iE][ix][iy];
                if(nmin>this->accum_E[iE][ix][iy]) nmin=this->accum_E[iE][ix][iy];
            }
        }
        printf("[INFO] Range of back-scattered electron distribution at energy %.5f: %d, %d\n", this->Ebins[iE], nmin, nmax);
    }
}

void DKD_MC::compute(int &count_bse, int &count_e, int ne, int id)
{
    printf("[INFO] Starting computation %d of spatial and energy distributions of back-scattered electrons...\n", id+1);
    double dirx=cos(omega*DEG_TO_RAD)*sin(sigma*DEG_TO_RAD);
    double diry=sin(omega*DEG_TO_RAD)*sin(sigma*DEG_TO_RAD);
    double dirz=cos(sigma*DEG_TO_RAD);
    int imp=nump/2, imz=numpz/2;
    rng->seed(id);
    for(int ie=0;ie<ne;ie++){
        //Set the initial energy, coordinate, and direction (cosine) for this incident electron
        double E0=EkeV, E=EkeV;
        double dir0[3]={dirx, diry, dirz}, dir[3];
        vector_copy(dir, dir0);
        double xyz[3]={0.0, 0.0, 0.0};
        double alpha=0.0, step=0.0;
        update_free_path(step, alpha, E, E0);
        update_incident_coordinate(xyz, dir, step);
        update_incident_energy(E, step);

        int jt=0; double R;
        while(jt<SCATTERING_EVENT_NUMBER){
            update_free_path(step, alpha, E, E0);
            update_incident_direction(dir, alpha);
            update_incident_coordinate(xyz, dir, step);
            update_incident_energy(E, step);
            if(xyz[2]<=0.0){//Determine whether the electron exit the crystal
                double dxy[2]; int ierr;
                compute_Lambert_Projection(dxy, ierr, dir);  //Coordinate in the Lambert projection
                if(ierr==1){
                    printf("[ERROR] Zero electron direction in the Monte Carlo simulation");
                    exit(1);
                }
                dxy[0]*=(double)imp; dxy[1]*=(double)imp;
                int ipx=int(round(dxy[1])), ipy=int(round(-dxy[0]));
                if(abs(ipx)<=imp&&abs(ipy)<=imp){
                    if(E>Emin){
                        int iE=floor((E-Emin)/Ebin);
                        int iz=round(fabs(xyz[2]/dir[2])/depthstep);
                        count_bse++;
                        accum_E[iE][ipx+imp][ipy+imp]++;
                        if((iz>=0)&&(iz<numzbin)){
                            ipx=int(round(ipx/10.0));
                            ipy=int(round(ipy/10.0));
                            accum_z[iE][iz][ipx+imz][ipy+imz]++;
                        }
                    }
                }
                break;
            }
            jt++;
        }
    }
    count_e+=ne;
    printf("[INFO] Completed incident electrons %d of %d\n", count_e, num_e);
    printf("[INFO] Back-scattered electrons hits = %d\n", count_bse);
    printf("[INFO] Ending computation %d of spatial and energy distributions of back-scattered electrons\n", id+1);
}

void DKD_MC::compute_depth_distribution()
{
    callocate(&depths, numEbin, 0.0);
    izmax=0;
    for(int iE=0;iE<numEbin;iE++){
        for(int ix=0;ix<numpz;ix++){
            for(int iy=0;iy<numpz;iy++){
                int sum_z=0;
                for(int i=0;i<numzbin;i++){
                    sum_z+=accum_z[iE][i][ix][iy];
                }
                int iz=1, isum_z=accum_z[iE][0][ix][iy];
                while(isum_z<0.99*double(sum_z)){
                    isum_z+=accum_z[iE][iz][ix][iy];
                    iz++;
                }
                if(iz>izmax) izmax=iz;
            }
        }
        depths[iE]=double(izmax)*depthstep;
    }
}

void DKD_MC::update_free_path(double &step, double &alpha, double E, double E0)
{
    alpha=3.4e-3*pow(ave_Z, 0.66666667)/E;
    double R=rng->random();
    double energy=E*(E0+1022.0)/ave_Z/(E+511.0);
    double area_inv=CONST_IN_AREA_INV*energy*energy*(alpha*(1.0+alpha));
    double lambda=1.0e7*ave_M/AVOGADRO_CONSTANT/density*area_inv;
    step=-lambda*log(R);
}

void DKD_MC::update_incident_coordinate(double xyz[3], double dir[3], double step)
{
    double sdir[3];
    vector_constant(sdir, step, dir);
    vector_plus(xyz, xyz, sdir);
}

void DKD_MC::update_incident_energy(double &E, double step)
{
    double J=(9.76*ave_Z+58.5*pow(ave_Z, -0.19))*1.0e-3;
    double dE=0.00785*density*ave_Z/ave_M/E*log(1.166*E/J+0.9911); dE=dE*step;
    E=E-dE;
}

//Find the deflection and azimuthal angle by the scattering event, Advance the direction
void DKD_MC::update_incident_direction(double dir[3], double alpha)
{
    double R=rng->random();
    double cphi=1.0-2.0*alpha*R/(1.0+alpha-R);
    double sphi=sin(acos(cphi));
    R=rng->random();
    double cpsi=cos(TWO_PI*R), spsi=sin(TWO_PI*R);
    if(fabs(dir[2])>=0.99999){
        dir[0]=sphi*cpsi; dir[1]=sphi*spsi; dir[2]=(dir[2]/fabs(dir[2]))*cphi;
    }
    else{
        double cdir[3]; vector_copy(cdir, dir);
        double dsq=sqrt(1.0-cdir[2]*cdir[2]);
        double dsqi=1.0/dsq;
        dir[0]=sphi*(cdir[0]*cdir[2]*cpsi-cdir[1]*spsi)*dsqi+cdir[0]*cphi;
        dir[1]=sphi*(cdir[1]*cdir[2]*cpsi+cdir[0]*spsi)*dsqi+cdir[1]*cphi;
        dir[2]=-sphi*cpsi*dsq+cdir[2]*cphi;
    }
}

// DKD_MC::DKD_MC(const char *hdf5_path)
// {
//     size_t size1, size2, size3, size4;
//     double ***aaccum_E, ****aacum_z;
//     HDF5 hdf;
//     hdf.open(hdf5_path);
//     hdf.read("/MonteCarlo/SampleRotateAngle", omega);
//     hdf.read("/MonteCarlo/SampleTiltAngle", sigma);
//     hdf.read("/MonteCarlo/AcceleratingVoltage", EkeV);
//     hdf.read("/MonteCarlo/MinimumEnergy", Emin);
//     hdf.read("/MonteCarlo/EnergyBinSize", Ebin);
//     hdf.read("/MonteCarlo/MaximumDepth", depthmax);
//     hdf.read("/MonteCarlo/DepthStepSize", depthstep);
//     hdf.read("/MonteCarlo/ElectronNumber", num_e);
//     hdf.read("/MonteCarlo/EnergyPixelNumber", nump);
//     hdf.read("/MonteCarlo/EnergyBinNumber", numEbin);
//     hdf.read("/MonteCarlo/DepthStepNumber", numzbin);
//     hdf.read("/MonteCarlo/DepthPixelNumber", numpz);
//     hdf.read("/MonteCarlo/MaximumDepthIndex", izmax);
//     hdf.read_array("/MonteCarlo/EnergyDistribution", &Ebins, size1);
//     hdf.read_array("/MonteCarlo/DepthDistribution", &depths, size1);
//     hdf.read_array_3d("/MonteCarlo/EnergyProjectionArray", &aaccum_E, size1, size2, size3);
//     hdf.read_array_4d("/MonteCarlo/DepthProjectionArray", &aacum_z, size1, size2, size3, size4);
//     // hdf.read("/NMLparameters/MCCLNameList/omega", omega);
//     // hdf.read("/NMLparameters/MCCLNameList/sig", sigma);
//     // hdf.read("/NMLparameters/MCCLNameList/EkeV", EkeV);
//     // hdf.read("/NMLparameters/MCCLNameList/Emin", Emin);
//     // hdf.read("/NMLparameters/MCCLNameList/Ebinsize", Ebin);
//     // hdf.read("/NMLparameters/MCCLNameList/depthmax", depthmax);
//     // hdf.read("/NMLparameters/MCCLNameList/depthstep", depthstep);
//     // hdf.read("/EMData/MCOpenCL/totnum_el", num_e);
//     // hdf.read("/NMLparameters/MCCLNameList/numsx", nump);
//     // hdf.read("/EMData/MCOpenCL/numEbins", numEbin);
//     // hdf.read("/EMData/MCOpenCL/numzbins", numzbin);
//     // hdf.read_array_3d("/EMData/MCOpenCL/accum_e", &aaccum_E, size1, size2, size3);
//     // hdf.read_array_4d("/EMData/MCOpenCL/accum_z", &aacum_z, size1, size2, size3, size4);
//     // numpz=nump/10+1;
//     mallocate_3d(&accum_E, numEbin, nump, nump);
//     mallocate_4d(&accum_z, numEbin, numzbin, numpz, numpz);
//     for(int i=0;i<numEbin;i++){
//         for(int j=0;j<nump;j++){
//             for(int k=0;k<nump;k++){
//                 accum_E[i][j][k]=aaccum_E[j][k][i];
//             }
//         }
//     }
//     for(int i=0;i<numEbin;i++){
//         for(int j=0;j<numzbin;j++){
//             for(int k=0;k<numpz;k++){
//                 for(int n=0;n<numpz;n++){
//                     accum_z[i][j][k][n]=aacum_z[k][n][j][i];
//                 }
//             }
//         }
//     }
//     deallocate_3d(aaccum_E, nump, nump);
//     deallocate_4d(aacum_z, numpz, numpz, numEbin);
// }

DKD_MC::~DKD_MC()
{
    deallocate_3d(accum_E, numEbin, nump);
    deallocate_4d(accum_z, numEbin, numzbin, numpz);
}

// void DKD_MC::hdf5(const char *hdf5_path)
// {
//     int ***aaccum_E, ****aacum_z;
//     mallocate_3d(&aaccum_E, nump, nump, numEbin);
//     mallocate_4d(&aacum_z, numpz, numpz, numzbin, numEbin);
//     for(int i=0;i<numEbin;i++){
//         for(int j=0;j<nump;j++){
//             for(int k=0;k<nump;k++){
//                 aaccum_E[j][k][i]=accum_E[i][j][k];
//             }
//         }
//     }
//     for(int i=0;i<numEbin;i++){
//         for(int j=0;j<numzbin;j++){
//             for(int k=0;k<numpz;k++){
//                 for(int n=0;n<numpz;n++){
//                     aacum_z[k][n][j][i]=accum_z[i][j][k][n];
//                 }
//             }
//         }
//     }
//     HDF5 hdf;
//     hdf.open(hdf5_path);
//     hdf.write_group("/MonteCarlo");
//     hdf.write("/MonteCarlo/SampleRotateAngle", omega);
//     hdf.write("/MonteCarlo/SampleTiltAngle", sigma);
//     hdf.write("/MonteCarlo/AcceleratingVoltage", EkeV);
//     hdf.write("/MonteCarlo/MinimumEnergy", Emin);
//     hdf.write("/MonteCarlo/EnergyBinSize", Ebin);
//     hdf.write("/MonteCarlo/MaximumDepth", depthmax);
//     hdf.write("/MonteCarlo/DepthStepSize", depthstep);
//     hdf.write("/MonteCarlo/ElectronNumber", num_e);
//     hdf.write("/MonteCarlo/EnergyPixelNumber", nump);
//     hdf.write("/MonteCarlo/EnergyBinNumber", numEbin);
//     hdf.write("/MonteCarlo/DepthStepNumber", numzbin);
//     hdf.write("/MonteCarlo/DepthPixelNumber", numpz);
//     hdf.write("/MonteCarlo/MaximumDepthIndex", izmax);
//     hdf.write_array("/MonteCarlo/EnergyDistribution", Ebins, numEbin);
//     hdf.write_array("/MonteCarlo/DepthDistribution", depths, numEbin);
//     hdf.write_array_3d("/MonteCarlo/EnergyProjectionArray", aaccum_E, nump, nump, numEbin);
//     hdf.write_array_4d("/MonteCarlo/DepthProjectionArray", aacum_z, numpz, numpz, numzbin, numEbin);
//     hdf.close();
//     deallocate_3d(aaccum_E, nump, nump);
//     deallocate_4d(aacum_z, numpz, numpz, numzbin);
//     printf("[INFO] Monte-carlo data stored in %s.\n", hdf5_path);
// }

void DKD_MC::img(char *img_path, double dimension, int resolution)
{
    char name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
    split_path(name, ext, img_path);
    if(0!=strcmp(ext, ".png")){
        printf("[ERROR] Unrecognized extension %s.", ext);
        exit(1);
    }
    char png_path[PATH_CHAR_NUMBER]; 
    char exts[3][EXT_CHAR_NUMBER]; strcpy(exts[0], ".DKD_MC."); strcpy(exts[2], ext);
    for(int i=0;i<numEbin;i++){
        int_to_str(exts[1], i);
        strcpy(png_path, name); merge_path(png_path, exts, 3);
        double E=Emin+i*Ebin;
        for(int j=0;j<nump;j++){
            for(int k=0;k<nump;k++){
                if(accum_E[i][j][k]>1000){
                    printf("%d %d\n", j, k);
                }
            }
        }
        image_array(png_path, accum_E[i], nump, nump, dimension, dimension, resolution);
        printf("[INFO] Image data for energy %.5f stored in %s.\n", E, png_path);
    }
}