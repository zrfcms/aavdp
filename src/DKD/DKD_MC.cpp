#include "DKD_MC.h"

DKD_MC::DKD_MC(const char *hdf5_path, double omega, double sigma, double Emax, double Emin, double Ebin, double zmax, double zstep, int num_e, int nump)
{
    this->omega=omega; this->sigma=sigma;
    this->EkeV=Emax; this->Ehistmin=Emin; this->Ebin=Ebin;
    this->depthmax=zmax; this->depthstep=zstep;
    this->num_e=num_e; this->nump=nump;

    this->numEbin=int((EkeV-Ehistmin)/Ebin)+1; 
    this->numzbin=int(depthmax/depthstep)+1;
    this->numz=(nump-1)/10+1;

    CELL cell(hdf5_path);
    this->ave_M=cell.ave_M; 
    this->ave_Z=cell.ave_Z;
    this->density=cell.density;

    callocate_3d(&this->accum_E, this->numEbin, this->nump, this->nump, 0);
    callocate_4d(&this->accum_z, this->numEbin, this->numzbin, this->numz, this->numz, 0);
    int count_E=0;
    for(int i=0;i<this->multiplier;i++){
        compute(count_E, this->prime_seed+i+1);
    }
}

void DKD_MC::compute(int &count, int seed)
{
    printf("[INFO] Starting computation of spatial and energy distributions of back-scattered electrons for seed %d...\n", seed);
    int imp=(nump-1)/2, imz=(numz-1)/2;
    double Emin=Ehistmin-Ebin/2.0;
    double dirx=cos((90.0-sigma)*DEG_TO_RAD), dirz=-sin((90.0-sigma)*DEG_TO_RAD);
    double tano=tan(omega*DEG_TO_RAD);
    rng.seed(seed);
    int count_e=0;
    for(int ie=0;ie<num_e;ie++){
        //Set the initial energy, coordinate, and direction (cosine) for this incident electron
        double E0=EkeV;
        double dir[3]={dirx, 0.0, dirz};
        double xyz[3]={0.0, 0.0, 0.0};
        double alpha=0.0, step=0.0;
        update_free_path(step, alpha, E0);
        update_incident_coordinate(xyz, dir, step);

        int jt=0; double R;
        while(jt<SCATTERING_EVENT_NUMBER){
            update_incident_energy(E0, step);
            if(E0<0.0) break; //Exit if the energy becomes low enough
            double dirt[3], xyzt[3];
            vector_copy(dirt, dir); vector_copy(xyzt, xyz);
            update_incident_direction(dirt, alpha);
            update_free_path(step, alpha, E0);
            update_incident_coordinate(xyzt, dirt, step);

            double zmax=xyzt[1]*tano;
            if(xyzt[2]>zmax){//Determine whether the electron exit the crystal
                double dxy[2]; int ierr;
                compute_square_Lambert(dxy, ierr, dirt);  //Coordinate in the Lambert projection
                dxy[0]*=(double)imp; dxy[1]*=(double)imp;
                int ipx=int(round(dxy[1])), ipy=int(round(-dxy[0]));
                if(fabs(ipx)<=imp&&fabs(ipy)<=imp){
                    if(E0>Emin){
                        int iE=round((E0-Ehistmin)/Ebin);
                        int iz=round(0.1*fabs(xyz[2]/dirt[2])/depthstep);
                        accum_E[iE][ipx+imp][ipy+imp]++;
                        count++;
                        if((iz>=0)&&(iz<numzbin)){
                            ipx=round(ipx/10.0);
                            ipy=round(ipy/10.0);
                            accum_z[iE][iz][ipx+imz][ipy+imz]++;
                        }
                    }
                }
                break;
            }
            vector_copy(dir, dirt); vector_copy(xyz, xyzt);
            jt++;
        }
        count_e++;
        if(0==count_e%1000){
            printf("[INFO] Completed incident electrons %d of %d\n", count_e, num_e);
            printf("[INFO] Back-scattered electrons hits = %d\n", count);
        }
    }
    printf("[INFO] Ending computation of spatial and energy distributions of back-scattered electrons for seed %d\n", seed);
}

//Set the free path for the energy and scale it by a random number
void DKD_MC::update_free_path(double &step, double &alpha, double E)
{
    alpha=3.4e-3*pow(ave_Z, 2.0/3.0)/E;
    double R=rng.uniform();
    double energy=E*(E+1024.0)/ave_Z/(E+511.0);
    double area_inv=CONST_IN_AREA_INV*energy*energy*(alpha*(1.0+alpha));
    double lambda=ave_M/AVOGADRO_CONSTANT/density*area_inv;
    step=-lambda*log(R);
}

void DKD_MC::update_incident_energy(double &E, double step)
{
    double J=(9.76*ave_Z+58.5/pow(ave_Z, 0.19))*1.0e-3/1.166; J=1.0/J;
    double dE=78500.0*density*ave_Z/ave_M*log(E*J+1.0)/E; dE=dE*step;
    E=E-dE;
}

//Find the deflection and azimuthal angle by the scattering event, Advance the direction
void DKD_MC::update_incident_direction(double dir[3], double alpha)
{
    double R=rng.uniform();
    double cphi=1.0-2.0*alpha*R/(1.0+alpha-R);
    double sphi=sin(acos(cphi));
    R=rng.uniform();
    double cpsi=cos(TWO_PI*R), spsi=sin(TWO_PI*R);
    if(fabs(dir[2])>0.99999){
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
    vector_normalize(dir, dir);
}

void DKD_MC::update_incident_coordinate(double xyz[3], double dir[3], double step)
{
    vector_constant(xyz, CEN_TO_ANG*step, dir);
}

DKD_MC::DKD_MC(const char *hdf5_path)
{
    size_t size1, size2, size3, size4;
    double ***aaccum_E, ****aacum_z;
    HDF5 hdf;
    hdf.open(hdf5_path);
    hdf.read("/MonteCarlo/SampleRotateAngle", omega);
    hdf.read("/MonteCarlo/SampleTiltAngle", sigma);
    hdf.read("/MonteCarlo/AcceleratingVoltage", EkeV);
    hdf.read("/MonteCarlo/MinimumEnergy", Ehistmin);
    hdf.read("/MonteCarlo/EnergyBinSize", Ebin);
    hdf.read("/MonteCarlo/MaximumDepth", depthmax);
    hdf.read("/MonteCarlo/DepthStepSize", depthstep);
    hdf.read("/MonteCarlo/ElectronNumber", num_e);
    hdf.read("/MonteCarlo/EnergyPixelNumber", nump);
    hdf.read("/MonteCarlo/EnergyBinNumber", numEbin);
    hdf.read("/MonteCarlo/DepthStepNumber", numzbin);
    hdf.read("/MonteCarlo/DepthPixelNumber", numz);
    hdf.read("/MonteCarlo/Multiplier", multiplier);
    hdf.read_array_3d("/MonteCarlo/EnergyProjectionArray", &aaccum_E, size1, size2, size3);
    hdf.read_array_4d("/MonteCarlo/DepthProjectionArray", &aacum_z, size1, size2, size3, size4);
    // hdf.read("/NMLparameters/MCCLNameList/omega", omega);
    // hdf.read("/NMLparameters/MCCLNameList/sig", sigma);
    // hdf.read("/NMLparameters/MCCLNameList/EkeV", EkeV);
    // hdf.read("/NMLparameters/MCCLNameList/Ehistmin", Ehistmin);
    // hdf.read("/NMLparameters/MCCLNameList/Ebinsize", Ebin);
    // hdf.read("/NMLparameters/MCCLNameList/depthmax", depthmax);
    // hdf.read("/NMLparameters/MCCLNameList/depthstep", depthstep);
    // hdf.read("/EMData/MCOpenCL/totnum_el", num_e);
    // hdf.read("/NMLparameters/MCCLNameList/multiplier", multiplier);
    // hdf.read("/NMLparameters/MCCLNameList/numsx", nump);
    // hdf.read("/EMData/MCOpenCL/numEbins", numEbin);
    // hdf.read("/EMData/MCOpenCL/numzbins", numzbin);
    // hdf.read_array_3d("/EMData/MCOpenCL/accum_e", &aaccum_E, size1, size2, size3);
    // hdf.read_array_4d("/EMData/MCOpenCL/accum_z", &aacum_z, size1, size2, size3, size4);
    mallocate_3d(&accum_E, numEbin, nump, nump);
    mallocate_4d(&accum_z, numEbin, numzbin, numz, numz);
    for(int i=0;i<numEbin;i++){
        for(int j=0;j<nump;j++){
            for(int k=0;k<nump;k++){
                accum_E[i][j][k]=aaccum_E[j][k][i];
            }
        }
    }
    for(int i=0;i<numEbin;i++){
        for(int j=0;j<numzbin;j++){
            for(int k=0;k<numz;k++){
                for(int n=0;n<numz;n++){
                    accum_z[i][j][k][n]=aacum_z[k][n][j][i];
                }
            }
        }
    }
    deallocate_3d(aaccum_E, nump, nump);
    deallocate_4d(aacum_z, numz, numz, numEbin);
}

DKD_MC::~DKD_MC()
{
    deallocate_3d(accum_E, numEbin, nump);
    deallocate_4d(accum_z, numEbin, numzbin, numz);
}

void DKD_MC::hdf5(const char *hdf5_path)
{
    int ***aaccum_E, ****aacum_z;
    mallocate_3d(&aaccum_E, nump, nump, numEbin);
    mallocate_4d(&aacum_z, numz, numz, numzbin, numEbin);
    for(int i=0;i<numEbin;i++){
        for(int j=0;j<nump;j++){
            for(int k=0;k<nump;k++){
                aaccum_E[j][k][i]=accum_E[i][j][k];
            }
        }
    }
    for(int i=0;i<numEbin;i++){
        for(int j=0;j<numzbin;j++){
            for(int k=0;k<numz;k++){
                for(int n=0;n<numz;n++){
                    aacum_z[k][n][j][i]=accum_z[i][j][k][n];
                }
            }
        }
    }
    HDF5 hdf;
    hdf.open(hdf5_path);
    hdf.write_group("/MonteCarlo");
    hdf.write("/MonteCarlo/SampleRotateAngle", omega);
    hdf.write("/MonteCarlo/SampleTiltAngle", sigma);
    hdf.write("/MonteCarlo/AcceleratingVoltage", EkeV);
    hdf.write("/MonteCarlo/MinimumEnergy", Ehistmin);
    hdf.write("/MonteCarlo/EnergyBinSize", Ebin);
    hdf.write("/MonteCarlo/MaximumDepth", depthmax);
    hdf.write("/MonteCarlo/DepthStepSize", depthstep);
    hdf.write("/MonteCarlo/ElectronNumber", num_e);
    hdf.write("/MonteCarlo/Multiplier", multiplier);
    hdf.write("/MonteCarlo/EnergyPixelNumber", nump);
    hdf.write("/MonteCarlo/EnergyBinNumber", numEbin);
    hdf.write("/MonteCarlo/DepthStepNumber", numzbin);
    hdf.write("/MonteCarlo/DepthPixelNumber", numz);
    hdf.write_array_3d("/MonteCarlo/EnergyProjectionArray", aaccum_E, nump, nump, numEbin);
    hdf.write_array_4d("/MonteCarlo/DepthProjectionArray", aacum_z, numz, numz, numzbin, numEbin);
    hdf.close();
    deallocate_3d(aaccum_E, nump, nump);
    deallocate_4d(aacum_z, numz, numz, numzbin);
    printf("Monte-carlo data stored in %s.\n", hdf5_path);
}

void DKD_MC::img(const char *img_path, double dimension, int resolution)
{
    char name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
    split_path(name, ext, img_path);
    if(0!=strcmp(ext, ".png")){
        printf("[ERROR] Unrecognized extension %s.", ext);
        exit(EXIT_FAILURE);
    }
    char png_path[PATH_CHAR_NUMBER]; 
    char exts[3][EXT_CHAR_NUMBER]; strcpy(exts[0], ".DKD_MC."); strcpy(exts[2], ext);
    for(int i=0;i<numEbin;i++){
        int_to_str(exts[1], i);
        strcpy(png_path, name); merge_path(png_path, exts, 3);
        double E=Ehistmin+i*Ebin;
        for(int j=0;j<nump;j++){
            for(int k=0;k<nump;k++){
                if(accum_E[i][j][k]>1000){
                    printf("%d %d\n", j, k);
                }
            }
        }
        image_array(png_path, accum_E[i], nump, nump, dimension, dimension, resolution);
        printf("Image data for energy %.5f stored in %s.\n", E, png_path);
    }
}


void DKD_MC::set_energies_and_depths()
{
    mallocate(&Ebins, numEbin);
    mallocate(&depths, numEbin);
    for(int i=0;i<numEbin;i++){
        Ebins[i]=Ehistmin+(double)i*Ebin;
    }
    izmax=0;
    for(int iE=0;iE<numEbin;iE++){
        for(int ix=0;ix<numz;ix++){
            for(int iy=0;iy<numz;iy++){
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