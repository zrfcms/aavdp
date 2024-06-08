#include "DKD_MC.h"

RNG::RNG(const char *rand_path, int nseed)
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
        exit(EXIT_FAILURE);
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

DKD_MC::DKD_MC(const char *hdf5_path, double omega, double sigma, double Emax, double Emin, double Ebin, double zmax, double zstep, int num_e, int nump)
{
    this->omega=omega; this->sigma=sigma;
    this->EkeV=Emax; this->Ehistmin=Emin; this->Ebin=Ebin;
    this->depthmax=zmax; this->depthstep=zstep;
    this->num_e=num_e; this->nump=nump;

    this->numEbin=int((EkeV-Ehistmin)/Ebin)+1;
    callocate(&this->Ebins, this->numEbin, 0.0);
    for(int i=0;i<this->numEbin;i++){
        this->Ebins[i]=this->Ehistmin+(double)i*this->Ebin;
    }
    this->numzbin=int(depthmax/depthstep)+1;
    this->numpz=(nump-1)/10+1;

    CELL cell(hdf5_path);
    this->ave_M=cell.ave_M; 
    this->ave_Z=cell.ave_Z;
    this->density=cell.density;

    callocate_3d(&this->accum_E, this->numEbin, this->nump, this->nump, 0);
    callocate_4d(&this->accum_z, this->numEbin, this->numzbin, this->numpz, this->numpz, 0);
    rng=new RNG(rand_path, nseed);
    int ebin=int(ceil(double(num_e)/double(nseed)));
    int count_bse=0, count_e=0;
    for(int i=0;i<this->nseed;i++){
        if((num_e-i*ebin)<ebin) ebin=num_e-i*ebin;
        compute(count_e, ebin, i);
    }
    printf("%d\n", count_e);
    count_bse=0;
    for(int iE=0;iE<this->numEbin;iE++){
        int nmax=0, nmin=1e8, x, y;
        for(int ix=0;ix<this->nump;ix++){
            for(int iy=0;iy<this->nump;iy++){
                if(nmax<this->accum_E[iE][ix][iy]) nmax=this->accum_E[iE][ix][iy];
                if(nmin>this->accum_E[iE][ix][iy]) nmin=this->accum_E[iE][ix][iy];
                count_bse+=this->accum_E[iE][ix][iy];
            }
        }
        printf("[INFO] Range of back-scattered electron %d distribution at energy %.5f: %d, %d\n", count_bse, this->Ebins[iE], nmin, nmax);
    }
}

void DKD_MC::compute(int &count_e, int ne, int id)
{
    rng->seed(id);
    int imp=nump/2, imz=numpz/2;
    double dir_cos[3];    
    int counter1, counter2;
    double c_new[3], r_new[3];
    double E_new, alpha, de_ds, phi, psi, mfp,sig_eNA,step, dsq, dsqi, absc0z;
    double sig=sigma;
    double z=ave_Z, A=ave_M, rho=density;
    
    double J;    // refer to Monte Carlo simulation for Electron Microscopy and Microanalysis, David C. Joy
    J = (9.76*z + 58.5*pow(z,-0.19))*1e-3;
    
    // double r0[3] = {0.0, 0.0, 0.0};
    // double c0[3] = {cos(omega)*sin(sig), sin(omega)*sin(sig), cos(sig)};
    double escape_depth;


// Setting all values to -10. Any value other than -10 will denote a backscattered electron with the x and y component of the Lambert Projection

    int num_el=ne;
    double *Lamx, *Lamy, *depth, *energy;
    callocate(&Lamx, num_el, -10.0);
    callocate(&Lamy, num_el, -10.0);
    callocate(&depth, num_el, 10.0);
    callocate(&energy, num_el, 0.0);
	// for (int i = 0; i < num_el; ++i){
	// 	Lamx[num_el*id + i] = -10.0f;
	// 	Lamy[num_el*id + i] = -10.0f;
    //     	depth[num_el*id + i] = 10.0f;
    //     	energy[num_el*id + i] = 0.0f;
	// }

    
    
    double E=EkeV;
    double rand;
    double rtemp[3];
    for (int i = 0; i < num_el; ++i){
        //rand_seed = rando();
        //seed = rand_seed;
        // z11 = seeds[4*id];
        // z22 = seeds[4*id + 1];
        // z33 = seeds[4*id + 2];
        // z44 = seeds[4*id + 3];
        // retrnd = lfsr113_Bits(z11,z22,z33,z44);
        // seeds[4*id] = retrnd.z1;
        // seeds[4*id + 1] = retrnd.z2;
        // seeds[4*id + 2] = retrnd.z3;
        // seeds[4*id + 3] = retrnd.z4;
        // rand = fabs(retrnd.rand/RAND_MAX); //some random no. generator in gpu
        rand=rng->random();
        double r0[3] = {0.0, 0.0, 0.0};
        double c0[3] = {cos(omega*DEG_TO_RAD)*sin(sig*DEG_TO_RAD), sin(omega*DEG_TO_RAD)*sin(sig*DEG_TO_RAD), cos(sig*DEG_TO_RAD)};
        E_new = E;
        vector_copy(c_new, c0);
        escape_depth = 0.0;

        // alpha = (3.4E-3)*pow(z,0.66667)/E_new;
        // sig_eNA = (5.21 * 602.2)*z*z/E_new/E_new*4.0*PI/alpha/(1.0+alpha)*pow(E_new+511.0,2.0)/pow(E+1022.0,2.0);
        // mfp = A * 1.0e7/(rho*sig_eNA);
        // step = -mfp * log(rand);
        update_free_path(step, alpha, E_new, E);

        // vector_constant(rtemp, step, c_new);
        // vector_plus(r_new, r0, rtemp);
        // vector_copy(r0, r_new);
        update_incident_coordinate(r0, c_new, step);
        // de_ds = -0.00785*(z/(A*E_new)) * log(1.166*E_new/J + 0.9911);
        // E_new += step*rho*de_ds;
        update_incident_energy(E_new, step);

        counter1 = 0;   // This is used as a counter for the number of monte carlo steps to carry out for each electron. This is due to the lock step nature of the GPU code. We have arbitly set this to a 1000, though this is material dependent
        
        counter2 = 0;   // This counter is used to figure out if the electron has left the sample or not. Again, this is because of the lock step nature. All steps have to be executed on each thread irrespective of the fact that the electron may have actually left the sample
        int steps=300;
        while (counter1 < steps){
// inline code rather than function call
// Taken from book Monte Carlo simulation for Electron Microscopy and Microanalysis, David C. Joy

            // alpha = (3.4e-3)*pow(z,0.66667)/E_new;
            // sig_eNA = (5.21f * 602.2)*z*z/E_new/E_new*4.0*PI/alpha/(1.0+alpha)*pow(E_new+511.0,2.0)/pow(E+1022.0,2.0);
   	        // mfp = A * 1.0e7f/(rho*sig_eNA);

            // // z11 = seeds[4*id];
            // // z22 = seeds[4*id + 1];
            // // z33 = seeds[4*id + 2];
            // // z44 = seeds[4*id + 3];
            // // retrnd = lfsr113_Bits(z11,z22,z33,z44);
            // // seeds[4*id] = retrnd.z1;
            // // seeds[4*id + 1] = retrnd.z2;
            // // seeds[4*id + 2] = retrnd.z3;
            // // seeds[4*id + 3] = retrnd.z4;
            // // rand = fabs(retrnd.rand/RAND_MAX); //some random no. generator in gpu
            // rand=rng->random();
            // step = -mfp * log(rand);
            update_free_path(step, alpha, E_new, E);

// This is the Continuous Slowing Down approximation that we want to get rid of

            //de_ds = -0.00785*(z/(A*E_new)) * log(1.166*E_new/J + 0.9911);

            // z11 = seeds[4*id];
            // z22 = seeds[4*id + 1];
            // z33 = seeds[4*id + 2];
            // z44 = seeds[4*id + 3];
            // retrnd = lfsr113_Bits(z11,z22,z33,z44);
            // seeds[4*id] = retrnd.z1;
            // seeds[4*id + 1] = retrnd.z2;
            // seeds[4*id + 2] = retrnd.z3;
            // seeds[4*id + 3] = retrnd.z4;
            // rand = fabs(retrnd.rand/RAND_MAX);
            // rand=rng->random();
            // phi = acos(1.0 - ((2.0*alpha*rand)/(1.0 + alpha - rand)));

            // z11 = seeds[4*id];
            // z22 = seeds[4*id + 1];
            // z33 = seeds[4*id + 2];
            // z44 = seeds[4*id + 3];
            // retrnd = lfsr113_Bits(z11,z22,z33,z44);
            // seeds[4*id] = retrnd.z1;
            // seeds[4*id + 1] = retrnd.z2;
            // seeds[4*id + 2] = retrnd.z3;
            // seeds[4*id + 3] = retrnd.z4;
            // rand = fabs(retrnd.rand/RAND_MAX);
//             rand=rng->random();
//             psi = 2.0*PI*rand;
            
            
// // new direction cosines of the electrons after scattering event
//             if ((c0[2] >= 0.99999) || (c0[2] <= -0.99999) ){
//                 absc0z = fabs(c0[2]);
//                 c_new[0]=sin(phi) * cos(psi);
//                 c_new[1]=sin(phi) * sin(psi);
//                 c_new[2]=(c0[2]/absc0z)*cos(phi);
//             }else{
//                 dsq = sqrt(1.0-c0[2]*c0[2]);
//                 dsqi = 1.0/dsq;
//                 c_new[0]=sin(phi)*(c0[0]*c0[2]*cos(psi) - c0[1]*sin(psi))*dsqi + c0[0]*cos(phi);
//                 c_new[1]=sin(phi) * (c0[1] * c0[2] * cos(psi) + c0[0] * sin(psi)) * dsqi + c0[1] * cos(phi);
//                 c_new[2]=-sin(phi) * cos(psi) * dsq + c0[2] * cos(phi);
//             }
            update_incident_direction(c_new, alpha);

            if (fabs(c_new[2]) > 1.0e-5){
                escape_depth = r_new[2]/c_new[2];
            }

            // vector_constant(rtemp, step, c_new);
            // vector_plus(r_new, r0, rtemp);
            // vector_copy(r0, r_new);
            update_incident_coordinate(r0, c_new, step);
            vector_copy(c0, c_new);
            //E_new += step*rho*de_ds;
            update_incident_energy(E_new, step);
            if (r0[2] <= 0 && counter2 == 0){
                dir_cos[0] = c0[0];
                dir_cos[1] = c0[1];
                dir_cos[2] = c0[2];
                double ret[2];
                int ierr;
                if(dir_cos[0] != 0.0 && dir_cos[1] != 0.0 &&dir_cos[2] != 0.0){
                    LambertSphereToPlane(ret, dir_cos);
                }
                count_e++;
                Lamx[i]=ret[0]; Lamy[i]=ret[1];
                depth[i]=escape_depth;
                // Lamx[num_el*id + i] = ret.x;
                // Lamy[num_el*id + i] = ret.y;
                depth[i] = escape_depth;
                energy[i] = E_new;
                counter2 = 1;
                int ipx=int(round(ret[1]*imp)), ipy=int(round(-ret[0]*imp));
                if(abs(ipx)>imp||abs(ipy)>imp) printf("%d %d\n", ipx, ipy);
                if(abs(ipx)<=imp&&abs(ipy)<=imp){
                    if(E_new>Ehistmin){
                        int iE=round((E_new-Ehistmin)/Ebin);
                        int iz=round(fabs(escape_depth)/depthstep);
                        accum_E[iE][ipx+imp][ipy+imp]++;
                        //if(ipx+imp==54&&ipy+imp==1) printf("%.5f %.5f %.5f %.5f %.5f\n", dxy[0], dxy[1], dir[0], dir[1], dir[2]);
                        if((iz>=0)&&(iz<numzbin)){
                            ipx=int(round(ipx/10.0));
                            ipy=int(round(ipy/10.0));
                            accum_z[iE][iz][ipx+imz][ipy+imz]++;
                        }
                    }
                }
            }
			
            counter1++ ;
            
        }
        
    }
}

void DKD_MC::compute(int &count_bse, int &count_e, int ne, int id)
{
    //printf("[INFO] Starting computation %d of spatial and energy distributions of back-scattered electrons...\n", id+1);
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
                compute_sphere_Lambert(dxy, ierr, dir);  //Coordinate in the Lambert projection
                if(ierr==1){
                    printf("[ERROR] Zero direction");
                    exit(EXIT_FAILURE);
                }
                dxy[0]*=(double)imp; dxy[1]*=(double)imp;
                int ipx=int(round(dxy[1])), ipy=int(round(-dxy[0]));
                if(abs(ipx)<=imp&&abs(ipy)<=imp){
                    count_e++;
                    if(E>Ehistmin){
                        int iE=round((E-Ehistmin)/Ebin);
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
    // printf("[INFO] Completed incident electrons %d of %d\n", count_e, num_e);
    // printf("[INFO] Back-scattered electrons hits = %d\n", count_bse);
    // printf("[INFO] Ending computation %d of spatial and energy distributions of back-scattered electrons\n", id+1);
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

//Set the free path for the energy and scale it by a random number
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
    hdf.read("/MonteCarlo/DepthPixelNumber", numpz);
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
    // hdf.read("/NMLparameters/MCCLNameList/numsx", nump);
    // hdf.read("/EMData/MCOpenCL/numEbins", numEbin);
    // hdf.read("/EMData/MCOpenCL/numzbins", numzbin);
    // hdf.read_array_3d("/EMData/MCOpenCL/accum_e", &aaccum_E, size1, size2, size3);
    // hdf.read_array_4d("/EMData/MCOpenCL/accum_z", &aacum_z, size1, size2, size3, size4);
    // numpz=nump/10+1;
    mallocate_3d(&accum_E, numEbin, nump, nump);
    mallocate_4d(&accum_z, numEbin, numzbin, numpz, numpz);
    for(int i=0;i<numEbin;i++){
        for(int j=0;j<nump;j++){
            for(int k=0;k<nump;k++){
                accum_E[i][j][k]=aaccum_E[j][k][i];
            }
        }
    }
    for(int i=0;i<numEbin;i++){
        for(int j=0;j<numzbin;j++){
            for(int k=0;k<numpz;k++){
                for(int n=0;n<numpz;n++){
                    accum_z[i][j][k][n]=aacum_z[k][n][j][i];
                }
            }
        }
    }
    deallocate_3d(aaccum_E, nump, nump);
    deallocate_4d(aacum_z, numpz, numpz, numEbin);
}

DKD_MC::~DKD_MC()
{
    deallocate_3d(accum_E, numEbin, nump);
    deallocate_4d(accum_z, numEbin, numzbin, numpz);
}

void DKD_MC::hdf5(const char *hdf5_path)
{
    int ***aaccum_E, ****aacum_z;
    mallocate_3d(&aaccum_E, nump, nump, numEbin);
    mallocate_4d(&aacum_z, numpz, numpz, numzbin, numEbin);
    for(int i=0;i<numEbin;i++){
        for(int j=0;j<nump;j++){
            for(int k=0;k<nump;k++){
                aaccum_E[j][k][i]=accum_E[i][j][k];
            }
        }
    }
    for(int i=0;i<numEbin;i++){
        for(int j=0;j<numzbin;j++){
            for(int k=0;k<numpz;k++){
                for(int n=0;n<numpz;n++){
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
    hdf.write("/MonteCarlo/EnergyPixelNumber", nump);
    hdf.write("/MonteCarlo/EnergyBinNumber", numEbin);
    hdf.write("/MonteCarlo/DepthStepNumber", numzbin);
    hdf.write("/MonteCarlo/DepthPixelNumber", numpz);
    hdf.write_array_3d("/MonteCarlo/EnergyProjectionArray", aaccum_E, nump, nump, numEbin);
    hdf.write_array_4d("/MonteCarlo/DepthProjectionArray", aacum_z, numpz, numpz, numzbin, numEbin);
    hdf.close();
    deallocate_3d(aaccum_E, nump, nump);
    deallocate_4d(aacum_z, numpz, numpz, numzbin);
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

}

// void DKD_MC::compute(int &count, int seed)
// {
//     printf("[INFO] Starting computation of spatial and energy distributions of back-scattered electrons for seed %d...\n", seed);
//     int imp=(nump-1)/2, imz=(numpz-1)/2;
//     double Emin=Ehistmin-Ebin/2.0;
//     double dirx=cos((90.0-sigma)*DEG_TO_RAD), dirz=-sin((90.0-sigma)*DEG_TO_RAD);
//     double tano=tan(omega*DEG_TO_RAD);
//     rng->seed(seed);
//     int count_e=0;
//     for(int ie=0;ie<num_e;ie++){
//         //Set the initial energy, coordinate, and direction (cosine) for this incident electron
//         double E0=EkeV;
//         double dir[3]={dirx, 0.0, dirz};
//         double xyz[3]={0.0, 0.0, 0.0};
//         double alpha=0.0, step=0.0;
//         update_free_path(step, alpha, E0);
//         update_incident_coordinate(xyz, dir, step);

//         int jt=0; double R;
//         while(jt<SCATTERING_EVENT_NUMBER){
//             update_incident_energy(E0, step);
//             if(E0<0.0) break; //Exit if the energy becomes low enough
//             double dirt[3], xyzt[3];
//             vector_copy(dirt, dir); vector_copy(xyzt, xyz);
//             update_incident_direction(dirt, alpha);
//             update_free_path(step, alpha, E0);
//             update_incident_coordinate(xyzt, dirt, step);

//             double zmax=xyzt[1]*tano;
//             if(xyzt[2]>zmax){//Determine whether the electron exit the crystal
//                 double dxy[2]; int ierr;
//                 compute_square_Lambert(dxy, ierr, dirt);  //Coordinate in the Lambert projection
//                 dxy[0]*=(double)imp; dxy[1]*=(double)imp;
//                 int ipx=int(round(dxy[1])), ipy=int(round(-dxy[0]));
//                 if(fabs(ipx)<=imp&&fabs(ipy)<=imp){
//                     if(E0>Emin){
//                         int iE=round((E0-Ehistmin)/Ebin);
//                         int iz=round(0.1*fabs(xyz[2]/dirt[2])/depthstep);
//                         accum_E[iE][ipx+imp][ipy+imp]++;
//                         count++;
//                         if((iz>=0)&&(iz<numzbin)){
//                             ipx=round(ipx/10.0);
//                             ipy=round(ipy/10.0);
//                             accum_z[iE][iz][ipx+imz][ipy+imz]++;
//                         }
//                     }
//                 }
//                 break;
//             }
//             vector_copy(dir, dirt); vector_copy(xyz, xyzt);
//             jt++;
//         }
//         count_e++;
//         if(0==count_e%1000){
//             printf("[INFO] Completed incident electrons %d of %d\n", count_e, num_e);
//             printf("[INFO] Back-scattered electrons hits = %d\n", count);
//         }
//     }
//     printf("[INFO] Ending computation of spatial and energy distributions of back-scattered electrons for seed %d\n", seed);
// }