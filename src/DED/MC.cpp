#include "MC.h"

RNG::RNG(char *seed_path, int nseed)
{
    this->nseed=nseed;
    callocate_2d(&this->default_seeds, this->nseed, this->ns, 0);
    int *seeds;
    int seed_num=this->nseed*this->ns;
    FILE *fp;
    fp=fopen(seed_path, "rb");
    mallocate(&seeds, seed_num);
    fread(seeds, sizeof(int), seed_num, fp);
    fclose(fp);
    for(int i=0;i<this->nseed;i++){
        for(int j=0;j<this->ns;j++){
            this->default_seeds[i][j]=seeds[i*this->ns+j];
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

MC::MC(CELL *cell, char *seed_path, double omega, double sigma, double Emax, double Emin, double zmax, double zstep, int nume, int nump)
{
    ave_M=cell->ave_M; 
    ave_Z=cell->ave_Z;
    density=cell->density;

    numpE=nump;
    if(nump%2==0) numpE=nump+1;
    numzbin=int(zmax/zstep)+1;
    numpz=(numpE-1)/10+1;
    callocate_2d(&accum_E, numpE, numpE, 0);
    callocate_3d(&accum_z, numzbin, numpz, numpz, 0);

    dz=zstep*0.1; EkeV=Emax;
    int ebin=10000;
    int nbatch=int(ceil(double(nume)/double(ebin)));
    int count_bse=0, count_e=0;
    double dir0[3]={cos(omega*DEG_TO_RAD)*sin(sigma*DEG_TO_RAD), sin(omega*DEG_TO_RAD)*sin(sigma*DEG_TO_RAD), cos(sigma*DEG_TO_RAD)}; 
    rng=new RNG(seed_path, nbatch);
    printf("[INFO] Starting computation of spatial and energy distributions of back-scattered electrons...\n");
    for(int i=0;i<nbatch;i++){
        if((nume-i*ebin)<ebin) ebin=nume-i*ebin;
        rng->seed(i);
        compute(count_bse, ebin, dir0, Emax, Emin);
        count_e+=ebin;
        printf("[INFO] Completed incident electrons %d of %d\n", count_e, nume);
        printf("[INFO] Back-scattered electrons hits = %d\n", count_bse);
    }
    izmax=0;
    for(int ix=0;ix<numpz;ix++){
        for(int iy=0;iy<numpz;iy++){
            int sum_z=0;
            for(int i=0;i<numzbin;i++){
                sum_z+=accum_z[i][ix][iy];
            }
            int iz=1, isum_z=accum_z[0][ix][iy];
            while(isum_z<0.99*double(sum_z)){
                isum_z+=accum_z[iz][ix][iy];
                iz++;
            }
            if(iz>izmax) izmax=iz;
        }
    }
    z=double(izmax)*dz;
    callocate(&weight, izmax, 0.0);
    for(int iz=0;iz<izmax;iz++){
        int sum_z=0;
        for(int ix=0;ix<numpz;ix++){
            for(int iy=0;iy<numpz;iy++){
                sum_z+=accum_z[iz][ix][iy];
            }
        }
        weight[iz]=double(sum_z)/double(nume)*exp(TWO_PI*double(iz)*dz/cell->fouri0.sigp);
    }
    printf("[INFO] Ending computation of spatial and energy distributions of back-scattered electrons\n");
}

MC::~MC()
{
    deallocate_2d(accum_E, numpE);
    deallocate_3d(accum_z, numzbin, numpz);
}

void MC::compute(int &count_bse, int ne, double dir0[3], double E0, double Ex)
{
    int imp=numpE/2, imz=numpz/2;
    for(int ie=0;ie<ne;ie++){
        //Set the initial energy, coordinate, and direction (cosine) for this incident electron
        double E=E0;
        double dir[3]; vector_copy(dir, dir0);
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
                double dxyz[3]; vector_normalize(dxyz, dir);
                compute_square_Lambert(dxy, ierr, dxyz);  //Coordinate in the Lambert projection
                if(ierr==1){
                    printf("[ERROR] Zero electron direction in the Monte Carlo simulation");
                    exit(1);
                }
                dxy[0]*=(double)imp; dxy[1]*=(double)imp;
                int ipx=int(round(dxy[1])), ipy=int(round(-dxy[0]));
                if(abs(ipx)<=imp&&abs(ipy)<=imp){
                    if(E>Ex){
                        int iz=round(fabs(xyz[2]/dir[2])/dz);
                        count_bse++;
                        accum_E[ipx+imp][ipy+imp]++;
                        if((iz>=0)&&(iz<numzbin)){
                            ipx=int(round(ipx/10.0));
                            ipy=int(round(ipy/10.0));
                            accum_z[iz][ipx+imz][ipy+imz]++;
                        }
                    }
                }
                break;
            }
            jt++;
        }
    }
}

void MC::update_free_path(double &step, double &alpha, double E, double E0)
{
    alpha=3.4e-3*pow(ave_Z, 0.66666667)/E;
    double R=rng->random();
    double energy=E*(E0+1022.0)/ave_Z/(E+511.0);
    double area_inv=CONST_IN_AREA_INV*energy*energy*(alpha*(1.0+alpha));
    double lambda=1.0e7*ave_M/AVOGADRO_CONSTANT/density*area_inv;
    step=-lambda*log(R);
}

void MC::update_incident_coordinate(double xyz[3], double dir[3], double step)
{
    double sdir[3];
    vector_constant(sdir, step, dir);
    vector_plus(xyz, xyz, sdir);
}

void MC::update_incident_energy(double &E, double step)
{
    double J=(9.76*ave_Z+58.5*pow(ave_Z, -0.19))*1.0e-3;
    double dE=0.00785*density*ave_Z/ave_M/E*log(1.166*E/J+0.9911); dE=dE*step;
    E=E-dE;
}

//Find the deflection and azimuthal angle by the scattering event, Advance the direction
void MC::update_incident_direction(double dir[3], double alpha)
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