#include "RDF.h"

RDF::RDF(const char *model_path, double rmax, int nbin, bool is_partial)
{
    printf("[INFO] Starting computation of radial distribution function...\n");
	QB_tools QB;
	QB_init(&QB);
	QB_read_file(&QB, model_path);

	set_volume(QB.mat);
	numrbin=nbin;
	double rbin=rmax/(double)nbin;
	dr=rbin;
	mallocate(&rij, numrbin);
	for(int i=0;i<numrbin;i++){
		rij[i]=((double)i+0.5)*rbin;
	}

	if(is_partial){
		numij=QB.TypeNumber*(QB.TypeNumber+1)/2+1;
	}else{
		numij=1;
	}

	int natom=QB.TotalNumber;
	QB_pbc(&QB, rmax); QB.MCN=QB.TotalNumber;
	QB_network_init(&QB, rmax);
	callocate(&rhoij, numij, 0.0);
	mallocate_2d(&ij, numij, 2);
	mallocate_2d(&gij, numij, numrbin);
	int countij=0;
	ij[countij][0]=0; ij[countij][1]=0; 
	compute(&QB, natom, rmax, rbin, countij); countij++;
	if(countij<numij){
		for(int i=1;i<=QB.TypeNumber;i++){
			for(int j=i;j<=QB.TypeNumber;j++){
				ij[countij][0]=i; ij[countij][1]=j; 
				compute(&QB, natom, rmax, rbin, i, j, countij); countij++;
			}
		}
	}
	QB_network_free(&QB);
	QB_pbc_clean(&QB);
	QB_free_atom(&QB);
	printf("[INFO] Ending computation of radial distribution function\n");
}

RDF::~RDF()
{
	if(0!=numij){
		deallocate(rhoij);
		deallocate(rij);
		deallocate_2d(ij, numij);
		deallocate_2d(gij, numij);
	}
}

void RDF::set_volume(double mat[3][3])
{
	double odd=0.0, even=0.0;
	for(int i=0;i<3;i++){
        odd+=mat[0][i]*mat[1][(i+1)%3]*mat[2][(i+2)%3];
    }
	for(int i=0;i<3;i++){
        even+=mat[0][i]*mat[1][(i+2)%3]*mat[2][(i+1)%3];
    }
	vol=fabs(odd-even);
}

void RDF::compute(QB_tools *QB, int natom, double rmax, double rbin, int pairij_id)
{
    printf("[INFO] Starting computation of radial distribution function of pair *-*...\n");
	int *natombin; callocate(&natombin, numrbin, 0);
	for(int i=0;i<natom;i++){
		QB_list_build(QB, i, rmax);
		double posi[3]={QB->atom[i].x, QB->atom[i].y, QB->atom[i].z};
		for(int j=0;j<QB->neb.num;j++){
			int    id=QB->neb.id[j];
			double posj[3]={QB->atom[id].x, QB->atom[id].y, QB->atom[id].z};
			double dist[3];
			vector_difference(dist, posi, posj);
			id=int(vector_length(dist)/rbin);
			if(id>=0&&id<numrbin){
				natombin[id]++;
			}
		}
		QB_list_clear(QB);		
	}

	rhoij[pairij_id]=(double)natom/vol;
	double constg=THREE_QUARTER_PIINV/pow(rbin, 3)/rhoij[pairij_id]/(double)natom;
	double gmax=0.0, gmin=1.0e8;
	for(int i=0;i<numrbin;i++){
		gij[0][i]=constg*natombin[i]/double(pow(i+1, 3)-pow(i, 3));
		if(gmax<gij[0][i]) gmax=gij[pairij_id][i];
		if(gmin>gij[0][i]) gmin=gij[pairij_id][i];
	}
	deallocate(natombin);
    printf("[INFO] Ending computation of radial distribution function of pair *-*\n");
	printf("[INFO] Density of * atoms [in Angstrom-3]: %.8f\n", rhoij[pairij_id]);
	printf("[INFO] Range of radial distribution function of pair *-*: %.8f %.8f\n", gmin, gmax);
}

void RDF::compute(QB_tools *QB, int natom, double rmax, double rbin, int typei, int typej, int pairij_id)
{
    printf("[INFO] Starting computation of radial distribution function of pair %d-%d...\n", typei, typej);
	int natomi=0, natomj=0;
	int *natombin; callocate(&natombin, numrbin, 0);
	for(int i=0;i<natom;i++){
		if(QB->atom[i].type==typei){
			natomi++;
			QB_list_build(QB, i, rmax);
			double posi[3]={QB->atom[i].x, QB->atom[i].y, QB->atom[i].z};
			for(int j=0;j<QB->neb.num;j++){
				int    id=QB->neb.id[j];
				if(QB->atom[id].type==typej){
					double posj[3]={QB->atom[id].x, QB->atom[id].y, QB->atom[id].z};
					double dist[3];
					vector_difference(dist, posi, posj);
					id=vector_length(dist)/rbin;
					if(id>=0&&id<numrbin){
						natombin[id]++;
					}
				}
			}
			QB_list_clear(QB);
		}else if(QB->atom[i].type==typej){
			natomj++;
		}	
	}
	if(typei==typej) natomj=natomi;

	rhoij[pairij_id]=natomj/vol;
	double constg=THREE_QUARTER_PIINV/pow(rbin, 3)/rhoij[pairij_id]/(double)natomi;
	double gmax=0.0, gmin=1.0e8;
	for(int i=0;i<numrbin;i++){
		gij[pairij_id][i]=constg*natombin[i]/double(pow(i+1, 3)-pow(i, 3));
		if(gmax<gij[pairij_id][i]) gmax=gij[pairij_id][i];
		if(gmin>gij[pairij_id][i]) gmin=gij[pairij_id][i];
	}
	deallocate(natombin);
    printf("[INFO] Ending computation of radial distribution function of pair %d-%d\n", typei, typej);
	printf("[INFO] Density of j atoms of type %d [in Angstrom-3]: %.8f\n", typej, rhoij[pairij_id]);
	printf("[INFO] Range of radial distribution function of pair %d-%d: %.8f %.8f\n", typei, typej, gmin, gmax);
}

void RDF::rdf(const char *rdf_path)
{
	FILE* fp=fopen(rdf_path,"w");
	fprintf(fp,"# Radial Distribution Function (%d data points)\n", numrbin);
	fprintf(fp,"# r\tg(r)\n");
	if(numij>1){
		fprintf(fp,"# *-*\t");
		for(int i=1;i<numij;i++){
			fprintf(fp,"%d-%d\t", ij[i][0], ij[i][1]);
		}
		fprintf(fp,"\n");
	}
	for(int i=0;i<numrbin;i++){
		fprintf(fp,"%lf\t", rij[i]);
		for(int j=0;j<numij;j++){
			fprintf(fp,"%lf\t", gij[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
    printf("[INFO] Visualized data for radial distribution function stored in %s.\n", rdf_path);
}

// void RDF::nml(const char *nml_path)
// {
// 	FILE* f=fopen(nml_path, "w");
// 	fprintf(f,"distance %lf\n", rbin);
// 	fprintf(f,"density0 %lf\n", rho0ij[0]);
// 	fprintf(f,"radial_distribution_function %d\n", numrbin);
// 	for(int i=0;i<numrbin;i++){
// 		fprintf(f,"%lf\t%lf\n", rij[i], gij[0][i]);
// 	}
// 	fclose(f);
//     printf("[INFO] Information for static structure factor stored in %s.\n", nml_path);
// }