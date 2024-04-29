#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "main.h"

int main(int argc, char* argv[])
{
    int    i=1;
    if(0==strcmp(argv[i], "--xrd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER];
        no_extension(name, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double lambda=1.54056;
        double c[3]={0.1, 0.1, 0.1};
        double min2Theta=0.0, max2Theta=180.0;
        bool   is_lorentz_flag=false;
        bool   is_spacing_auto=false;
        bool   is_Bragg_flag=false;
        char   xrd_path[PATH_CHAR_NUMBER]; strcpy(xrd_path, name);
        while(i<argc){
            if(0==strcmp(argv[i], "-e")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    strcpy(types[count++], argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-l")){
                i++;
                lambda=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c[0]=be_double(argv[i++]);
                c[1]=be_double(argv[i++]);
                c[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-2t")){
                i++;
                min2Theta=be_double(argv[i++]);
                max2Theta=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-lp")){
                i++;
                is_lorentz_flag=true;
                continue;
            }else if(0==strcmp(argv[i], "-B")){
                i++;
                is_Bragg_flag=true;
                continue;
            }else if(0==strcmp(argv[i], "-auto")){
                i++;
                is_spacing_auto=true;
                continue;
            }
        }
        strcat(xrd_path, ".xrd");
        MODEL model(input_path, types, lambda);
        if(is_Bragg_flag){
            XRD xrd(&model, is_lorentz_flag);
            xrd.xrd(xrd_path);
        }else{
            XRD xrd(&model, c, min2Theta, max2Theta, is_lorentz_flag, is_spacing_auto);
            xrd.xrd(xrd_path);
        }
    }else if(0==strcmp(argv[i], "--rdf")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER];
        no_extension(name, input_path);
        double rmax=5.0;
        int    nbin=200;
        bool   is_partial_flag=false;
        bool   is_nml_flag=false;
        char   rdf_path[PATH_CHAR_NUMBER]; strcpy(rdf_path, name);
        char   nml_path[PATH_CHAR_NUMBER]; strcpy(nml_path, name);
        while(i<argc){
            if(0==strcmp(argv[i], "-r")){
                i++;
                rmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-n")){
                i++;
                nbin=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-partial")){
                i++;
                is_partial_flag=true;
                strcat(rdf_path, ".partial");
                continue;
            }else if(0==strcmp(argv[i], "-nml")){
                i++;
                is_nml_flag=true;
                continue;
            }
        }
        strcat(rdf_path, ".rdf"); strcat(nml_path, ".nml");
        RDF rdf(input_path, rmax, nbin, is_partial_flag);
        rdf.rdf(rdf_path);
        if(is_nml_flag) rdf.nml(nml_path);
    }else if(0==strcmp(argv[i], "--ssf")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER];
        no_extension(name, input_path);
        double qmax=5.0;
        int    nbin=200;
        char   ssf_path[PATH_CHAR_NUMBER]; strcpy(ssf_path, name); 
        while(i<argc){
            if(0==strcmp(argv[i], "-q")){
                i++;
                qmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-n")){
                i++;
                nbin=(int)be_double(argv[i++]);
                continue;
            }
        }
        strcat(ssf_path, ".ssf");
        SSF ssf(input_path, qmax, nbin);
        ssf.ssf(ssf_path);
    }else if(0==strcmp(argv[i], "--sed")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER];
        no_extension(name, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double lambda=0.0251;
        double Kmax=1.70;
        double c[3]={0.1, 0.1, 0.1};
        int    zone[3]={0, 0, 0};
        double thickness=0.1;
        bool   is_zone_existing=false;
        bool   is_kkd_output=false;
        while(i<argc){
            if(0==strcmp(argv[i], "-e")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    strcpy(types[count++], argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-l")){
                i++;
                lambda=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-r")){
                i++;
                Kmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c[0]=be_double(argv[i++]);
                c[1]=be_double(argv[i++]);
                c[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zone[0]=(int)be_double(argv[i++]);
                zone[1]=(int)be_double(argv[i++]);
                zone[2]=(int)be_double(argv[i++]);
                is_zone_existing=true;
                continue;
            }else if(0==strcmp(argv[i], "-t")){
                i++;
                thickness=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-kkd")){
                i++;
                is_kkd_output=true;
                continue;
            }
        }
        char   vtk_path[PATH_CHAR_NUMBER]; strcpy(vtk_path, name); 
        char   *u; to_string(zone[0], u);
        char   *v; to_string(zone[1], v);
        char   *w; to_string(zone[2], w);
        char   uvw[]="."; strcat(uvw, u); strcat(uvw, v); strcat(uvw, w);
        strcat(vtk_path, uvw); strcat(vtk_path, ".vtk");
        EMODEL model(input_path, types, lambda);
        if(is_zone_existing){
            SED sed(&model, c, zone, thickness, Kmax);
            sed.vtk(vtk_path);
        }else{
            SED sed(&model, c, Kmax);
            sed.vtk(vtk_path);
            if(is_kkd_output){
                char   kkd_path[PATH_CHAR_NUMBER]; strcpy(kkd_path, name);
                strcat(kkd_path, uvw); strcat(kkd_path, ".kkd");
                sed.kkd(kkd_path);
            }
        }
    }else if(0==strcmp(argv[i], "--kkd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER];
        no_extension(name, input_path);
        int zone[3]={0};
        double threshold=0.001;
        double dist=3.0;
        int    dpi=512;
        int    width=6;
        int    height=4;
        char   background='w';
        bool   is_zone_existing=false;
        while(i<argc){
            if(0==strcmp(argv[i], "-z")){
                i++;
                zone[0]=(int)be_double(argv[i++]);
                zone[1]=(int)be_double(argv[i++]);
                zone[2]=(int)be_double(argv[i++]);
                is_zone_existing=true;
                continue;
            }else if(0==strcmp(argv[i], "-thr")){
                i++;
                threshold=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-d")){
                i++;
                dist=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-dpi")){
                i++;
                dpi=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-w")){
                i++;
                width=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-h")){
                i++;
                height=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-b")){
                i++;
                background=(argv[i++])[0];
            }
        }
        char   img_path[PATH_CHAR_NUMBER]; strcpy(img_path, name);
        char   *u; to_string(zone[0], u);
        char   *v; to_string(zone[1], v);
        char   *w; to_string(zone[2], w);
        char   uvw[]="."; strcat(uvw, u); strcat(uvw, v); strcat(uvw, w);
        strcat(img_path, uvw); strcat(img_path, ".png");
        if(is_zone_existing){
            KKD kkd(input_path, zone, threshold, dist, dpi, width, height);
            kkd.img(img_path, background);
        }else{
            KKD kkd(input_path, threshold, dist, dpi, width, height);
            kkd.img(img_path, background);
        }
    }
    else{
        printf("[ERROR] Unrecognized command %s.\n", argv[i]);
    }
}
