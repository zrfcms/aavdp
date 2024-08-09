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
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double lambda=1.54184;
        double min2Theta=0.0, max2Theta=180.0;
        int    lp_type=3;
        double c[3]={0.0};
        bool   is_spacing_auto=true;
        double pk_param=-1;
        double D=400000.0;
        double dt=0.02;
        char   xrd_path[PATH_CHAR_NUMBER]; strcpy(xrd_path, name); strcat(xrd_path, ".xrd");
        char   png_path[PATH_CHAR_NUMBER]; strcpy(png_path, name); strcat(png_path, ".png");
        while(i<argc){
            if(0==strcmp(argv[i], "-e")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    strcpy(types[count++], argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-dw")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    DWs[count++]=be_double(argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-l")){
                i++;
                lambda=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-2t")){
                i++;
                min2Theta=be_double(argv[i++]);
                max2Theta=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-lp")){
                i++;
                lp_type=int(be_double(argv[i++]));
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c[0]=be_double(argv[i++]);
                c[1]=be_double(argv[i++]);
                c[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-auto")){
                i++;
                is_spacing_auto=true;
                if(fabs(c[0])<ZERO_LIMIT&&fabs(c[1])<ZERO_LIMIT&&fabs(c[2])<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=1.0;
                }
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(fabs(c[0])<ZERO_LIMIT&&fabs(c[1])<ZERO_LIMIT&&fabs(c[2])<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=0.1;
                }
                continue;
            }else if(0==strcmp(argv[i], "-pk")){
                i++;
                pk_param=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-D")){
                i++;
                D=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-dt")){
                i++;
                dt=be_double(argv[i++]);
                continue;
            }else{
                i++;
                continue;
            }
        }
        XMODEL model(input_path, types, DWs, lambda);
        XRD xrd(&model, min2Theta, max2Theta, lp_type, c, is_spacing_auto);
        if(pk_param<0.0){
            xrd.xrd(xrd_path, png_path);
        }else{
            xrd.xrd(xrd_path, png_path, pk_param, lambda, D, dt);
        }
    }else if(0==strcmp(argv[i], "--ned")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double lambda=1.54184;
        double min2Theta=0.0, max2Theta=180.0;
        int    lp_type=2;
        double c[3]={0.0};
        bool   is_spacing_auto=true;
        double pk_param=-1;
        double D=400000.0;;
        double dt=0.02;
        char   ned_path[PATH_CHAR_NUMBER]; strcpy(ned_path, name); strcat(ned_path, ".ned");
        char   png_path[PATH_CHAR_NUMBER]; strcpy(png_path, name); strcat(png_path, ".png");
        while(i<argc){
            if(0==strcmp(argv[i], "-e")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    strcpy(types[count++], argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-dw")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    DWs[count++]=be_double(argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-l")){
                i++;
                lambda=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-2t")){
                i++;
                min2Theta=be_double(argv[i++]);
                max2Theta=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-lp")){
                i++;
                lp_type=int(be_double(argv[i++]));
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c[0]=be_double(argv[i++]);
                c[1]=be_double(argv[i++]);
                c[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-auto")){
                i++;
                is_spacing_auto=true;
                if(fabs(c[0])<ZERO_LIMIT&&fabs(c[1])<ZERO_LIMIT&&fabs(c[2])<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=1.0;
                }
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(fabs(c[0])<ZERO_LIMIT&&fabs(c[1])<ZERO_LIMIT&&fabs(c[2])<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=0.1;
                }
                continue;
            }else if(0==strcmp(argv[i], "-pk")){
                i++;
                pk_param=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-D")){
                i++;
                D=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-dt")){
                i++;
                dt=be_double(argv[i++]);
                continue;
            }else{
                i++;
                continue;
            }
        }
        NMODEL model(input_path, types, DWs, lambda);
        NED ned(&model, min2Theta, max2Theta, lp_type, c, is_spacing_auto);
        if(pk_param<0.0){
            ned.ned(ned_path, png_path);
        }else{
            ned.ned(ned_path, png_path, pk_param, lambda, D, dt);
        }
    }else if(0==strcmp(argv[i], "--sed")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double lambda=0.0251;
        double Kmax=1.70;
        int    zone[3]={0};
        int    xaxs[3]={0}, yaxs[3]={0};
        double c[3]={0.0, 0.0, 0.0};
        double thickness=0.1;
        double threshold=0.001;
        bool   is_spacing_auto=true;
        bool   is_vtk_output=false;
        char   sed_path[PATH_CHAR_NUMBER]; 
        char   png_path[PATH_CHAR_NUMBER];
        char   vtk_path[PATH_CHAR_NUMBER]; 
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
            }else if(0==strcmp(argv[i], "-q")){
                i++;
                Kmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zone[0]=(int)be_double(argv[i++]);
                zone[1]=(int)be_double(argv[i++]);
                zone[2]=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-x")){
                i++;
                xaxs[0]=(int)be_double(argv[i++]);
                xaxs[1]=(int)be_double(argv[i++]);
                xaxs[2]=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-y")){
                i++;
                yaxs[0]=(int)be_double(argv[i++]);
                yaxs[1]=(int)be_double(argv[i++]);
                yaxs[2]=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-t")){
                i++;
                thickness=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c[0]=be_double(argv[i++]);
                c[1]=be_double(argv[i++]);
                c[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-auto")){
                i++;
                is_spacing_auto=true;
                if(fabs(c[0])<ZERO_LIMIT&&fabs(c[1])<ZERO_LIMIT&&fabs(c[2])<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=1.0;
                }
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(fabs(c[0])<ZERO_LIMIT&&fabs(c[1])<ZERO_LIMIT&&fabs(c[2])<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=0.1;
                }
                continue;
            }else if(0==strcmp(argv[i], "-thr")){
                i++;
                threshold=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-vtk")){
                i++;
                is_vtk_output=true;
                continue;
            }else{
                i++;
                continue;
            }
        }
        EMODEL model(input_path, types, DWs, lambda);
        if((!model.is_orthogonal)&&is_vtk_output){
            printf("[WARNING] Unable to create vtk file since the input model is not orthogonal");
            is_vtk_output=false;
        }
        if(zone[0]==0&&zone[1]==0&&zone[2]==0){
            strcpy(sed_path, name); strcat(sed_path, ".sed");
            strcpy(vtk_path, name); strcat(vtk_path, ".vtk");
            SED sed(&model, Kmax, c, is_spacing_auto);
            sed.sed(sed_path, threshold);
            if(is_vtk_output) sed.vtk(vtk_path, threshold);
        }else{
            SED sed(&model, zone, thickness, Kmax, c, is_spacing_auto);
            char uvw[4][EXT_CHAR_NUMBER]; strcpy(uvw[0], ".");
            int_to_str(uvw[1], zone[0]); int_to_str(uvw[2], zone[1]); int_to_str(uvw[3], zone[2]);
            strcpy(sed_path, name); merge_path(sed_path, uvw, 4); strcat(sed_path, ".sed");
            strcpy(vtk_path, name); merge_path(vtk_path, uvw, 4); strcat(vtk_path, ".sed");
            strcpy(png_path, name); merge_path(png_path, uvw, 4); strcat(png_path, ".png");
            sed.sed(sed_path, png_path, xaxs, yaxs, zone, threshold);
            if(is_vtk_output) sed.vtk(vtk_path, threshold);
        }
    }else if(0==strcmp(argv[i], "--kkd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double lambda=0.0251;
        double Kmax=1.70;
        double c[3]={0.0};
        double threshold=0.001;
        bool   is_spacing_auto=true;
        int    zone[3]={0};
        double dist=3.0;
        int    dpi=300;
        int    width=6;
        int    height=4;
        char   background='w';
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
            }else if(0==strcmp(argv[i], "-q")){
                i++;
                Kmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c[0]=be_double(argv[i++]);
                c[1]=be_double(argv[i++]);
                c[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-auto")){
                i++;
                is_spacing_auto=true;
                if(fabs(c[0])<ZERO_LIMIT&&fabs(c[1])<ZERO_LIMIT&&fabs(c[2])<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=1.0;
                }
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(fabs(c[0])<ZERO_LIMIT&&fabs(c[1])<ZERO_LIMIT&&fabs(c[2])<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=0.1;
                }
                continue;
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zone[0]=(int)be_double(argv[i++]);
                zone[1]=(int)be_double(argv[i++]);
                zone[2]=(int)be_double(argv[i++]);
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
        if(0==strcmp(ext, ".restart")){
            char   img_path[PATH_CHAR_NUMBER]; 
            if(0==zone[0]&&0==zone[1]&&0==zone[2]){
                strcpy(img_path, name); strcat(img_path, ".png");
                KKD kkd(input_path, threshold, dist, width, dpi);
                kkd.img(img_path, background);
            }else{
                strcpy(img_path, name);
                char uvw[4][EXT_CHAR_NUMBER]; strcpy(uvw[0], ".");
                int_to_str(uvw[1], zone[0]); int_to_str(uvw[2], zone[1]); int_to_str(uvw[3], zone[2]);
                strcpy(img_path, name); merge_path(img_path, uvw, 4); strcat(img_path, ".png");
                KKD kkd(input_path, zone, threshold, dist, width, height, dpi);
                kkd.img(img_path, background);
            }
        }else{
            char   restart_path[PATH_CHAR_NUMBER];
            strcpy(restart_path, name); strcat(restart_path, ".restart");
            EMODEL model(input_path, types, DWs, lambda);
            SED sed(&model, Kmax, c, is_spacing_auto);
            sed.restart(restart_path);
        }
    }else if(0==strcmp(argv[i], "--dkd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double omega=0.0, sigma=75.7;
        double EkeV=20.0, Emin=19.0;
        double zmax=100.0, zstep=1.0;
        int    num_e=20000;
        int    nump=501;
        double dmin=0.1;
        double c1=4.0, c2=8.0, c3=50.0, c_sg=1.0;
        double img_L=6.0;
        int    img_dpi=512;
        char   img_path[PATH_CHAR_NUMBER]; strcpy(img_path, name); strcat(img_path, ".png");
        while(i<argc){
            if(0==strcmp(argv[i], "-e")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    strcpy(types[count++], argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-dw")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    DWs[count++]=be_double(argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-E")){
                i++;
                EkeV=be_double(argv[i++]);
                Emin=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zmax=be_double(argv[i++]);
                zstep=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-ne")){
                i++;
                num_e=(int)be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-np")){
                i++;
                nump=(int)be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-d")){
                i++;
                dmin=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c1=be_double(argv[i++]);
                c2=be_double(argv[i++]);
                c3=be_double(argv[i++]);
                c_sg=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-L")){
                i++;
                img_L=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-dpi")){
                i++;
                img_dpi=(int)be_double(argv[i++]);
            }else{
                i++;
                continue;
            }
        }
        CELL cell(input_path, types, DWs);
        DKD_MC mc(&cell, omega, sigma, EkeV, Emin, EkeV-Emin, zmax, zstep, num_e, nump);
        mc.img(img_path, img_L, img_dpi);
        DKD dkd(&mc, &cell, dmin, c1, c2, c3, c_sg);
        dkd.img(img_path, img_L, img_dpi);
    }else if(0==strcmp(argv[i], "--rdf")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        double rmax=5.0;
        int    nbin=200;
        bool   is_partial=false;
        char   rdf_path[PATH_CHAR_NUMBER]; strcpy(rdf_path, name);
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
                is_partial=true;
                continue;
            }else{
                i++;
                continue;
            }
        }
        strcat(rdf_path, ".rdf");
        RDF rdf(input_path, rmax, nbin, is_partial);
        rdf.rdf(rdf_path);
    }else if(0==strcmp(argv[i], "--ssf")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        double qmax=5.0;
        int    nbin=200;
        bool   is_partial=false;
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
            }else if(0==strcmp(argv[i], "-partial")){
                i++;
                is_partial=true;
                continue;
            }else{
                i++;
                continue;
            }
        }
        strcat(ssf_path, ".ssf");
        SSF ssf(input_path, qmax, nbin, is_partial);
        ssf.ssf(ssf_path);
    }else{
        printf("[ERROR] Unrecognized command %s.\n", argv[i]);
    }
}
