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
        char   xrd_path[PATH_CHAR_NUMBER]; 
        char   png_path[PATH_CHAR_NUMBER];
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
            }else if(0==strcmp(argv[i], "-d")){
                i++;
                D=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-dt")){
                i++;
                dt=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(name, argv[i++]);
                continue;
            }else{
                i++;
                continue;
            }
        }
        strcpy(xrd_path, name); strcat(xrd_path, ".xrd");
        strcpy(png_path, name); strcat(png_path, ".xrd"); strcat(png_path, ".png");
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
        char   ned_path[PATH_CHAR_NUMBER]; 
        char   png_path[PATH_CHAR_NUMBER]; 
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
            }else if(0==strcmp(argv[i], "-d")){
                i++;
                D=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-dt")){
                i++;
                dt=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(name, argv[i++]);
                continue;
            }else{
                i++;
                continue;
            }
        }
        strcpy(ned_path, name); strcat(ned_path, ".ned");
        strcpy(png_path, name); strcat(png_path, ".ned"); strcat(png_path, ".png");
        NMODEL model(input_path, types, DWs, lambda);
        NED ned(&model, min2Theta, max2Theta, lp_type, c, is_spacing_auto);
        if(pk_param<0.0){
            ned.ned(ned_path, png_path);
        }else{
            ned.ned(ned_path, png_path, pk_param, lambda, D, dt);
        }
    }else if(0==strcmp(argv[i], "--ked")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double lambda=0.0251;
        double Kmax=1.70;
        int    zaxs[3]={0};
        int    xaxs[3]={0}, yaxs[3]={0};
        double c[3]={0.0, 0.0, 0.0};
        double thickness=0.1;
        double threshold=0.001;
        bool   is_spacing_auto=true;
        bool   is_vtk_output=false;
        char   ked_path[PATH_CHAR_NUMBER]; 
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
                zaxs[0]=(int)be_double(argv[i++]);
                zaxs[1]=(int)be_double(argv[i++]);
                zaxs[2]=(int)be_double(argv[i++]);
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
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(name, argv[i++]);
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
        if(zaxs[0]==0&&zaxs[1]==0&&zaxs[2]==0){
            strcpy(ked_path, name); strcat(ked_path, ".ked");
            strcpy(vtk_path, name); strcat(vtk_path, ".vtk");
            KED ked(&model, Kmax, c, is_spacing_auto);
            ked.ked(ked_path, threshold);
            if(is_vtk_output) ked.vtk(vtk_path, threshold);
        }else{
            KED ked(&model, zaxs, thickness, Kmax, c, is_spacing_auto);
            char uvw[4][EXT_CHAR_NUMBER]; strcpy(uvw[0], ".");
            int_to_str(uvw[1], zaxs[0]); int_to_str(uvw[2], zaxs[1]); int_to_str(uvw[3], zaxs[2]);
            strcpy(ked_path, name); merge_path(ked_path, uvw, 4); strcat(ked_path, ".ked");
            strcpy(vtk_path, name); merge_path(vtk_path, uvw, 4); strcat(vtk_path, ".vtk");
            strcpy(png_path, name); merge_path(png_path, uvw, 4); strcat(png_path, ".ked.png");
            ked.ked(ked_path, png_path, xaxs, yaxs, zaxs, threshold);
            if(is_vtk_output) ked.vtk(vtk_path, threshold);
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
        int    zaxs[3]={0};
        double dist=3.0;
        int    dpi=300;
        double width=3;
        double height=3;
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
                zaxs[0]=(int)be_double(argv[i++]);
                zaxs[1]=(int)be_double(argv[i++]);
                zaxs[2]=(int)be_double(argv[i++]);
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
            if(0==zaxs[0]&&0==zaxs[1]&&0==zaxs[2]){
                strcpy(img_path, name); strcat(img_path, ".kkd.png");
                KKD kkd(input_path, threshold, dist, width, dpi);
                kkd.img(img_path, background);
            }else{
                strcpy(img_path, name);
                char uvw[4][EXT_CHAR_NUMBER]; strcpy(uvw[0], ".");
                int_to_str(uvw[1], zaxs[0]); int_to_str(uvw[2], zaxs[1]); int_to_str(uvw[3], zaxs[2]);
                strcpy(img_path, name); merge_path(img_path, uvw, 4); strcat(img_path, ".kkd.png");
                KKD kkd(input_path, zaxs, threshold, dist, width, height, dpi);
                kkd.img(img_path, background);
            }
        }else{
            char   restart_path[PATH_CHAR_NUMBER];
            strcpy(restart_path, name); strcat(restart_path, ".restart");
            EMODEL model(input_path, types, DWs, lambda);
            KED ked(&model, Kmax, c, is_spacing_auto);
            ked.restart(restart_path);
        }
    }else if(0==strcmp(argv[i], "--ded")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        for(int i=0;i<TYPE_INPUT_NUMBER;i++) DWs[i]=1.0e-6;
        int    zone[3]={0, 0, 1};
        int    fnorm[3]={0, 0, 1};
        double thickness=10.0;
        double dmin=0.1;
        double voltage=200.0;
        double c1=4.0, c2=8.0, c3=50.0, c_sg=0.2;
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
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zone[0]=(int)be_double(argv[i++]);
                zone[1]=(int)be_double(argv[i++]);
                zone[2]=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-n")){
                i++;
                fnorm[0]=(int)be_double(argv[i++]);
                fnorm[1]=(int)be_double(argv[i++]);
                fnorm[2]=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-t")){
                i++;
                thickness=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-d")){
                i++;
                dmin=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-v")){
                i++;
                voltage=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c1=be_double(argv[i++]);
                c2=be_double(argv[i++]);
                c3=be_double(argv[i++]);
                c_sg=be_double(argv[i++]);
            }else{
                i++;
                continue;
            }
        }
        char   ded_path[PATH_CHAR_NUMBER]; strcpy(ded_path, name);
        char uvw[4][EXT_CHAR_NUMBER]; strcpy(uvw[0], ".");
        int_to_str(uvw[1], zone[0]); int_to_str(uvw[2], zone[1]); int_to_str(uvw[3], zone[2]);
        strcpy(ded_path, name); merge_path(ded_path, uvw, 4); strcat(ded_path, ".ded");
        CELL cell(input_path, types, DWs);
        DED_BETHE bethe;
        bethe.c1=c1; bethe.c2=c2; bethe.c3=c3; bethe.c_sg=c_sg;
        DED ded(&cell, &bethe, zone, fnorm, voltage, thickness, dmin);
        ded.ded(ded_path);
    }else if(0==strcmp(argv[i], "--dkd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        for(int i=0;i<TYPE_INPUT_NUMBER;i++){
            DWs[i]=1.0e-6;
        }
        double omega=0.0, sigma=75.7;
        double EkeV=20.0, Eexit=19.0;
        double zmax=100.0, zstep=1.0;
        int    num_e=20000;
        int    nump=501;
        double dmin=0.1;
        double c1=4.0, c2=8.0, c_rg=50.0, c_sg=1.0;
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
            }else if(0==strcmp(argv[i], "-en")){
                i++;
                EkeV=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-ex")){
                i++;
                Eexit=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-z0")){
                i++;
                zmax=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-dz")){
                i++;
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
            }else if(0==strcmp(argv[i], "-rg")){
                i++;
                c_rg=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-sg")){
                i++;
                c_sg=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-w")){
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
        DKD_MC mc(&cell, omega, sigma, EkeV, Eexit, EkeV-Eexit, zmax, zstep, num_e, nump);
        mc.img(img_path, img_L, img_dpi);
        DKD dkd(&mc, &cell, dmin, c1, c2, c_rg, c_sg);
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
