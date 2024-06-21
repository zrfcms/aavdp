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
        double lambda=1.54180;
        double c[3]={0.0, 0.0, 0.0};
        double min2Theta=0.0, max2Theta=180.0;
        int    nbin=0;
        bool   is_lorentz=true;
        bool   is_spacing_auto=true;
        char   xrd_path[PATH_CHAR_NUMBER]; strcpy(xrd_path, name); strcat(xrd_path, ".xrd");
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
            }else if(0==strcmp(argv[i], "-n")){
                i++;
                nbin=(int)be_double(argv[i++]);
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
                if(fabs(c[0]-0.0)<ZERO_LIMIT&&fabs(c[1]-0.0)<ZERO_LIMIT&&fabs(c[2]-0.0)<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=1.0;
                }
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(fabs(c[0]-0.0)<ZERO_LIMIT&&fabs(c[1]-0.0)<ZERO_LIMIT&&fabs(c[2]-0.0)<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=0.1;
                }
                continue;
            }else if(0==strcmp(argv[i], "-without_lorentz_factor")){
                i++;
                is_lorentz=false;
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(xrd_path, argv[i++]);
            }else{
                i++;
                continue;
            }
        }
        XMODEL model(input_path, types, DWs, lambda);
        XRD xrd(&model, min2Theta, max2Theta, c, is_spacing_auto, is_lorentz);
        xrd.xrd(xrd_path, nbin);
    }else if(0==strcmp(argv[i], "--snd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double lambda=1.54180;
        double c[3]={1.0, 1.0, 1.0};
        double min2Theta=0.0, max2Theta=180.0;
        int    nbin=0;
        bool   is_lorentz=true;
        bool   is_spacing_auto=true;
        char   snd_path[PATH_CHAR_NUMBER]; strcpy(snd_path, name); strcat(snd_path, ".snd");
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
            }else if(0==strcmp(argv[i], "-n")){
                i++;
                nbin=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c[0]=be_double(argv[i++]);
                c[1]=be_double(argv[i++]);
                c[2]=be_double(argv[i++]);
                is_spacing_auto=false;
                continue;
            }else if(0==strcmp(argv[i], "-auto")){
                i++;
                is_spacing_auto=true;
                if(fabs(c[0]-0.0)<ZERO_LIMIT&&fabs(c[1]-0.0)<ZERO_LIMIT&&fabs(c[2]-0.0)<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=1.0;
                }
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(fabs(c[0]-0.0)<ZERO_LIMIT&&fabs(c[1]-0.0)<ZERO_LIMIT&&fabs(c[2]-0.0)<ZERO_LIMIT){
                    c[0]=c[1]=c[2]=0.1;
                }
                continue;
            }else if(0==strcmp(argv[i], "-without_lorentz_factor")){
                i++;
                is_lorentz=false;
                continue;
            }else{
                i++;
                continue;
            }
        }
        NMODEL model(input_path, types, DWs, lambda);
        SND snd(&model, min2Theta, max2Theta, c, is_spacing_auto, is_lorentz);
        snd.snd(snd_path, nbin);
    }else if(0==strcmp(argv[i], "--sed")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        if(0==strcmp(ext, ".h5")||0==strcmp(ext, ".hdf5")){
            int    xaxis[3]={0, 0, 0};
            int    zaxis[3]={0, 0, 1};
            double thickness=0.1;
            double threshold=0.001;
            char   sed_path[PATH_CHAR_NUMBER]; strcpy(sed_path, name);
            while(i<argc){
                if(0==strcmp(argv[i], "-x")){
                    i++;
                    xaxis[0]=(int)be_double(argv[i++]);
                    xaxis[1]=(int)be_double(argv[i++]);
                    xaxis[2]=(int)be_double(argv[i++]);
                    continue;
                }else if(0==strcmp(argv[i], "-z")){
                    i++;
                    zaxis[0]=(int)be_double(argv[i++]);
                    zaxis[1]=(int)be_double(argv[i++]);
                    zaxis[2]=(int)be_double(argv[i++]);
                    continue;
                }else if(0==strcmp(argv[i], "-t")){
                    i++;
                    thickness=be_double(argv[i++]);
                    continue;
                }else if(0==strcmp(argv[i], "-thr")){
                    i++;
                    threshold=be_double(argv[i++]);
                    continue;
                }else{
                    i++;
                    continue;
                }
            }
            char uvw[4][EXT_CHAR_NUMBER]; strcpy(uvw[0], ".");
            int_to_str(uvw[1], zaxis[0]); int_to_str(uvw[2], zaxis[1]); int_to_str(uvw[3], zaxis[2]);
            merge_path(sed_path, uvw, 4); strcat(sed_path, ".sed");
            SED sed(input_path);
            sed.sed(sed_path, xaxis, zaxis, thickness, threshold);
        }else{
            char   types[TYPE_INPUT_NUMBER][10]={0};
            double DWs[TYPE_INPUT_NUMBER]={0.0};
            double lambda=0.0251;
            double Kmax=1.70;
            double c[3]={0.0, 0.0, 0.0};
            int    zone[3]={0, 0, 0};
            double thickness=0.1;
            bool   is_spacing_auto=true;
            bool   is_vtk_output=false;
            char   hdf5_path[PATH_CHAR_NUMBER]; strcpy(hdf5_path, name);
            char   vtk_path[PATH_CHAR_NUMBER]; strcpy(vtk_path, name);
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
                    if(fabs(c[0]-0.0)<ZERO_LIMIT&&fabs(c[1]-0.0)<ZERO_LIMIT&&fabs(c[2]-0.0)<ZERO_LIMIT){
                        c[0]=c[1]=c[2]=1.0;
                    }
                    continue;
                }else if(0==strcmp(argv[i], "-manu")){
                    i++;
                    is_spacing_auto=false;
                    if(fabs(c[0]-0.0)<ZERO_LIMIT&&fabs(c[1]-0.0)<ZERO_LIMIT&&fabs(c[2]-0.0)<ZERO_LIMIT){
                        c[0]=c[1]=c[2]=0.1;
                    }
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
            char uvw[4][EXT_CHAR_NUMBER]; strcpy(uvw[0], ".");
            int_to_str(uvw[1], zone[0]); int_to_str(uvw[2], zone[1]); int_to_str(uvw[3], zone[2]);
            merge_path(hdf5_path, uvw, 4); strcat(hdf5_path, ".h5");
            merge_path(vtk_path, uvw, 4); strcat(vtk_path, ".vtk");
            EMODEL model(input_path, types, DWs, lambda);
            if((!model.is_orthogonal)&&is_vtk_output){
                printf("[WARNING] Unable to create vtk file since the input model is not orthogonal");
                is_vtk_output=false;
            }
            if(0==zone[0]&&0==zone[1]&&0==zone[2]){
                SED sed(&model, Kmax, c, is_spacing_auto);
                sed.hdf5(hdf5_path);
                if(is_vtk_output) sed.vtk(vtk_path);
            }else{
                SED sed(&model, zone, thickness, Kmax, c, is_spacing_auto);
                sed.hdf5(hdf5_path);
                if(is_vtk_output) sed.vtk(vtk_path);
            }
        }
    }else if(0==strcmp(argv[i], "--kkd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        int    zone[3]={0};
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
        char   exts[5][EXT_CHAR_NUMBER]; strcpy(exts[0], "."); strcpy(exts[4], ".png");
        int_to_str(exts[1], zone[0]); int_to_str(exts[2], zone[1]); int_to_str(exts[3], zone[2]);
        merge_path(img_path, exts, 5);
        if(is_zone_existing){
            KKD kkd(input_path, zone, threshold, dist, width, height, dpi);
            kkd.img(img_path, background);
        }else{
            KKD kkd(input_path, threshold, dist, width, dpi);
            kkd.img(img_path, background);
        }
    }else if(0==strcmp(argv[i], "--dkd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);
        if(0==strcmp(ext, ".h5")||0==strcmp(ext, ".hdf5")){
            if(0==strcmp(argv[i], "-m")){
                i++;
                char   img_path[PATH_CHAR_NUMBER]; strcpy(img_path, name);
                double omega=0.0, sigma=75.7;
                double Emax=20.0, Emin=10.0, Ebin=1.0;
                double zmax=100.0, zstep=1.0;
                int    num_e=10000;
                int    nump=501;
                double img_L=6.0;
                int    img_dpi=512;
                char   hdf5_path[PATH_CHAR_NUMBER];
                bool   is_hdf5_flag=false;
                while(i<argc){
                    if(0==strcmp(argv[i], "-t")){
                        i++;
                        omega=be_double(argv[i++]);
                        sigma=be_double(argv[i++]);
                    }else if(0==strcmp(argv[i], "-E")){
                        i++;
                        Emax=be_double(argv[i++]);
                        Emin=be_double(argv[i++]);
                        Ebin=be_double(argv[i++]);
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
                    }else if(0==strcmp(argv[i], "-L")){
                        i++;
                        img_L=be_double(argv[i++]);
                    }else if(0==strcmp(argv[i], "-dpi")){
                        i++;
                        img_dpi=(int)be_double(argv[i++]);
                    }else if(0==strcmp(argv[i], "-h5")){
                        i++;
                        strcpy(hdf5_path, argv[i++]);
                        is_hdf5_flag=true;
                    }else{
                        i++;
                        continue;
                    }
                }
                strcat(img_path, ".png");
                if(is_hdf5_flag){
                    DKD_MC mc(hdf5_path);
                    mc.hdf5(input_path);
                    mc.img(img_path, img_L, img_dpi);
                }else{
                    DKD_MC mc(input_path, omega, sigma, Emax, Emin, Ebin, zmax, zstep, num_e, nump);
                    mc.hdf5(input_path);
                    mc.img(img_path, img_L, img_dpi);
                }
            }else if(0==strcmp(argv[i], "-k")){
                i++;
                char   img_path[PATH_CHAR_NUMBER]; strcpy(img_path, name);
                double dmin=0.05;
                double c1=4.0, c2=8.0, c3=50.0, c_sg=1.0;
                double img_L=6.0;
                int    img_dpi=512;
                while(i<argc){
                    if(0==strcmp(argv[i], "-d")){
                        i++;
                        dmin=be_double(argv[i++]);
                    }else if(0==strcmp(argv[i], "-b")){
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
                strcat(img_path, ".png");
                DKD dkd(input_path, dmin, c1, c2, c3, c_sg);
                dkd.hdf5(input_path);
                dkd.img(img_path, img_L, img_dpi);
            }else if(0==strcmp(argv[i], "-g")){
                i++;
                char   img_path[PATH_CHAR_NUMBER]; strcpy(img_path, name);
                double img_L=6.0;
                int    img_dpi=512;
                while(i<argc){
                    if(0==strcmp(argv[i], "-L")){
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
                strcat(img_path, ".png");
                DKD_MC mc(input_path);
                mc.img(img_path, img_L, img_dpi);
                DKD dkd(input_path);
                dkd.img(img_path, img_L, img_dpi);
            }
        }else{
            char   hdf5_path[PATH_CHAR_NUMBER]; strcpy(hdf5_path, name);
            char   types[TYPE_INPUT_NUMBER][10]={0};
            double DWs[TYPE_INPUT_NUMBER]={0.0};
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
                }else{
                    i++;
                    continue;
                }
            }
            strcat(hdf5_path, ".h5");
            CELL cell(input_path, types, DWs);
            cell.hdf5(hdf5_path);
        }
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
