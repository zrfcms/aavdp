#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "main.h"

int main(int argc, char* argv[])
{
    int    i=1;
    if(argc<2){
        print_version();
    }else if(argc<3){
        if(0==strcmp(argv[i], "-h")||0==strcmp(argv[i], "--help")){
            print_help();
        }else if(0==strcmp(argv[i], "-v")||0==strcmp(argv[i], "--version")){
            print_version();
        }else if(0==strcmp(argv[i], "--xrd")||0==strcmp(argv[i], "--ned")||0==strcmp(argv[i], "--ked")||0==strcmp(argv[i], "--kkd")||0==strcmp(argv[i], "--ded")||0==strcmp(argv[i], "--dkd")||0==strcmp(argv[i], "--rdf")||0==strcmp(argv[i], "--ssf")){
            printf("[ERROR] inputfile is required for mode '%s'\n", argv[i]);
        }else{
            printf("AAVDP: unrecognized mode '%s'\n", argv[i]);
            printf("Try 'AAVDP -h' for more information'\n");
        }
    }else if(0==strcmp(argv[i], "--xrd")||0==strcmp(argv[i], "--ned")){
        char   *mode; mode=argv[i]+2;
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   output_path[PATH_CHAR_NUMBER]; strcpy(output_path, "");
        is_path_accessible(input_path);

        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double lambda=1.54184;
        double min2Theta=0.0;
        double max2Theta=90.0;
        int    lp_type=3;
        if(0==strcmp(mode, "ned")) lp_type=2;
        double c[3]={1.0, 1.0, 1.0};
        bool   is_spacing_auto=true;
        bool   is_c_default=true;
        double threshold=0.001;

        double mixing_param=0.0;
        double grain_diameter=500.0;
        double d2t=0.01;
        bool   is_scherrer=false;
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
                is_c_default=false;
                continue;
            }else if(0==strcmp(argv[i], "-auto")){
                i++;
                is_spacing_auto=true;
                if(is_c_default) c[0]=c[1]=c[2]=1.0;
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(is_c_default) c[0]=c[1]=c[2]=0.1;
                continue;
            }else if(0==strcmp(argv[i], "-thr")){
                i++;
                threshold=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(output_path, argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-scherrer")){
                i++;
                is_scherrer=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-scherrer_m")){
                        i++;
                        mixing_param=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-scherrer_d")){
                        i++;
                        grain_diameter=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-scherrer_d2t")){
                        i++;
                        d2t=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-scherrer")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h' for more information");
                exit(1);
            }
        }
        char   xrd_path[PATH_CHAR_NUMBER];
        if(0==strcmp(output_path, "")){
            strcpy(xrd_path, "./AAVDP."); strcat(xrd_path, mode);
        }else{
            strcpy(xrd_path, output_path);
        }
        if(0==strcmp(mode, "xrd")){
            XMODEL model(input_path, types, DWs, lambda);
            model.update_lorentzP_type(lp_type);
            XRD xrd(&model, min2Theta, max2Theta, threshold, c, is_spacing_auto);
            if(is_scherrer){
                xrd.xrd(xrd_path, mixing_param, lambda, grain_diameter, d2t);
            }else{
                xrd.xrd(xrd_path);
            }
        }else{
            NMODEL model(input_path, types, DWs, lambda);
            model.update_lorentzP_type(lp_type);
            XRD xrd(&model, min2Theta, max2Theta, threshold, c, is_spacing_auto);
            if(is_scherrer){
                xrd.xrd(xrd_path, mixing_param, lambda, grain_diameter, d2t);
            }else{
                xrd.xrd(xrd_path);
            } 
        }
    }else if(0==strcmp(argv[i], "--ked")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   output_path[PATH_CHAR_NUMBER]; strcpy(output_path, "");
        is_path_accessible(input_path);

        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double voltage=200.0;
        double qmax=1.0;
        double thickness=0.1;
        int    zaxis[3]={0, 0, 1};
        double c[3]={1.0, 1.0, 1.0};
        bool   is_spacing_auto=true;
        bool   is_c_default=true;
        double threshold=0.001;

        double xaxis[3]={1.0, 0.0, 0.0};
        double yaxis[3]={0.0, 1.0, 0.0};
        bool   is_rotate=false;

        double sigma=0.01;
        double dx=0.005;
        bool   is_gauss=false;
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
                voltage=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-q")){
                i++;
                qmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zaxis[0]=be_double(argv[i++]);
                zaxis[1]=be_double(argv[i++]);
                zaxis[2]=be_double(argv[i++]);
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
                if(is_c_default) c[0]=c[1]=c[2]=1.0;
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(is_c_default) c[0]=c[1]=c[2]=0.1;
                continue;
            }else if(0==strcmp(argv[i], "-thr")){
                i++;
                threshold=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(output_path, argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-rotate")){
                i++;
                is_rotate=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-rotate_x")){
                        i++;
                        xaxis[0]=be_double(argv[i++]);
                        xaxis[1]=be_double(argv[i++]);
                        xaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-rotate_y")){
                        i++;
                        yaxis[0]=be_double(argv[i++]);
                        yaxis[1]=be_double(argv[i++]);
                        yaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-rotate")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "-gauss")){
                i++;
                is_gauss=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-gauss_sig")){
                        i++;
                        sigma=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-gauss_dx")){
                        i++;
                        dx=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-gauss")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h' for more information");
                exit(1);
            }
        }
        char ked_path[PATH_CHAR_NUMBER];
        if(0==strcmp(output_path, "")){
            strcpy(ked_path, "./AAVDP.ked");
        }else{
            strcpy(ked_path, output_path);
        }
        VMODEL model(input_path, types, DWs, voltage);
        KED ked(&model, zaxis, thickness, qmax, threshold, c, is_spacing_auto);
        if(is_rotate) ked.rotate(xaxis, yaxis);
        if(is_gauss){
            ked.ked(ked_path, sigma, dx);
        }else{
            ked.ked(ked_path);
        }
    }else if(0==strcmp(argv[i], "--kkd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   output_path[PATH_CHAR_NUMBER]; strcpy(output_path, "");
        is_path_accessible(input_path);
        char   name[PATH_CHAR_NUMBER], ext[EXT_CHAR_NUMBER];
        split_path(name, ext, input_path);

        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double voltage=200.0;
        double qmax=1.0;
        double c[3]={1.0, 1.0, 1.0};
        bool   is_spacing_auto=true;
        bool   is_c_default=true;
        double threshold=0.001;

        double zaxis[3]={0.0, 0.0, 1.0};
        double ratiox=1.0, ratioy=1.0;
        double thickness=0.2;
        int    npx=501, npy=501;
        char   projection[10]="stereo";

        double xaxis[3]={1.0, 0.0, 0.0};
        double yaxis[3]={0.0, 1.0, 0.0};
        bool   is_rotate=false;

        double imax=1.0e3, imin=0.0;
        char   background='b';
        bool   is_scale=false;

        bool   is_ked3=false;
        while(i<argc){
            if(0==strcmp(argv[i], "-e")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    strcpy(types[count++], argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-en")){
                i++;
                voltage=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-q")){
                i++;
                qmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-c")){
                i++;
                c[0]=be_double(argv[i++]);
                c[1]=be_double(argv[i++]);
                c[2]=be_double(argv[i++]);
                is_c_default=false;
                continue;
            }else if(0==strcmp(argv[i], "-auto")){
                i++;
                is_spacing_auto=true;
                if(is_c_default) c[0]=c[1]=c[2]=1.0;
                continue;
            }else if(0==strcmp(argv[i], "-manu")){
                i++;
                is_spacing_auto=false;
                if(is_c_default) c[0]=c[1]=c[2]=0.1;
                continue;
            }else if(0==strcmp(argv[i], "-thr")){
                i++;
                threshold=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zaxis[0]=be_double(argv[i++]);
                zaxis[1]=be_double(argv[i++]);
                zaxis[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-rx")){
                i++;
                ratiox=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-ry")){
                i++;
                ratioy=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-t")){
                i++;
                thickness=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-px")){
                i++;
                npx=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-py")){
                i++;
                npy=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-projection")){
                i++;
                strcpy(projection, argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-background")){
                i++;
                background=(argv[i++])[0];
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(output_path, argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-rotate")){
                i++;
                is_rotate=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-rotate_x")){
                        i++;
                        xaxis[0]=be_double(argv[i++]);
                        xaxis[1]=be_double(argv[i++]);
                        xaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-rotate_y")){
                        i++;
                        yaxis[0]=be_double(argv[i++]);
                        yaxis[1]=be_double(argv[i++]);
                        yaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-rotate")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "-scale")){
                i++;
                is_scale=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-scale_i")){
                        i++;
                        imin=be_double(argv[i++]);
                        imax=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-scale")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "-ked3")){
                i++;
                is_ked3=true;
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h' for more information");
                exit(1);
            }
        }
        if(is_ked3){
            char ked3_path[PATH_CHAR_NUMBER];
            if(0==strcmp(output_path, "")){
                strcpy(ked3_path, "./AAVDP.ked3");
            }else{
                strcpy(ked3_path, output_path);
            }
            VMODEL model(input_path, types, DWs, voltage);
            KED ked(&model, qmax, threshold, c, is_spacing_auto);
            ked.ked3(ked3_path);
        }else{
            char kkd_path[PATH_CHAR_NUMBER];
            if(0==strcmp(output_path, "")){
                strcpy(kkd_path, "./AAVDP.kkd");
            }else{
                strcpy(kkd_path, output_path);
            }
            if((!is_rotate)&&(0!=zaxis[0]||0!=zaxis[1])){
                xaxis[0]=-zaxis[1]; xaxis[1]=zaxis[0]; xaxis[2]=0.0;
                yaxis[0]=zaxis[2]*zaxis[0]; yaxis[1]=zaxis[2]*zaxis[1]; yaxis[2]=-1.0*(zaxis[0]*zaxis[0]+zaxis[1]*zaxis[1]);
            }
            if(0==npx%2) npx++;
            if(0==npy%2) npy++;
            if(0==strcmp(ext, ".ked3")){
                KED ked(input_path);
                KKD kkd(&ked, xaxis, yaxis, zaxis, thickness, ratiox, ratioy, npx, npy, projection);
                if(is_scale){
                    kkd.kkd(kkd_path, imax, imin, background);
                }else{
                    kkd.kkd(kkd_path, background);
                }
            }else{
                VMODEL model(input_path, types, DWs, voltage);
                KED ked(&model, qmax, threshold, c, is_spacing_auto);
                KKD kkd(&ked, xaxis, yaxis, zaxis, thickness, ratiox, ratioy, npx, npy, projection);
                if(is_scale){
                    kkd.kkd(kkd_path, imax, imin, background);
                }else{
                    kkd.kkd(kkd_path, background);
                }
            }
        }
    }else if(0==strcmp(argv[i], "--ded")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   output_path[PATH_CHAR_NUMBER]; strcpy(output_path, "");
        is_path_accessible(input_path);
        
        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double voltage=200.0;
        double qmax=1.0;
        int    zone[3]={0, 0, 1};
        int    fnorm[3]={0, 0, 1};
        double fthick=100.0;
        double c1=4.0, c2=8.0, c3=50.0, c_sg=1.0;
        double threshold=0.001;

        double xaxis[3]={1.0, 0.0, 0.0};
        double yaxis[3]={0.0, 1.0, 0.0};
        bool   is_rotate=false;

        double sigma=0.01;
        double dx=0.005;
        bool   is_gauss=false;
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
                voltage=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-q")){
                i++;
                qmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zone[0]=(int)be_double(argv[i++]);
                zone[1]=(int)be_double(argv[i++]);
                zone[2]=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-fn")){
                i++;
                fnorm[0]=(int)be_double(argv[i++]);
                fnorm[1]=(int)be_double(argv[i++]);
                fnorm[2]=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-ft")){
                i++;
                fthick=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-bethe")){
                i++;
                c1=be_double(argv[i++]);
                c2=be_double(argv[i++]);
                c3=be_double(argv[i++]);
                c_sg=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-thr")){
                i++;
                threshold=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(output_path, argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-rotate")){
                i++;
                is_rotate=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-rotate_x")){
                        i++;
                        xaxis[0]=be_double(argv[i++]);
                        xaxis[1]=be_double(argv[i++]);
                        xaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-rotate_y")){
                        i++;
                        yaxis[0]=be_double(argv[i++]);
                        yaxis[1]=be_double(argv[i++]);
                        yaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-rotate")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "-gauss")){
                i++;
                is_gauss=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-gauss_sig")){
                        i++;
                        sigma=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-gauss_dx")){
                        i++;
                        dx=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-gauss")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h' for more information");
                exit(1);
            }
        }
        char   ded_path[PATH_CHAR_NUMBER];
        if(0==strcmp(output_path, "")){
            strcpy(ded_path, "./AAVDP.ded");
        }else{
            strcpy(ded_path, output_path);
        }
        CELL cell(input_path, types, DWs, voltage);
        BETHE bethe;
        bethe.c1=c1; bethe.c2=c2; bethe.c3=c3; bethe.c_sg=c_sg;
        DED ded(&cell, &bethe, zone, fnorm, fthick, voltage, qmax, threshold);
        if(is_rotate) ded.rotate(xaxis, yaxis);
        if(is_gauss){
            ded.ded(ded_path, sigma, dx);
        }else{
            ded.ded(ded_path);
        }
    }else if(0==strcmp(argv[i], "--dkd")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   output_path[PATH_CHAR_NUMBER]; strcpy(output_path, "");
        char   mc_path[PATH_CHAR_NUMBER]; strcpy(mc_path, "");
        is_path_accessible(input_path);

        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        double voltage=200.0;
        double qmax=1.0;
        double fthick=100.0;
        double c1=4.0, c2=8.0, c3=50.0, c_sg=1.0;

        double zaxis[3]={0.0, 0.0, 1.0};
        double ratiox=1.0, ratioy=1.0;
        int    numpx=501, numpy=501;
        char   projection[10]="stereo";
        char   background='b';

        double omega=0.0, sigma=75.7;
        double Eexit=190;
        double dthick=10.0;
        int    nume=20000;
        int    nump=501;
        char   seed_path[PATH_CHAR_NUMBER]="./RandomSeeds.data";
        bool   is_monte=false;

        double xaxis[3]={1.0, 0.0, 0.0};
        double yaxis[3]={0.0, 1.0, 0.0};
        bool   is_rotate=false;

        double imax=1.0e3, imin=0.0;
        bool   is_scale=false;
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
                voltage=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-q")){
                i++;
                qmax=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-ft")){
                i++;
                fthick=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-bethe")){
                i++;
                c1=be_double(argv[i++]);
                c2=be_double(argv[i++]);
                c3=be_double(argv[i++]);
                c_sg=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zaxis[0]=be_double(argv[i++]);
                zaxis[1]=be_double(argv[i++]);
                zaxis[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-rx")){
                i++;
                ratiox=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-ry")){
                i++;
                ratioy=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-px")){
                i++;
                numpx=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-py")){
                i++;
                numpy=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-projection")){
                i++;
                strcpy(projection, argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-background")){
                i++;
                background=(argv[i++])[0];
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(output_path, argv[i++]);
            }else if(0==strcmp(argv[i], "-monte")){
                i++;
                is_monte=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-monte_rd")){
                        i++;
                        omega=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-monte_td")){
                        i++;
                        sigma=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-monte_ex")){
                        i++;
                        Eexit=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-monte_dt")){
                        i++;
                        dthick=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-monte_ne")){
                        i++;
                        nume=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-monte_np")){
                        i++;
                        nump=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-monte_seed")){
                        i++;
                        strcpy(seed_path, argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-monte_o")){
                        i++;
                        strcpy(mc_path, argv[i++]);
                    }else if(not_in_switch(argv[i], "-monte")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "-rotate")){
                i++;
                is_rotate=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-rotate_x")){
                        i++;
                        xaxis[0]=be_double(argv[i++]);
                        xaxis[1]=be_double(argv[i++]);
                        xaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-rotate_y")){
                        i++;
                        yaxis[0]=be_double(argv[i++]);
                        yaxis[1]=be_double(argv[i++]);
                        yaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-rotate")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "-scale")){
                i++;
                is_scale=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-scale_i")){
                        i++;
                        imin=be_double(argv[i++]);
                        imax=be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-scale")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h' for more information");
                exit(1);
            }
        }
        char   dkd_path[PATH_CHAR_NUMBER];
        if(0==strcmp(output_path, "")){
            strcpy(dkd_path, "./AAVDP.dkd");
        }else{
            strcpy(dkd_path, output_path);
        }
        if((!is_rotate)&&(0!=zaxis[0]||0!=zaxis[1])){
            xaxis[0]=-zaxis[1]; xaxis[1]=zaxis[0]; xaxis[2]=0.0;
            yaxis[0]=zaxis[2]*zaxis[0]; yaxis[1]=zaxis[2]*zaxis[1]; yaxis[2]=-1.0*(zaxis[0]*zaxis[0]+zaxis[1]*zaxis[1]);
        }
        if(0==numpx%2) numpx++;
        if(0==numpy%2) numpy++;
        if(0==nump%2)  nump++;
        CELL cell(input_path, types, DWs, voltage);
        BETHE bethe;
        bethe.c1=c1; bethe.c2=c2; bethe.c3=c3; bethe.c_sg=c_sg;
        if(is_monte){
            if(0==strcmp(mc_path, "")){
                strcpy(mc_path, "./AAVDP.mc");
            }
            MC mc(&cell, seed_path, omega, sigma, voltage, Eexit, fthick, dthick, nume, nump);
            mc.mc(mc_path, 'b');
            // DKD dkd(&cell, &mc, &bethe, Kmax, projection);
            DKD dkd(&cell, &mc, &bethe, xaxis, yaxis, zaxis, ratiox, ratioy, qmax, numpx, numpy, projection);
            if(is_scale){
                dkd.dkd(dkd_path, imax, imin, background);
            }else{
                dkd.dkd(dkd_path, background);
            }
        }else{
            // DKD dkd(&cell, &bethe, voltage, fthick, Kmax, nump, projection);
            DKD dkd(&cell, &bethe, xaxis, yaxis, zaxis, ratiox, ratioy, voltage, fthick, qmax, numpx, numpy, projection);
            if(is_scale){
                dkd.dkd(dkd_path, imax, imin, background);
            }else{
                dkd.dkd(dkd_path, background);
            }
        }
    }else if(0==strcmp(argv[i], "--rdf")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   output_path[PATH_CHAR_NUMBER]; strcpy(output_path, "");
        is_path_accessible(input_path);

        double rmax=5.0;
        int    nbin=100;
        bool   is_partial=false;
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
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(output_path, argv[i++]);
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h' for more information");
                exit(1);
            }
        }
        char   rdf_path[PATH_CHAR_NUMBER];
        if(0==strcmp(output_path, "")){
            strcpy(rdf_path, "./AAVDP.rdf");
        }else{
            strcpy(rdf_path, output_path);
        }
        RDF rdf(input_path, rmax, nbin, is_partial);
        rdf.rdf(rdf_path);
    }else if(0==strcmp(argv[i], "--ssf")){
        i++;
        char   input_path[PATH_CHAR_NUMBER]; strcpy(input_path, argv[i++]);
        char   output_path[PATH_CHAR_NUMBER]; strcpy(output_path, "");
        is_path_accessible(input_path);

        double qmax=5.0;
        int    nqbin=100;
        bool   is_partial=false;

        double rmax=5.0;
        int    nrbin=100;
        bool   is_rdf=false;
        while(i<argc){
            if(0==strcmp(argv[i], "-q")){
                i++;
                qmax=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-n")){
                i++;
                nqbin=(int)be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-partial")){
                i++;
                is_partial=true;
                continue;
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(output_path, argv[i++]);
            }else if(0==strcmp(argv[i], "-rdf")){
                i++;
                is_rdf=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-rdf_r")){
                        i++;
                        rmax=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-rdf_n")){
                        i++;
                        nrbin=(int)be_double(argv[i++]);
                        continue;
                    }else if(not_in_switch(argv[i], "-rdf")){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h' for more information");
                exit(1);
            }
        }
        char   ssf_path[PATH_CHAR_NUMBER];
        if(0==strcmp(output_path, "")){
            strcpy(ssf_path, "./AAVDP.ssf");
        }else{
            strcpy(ssf_path, output_path);
        }
        if(is_rdf){
            RDF rdf(input_path, rmax, nrbin, is_partial);
            SSF ssf(&rdf, qmax, nqbin, is_partial);
            ssf.ssf(ssf_path);
        }else{
            SSF ssf(input_path, qmax, nqbin, is_partial);
            ssf.ssf(ssf_path);
        }
    }else{
        printf("AAVDP: unrecognized mode '%s'\n", argv[i]);
        printf("Try 'AAVDP -h' for more information'\n");
    }
}
