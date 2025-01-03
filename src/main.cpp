#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "main.h"

int main(int argc, char* argv[])
{
    int    i=1;
    if(argc<2){
        printf("AAVDP: unrecognized mode ''\n");
        printf("Try 'AAVDP -h' or 'AAVDP --help' for more information'\n");
    }
    if(0==strcmp(argv[i], "--xrd")||0==strcmp(argv[i], "--ned")){
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
            }else if(0==strcmp(argv[i], "--scherrer")){
                i++;
                is_scherrer=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-m")){
                        i++;
                        mixing_param=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-d")){
                        i++;
                        grain_diameter=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-d2t")){
                        i++;
                        d2t=be_double(argv[i++]);
                        continue;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h %s' for more information", mode);
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h %s' for more information", mode);
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
            }else if(0==strcmp(argv[i], "--rotate")){
                i++;
                is_rotate=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-x")){
                        i++;
                        xaxis[0]=be_double(argv[i++]);
                        xaxis[1]=be_double(argv[i++]);
                        xaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-y")){
                        i++;
                        yaxis[0]=be_double(argv[i++]);
                        yaxis[1]=be_double(argv[i++]);
                        yaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h ked' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "--gauss")){
                i++;
                is_gauss=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-sig")){
                        i++;
                        sigma=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-dx")){
                        i++;
                        dx=be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h ked' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h ked' for more information");
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
            }else if(0==strcmp(argv[i], "--rotate")){
                i++;
                is_rotate=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-x")){
                        i++;
                        xaxis[0]=be_double(argv[i++]);
                        xaxis[1]=be_double(argv[i++]);
                        xaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-y")){
                        i++;
                        yaxis[0]=be_double(argv[i++]);
                        yaxis[1]=be_double(argv[i++]);
                        yaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h kkd' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "--scale")){
                i++;
                is_scale=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-i")){
                        i++;
                        imin=be_double(argv[i++]);
                        imax=be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h kkd' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "--ked3")){
                i++;
                is_ked3=true;
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h kkd' for more information");
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
            }else if(0==strcmp(argv[i], "--rotate")){
                i++;
                is_rotate=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-x")){
                        i++;
                        xaxis[0]=be_double(argv[i++]);
                        xaxis[1]=be_double(argv[i++]);
                        xaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-y")){
                        i++;
                        yaxis[0]=be_double(argv[i++]);
                        yaxis[1]=be_double(argv[i++]);
                        yaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h ded' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "--gauss")){
                i++;
                is_gauss=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-sig")){
                        i++;
                        sigma=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-dx")){
                        i++;
                        dx=be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h ded' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h ded' for more information");
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
            }else if(0==strcmp(argv[i], "--monte")){
                i++;
                is_monte=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-rd")){
                        i++;
                        omega=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-td")){
                        i++;
                        sigma=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-ex")){
                        i++;
                        Eexit=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-dt")){
                        i++;
                        dthick=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-ne")){
                        i++;
                        nume=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-np")){
                        i++;
                        nump=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-seed")){
                        i++;
                        strcpy(seed_path, argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-o")){
                        i++;
                        strcpy(mc_path, argv[i++]);
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h dkd' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "--rotate")){
                i++;
                is_rotate=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-x")){
                        i++;
                        xaxis[0]=be_double(argv[i++]);
                        xaxis[1]=be_double(argv[i++]);
                        xaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-y")){
                        i++;
                        yaxis[0]=be_double(argv[i++]);
                        yaxis[1]=be_double(argv[i++]);
                        yaxis[2]=be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h dkd' for more information");
                        exit(1);
                    }
                }
            }else if(0==strcmp(argv[i], "--scale")){
                i++;
                is_scale=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-i")){
                        i++;
                        imin=be_double(argv[i++]);
                        imax=be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h dkd' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h dkd' for more information");
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
                printf("Try 'AAVDP -h rdf' for more information");
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
            }else if(0==strcmp(argv[i], "--rdf")){
                i++;
                is_rdf=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-r")){
                        i++;
                        rmax=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-n")){
                        i++;
                        nrbin=(int)be_double(argv[i++]);
                        continue;
                    }else if(is_mode(argv[i])){
                        break;
                    }else{
                        printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                        printf("Try 'AAVDP -h ssf' for more information");
                        exit(1);
                    }
                }
            }else{
                printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                printf("Try 'AAVDP -h ssf' for more information");
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
    }else if(0==strcmp(argv[i], "-h")||0==strcmp(argv[i], "--help")){
        i++;
        if(i==argc){
            printf("The syntax format and rules for AAVDP:\n");
            printf("    AAVDP <--mode> (inputfile) <-parameter> [value]   \n");
            printf("    AAVDP --xrd     # X-ray diffraction               \n");
            printf("    AAVDP --ned     # Neutron diffraction             \n");
            printf("    AAVDP --ked     # Kinematic electron diffraction\n");
            printf("    AAVDP --ded     # Dynamical electron diffraction  \n");
            printf("    AAVDP --kkd     # Kinematic Kikuchi diffraction \n");
            printf("    AAVDP --dkd     # Dynamical Kikuchi diffraction   \n");
            printf("    AAVDP --rdf     # Radial distribution function    \n");
            printf("    AAVDP --ssf     # Static structure factor         \n");
            printf("Try 'AAVDP -h [mode]' for more information, e.g., 'AAVDP -h xrd'\n");
        }else{
            if(0==strcmp(argv[i], "xrd")){
                printf("The syntax format and rules for xrd in AAVDP:\n");
                printf("    AAVDP --xrd [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -l [lambda] -2t [2theta_min] [2theta_max] -lp [Lorentz_polarization_type] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -o [outputfile] --scherrer -m [mixing_parameter] -d [grain_diameter] -d2t [2theta_bin]\n");
                printf("Check /man/manual.pdf for more information\n");
            }else if(0==strcmp(argv[i], "ned")){
                printf("The syntax format and rules for ned in AAVDP:\n");
                printf("    AAVDP --ned [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -l [lambda] -2t [2theta_min] [2theta_max] -lp [Lorentz_type] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -o [outputfile] --scherrer -m [mixing_parameter] -d [grain_diameter] -d2t [2theta_bin]\n");
                printf("Check /man/manual.pdf for more information\n");
            }else if(0==strcmp(argv[i], "ked")){
                printf("The syntax format and rules for ked in AAVDP:\n");
                printf("    AAVDP --ked [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -en [voltage] -z [z_1] [z_2] [z_3] -q [qmax] -t [intersection_thickness] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -o [output_file] --gauss -sig [standard_deviation] -dx [qx_bin] --rotate -x [x_1] [x_2] [x_3] -y [y_1] [y_2] [y_3]\n");
                printf("Check /man/manual.pdf for more information\n");
            }else if(0==strcmp(argv[i], "kkd")){
                printf("The syntax format and rules for ked in AAVDP:\n");
                printf("    AAVDP --kkd [inputfile] -e [type1] [type2] … [typeN] -en [voltage] -q [qmax] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -z [z_1] [z_2] [z_3] -rx [projection_ratiox] -ry [projection_ratioy] -t [Kikuchi_line_thickness] -px [numpx] -py [numpy] -background [background_color] -o [outputfile] --rotate -x [x_1] [x_2] [x_3] -y [y_1] [y_2] [y_3] --scale -i [intensity_min] [intensity_max]\n");
                printf("    or\n");
                printf("    AAVDP --kkd [inputfile] -e [type1] [type2] … [typeN] -en [voltage] -q [qmax] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -o [ked3file] --ked3\n");
                printf("    AAVDP --kkd [ked3file] -z [z_1] [z_2] [z_3] -rx [projection_ratiox] -ry [projection_ratioy] -t [Kikuchi_line_thickness] -px [numpx] -py [numpy] -background [background_color] -o [outputfile] --rotate -x [x_1] [x_2] [x_3] -y [y_1] [y_2] [y_3] --scale -i [intensity_min] [intensity_max]\n");
                printf("Check /man/manual.pdf for more information\n");
            }else if(0==strcmp(argv[i], "ded")){
                printf("The syntax format and rules for ded in AAVDP:\n");
                printf("    AAVDP --ded [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] ... -en [voltage] -z [z_1] [z_2] [z_3] -q [qmax] -fn [n_1] [n_2] [n_3] -ft [foil_thickness] -bethe [rg_c1] [rg_c2] [rg_c3] [sg_c] -thr [threshold] -o [outputfile] --gauss -sig [standard_deviation] -dx [qx_bin] --rotate -x [x_1] [x_2] [x_3] -y [y_1] [y_2] [y_3]\n");
                printf("Check /man/manual.pdf for more information\n");
            }else if(0==strcmp(argv[i], "dkd")){
                printf("The syntax format and rules for dkd in AAVDP:\n");
                printf("    AAVDP --dkd [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -en [voltage] -q [qmax] -ft [foil_thickness] -bethe [rg_c1] [rg_c2] [rg_c3] [sg_c] -z [z_1] [z_2] [z_3] -rx [projection_ratiox] -ry [projection_ratioy] -px [numpx] -py [numpy] -background [background_color] -o [outputfile] --monte -rd [rotation_angle] -td [tilt_angle] -ex [energy_exit] -dt [foil_depth_bin] -ne [nume] -np [nump] -seed [randomfile] -o [outputfile] --rotate -x [x_1] [x_2] [x_3] -y [y_1] [y_2] [y_3] --scale -i [intensity_min] [intensity_max]\n");
                printf("Check /man/manual.pdf for more information\n");
            }else if(0==strcmp(argv[i], "rdf")){
                printf("The syntax format and rules for dkd in AAVDP:\n");
                printf("    AAVDP --rdf [inputfile] -r [rmax] -n [nbin] -partial -o [outputfile]\n");
                printf("Check /man/manual.pdf for more information\n");
            }else if(0==strcmp(argv[i], "ssf")){
                printf("The syntax format and rules for dkd in AAVDP:\n");
                printf("    AAVDP --ssf [inputfile] -q [qmax] -n [nbin] -partial -o [outputfile] --rdf -r [rmax] -n [nbin]\n");
                printf("Check /man/manual.pdf for more information\n");
            }else{
                printf("AAVDP: unrecognized mode '%s'\n", argv[i]);
                printf("Try 'AAVDP -h' or 'AAVDP --help' for more information'\n");
            }
        }
    }else if(0==strcmp(argv[i], "-v")||0==strcmp(argv[i], "--version")){
        printf("***------------------------------------------------AAVDP Version 0.0.4 (2024.10.30)-----------------------------------------------***\n");
        printf("*** An integrated command-line program for automatic analysis of virtual diffraction patterns of artificial atomistic structures. ***\n");
        printf("***                           Copyright[c] 2022-2024, Beihang University by Zhang Yan and Ruifeng Zhang.                          ***\n");
        printf("***                                     Please send bugs and suggestions to zrfcms@buaa.edu.cn                                    ***\n");
        printf("***-------------------------------------------------------------------------------------------------------------------------------***\n");
    }else{
        printf("AAVDP: unrecognized mode '%s'\n", argv[i]);
        printf("Try 'AAVDP -h' or 'AAVDP --help' for more information'\n");
    }
}
