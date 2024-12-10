#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "main.h"

int main(int argc, char* argv[])
{
    int    i=1;
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
        double max2Theta=180.0;
        int    lp_type=3;
        if(0==strcmp(mode, "ned")) lp_type=2;
        double c[3]={1.0, 1.0, 1.0};
        bool   is_spacing_auto=true;
        bool   is_c_default=true;
        double threshold=0.001;

        double mixing_param=0.0;
        double diameter=500.0;
        double FWHM=1.0;
        double d2t=0.01;
        bool   is_scherrer=false;
        bool   is_smear=false;
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
                        diameter=be_double(argv[i++]);
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
            }else if(0==strcmp(argv[i], "--smear")){
                i++;
                is_smear=true;
                while(i<argc){
                    if(0==strcmp(argv[i], "-m")){
                        i++;
                        mixing_param=be_double(argv[i++]);
                        continue;
                    }else if(0==strcmp(argv[i], "-w")){
                        i++;
                        FWHM=be_double(argv[i++]);
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
                xrd.xrd(xrd_path, mixing_param, lambda, diameter, d2t);
            }else if(is_smear){
                xrd.xrd(xrd_path, mixing_param, FWHM, d2t);
            }else{
                xrd.xrd(xrd_path);
            }
        }else{
            NMODEL model(input_path, types, DWs, lambda);
            model.update_lorentzP_type(lp_type);
            XRD xrd(&model, min2Theta, max2Theta, threshold, c, is_spacing_auto);
            if(is_scherrer){
                xrd.xrd(xrd_path, mixing_param, lambda, diameter, d2t);
            }else if(is_smear){
                xrd.xrd(xrd_path, mixing_param, FWHM, d2t);
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
        double Kmax=1.0;
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
                Kmax=be_double(argv[i++]);
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
        KED ked(&model, zaxis, thickness, Kmax, threshold, c, is_spacing_auto);
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
        double Kmax=1.0;
        double c[3]={1.0, 1.0, 1.0};
        bool   is_spacing_auto=true;
        bool   is_c_default=true;
        double threshold=0.001;

        double xaxis[3]={1.0, 0.0, 0.0};
        double yaxis[3]={0.0, 1.0, 0.0};
        double zaxis[3]={0.0, 0.0, 1.0};
        double ratiox=1.0, ratioy=1.0;
        double thickness=0.2;
        int    npx=500, npy=500;
        bool   is_stereo_proj=true;

        double vmax=1.0e6, vmin=0.0;
        char   background='b';
        bool   is_scale=false;
        if(0!=strcmp(ext, ".ked3")){
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
                    Kmax=be_double(argv[i++]);
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
                }else{
                    printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                    printf("Try 'AAVDP -h kkd' for more information");
                    exit(1);
                }
            }
            char ked3_path[PATH_CHAR_NUMBER];
            if(0==strcmp(output_path, "")){
                strcpy(ked3_path, "./AAVDP.ked3");
            }else{
                strcpy(ked3_path, output_path);
            }
            VMODEL model(input_path, types, DWs, voltage);
            KED ked(&model, Kmax, threshold, c, is_spacing_auto);
            ked.ked3(ked3_path);
        }else{
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
                }else if(0==strcmp(argv[i], "-stereo")){
                    i++;
                    is_stereo_proj=true;
                    continue;
                }else if(0==strcmp(argv[i], "-ortho")){
                    i++;
                    is_stereo_proj=false;
                    continue;
                }else if(0==strcmp(argv[i], "-background")){
                    i++;
                    background=(argv[i++])[0];
                    continue;
                }else if(0==strcmp(argv[i], "-thr")){
                    i++;
                    threshold=be_double(argv[i++]);
                    continue;
                }else if(0==strcmp(argv[i], "-o")){
                    i++;
                    strcpy(output_path, argv[i++]);
                    continue;
                }else if(0==strcmp(argv[i], "--scale")){
                    i++;
                    is_scale=true;
                    while(i<argc){
                        if(0==strcmp(argv[i], "-max")){
                            i++;
                            vmax=be_double(argv[i++]);
                            continue;
                        }else if(0==strcmp(argv[i], "-min")){
                            i++;
                            vmin=be_double(argv[i++]);
                            continue;
                        }else if(is_mode(argv[i])){
                            break;
                        }else{
                            printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                            printf("Try 'AAVDP -h kkd' for more information");
                            exit(1);
                        }
                    }
                }else{
                    printf("AAVDP: unrecognized option '%s'\n", argv[i]);
                    printf("Try 'AAVDP -h kkd' for more information");
                    exit(1);
                }
            }
            char kkd_path[PATH_CHAR_NUMBER];
            if(0==strcmp(output_path, "")){
                strcpy(kkd_path, "./AAVDP.kkd");
            }else{
                strcpy(kkd_path, output_path);
            }
            KKD kkd(input_path, xaxis, yaxis, zaxis, thickness, threshold, ratiox, ratioy, npx, npy, is_stereo_proj);
            if(is_scale){
                kkd.kkd(kkd_path, vmax, vmin, background);
            }else{
                kkd.kkd(kkd_path, background);
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
        double Kmax=1.0;
        int    zone[3]={0, 0, 1};
        int    fnorm[3]={0, 0, 1};
        double fthick=100.0;
        double c1=8.0, c2=50.0, c3=100.0, c_sg=1.0;
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
                Kmax=be_double(argv[i++]);
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
        DED ded(&cell, &bethe, zone, fnorm, fthick, voltage, Kmax, threshold);
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

        char   types[TYPE_INPUT_NUMBER][10]={0};
        double DWs[TYPE_INPUT_NUMBER]={0.0};
        for(int i=0;i<TYPE_INPUT_NUMBER;i++) DWs[i]=1.0e-6;
        double omega=0.0, sigma=75.7;
        double Emax=20.0, Emin=19.0;
        double zmax=100.0, zstep=1.0;
        int    nume=20000;
        int    nump=501;
        double dmin=0.1;
        double c1=4.0, c2=8.0, c3=50.0, c_sg=1.0;
        double img_L=6.0;
        int    img_dpi=512;
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
                Emax=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-ex")){
                i++;
                Emin=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-z0")){
                i++;
                zmax=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-dz")){
                i++;
                zstep=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-ne")){
                i++;
                nume=(int)be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-np")){
                i++;
                nump=(int)be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-dm")){
                i++;
                dmin=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-bethe")){
                i++;
                c1=be_double(argv[i++]);
                c2=be_double(argv[i++]);
                c3=be_double(argv[i++]);
                c_sg=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-w")){
                i++;
                img_L=be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-dpi")){
                i++;
                img_dpi=(int)be_double(argv[i++]);
            }else if(0==strcmp(argv[i], "-o")){
                i++;
                strcpy(output_path, argv[i++]);
            }else{
                i++;
                continue;
            }
        }
        char   mc_path[PATH_CHAR_NUMBER];
        char   LPNH_path[PATH_CHAR_NUMBER], LPSH_path[PATH_CHAR_NUMBER], SPNH_path[PATH_CHAR_NUMBER], SPSH_path[PATH_CHAR_NUMBER];
        if(0==strcmp(output_path, "")){
            strcpy(mc_path, "./AAVDP.ded.mc.png");
            strcpy(LPNH_path, "./AAVDP.ded.LPNH.png"); strcpy(LPSH_path, "./AAVDP.ded.LPSH.png"); strcpy(SPNH_path, "./AAVDP.ded.SPNH.png"); strcpy(SPSH_path, "./AAVDP.ded.SPSH.png");
        }else{
            strcpy(mc_path, output_path); strcat(mc_path, ".mc.png");
            strcpy(LPNH_path, output_path); strcat(LPNH_path, ".LPNH.png");
            strcpy(LPSH_path, output_path); strcat(LPSH_path, ".LPSH.png");
            strcpy(SPNH_path, output_path); strcat(SPNH_path, ".SPNH.png");
            strcpy(SPSH_path, output_path); strcat(SPSH_path, ".SPSH.png");
        }
        BETHE bethe;
        bethe.c1=c1; bethe.c2=c2; bethe.c3=c3; bethe.c_sg=c_sg;
        CELL cell(input_path, types, DWs, Emax);
        DED_MC mc(&cell, omega, sigma, Emax, Emin, zmax, zstep, nume, nump);
        mc.img(mc_path, img_L, img_dpi);
        DKD dkd(&cell, &mc, &bethe, Emax, dmin, nump);
        dkd.img(LPNH_path, LPSH_path, SPNH_path, SPSH_path, img_L, img_dpi);
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
        char   png_path[PATH_CHAR_NUMBER]; strcpy(png_path, name); 
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
        RDF rdf(input_path, 7.5, 150, is_partial);
        SSF ssf(&rdf, qmax, nbin, is_partial);
        ssf.ssf(ssf_path, png_path);
    }else if(0==strcmp(argv[i], "-h")||0==strcmp(argv[i], "--help")){
        i++;
        if(i==argc){
            printf("The syntax format and rules for AAVDP:\n");
            printf("    AAVDP <--mode> (inputfile) <-parameter> [value]   \n");
            printf("    AAVDP --xrd     # X-ray diffraction               \n");
            printf("    AAVDP --ned     # Neutron diffraction             \n");
            printf("    AAVDP --ked     # Kinematical electron diffraction\n");
            printf("    AAVDP --ded     # Dynamical electron diffraction  \n");
            printf("    AAVDP --kkd     # Kinematical Kikuchi diffraction \n");
            printf("    AAVDP --dkd     # Dynamical Kikuchi diffraction   \n");
            printf("    AAVDP --rdf     # Radial distribution function    \n");
            printf("    AAVDP --ssf     # Static structure factor         \n");
            printf("Try 'AAVDP -h [mode]' for more information, e.g., 'AAVDP -h xrd'\n");
        }else{
            if(0==strcmp(argv[i], "xrd")){
                printf("The syntax format and rules for xrd in AAVDP:\n");
                printf("    AAVDP --xrd [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -l [lambda] -2t [2theta_min] [2theta_max] -lp [lorentz_type] -c [c1] [c2] [c3] -auto/manu -o [outputfile] --smear -m [mixing_parameter] -d [scherrer_diameter] -d2t [2theta_bin]\n");
                printf("    e.g., AAVDP --xrd Ni.vasp              # The prerequisite is that atomic types have been stored in Ni.vasp\n");
                printf("    e.g., AAVDP --xrd Ni.lmp -e Ni\n");
                printf("    e.g., AAVDP --xrd Ni.lmp -e Ni --smear\n");
                printf("Check /man/manual.pdf for more information\n");
            }else if(0==strcmp(argv[i], "ned")){
                printf("The syntax format and rules for ned in AAVDP:\n");
                printf("    AAVDP --ned [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -l [lambda] -2t [2theta_min] [2theta_max] -lp [lorentz_type] -c [c1] [c2] [c3] -auto/manu -o [outputfile] --smear -m [mixing_parameter] -d [scherrer_diameter] -d2t [2theta_bin]\n");
                printf("    e.g., AAVDP --ned Ni.vasp              # The prerequisite is that atomic types have been stored in Ni.vasp\n");
                printf("    e.g., AAVDP --ned Ni.lmp -e Ni\n");
                printf("    e.g., AAVDP --ned Ni.lmp -e Ni --smear\n");
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
