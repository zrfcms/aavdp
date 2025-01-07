#ifndef __AAVDP_MAIN__
#define __AAVDP_MAIN__
#include <unistd.h>
#include "./MODEL/MODEL.h"
#include "./XRD/XRD.h"
#include "./KED/KED.h"
#include "./DED/DED.h"
#include "./RDF/RDF.h"

#define TYPE_INPUT_NUMBER 10
#define PATH_CHAR_NUMBER 100
#define EXT_CHAR_NUMBER 20

extern void print_version();
extern void print_help();
extern bool is_path_accessible(char path[]);
extern double be_double(char str[]);
extern bool is_parameter(char str[]);
extern bool not_in_switch(char option[], char swich[]);
extern void merge_path(char file_path[], char exts[][EXT_CHAR_NUMBER], int num);
extern void split_path(char name[], char ext[], char file_path[]);

void print_version()
{
    printf("AAVDP Version 1.0.0 (2025.01.04)\n"
           "An integrated command-line program for automatic analysis of virtual diffraction patterns of artificial atomistic structures.\n"
           "Copyright[c] 2025-2027, Beihang University by Yan Zhang and Ruifeng Zhang.\n"
           "Please send bugs and suggestions to zrfcms@buaa.edu.cn.\n");
}

void print_help()
{
    printf("The syntax format and rules for AAVDP:\n"
           "    AAVDP <--mode> (inputfile) <-parameter> [value]   \n"
           "    AAVDP --xrd     # X-ray diffraction               \n"
           "    AAVDP --ned     # Neutron diffraction             \n"
           "    AAVDP --ked     # Kinematic electron diffraction\n"
           "    AAVDP --ded     # Dynamical electron diffraction  \n"
           "    AAVDP --kkd     # Kinematic Kikuchi diffraction \n"
           "    AAVDP --dkd     # Dynamical Kikuchi diffraction   \n"
           "    AAVDP --rdf     # Radial distribution function    \n"
           "    AAVDP --ssf     # Static structure factor         \n"
           "\n"
           "The syntax format and rules for xrd in AAVDP:\n"
           "    AAVDP --xrd [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -l [lambda] -2t [2theta_min] [2theta_max] -lp [Lorentz_polarization_type] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -o [outputfile] -scherrer -scherrer_m [mixing_parameter] -scherrer_d [grain_diameter] -scherrer_d2t [2theta_bin]\n"
           "\n"
           "The syntax format and rules for ned in AAVDP:\n"
           "    AAVDP --ned [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -l [lambda] -2t [2theta_min] [2theta_max] -lp [Lorentz_type] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -o [outputfile] -scherrer -scherrer_m [mixing_parameter] -scherrer_d [grain_diameter] -scherrer_d2t [2theta_bin]\n"
           "\n"
           "The syntax format and rules for ked in AAVDP:\n"
           "    AAVDP --ked [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -en [voltage] -z [z_1] [z_2] [z_3] -q [qmax] -t [intersection_thickness] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -o [output_file] -gauss -gauss_sig [standard_deviation] -gauss_dx [qx_bin] -rotate -rotate_x [x_1] [x_2] [x_3] -rotate_y [y_1] [y_2] [y_3]\n"
           "\n"
           "The syntax format and rules for kkd in AAVDP:\n"
           "    AAVDP --kkd [inputfile] -e [type1] [type2] … [typeN] -en [voltage] -q [qmax] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -z [z_1] [z_2] [z_3] -rx [projection_ratiox] -ry [projection_ratioy] -t [Kikuchi_line_thickness] -px [numpx] -py [numpy] -background [background_color] -o [outputfile] -rotate -rotate_x [x_1] [x_2] [x_3] -rotate_y [y_1] [y_2] [y_3] -scale -scale_i [intensity_min] [intensity_max]\n"
           "    or\n"
           "    AAVDP --kkd [inputfile] -e [type1] [type2] … [typeN] -en [voltage] -q [qmax] -c [c1] [c2] [c3] -auto/manu -thr [threshold] -o [ked3file] -ked3\n"
           "    AAVDP --kkd [ked3file] -z [z_1] [z_2] [z_3] -rx [projection_ratiox] -ry [projection_ratioy] -t [Kikuchi_line_thickness] -px [numpx] -py [numpy] -background [background_color] -o [outputfile] -rotate -rotate_x [x_1] [x_2] [x_3] -rotate_y [y_1] [y_2] [y_3] -scale -scale_i [intensity_min] [intensity_max]\n"
           "\n"
           "The syntax format and rules for ded in AAVDP:\n"
           "    AAVDP --ded [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] ... -en [voltage] -z [z_1] [z_2] [z_3] -q [qmax] -fn [n_1] [n_2] [n_3] -ft [foil_thickness] -bethe [rg_c1] [rg_c2] [rg_c3] [sg_c] -thr [threshold] -o [outputfile] -gauss -gauss_sig [standard_deviation] -gauss_dx [qx_bin] -rotate -rotate_x [x_1] [x_2] [x_3] -rotate_y [y_1] [y_2] [y_3]\n"
           "\n"
           "The syntax format and rules for dkd in AAVDP:\n"
           "    AAVDP --dkd [inputfile] -e [type1] [type2] … [typeN] -dw [DW1] [DW2] … [DWN] -en [voltage] -q [qmax] -ft [foil_thickness] -bethe [rg_c1] [rg_c2] [rg_c3] [sg_c] -z [z_1] [z_2] [z_3] -rx [projection_ratiox] -ry [projection_ratioy] -px [numpx] -py [numpy] -background [background_color] -o [outputfile] -monte -monte_rd [rotation_angle] -monte_td [tilt_angle] -monte_ex [energy_exit] -monte_dt [foil_depth_bin] -monte_ne [nume] -monte_np [nump] -monte_seed [randomfile] -o [outputfile] -rotate -rotate_x [x_1] [x_2] [x_3] -rotate_y [y_1] [y_2] [y_3] -scale -scale_i [intensity_min] [intensity_max]\n"
           "\n"
           "The syntax format and rules for dkd in AAVDP:\n"
           "    AAVDP --rdf [inputfile] -r [rmax] -n [nbin] -partial -o [outputfile]\n"
           "\n"
           "The syntax format and rules for dkd in AAVDP:\n"
           "    AAVDP --ssf [inputfile] -q [qmax] -n [nbin] -partial -o [outputfile] -rdf -rdf_r [rmax] -rdf_n [nbin]\n"
           "Check /man/manual.pdf for more information\n");
}

bool is_path_accessible(char path[])
{
    if(access(path, F_OK)==-1){
        printf("[ERROR] Path %s does not exist", path);
        exit(1);
    }else if(access(path, R_OK)==-1){
        printf("[ERROR] Path %s is not readable", path);
        exit(1);
    }
    return true;
}

double be_double(char str[])
{
    char *err;
    double value=strtod(str, &err);
    if(strlen(err)!='\0')
    {
        printf("[ERROR] Input parameter %s should be float or int type", str);
        exit(1);
    }
    return value;
}

bool is_parameter(char str[])
{
    if(0==strncmp(str, "-", 1)){
        return false;
    }
    return true;
}

bool not_in_switch(char option[], char swich[])
{
    if((!is_parameter(option))&&strstr(option, swich)==nullptr){
        return true;
    }
    return false;
}

void merge_path(char file_path[], char exts[][EXT_CHAR_NUMBER], int num)
{
    for(int i=0;i<num;i++){
        strcat(file_path, exts[i]);
    }
}

void split_path(char name[], char ext[], char file_path[])
{
    char *ch=strrchr(file_path,'.');
    strcpy(name, file_path);
    strcpy(ext, ch);
    if(ch!=nullptr){
        int pos=ch-file_path;
        name[pos]='\0';
    }
}

#endif