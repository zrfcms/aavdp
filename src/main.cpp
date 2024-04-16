#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "main.h"

double be_double(char str[])
{
    char *err;
    double value=strtod(str, &err);
    if(strlen(err)!='\0')
    {
        printf("Error! Input parameter %s should be float or int type", str);
        exit(EXIT_FAILURE);
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

void extension_delete(char c_str[], char str[])
{
    strcpy(c_str, str);
    int pos=strrchr(c_str,'.')-c_str;
    c_str[pos]='\0';
}

int main(int argc, char* argv[])
{
    int i=1;
    if(0==strcmp(argv[i], "-kkd")){
        i++;
        char   sed_path[PATH_CHAR_NUMBER];
        strcpy(sed_path, argv[i++]);
        SED    sed(sed_path)
        return 0;
    }else if(0==strcmp(argv[i], "-sed")){
        i++;
        char   model_path[PATH_CHAR_NUMBER];
        strcpy(model_path, argv[i++]);
        char   types[10][10]={0};
        double Kmagnitude_max=1.70;
        double lambda=0.0251;
        double spacing[3]={1.0, 1.0, 1.0};
        int    zone[3]={0, 0, 0};
        double intercept_radius=1.0/lambda;
        char   name[PATH_CHAR_NUMBER];
        extension_delete(name, model_path);
        char   vtk_path[PATH_CHAR_NUMBER]; strcpy(vtk_path, name);
        strcat(vtk_path, ".sed.vtk");
        char   text_path[PATH_CHAR_NUMBER]; strcpy(text_path, name);
        strcat(text_path, ".sed.vtk");
        char   log_path[PATH_CHAR_NUMBER]; strcpy(log_path, name);
        strcat(log_path, ".sed.vtk");
        bool   is_spacing_ratio=false;
        bool   is_zone_existing=false;
        while(i<argc){
            if(0==strcmp(argv[i], "-e")){
                i++;
                int count=0;
                while(i<argc&&is_parameter(argv[i])){
                    strcpy(types[count++], argv[i++]);
                }
                continue;
            }else if(0==strcmp(argv[i], "-Kmax")){
                i++;
                Kmagnitude_max=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-lambda")){
                i++;
                lambda=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-s")){
                i++;
                spacing[0]=be_double(argv[i++]);
                spacing[1]=be_double(argv[i++]);
                spacing[2]=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-sr")){
                i++;
                spacing[0]=be_double(argv[i++]);
                spacing[1]=be_double(argv[i++]);
                spacing[2]=be_double(argv[i++]);
                is_spacing_ratio=true;
                continue;
            }else if(0==strcmp(argv[i], "-z")){
                i++;
                zone[0]=(int)be_double(argv[i++]);
                zone[1]=(int)be_double(argv[i++]);
                zone[2]=(int)be_double(argv[i++]);
                intercept_radius=1.0/lambda;
                is_zone_existing=true;
                continue;
            }else if(0==strcmp(argv[i], "-dr")){
                i++;
                intercept_radius=be_double(argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-vtk")){
                i++;
                strcpy(vtk_path, argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-text")){
                i++;
                strcpy(text_path, argv[i++]);
                continue;
            }else if(0==strcmp(argv[i], "-log")){
                i++;
                strcpy(log_path, argv[i++]);
                continue;
            }
        }
        EMODEL model(model_path, types, lambda);
        if(is_spacing_ratio){
            model.compute_reciprocal_spacing(spacing, spacing);
        }
        if(is_zone_existing){
            SED sed(&model, spacing, zone, intercept_radius, Kmagnitude_max);
            sed.result(vtk_path, text_path);
        }else{
            SED sed(&model, spacing, Kmagnitude_max);
            sed.result(vtk_path, text_path);
        }
    }else{
        printf("Error! Unrecognized command %s.\n", argv[1]);
    }
}
