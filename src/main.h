#ifndef __AAVDP_MAIN__
#define __AAVDP_MAIN__
#include <unistd.h>
#include "./MODEL/MODEL.h"
#include "./XRD/XRD.h"
#include "./KED/KED.h"
#include "./DED/DED.h"
#include "./KKD/KKD.h"
#include "./RDF/RDF.h"
#include "./SSF/SSF.h"

#define TYPE_INPUT_NUMBER 10

extern bool is_path_accessible(char path[]);
extern double be_double(char str[]);
extern bool is_parameter(char str[]);

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

#endif