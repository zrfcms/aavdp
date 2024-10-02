#ifndef __AAVDP_MAIN__
#define __AAVDP_MAIN__
#include "./MODEL/MODEL.h"
#include "./XRD/XRD.h"
#include "./NED/NED.h"
#include "./KED/KED.h"
#include "./DED/DED.h"
#include "./KKD/KKD.h"
#include "./DKD/DKD.h"
#include "./RDF/RDF.h"
#include "./SSF/SSF.h"

#define TYPE_INPUT_NUMBER 10
#define ZERO_LIMIT 1.0e-4

extern double be_double(char str[]);
extern bool is_parameter(char str[]);

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