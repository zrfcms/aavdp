#include "MATH.h"

void matrix_multiply(double mat[3][3], double mat1[3][3], double mat2[3][3])
{
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            mat[i][j]=0.0;
            for(int k=0;k<3;k++){
                mat[i][j]+=mat1[i][k]*mat2[k][j];
            }
        }
    }
}

void vector_rotate(double r_v[3], double v[3], double R[3][3])
{
    double c_v[3]={v[0], v[1], v[2]};
    for(int i=0;i<3;i++){
        r_v[i]=0.0;
        for(int j=0;j<3;j++){
            r_v[i]+=R[i][j]*c_v[j];
        }
    }
}

double vector_dot(double v1[3], double v2[3])
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

void vector_cross(double v[3], double v1[3], double v2[3])
{
    v[0]=v1[1]*v2[2]-v1[2]*v2[1];
    v[1]=v1[2]*v2[0]-v1[0]*v2[2];
    v[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

double vector_length(double v[3])
{
    return sqrt(vector_dot(v, v));
}

void vector_differance(double v[3], double v1[3], double v2[3])
{
    v[0]=v1[0]-v2[0]; v[1]=v1[1]-v2[1]; v[2]=v1[2]-v2[2];
}

void vector_normalize(double n_v[3], double v[3])
{
    double len_v=vector_length(v);
    if(0.0!=len_v){
        n_v[0]=v[0]/len_v; n_v[1]=v[1]/len_v; n_v[2]=v[2]/len_v;
    }else{
        n_v[0]=0.0; n_v[1]=0.0; n_v[2]=0.0; 
    }
}

QUATERNION quate_conjg(QUATERNION q)
{
    QUATERNION qc={q.c1, -q.c2, -q.c3, -q.c4};
    return qc;
}

QUATERNION quate_multi(QUATERNION q1, QUATERNION q2)
{
    QUATERNION qm;
    double eps=1.0;
    qm.c1=q1.c1*q2.c1-q1.c2*q2.c2-q1.c3*q2.c3-q1.c4*q2.c4;
    qm.c2=q1.c1*q2.c2+q1.c2*q2.c1+eps*(q1.c3*q2.c4-q1.c4*q2.c3);
    qm.c3=q1.c1*q2.c3+q1.c3*q2.c1+eps*(q1.c4*q2.c2-q1.c2*q2.c4);
    qm.c4=q1.c1*q2.c4+q1.c4*q2.c1+eps*(q1.c2*q2.c3-q1.c3*q2.c2); 
    return qm;
}

void quate_rotate(double r_v[3], double v[3], QUATERNION q)
{
    QUATERNION qv={0.0, v[0], v[1], v[2]};
    QUATERNION r_qv=quate_multi(q, quate_multi(qv, quate_conjg(q)));
    r_v[0]=r_qv.c2; r_v[1]=r_qv.c3; r_v[2]=r_qv.c4; 
}

QUATERNION quate_convert(double v[3], double angle)
{
    QUATERNION q;
    if(fabs(angle)<1e-12){
        q.c1=1.0; q.c2=q.c3=q.c4=0.0;
    }else{
        double c=cos(0.5*angle);
        double s=sin(0.5*angle);
        q.c1=c; q.c2=s*v[0]; q.c3=s*v[1]; q.c4=s*v[2];
    }
    return q;
}