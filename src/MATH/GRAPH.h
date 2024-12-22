#ifndef __AAVDP_GRAPH_H__
#define __AAVDP_GRAPH_H__
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "MATH.h"
#include "../include/png.h"

#define PI_SQRT_HALF 0.886226925452758 //sqrt(PI)/2
#define PI_HALF_SQRT 1.253314137315500 //sqrt(PI/2)
#define PI_INVERSE 0.318309886183791 //1/PI
#define SQRT_3_INVERSE 0.577350269190 //1/sqrt(3)

#define PI_PREE 0.7598356856515930 //3^(-1/4)
#define PI_PREF 1.3819765978853420 //sqrt(6/PI)
#define PI_PREA 0.5250375679043320 //3^(1/4)/sqrt(2*PI)
#define PI_PREG 1.55512030155621410 //2*sqrt(PI)/3^(3/4)
#define PI_PREC 0.9068996821171090 //PI/2*sqrt(3)
#define PI_PREB 1.0500751358086640 //3^(1/4)*sqrt(2/PI)
#define PI_PRED 2.0943951023931950 //2*PI/3
#define SQRT_3 1.732050807568877 //sqrt(3)

#define SQRT_HALF_3 0.866025403780 //sqrt(3)/2

// extern void compute_Lambert_Projection(double xy[2], int &ierr, double xyz[3]);
extern void compute_square_Lambert(double xy[2], int &ierr, double xyz[3]);
extern void compute_sphere_from_square_Lambert(double xyz[3], int &ierr, double xy[2]);
extern  int get_sextant(double x, double y);
extern void compute_hexagonal_Lambert(double xy[2], int &ierr, double xyz[3]);
extern void compute_sphere_from_hexagonal_Lambert(double xyz[3], int &ierr, double xy[2]);
extern void compute_sphere_from_stereographic_projection(double xyz[3], int &ierr, double xy[2]);
extern void compute_sphere_from_orthographic_projection(double xyz[3], int &ierr, double xy[2]);

extern void image_pixels(char* png_path, unsigned char *pixels, int nrow, int ncol);
extern void image_array(char* png_path, double **value, double vmax, double vmin, int nrow, int ncol, char background);

struct POINT
{
    double x=0.0, y=0.0;
    double black=1.0;
    char   style='c';
    double psize=5.0;
    bool   is_in_area_checked=true;
};

struct AXIS
{
    int    n_major_tick;
    double *major_ticks;
    double major_tick_psize=20.0;
    bool   is_tick_in=true;
    bool   is_visible=true;
    double black=1.0;
    double line_pwidth=5.0;
    bool   is_in_area_checked=false;
};

struct AREA
{
    double border_width_factor=0.1;
    double border_height_factor=0.1;
    int    pwidth, pheight;
    int    start_px, start_py;
    int    end_px, end_py;
    double spacing_x, spacing_y;
    double start_x, end_x;
    double start_y, end_y;
    double ratio_x, ratio_y;
};

struct IMAGE
{
    int    pwidth, pheight;
    double black=0.0;
    double **pixels=nullptr;
};

class GRAPH
{
public:
    char   style[10];
    IMAGE  image;
    AREA   area;
    AXIS   left, right, top, down;
    int    nump=0;
    POINT  *points;
    GRAPH(double width, double height, int resolution);
    ~GRAPH();
    void set_xlim(double xmin, double xmax);
    void set_ylim(double ymin, double ymax);
    void set_xticks(double *major_ticks, double n_major_tick);
    void set_yticks(double *major_ticks, double n_major_tick);
    void set_tick_in(bool is_tick_in);
    void set_top_visible(bool is_visible);
    void set_right_visible(bool is_visible);
    void scatter(double *x, double *y, double *value, int num, char marker_style='c', double marker_size=20.0);
    void line(double *x, double *y, int num, double line_width=5.0);
    void hist(double *x, double *y, int num, double line_width=5.0);
    void draw(char *png_path);
private:
    void get_pixel(int &i, int &j, double x, double y);
    void draw_marker(POINT *point);
    void draw_line(POINT *start_point, POINT *end_point);
    void draw_xaxis(AXIS *axis, double x1, double x2, double y);
    void draw_yaxis(AXIS *axis, double y1, double y2, double x);
    // inline void rgb_to_hsl(int hsl[3], int rgb[3]){
    //     double r=rgb[0]/255.0;
    //     double g=rgb[1]/255.0;
    //     double b=rgb[2]/255.0;
    //     double max_rgb=max3(r, g, b);
    //     double min_rgb=min3(r, g, b);

    //     double h, s, l;
    //     int    hsl[3];
    //     l=(max_rgb+min_rgb)/2.0;
    //     if(l<0.5){
    //         s=(max_rgb-min_rgb)/(max_rgb+min_rgb);
    //     }else{
    //         s=(max_rgb-min_rgb)/(2.0-max_rgb-min_rgb);
    //     }
    //     if(r==max_rgb){
    //         h=60.0*(g-b)/(max_rgb-min_rgb);
    //     }else if(g==max_rgb){
    //         h=60.0*((b-r)/(max_rgb-min_rgb)+2.0);
    //     }else{
    //         h=60.0*((r-g)/(max_rgb-min_rgb)+4.0);
    //     }
    //     hsl[0]=round(h*255.0);
    //     hsl[1]=round(s*255.0);
    //     hsl[2]=round(l*255.0);
    // };
    // inline void hsl_to_rgb(int rgb[3], int hsl[3]){
    //     double h=hsl[0]/255.0;
    //     double s=hsl[1]/255.0;
    //     double l=hsl[2]/255.0;
    //     int    f=int(h/60.0);
    //     double p=l*(1-s), q=l*(1-s*f);

    //     double r, g, b;
    //     int    rgb[3];
    //     if(f==0){
    //         r=q; g=p; b=p;
    //     }else if(f==1){
    //         r=p; g=q; b=p;
    //     }else if(f==2){
    //         r=p; g=q; b=p;
    //     }
    //     rgb[0]=round(r*255.0);
    //     rgb[1]=round(g*255.0);
    //     rgb[2]=round(b*255.0);
    // };
};

#endif