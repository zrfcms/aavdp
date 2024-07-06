#include "GRAPH.h"

void image_pixels(const char* png_path, unsigned char *pixels, int numpx, int numpy)
{	
    png_structp png_ptr;  
    png_infop info_ptr;  
    FILE *fp=fopen(png_path, "wb");  
    if(!fp){
        printf("[ERROR] Unable to create %s using fopen.\n", png_path);
        exit(EXIT_FAILURE);
    }
    png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);  
    if(png_ptr==NULL){  
        printf("[ERROR] Unable to create %s using png_create_write_struct.\n", png_path);
        fclose(fp);
        exit(EXIT_FAILURE);
    }  
    info_ptr=png_create_info_struct(png_ptr);
    if(info_ptr==NULL){  
        printf("[ERROR] Unable to create %s using png_create_info_struct.\n", png_path);
        png_destroy_write_struct(&png_ptr, NULL);  
        exit(EXIT_FAILURE);
    }  
    png_init_io(png_ptr, fp);  
    png_set_IHDR(png_ptr, info_ptr, numpx, numpy, 8, 
                 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE); 
    png_colorp palette=(png_colorp)png_malloc(png_ptr, PNG_MAX_PALETTE_LENGTH*sizeof(png_color));
    if(!palette){
        printf("[ERROR] Unable to create palette using png_malloc.\n");
        fclose(fp);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        exit(EXIT_FAILURE);
    }
    png_set_PLTE(png_ptr, info_ptr, palette, PNG_MAX_PALETTE_LENGTH);  
    png_write_info(png_ptr, info_ptr);  
    png_set_packing(png_ptr);
    png_bytepp rows=(png_bytepp)png_malloc(png_ptr, numpy*sizeof(png_bytep));
    for(int i=0;i<numpy;++i){
        rows[numpy-i-1]=(png_bytep)(pixels+(i)*numpx*3);
    }
    png_write_image(png_ptr, rows);  
    delete[] rows;
    png_write_end(png_ptr, info_ptr);  
    png_free(png_ptr, palette);  
    palette=NULL;  
    png_destroy_write_struct(&png_ptr, &info_ptr);  
    fclose(fp);   
}

GRAPH::GRAPH(double width, double height, int resolution)
{
    image.pwidth=round(width*resolution); 
    image.pheight=round(height*resolution);
    image.area_start_px=round(image.border_width_factor*image.pwidth);
    image.area_start_py=round(image.border_height_factor*image.pheight);
    image.area_pwidth=round((1.0-image.border_width_factor*2.0)*image.pwidth); 
    image.area_pheight=round((1.0-image.border_height_factor*2.0)*image.pheight);
    image.area_end_px=image.area_start_px+image.area_pwidth;
    image.area_end_py=image.area_start_py+image.area_pheight;
    image.pixels=new double*[image.pheight];
    for(int i=0;i<image.pheight;i++){
        image.pixels[i]=new double[image.pwidth];
    }
    for(int i=0;i<image.pheight;i++){
        for(int j=0;j<image.pwidth;j++){
            image.pixels[i][j]=image.black;
        }
    }
}

void GRAPH::get_pixel(int &i, int &j, double x, double y)
{
    i=image.area_start_px+round(image.area_pwidth*(x-xaxis.limit[0])/xaxis.spacing);
    j=image.area_start_py+round(image.area_pheight*(y-yaxis.limit[0])/yaxis.spacing);
}

void GRAPH::draw_marker(POINT *point)
{
    double sizex=point->psize/double(image.area_pwidth)*xaxis.spacing;
    double sizey=point->psize/double(image.area_pheight)*yaxis.spacing;
    double start_x=point->x-sizex/2.0, start_y=point->y-sizey/2.0;
    double end_x=point->x+sizex/2.0, end_y=point->y+sizey/2.0;
    int start_i, start_j, end_i, end_j;
    get_pixel(start_i, start_j, start_x, start_y);
    get_pixel(end_i, end_j, end_x, end_y);
    int dii=end_i-start_i, djj=end_j-start_j;
    int istep=dii>0?1:-1, jstep=djj>0?1:-1;
    switch(point->style)
    {
    case 's':
        for(int i=start_i-1;i<end_i;i+=istep){
            for(int j=start_j-1;j<end_j;j+=jstep){
                if(((!point->is_in_area_checked)&&i>=0&&i<image.pwidth&&j>=0&&j<image.pheight)||(i>=image.area_start_px-1&&i<image.area_end_px&&j>=image.area_start_py&&j<image.area_end_py)){
                    image.pixels[j][i]=point->black;
                }
            }
        }
        break;
    case 'c':
        {
            int num=0, sum_i=0, sum_j=0;
            for(int i=start_i-1;i<end_i;i+=istep){
                for(int j=start_j-1;j<end_j;j+=jstep){
                    num++;
                    sum_i+=i; sum_j+=j;
                }
            }
            double center_i=double(sum_i)/double(num);
            double center_j=double(sum_j)/double(num);
            double r2=point->psize*point->psize/4.0;
            for(int i=start_i-1;i<end_i;i+=istep){
                for(int j=start_j-1;j<end_j;j+=jstep){
                    if((i-center_i)*(i-center_i)+(j-center_j)*(j-center_j)<r2){
                        if(((!point->is_in_area_checked)&&i>=0&&i<image.pwidth&&j>=0&&j<image.pheight)||(i>=image.area_start_px-1&&i<image.area_end_px&&j>=image.area_start_py&&j<image.area_end_py)){
                            image.pixels[j][i]=point->black;
                        }
                    }
                }
            }
        }
        break;
    default:
        printf("[ERROR] Unrecognized marker style %c", point->style);
        exit(EXIT_FAILURE);
    }
}

void GRAPH::draw_line(POINT *start_point, POINT *end_point)
{
    if(fabs(start_point->x-end_point->x)<fabs(start_point->y-end_point->y)){
        double dy=yaxis.spacing/double(image.area_pheight);
        POINT  p={start_point->x, start_point->y, start_point->black, start_point->style, start_point->psize, start_point->is_in_area_checked};
        double k=(start_point->x-end_point->x)/(start_point->y-end_point->y);
        double b=end_point->x-k*end_point->y;
        if(start_point->y<end_point->y){
            for(p.y=start_point->y;p.y<end_point->y;p.y+=dy){
                p.x=k*p.y+b;
                draw_marker(&p);
            }
        }else{
            for(p.y=start_point->y;p.y>end_point->y;p.y-=dy){
                p.x=k*p.y+b;
                draw_marker(&p);
            } 
        }
        draw_marker(end_point);
    }else{
        double dx=xaxis.spacing/double(image.area_pwidth);
        POINT  p={start_point->x, start_point->y, start_point->black, start_point->style, start_point->psize, start_point->is_in_area_checked};
        double k=(start_point->y-end_point->y)/(start_point->x-end_point->x);
        double b=end_point->y-k*end_point->x;
        if(start_point->x<end_point->x){
            for(p.x=start_point->x;p.x<end_point->x;p.x+=dx){
                p.y=k*p.x+b;
                draw_marker(&p);
            }
        }else{
            for(p.x=start_point->x;p.x>end_point->x;p.x-=dx){
                p.y=k*p.x+b;
                draw_marker(&p);
            } 
        }
        draw_marker(end_point);
    }
}

void GRAPH::auto_xlim()
{
    double min_x, max_x;
    min_x=points[0].x;
    max_x=points[0].x;
    for(int i=1;i<nump;i++){
        if(min_x>points[i].x) min_x=points[i].x;
        if(max_x<points[i].x) max_x=points[i].x;
    }
    if(image.pwidth<image.pheight){
        xaxis.n_major_tick=round(double(image.pwidth)/double(image.pheight)*(yaxis.n_major_tick+1));
    }
    xaxis.major_tick_spacing=(max_x-min_x)/(xaxis.n_major_tick-1);
    xaxis.minor_tick_spacing=xaxis.major_tick_spacing/(xaxis.n_minor_tick+1);
    xaxis.limit[0]=min_x-xaxis.major_tick_spacing;
    xaxis.limit[1]=max_x+xaxis.major_tick_spacing;
    xaxis.spacing=xaxis.limit[1]-xaxis.limit[0];
}

void GRAPH::auto_ylim()
{
    double min_y, max_y;
    min_y=points[0].y;
    max_y=points[0].y;
    for(int i=1;i<nump;i++){
        if(min_y>points[i].y) min_y=points[i].y;
        if(max_y<points[i].y) max_y=points[i].y;
    }
    if(image.pwidth>image.pheight){
        yaxis.n_major_tick=round(double(image.pheight)/double(image.pwidth)*(xaxis.n_major_tick+1));
    }
    yaxis.major_tick_spacing=(max_y-min_y)/(yaxis.n_major_tick-1);
    yaxis.minor_tick_spacing=yaxis.major_tick_spacing/(yaxis.n_minor_tick+1);
    yaxis.limit[0]=min_y-yaxis.major_tick_spacing;
    yaxis.limit[1]=max_y+yaxis.major_tick_spacing;
    yaxis.spacing=yaxis.limit[1]-yaxis.limit[0];
}

void GRAPH::draw_xaxes()
{
    POINT p1, p2;
    p1.black=p2.black=xaxis.black;
    p1.style=p2.style='s';
    p1.psize=p2.psize=xaxis.line_pwidth;
    p1.is_in_area_checked=p2.is_in_area_checked=false;
    for(int i=0;i<2;i++){
        p1.x=xaxis.limit[0]; p1.y=yaxis.limit[i];
        p2.x=xaxis.limit[1]; p2.y=yaxis.limit[i];
        draw_line(&p1, &p2);
    }
    double x_major_tick_length=xaxis.major_tick_psize/double(image.area_pheight)*yaxis.spacing;
    double x_minor_tick_length=xaxis.minor_tick_psize/double(image.area_pheight)*yaxis.spacing;
    int    dr=xaxis.is_tick_in?1:-1;
    for(int j=0;j<xaxis.n_minor_tick;j++){
            p1.x=xaxis.limit[0]+(j+1)*xaxis.minor_tick_spacing; p1.y=yaxis.limit[0];
            p2.x=p1.x; p2.y=p1.y+x_minor_tick_length*dr;
            draw_line(&p1, &p2); 
    }
    for(int i=0;i<xaxis.n_major_tick;i++){
        p1.x=xaxis.limit[0]+(i+1)*xaxis.major_tick_spacing; p1.y=yaxis.limit[0];
        p2.x=p1.x; p2.y=p1.y+x_major_tick_length*dr;
        draw_line(&p1, &p2);
        for(int j=0;j<yaxis.n_minor_tick;j++){
            p1.x=p1.x+(j+1)*xaxis.minor_tick_spacing; p1.y=yaxis.limit[0];
            p2.x=p1.x; p2.y=p1.y+x_minor_tick_length*dr;
            draw_line(&p1, &p2);
        }
    }
}

void GRAPH::draw_yaxes()
{
    bool is_in_area_checked=false;
    POINT p1, p2;
    p1.black=p2.black=yaxis.black;
    p1.style=p2.style='s';
    p1.psize=p2.psize=yaxis.line_pwidth;
    p1.is_in_area_checked=p2.is_in_area_checked=false;
    for(int i=0;i<2;i++){
        p1.x=xaxis.limit[i]; p1.y=yaxis.limit[0];
        p2.x=xaxis.limit[i]; p2.y=yaxis.limit[1];
        draw_line(&p1, &p2);
    }
    double y_major_tick_length=yaxis.major_tick_psize/double(image.area_pwidth)*xaxis.spacing;
    double y_minor_tick_length=yaxis.minor_tick_psize/double(image.area_pwidth)*xaxis.spacing;
    int    dr=yaxis.is_tick_in?1:-1;
    for(int j=0;j<yaxis.n_minor_tick;j++){
        p1.x=xaxis.limit[0]; p1.y=yaxis.limit[0]+(j+1)*yaxis.minor_tick_spacing;
        p2.x=p1.x+y_minor_tick_length*dr; p2.y=p1.y;
        draw_line(&p1, &p2); 
    }
    for(int i=0;i<yaxis.n_major_tick;i++){
        p1.x=xaxis.limit[0]; p1.y=yaxis.limit[0]+(i+1)*yaxis.major_tick_spacing;
        p2.x=p1.x+y_major_tick_length*dr; p2.y=p1.y;
        draw_line(&p1, &p2); 
        for(int j=0;j<yaxis.n_minor_tick;j++){
            p1.x=xaxis.limit[0]; p1.y=p1.y+(j+1)*yaxis.minor_tick_spacing;
            p2.x=p1.x+y_minor_tick_length*dr; p2.y=p1.y;
            draw_line(&p1, &p2);
        }
    }
}

void GRAPH::draw(const char *png_path)
{
    int npixel=image.pwidth*image.pheight;
    double *wpixels=new double[npixel];
    for(int i=0;i<image.pheight;i++){
        for(int j=0;j<image.pwidth;j++){
            wpixels[i*image.pwidth+j]=image.pixels[i][j];
        }
    }
    unsigned char *pixels=new unsigned char[3*npixel];
    for(int i=0;i<npixel;i++){
        pixels[i*3]=pixels[i*3+1]=pixels[i*3+2]=round((1.0-wpixels[i])*255.0);
    }
    image_pixels(png_path, pixels, image.pwidth, image.pheight);
}

void GRAPH::scatter(const char *png_path, double *x, double *y, double *value, int num)
{
    nump=num;
    points=new POINT[num];
    double max_value=0.0;
    for(int i=0;i<num;i++){
        points[i].x=x[i];
        points[i].y=y[i];
        if(max_value<value[i]){
            max_value=value[i];
        }
    }
    auto_xlim();
    auto_ylim();
    draw_xaxes();
    draw_yaxes();
    for(int i=0;i<num;i++){
        points[i].black=value[i]/max_value;
        draw_marker(&points[i]);
    }
    draw(png_path);
}

void GRAPH::line(const char *png_path, double *x, double *y, int num)
{
    nump=num;
    points=new POINT[num];
    for(int i=0;i<num;i++){
        points[i].x=x[i];
        points[i].y=y[i];
    }
    auto_xlim();
    auto_ylim();
    draw_xaxes();
    draw_yaxes();
    for(int i=1;i<num;i++){
        draw_line(&points[i-1], &points[i]);
    }
    draw(png_path);
}