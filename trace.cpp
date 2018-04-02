#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"

//
// Global variables
//
extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];  

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;

// light 1 position and color
extern Point light1;
extern float light1_ambient[3];
extern float light1_diffuse[3];
extern float light1_specular[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern int shadow_on;
extern int step_max;

/////////////////////////////////////////////////////////////////////

/*********************************************************************
 * Phong illumination - you need to implement this!
 *********************************************************************/
RGB_float phong(Point q, Vector v, Vector n, Vector l, Vector r, Spheres *sph) {
	RGB_float color;
  float d = vec_len(l);  
  
  normalize(&v);  
  normalize(&n);
  normalize(&l);
  normalize(&r);

  float igar = global_ambient[0];
  float igag = global_ambient[1];
  float igab = global_ambient[2];
  float kga = sph->reflectance;
  
  float iar  = light1_ambient[0];
  float iag  = light1_ambient[1];
  float iab  = light1_ambient[2];
  float kar  = sph->mat_ambient[0];
  float kag  = sph->mat_ambient[1];
  float kab  = sph->mat_ambient[2];
  float idr  = light1_diffuse[0];
  float idg  = light1_diffuse[1];
  float idb  = light1_diffuse[2];
  float kdr  = sph->mat_diffuse[0];
  float kdg  = sph->mat_diffuse[1];
  float kdb  = sph->mat_diffuse[2];
  float isr  = light1_specular[0];
  float isg  = light1_specular[1];
  float isb  = light1_specular[2];
  float ksr  = sph->mat_specular[0];
  float ksg  = sph->mat_specular[1];
  float ksb  = sph->mat_specular[2];
  float N = sph->mat_shineness;
  
  float decay = 1 / (decay_a + decay_b*d + decay_c*(d*d));

  float diffuse_specular_r = idr*kdr*vec_dot(n, l) + isr*ksr*pow(fmax(vec_dot(r, v), 0), N); //rv^N or 0
  float diffuse_specular_g = idg*kdg*vec_dot(n, l) + isg*ksg*pow(fmax(vec_dot(r, v), 0), N); //rv^N or 0
  float diffuse_specular_b = idb*kdb*vec_dot(n, l) + isb*ksb*pow(fmax(vec_dot(r, v), 0), N); //rv^N or 0
  
  Spheres * shadow_sph = 0;
  if (shadow_on){
    Point shadow_hit;
    shadow_sph = intersect_scene(q, l, scene, &shadow_hit, 0);
  }
  if (shadow_sph != 0){
    color.r = igar*kga + iar*kar;
    color.g = igag*kga + iag*kag;
    color.b = igab*kga + iab*kab;
  }
  else {
    color.r = igar*kga + iar*kar + decay*diffuse_specular_r;
    color.g = igag*kga + iag*kag + decay*diffuse_specular_g;
    color.b = igab*kga + iab*kab + decay*diffuse_specular_b;
  }

	return color;
}

Vector calculate_r(Vector v, Vector n){
    float v_dot_n = vec_dot(v, n);
    Vector n_scaled = vec_scale(n, v_dot_n);
    Vector n_scaled_2 = vec_scale(n_scaled, 2.0);
    Vector return_vec = vec_minus(v, n_scaled_2);

    return return_vec;
}

/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace(Point initial_pos, Vector ray, int step) {
	RGB_float color = background_clr;
  RGB_float next_color = background_clr;

  Vector next_ray;
  Point hit;

  Spheres * sph = intersect_scene(initial_pos, ray, scene, &hit, 0);

  if (sph == NULL){ return color; }

  Vector v = get_vec(eye_pos, hit);
  normalize(&v);
  Vector n = sphere_normal(hit, sph);
  normalize(&n);
  Vector l = get_vec(hit,  light1);
  normalize(&l);
  Vector r = calculate_r(l, n);
  normalize(&r);
  next_ray = calculate_r(v, n);

  color = phong(hit, v, n, l, r, sph);

  if (step < step_max){
    next_color = recursive_ray_trace(hit, next_ray, ++step);
    
    // next_color didn't hit a sphere, set to nothing
    // if (next_color.r == background_clr.r){
    //   next_color.r = 0;
    //   next_color.g = 0;
    //   next_color.b = 0;
    // }

    color = clr_add(color, clr_scale(next_color, sph->reflectance));
  }

  return color;
}

/*********************************************************************
 * This function traverses all the pixels and cast rays. It calls the
 * recursive ray tracer and assign return color to frame
 *
 * You should not need to change it except for the call to the recursive
 * ray tracer. Feel free to change other parts of the function however,
 * if you must.
 *********************************************************************/
void ray_trace() {
  int i, j;
  float x_grid_size = image_width / float(win_width);
  float y_grid_size = image_height / float(win_height);
  float x_start = -0.5 * image_width;
  float y_start = -0.5 * image_height;
  RGB_float ret_color;
  Point cur_pixel_pos;
  Vector ray;
  Point * hit = (Point*)malloc(sizeof(Point*));
  Spheres * sph;

  // ray is cast through center of pixel
  cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
  cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
  cur_pixel_pos.z = image_plane;

  // for each pixel
  for (i=0; i<win_height; i++) {
    for (j=0; j<win_width; j++) {
      // initial ray from eye through pixel
      ray = get_vec(eye_pos, cur_pixel_pos);

      ret_color = recursive_ray_trace(eye_pos, ray, 0);

      // Parallel rays can be cast instead using below
      //
      // ray.x = ray.y = 0;
      // ray.z = -1.0;
      // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);

      frame[i][j][0] = GLfloat(ret_color.r);
      frame[i][j][1] = GLfloat(ret_color.g);
      frame[i][j][2] = GLfloat(ret_color.b);

      cur_pixel_pos.x += x_grid_size;
    }

    cur_pixel_pos.y += y_grid_size;
    cur_pixel_pos.x = x_start;
  }

  free(hit);
}
