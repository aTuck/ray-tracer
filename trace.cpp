#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"
#include "trace.h"

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

// command line arguments
extern int shadow_on;
extern int reflection_on;
extern int refraction_on;
extern int stochastic_on;
extern int chessboard_on;
extern int supersample_on;
extern int step_max;

// function declaration so stochastic/recursive can call each other
RGB_float recursive_ray_trace(Point initial_pos, Vector ray, int step);

/////////////////////////////////////////////////////////////////////



Vector reflect(Vector v, Vector n){
    float v_dot_n = vec_dot(v, n);
    Vector n_scaled = vec_scale(n, v_dot_n);
    Vector n_scaled_2 = vec_scale(n_scaled, 2.0);
    Vector return_vec = vec_minus(v, n_scaled_2);

    return return_vec;
}

Vector refract(Vector v, Vector n, Spheres * sph){
  // This function calculates refracted ray based on index of refraction
  //
  // Based off of this web page
  // https://www.cs.unc.edu/~rademach/xroads-RT/RTarticle.html

  float c1 = -vec_dot(n, v);
  Vector R1 = vec_plus(v, vec_scale(vec_scale(n, 2), c1));

  float n1 = 1;
  float n2 = sph->ior;
  float N = n1 / n2;

  float c2 = sqrt( 1 - N*N * (1 - c1*c1));
  float norm_scale = ((N * c1) - c2);

  Vector scaled_v = vec_scale(v, N);
  Vector scaled_n = vec_scale(n, norm_scale);

  Vector refracted_ray = vec_plus(scaled_v, scaled_n);
  normalize(&refracted_ray);

  return refracted_ray;
}

Vector generate_random_ray(){
  Vector random_ray;

  random_ray.x = (float)rand() / (float)(RAND_MAX/2) - 1;
  random_ray.y = (float)rand() / (float)(RAND_MAX/2) - 1;
  random_ray.z = (float)rand() / (float)(RAND_MAX/2) - 1;

  return random_ray;
}

void vec_print(Vector v){
  printf("%f, %f, %f\n", v.x, v.y, v.z);
}

void set_to_black(RGB_float * c){
  c->r = 0;
  c->g = 0;
  c->b = 0;
}

RGB_float check_for_negative_rgb(RGB_float c){
  if (c.r < 0){c.r = 0;}
  if (c.g < 0){c.g = 0;}
  if (c.b < 0){c.b = 0;}

  return c;
}

/***********************
 * Phong illumination 
 ***********************/
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

RGB_float stochastic_ray_trace(int step, Point hit, Spheres * sph){
  Vector random_ray;
  RGB_float temp_color = background_clr;
  RGB_float color = background_clr;
  
  for (int i = 0; i < MAX_STOCHASTIC; i++){
    do {
      random_ray = generate_random_ray();
    } while (vec_dot(random_ray, );
    temp_color = recursive_ray_trace(hit, random_ray, ++step);

    // combine colors
    color =  clr_add(color, temp_color);
  }

  // modify color by diffuse
  color.r = color.r * sph->mat_diffuse[0];
  color.g = color.g * sph->mat_diffuse[1];
  color.b = color.b * sph->mat_diffuse[2];

  // take average color
  color.r = color.r/MAX_STOCHASTIC;
  color.g = color.g/MAX_STOCHASTIC;
  color.b = color.b/MAX_STOCHASTIC;

  return color;
}

/**********************************
 * This is the recursive ray tracer
 ***********************************/
RGB_float recursive_ray_trace(Point initial_pos, Vector ray, int step) {
	RGB_float ret_color = background_clr;
  RGB_float phong_color = background_clr;
  RGB_float stochastic_color = background_clr;
  RGB_float reflected_color = background_clr;
  RGB_float refracted_color = background_clr;
  Vector reflected_ray;
  Vector refracted_ray;
  Vector random_ray;
  Point hit;

  // check for intersect
  normalize(&ray);
  Spheres * sph = intersect_scene(initial_pos, ray, scene, &hit, 0);

  // no intersect, return background color
  if (sph == NULL){ return background_clr; }
  
  // intersect exists, do phong shading
  Vector v = get_vec(eye_pos, hit);
  normalize(&v);
  Vector n;
  if (sph->shape == 0){
      n = sphere_normal(hit, sph);
  }
  else{
      n = {0.0, 1.0, 0.0};
  }
  normalize(&n);
  Vector l = get_vec(hit,  light1);
  normalize(&l);
  Vector r = reflect(l, n);
  normalize(&r);
  phong_color = phong(hit, v, n, l, r, sph);
  
  // initialize return color to phong shaded color
  ret_color = phong_color;

  if (step < step_max){
    int avg = 1;

    if (stochastic_on){
      stochastic_color = stochastic_ray_trace(step, hit, sph);

      ret_color = clr_add(ret_color, clr_scale(stochastic_color, sph->reflectance));
      avg++;
    }
    if (reflection_on){
      reflected_ray = reflect(v, n);
      reflected_color = recursive_ray_trace(hit, reflected_ray, ++step);

      ret_color = clr_add(ret_color, clr_scale(reflected_color, sph->reflectance));
      avg++;
    }
    if (refraction_on){
      // get rid of phong shading
      set_to_black(&ret_color);

      refracted_ray = refract(ray, n, sph);
      refracted_color = recursive_ray_trace(hit, refracted_ray, ++step);

      ret_color = clr_add(ret_color, refracted_color);
      avg++;
    }

    // take average color
    ret_color.r = ret_color.r/avg;
    ret_color.g = ret_color.g/avg;
    ret_color.b = ret_color.b/avg;

    // give back color to recursing functions
    return ret_color;
  }

  return ret_color;
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

      // do four corner rays and average
      if (supersample_on){
        // extra ray cast positions
        Point temp_pixel_pos;
        temp_pixel_pos.z = image_plane;
        
        // colors
        RGB_float topleft_clr;
        RGB_float topright_clr;
        RGB_float botleft_clr;
        RGB_float botright_clr;

        // distance from initial ray modifier (higher value is closer to center)
        int dist_from_center = 4;
        
        // top left
        temp_pixel_pos.x = cur_pixel_pos.x - (x_grid_size/dist_from_center);
        temp_pixel_pos.y = cur_pixel_pos.y - (y_grid_size/dist_from_center);
        ray = get_vec(eye_pos, temp_pixel_pos);
        topleft_clr = recursive_ray_trace(eye_pos, ray, 0);

        // top right
        temp_pixel_pos.x = cur_pixel_pos.x + (x_grid_size/dist_from_center);
        temp_pixel_pos.y = cur_pixel_pos.y - (y_grid_size/dist_from_center);
        ray = get_vec(eye_pos, temp_pixel_pos);
        topright_clr = recursive_ray_trace(eye_pos, ray, 0);

        // bottom left
        temp_pixel_pos.x = cur_pixel_pos.x - (x_grid_size/dist_from_center);
        temp_pixel_pos.y = cur_pixel_pos.y + (y_grid_size/dist_from_center);
        ray = get_vec(eye_pos, temp_pixel_pos);
        botleft_clr = recursive_ray_trace(eye_pos, ray, 0);
        
        // bottom right
        temp_pixel_pos.x = cur_pixel_pos.x + (x_grid_size/dist_from_center);
        temp_pixel_pos.y = cur_pixel_pos.y + (y_grid_size/dist_from_center);
        ray = get_vec(eye_pos, temp_pixel_pos);
        botright_clr = recursive_ray_trace(eye_pos, ray, 0);

        // combine
        ret_color = clr_add(ret_color, topleft_clr );
        ret_color = clr_add(ret_color, topright_clr);
        ret_color = clr_add(ret_color, botleft_clr );
        ret_color = clr_add(ret_color, botright_clr);

        // average
        ret_color.r = ret_color.r/5;
        ret_color.g = ret_color.g/5;
        ret_color.b = ret_color.b/5;
      }

      check_for_negative_rgb(ret_color);
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

// Parallel rays can be cast instead using below
//
// ray.x = ray.y = 0;
// ray.z = -1.0;
// ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);