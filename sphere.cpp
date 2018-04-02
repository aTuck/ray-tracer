#include <stdio.h>
#include "sphere.h"
#include <stdlib.h>
#include <math.h>

extern Spheres *scene;

/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection, 
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/
float intersect_sphere(Point u, Vector d, Spheres *sph, Point *hit) {
  // u (x0,y0,z0) is eye 
  // d (dx,dy,dz) is direction 
  // t (t1,t2) is parameter
  
  // initial point coords
  float x0 = u.x;
  float y0 = u.y;
  float z0 = u.z;

  // directional Vector coords
  float dx = d.x;
  float dy = d.y;
  float dz = d.z;

  // sphere values
  float xc =   0;
  float yc =   0;
  float zc =   0;
  float r  =   0;

  // sphere values
  xc = sph->center.x;
  yc = sph->center.y;
  zc = sph->center.z;
  r  =   sph->radius;

  // formula values
  float a = (dx*dx + dy*dy + dz*dz);
  float b = 2 * (dx*(x0 - xc) + dy*(y0 - yc) + dz*(z0 - zc));
  float c = (x0 - xc)*(x0 - xc) + (y0 - yc)*(y0 - yc) + (z0 - zc)*(z0 - zc) - r*r;

  float discrim = (b*b) - (4*a*c);
  if (discrim < 0){
    return -1;
  }
  else{
    float t1 = ( -b + sqrt(discrim) ) / (2 * a);
    float t2 = ( -b - sqrt(discrim) ) / (2 * a);

    if (t2 < t1){
      Vector scaled_vec = vec_scale(d, t2);
      hit->x = x0 + scaled_vec.x;
      hit->y = y0 + scaled_vec.y;
      hit->z = z0 + scaled_vec.z;

      return t2;
    }
  }

  return -1;
}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres* intersect_scene(Point u, Vector d, Spheres* sph, Point* hit, int n) {
  Spheres * temp = sph;
  Spheres * closest_sph = 0;
  Point temp_hit;
  
  float closest_param = -2;
  float param;
  float tol = 0.0001;

  normalize(&d);

  do {
    param = intersect_sphere(u, d, temp, &temp_hit);
    if ((param < closest_param && param != -1) || (closest_param == -2 && param != -1)){
      if (param > tol){
        closest_param = param;
        closest_sph = temp;
        hit->x = temp_hit.x;
        hit->y = temp_hit.y;
        hit->z = temp_hit.z;
      }
    }
  } while (temp = temp->next);

  if (closest_sph != 0){
    return closest_sph;
  }
  else {
    return NULL;
  }
}


/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine, 
		    float refl, int sindex) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;

  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}
