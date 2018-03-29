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
float intersect_sphere(Point o, Vector u, Spheres *sph, Point *hit) {
	return 0.0;
}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres* intersect_scene(Point u, Vector d, Spheres* sph, Point* hit, int n) {
  // flag to determine return value
  bool didIntersect = 0;
  Spheres * temp;
  temp = sph;

  // initial point coords
  float x0 = u.x;
  float y0 = u.y;
  float z0 = u.z;

  // directional Vector coords
  float dx = d.x;
  float dy = d.y;
  float dz = d.z;

  // parametric values
  float t1 =   0;
  float t2 =   0;

  // intermediate values
  float a =    0;
  float b =    0;
  float c =    0;

  // sphere values
  float xc =   0;
  float yc =   0;
  float zc =   0;
  float r  =   0;

  // u (x0,y0,z0) is eye 
  // d (dx,dy,dz) is direction 
  // t (t1,t2) is parameter
  //
  // x = x0 + dx * t
  // y = y0 + dy * t
  // z = z0 + dz * t
  //
  // iterating through spheres decrements index
  while (temp->index > 0){
    // printf("sph is 0x%d\n", &sph);
    didIntersect = 0;

    xc = temp->center.x;
    yc = temp->center.y;
    zc = temp->center.z;
    r  =   temp->radius;

    a = (dx*dx + dy*dy + dz*dz);
    b = 2 * (dx*(x0 - xc) + dy*(y0 - yc) + dz*(z0 - zc));
    c = (x0 - xc)*(x0 - xc) + (y0 - yc)*(y0 - yc) + (z0 - zc)*(z0 - zc) - r*r;

    float discrim = (b*b) - (4*a*c);
    if (discrim > 0){
      // printf("[2] intersects discrim: %f\n", discrim);
      didIntersect = 1;

      t1 = (-b + sqrt(discrim))/2*a;
      t2 = (-b - sqrt(discrim))/2*a;

      break;
    }
    else if (abs(discrim) <= 0.00001){
      // printf("[1] intersect discrim: %f\n", discrim);
      didIntersect = 1;

      t1 = -b/2*a;
      t2 =      0;

      break;
    }
    else{
      //printf("[0] no intersect\n");
      didIntersect = 0;
      
      t1 = 0;
      t2 = 0;
            
      // printf("didnt intersect with %i, go next\n", temp->index);
    }
    if (temp->next != NULL){
      temp = temp->next;
    }
    else{
      hit->x = 0;
      hit->y = 0;
      hit->z = 0;

      temp = NULL;
      return temp;
    }
  }

  float ix1 = x0 + dx * t1;
  float iy1 = y0 + dy * t1;
  float iz1 = z0 + dz * t1;

  hit->x = ix1;
  hit->y = iy1;
  hit->z = iz1;
  return temp;
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
