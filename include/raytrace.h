#ifndef _RAYTRACE_H_
#define _RAYTRACE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "json_parser.h"
#include "vector_math.h"
#include "pplib.h"
#include "illumination.h"

#define SHININESS 15
#define MAX_REC_LEVEL 10

// Struct to represent both position and direction
typedef struct Ray {
        V3 origin;
        V3 direction;
} Ray;

// Default background color
V3 background_color = {0, 0, 0};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Searches an array of objects and returns the index of the camera.
  Returns -1 if not camera is found.

  Parameters
  ==========
  OBJECT objects
         An OBJECT array containing different objects
  Returns
  =======
  size_t index
         The index of the camera object
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
size_t get_camera(OBJECT *objects) {
        int i = 0;
        while (i < MAX_OBJECTS && objects[i].type != 0) {
                if (objects[i].type == CAM) {
                        return i;
                }
                i++;
        }
        // no camera found in data
        return -1;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Sets the color value for a pixel in an image

  Parameters
  ==========
  double color
         The color used to set
  size_t row
         The pixel row
  size_t col
         The pixel column
  Image image
        The image the pixel comes from

  Returns
  =======
  void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
static inline void shade_pixel(double *color, size_t row, size_t col, Image *image) {
        image->data[row * image->width*4 + col*4] = (uint8_t)(MAX_COLOR_VAL * check_value(color[0]));
        image->data[row * image->width*4 + col*4+1] = (uint8_t)(MAX_COLOR_VAL* check_value(color[1]));
        image->data[row * image->width*4 + col*4+2]= (uint8_t)(MAX_COLOR_VAL* check_value(color[2]));
        image->data[row * image->width*4 + col*4+3]= 255;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Computes if ray intersects with plane, returns the distance of
  intersection if so.

  Parameters
  ==========
  Ray ray
      The ray used to compute the intersection with an object
  double pos
         The position of the plane
  V3 norm
     Plane surface normal

  Returns
  =======
  double distance
         The distance of intersection
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double plane_intersection(Ray *ray, double *Pos, V3 Norm){
        double alph,delta;
        normalize(Norm);

        alph = v3_dot(Norm, ray->direction);

        // the plane is parallel to the ray
        if (fabs(alph) <0.0001) {
                return -1;
        }

        V3 incident_vector;
        v3_subtract(Pos, ray->origin, incident_vector);
        delta = v3_dot(incident_vector, Norm);


        double t = delta/alph; // whcih means we check thea1 and thea2

        if (t<0.0) { // reflection, no intersection
                return -1;
        }
        return t; // return something, but not t , need to figure out it

}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Computes if ray intersects with sphere, returns the distance of
  intersection if so.

  Parameters
  ==========
  Ray ray
      The ray used to compute the intersection with an object
  V3 center
     The center of the sphere
  double r
         The radius of sphere

  Returns
  =======
  double distance
         The distance of intersection
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double sphere_intersection(Ray *ray, V3 C, double r) {
        double a, b;
        V3 vector_diff;
        v3_subtract(ray->origin, C, vector_diff);

        a = 2 * (ray->direction[0]*vector_diff[0] + ray->direction[1]*vector_diff[1] + ray->direction[2]*vector_diff[2]);
        b = sqr(vector_diff[0]) + sqr(vector_diff[1]) + sqr(vector_diff[2]) - sqr(r);

        // check that discriminant is <, =, or > 0
        double disc = sqr(a) - 4*b;
        double t; // solutions
        if (disc < 0) {
                return -1; // no solution
        }
        disc = sqrt(disc);
        t = (-a - disc) / 2.0;
        if (t < 0.0)
                t = (-a + disc) / 2.0;

        if (t < 0.0)
                return -1;
        return t;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Gets the location of intersection.

  Parameters
  ==========
  V3 intersection
     new vector with intersection position vector
  Ray ray
      Position used to compute intersection
  double t
         distance of intersection

  Returns
  =======
  void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void get_intersection(V3 intersection, Ray *ray, double t){
        intersection[0] = ray->origin[0] + ray->direction[0] * t;
        intersection[1] = ray->origin[1] + ray->direction[1] * t;
        intersection[2] = ray->origin[2] + ray->direction[2] * t;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Computes the object that intersects closest to some ray

  Parameters
  ==========
  Ray ray
      The ray to use for computing the intersection with an object
  size_t self_index
         index of object to find the solution for
  double max_distance
         Farthest intersection distance acceptable
  size_t* ret_index
         The variable to store the index of the best solution object
  double* ret_best_t
          The variable to store the distance of intersection with the best solution object

  Returns
  =======
  void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void get_best_solution(Ray *ray, int self_index, double max_distance, int *ret_index, double *ret_best_t) {
        int best_o = -1;
        double best_t = INFINITY;
        for (int i=0; objects[i].type != 0; i++) {
                // ignore self
                if (self_index == i) continue;

                // we need to run intersection test on each object
                double t = 0;
                switch(objects[i].type) {
                case 0:
                        printf("no object found\n");
                        break;
                case CAM:
                        break;
                case SPH:
                        t = sphere_intersection(ray, objects[i].sphere.position,
                                                objects[i].sphere.radius);
                        break;
                case PLAN:
                        t = plane_intersection(ray, objects[i].plane.position,
                                               objects[i].plane.normal);
                        break;
                default:
                        printf("No intersection\n");
                        exit(1);
                }
                if (max_distance != INFINITY && t > max_distance)
                        continue;
                if (t > 0 && t < best_t) {
                        best_t = t;
                        best_o = i;
                }
        }
        (*ret_index) = best_o;
        (*ret_best_t) = best_t;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Get the surface normal of object at intersection location

  Parameters
  ==========
  size_t object_index
         The index of the object in the object array
  V3 position
     The position to compute the normal for
  V3 normal
     The variable to store the final normal

  Returns
  =======
  void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void get_normal(int object_index, V3 position, V3 normal){
        if (objects[object_index].type == PLAN) {
                v3_copy(objects[object_index].plane.normal, normal);
        } else if (objects[object_index].type == SPH) {
                v3_subtract(position, objects[object_index].sphere.position, normal);
        } else {
                fprintf(stderr, "Error: normal_vector: Can't get normal vector for this project\n");
        }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Accessor for reflectivity

  Parameters
  ==========
  size_t object_index
         The index of the object in the object array

  Returns
  =======
  double reflectivity
         The reflectivity of the object
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double get_reflectivity(size_t object_index){
        if (objects[object_index].type == PLAN) {
                return objects[object_index].plane.reflect;
        }
        else if (objects[object_index].type == SPH) {
                return objects[object_index].sphere.reflect;
        }
        else {
                fprintf(stderr, "Error: get_reflectivity: Can't find reflect property for this project\n");
                return -1;
        }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  Accessor for refractivity

  Parameters
  ==========
  size_t object_index
         The index of the object in the object array

  Returns
  =======
  double refractivity
         The refractivity of the object
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double get_refractivity(int object_index) {
        if (objects[object_index].type == PLAN) {
                return objects[object_index].plane.refract;
        }
        else if (objects[object_index].type == SPH) {
                return objects[object_index].sphere.refract;
        }
        else {
                fprintf(stderr, "Error: get_reflectivity: Can't find refract property for this project\n");
                return -1;
        }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 Accessor for the index of refraction

 Parameters
 ==========
 size_t object_index
        The index of the object in the object array

 Returns
 =======
 double ior
        The index of refraction of the object
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double get_ior(int object_index){
        double ior;
        if (objects[object_index].type == PLAN) {
                ior = objects[object_index].plane.ior;
        }
        else if (objects[object_index].type == SPH) {
                ior = objects[object_index].sphere.ior;
        }
        else {
                fprintf(stderr, "Error: get_ior: Can't find ior property for this object\n");
                exit(1);
        }
        if (fabs(ior) < 0.0001)
                return 1;
        else
                return ior;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 Get the reflected vector of object at intersection location
 from some direction

 Parameters
 ==========
 size_t object_index
 The index of the object in the object array
 V3 position
 The position to compute the normal for
 V3 normal
 The variable to store the final reflected vector

 Returns
 =======
 void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void reflection_vector(V3 direction, V3 position, int object_index, V3 reflection) {
        V3 normal;
        get_normal(object_index, position, normal);
        normalize(normal);
        v3_reflect(direction, normal, reflection);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 Get the refraction vector of object at intersection location
 from some direction

 Parameters
 ==========
 V3 direction
    Direction of vector to refract
 V3 position
    position of vector to refract
 int object_index
     The index of the object which the ray is intersection, its
     properties define the final vector
 double out_ior
        The output index of refraction
 V3 refracted_vector
    The final vector, it is stored in this variable

 Returns
 =======
 void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void refraction_vector(V3 direction, V3 position, int object_index, double out_ior, V3 refracted_vector) {
        // initializations and variables setup
        V3 dir, pos;
        v3_copy(direction, dir);
        v3_copy(position, pos);
        normalize(dir);
        normalize(pos);

        double in_ior = get_ior(object_index);

        if (in_ior == out_ior) {
                in_ior = 1;
        }

        // define all the normals
        V3 normal, a, b;

        // find normal vector of current object
        get_normal(object_index, pos, normal);

        normalize(normal);

        // create coordinate frame with a and b, where b is tangent to the object intersection
        v3_cross(normal, dir, a);
        normalize(a);
        v3_cross(a, normal, b);

        // find transmission vector angle and direction
        double sin_theta = v3_dot(dir, b);
        double sin_phi = (out_ior / in_ior) * sin_theta;
        double cos_phi = sqrt(1 - sqr(sin_phi));

        v3_scale(normal, -1*cos_phi, normal);
        v3_scale(b, sin_phi, b);
        v3_add(normal, b, refracted_vector);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 Simple raycast function, computes intersection between a ray an an object,
 then computes the final color based on the object's properties

 Parameters
 ==========
 Ray* ray
      The ray used in ray casting
 int object_index
     The index of the object to compute the color for
 V3 position
     the position of intersection with the object
 LIGHT light
       The light used used to compute lighting
 double max_dist
        The max light distance
 V3 color
    Where the final color is stored

 Returns
 =======
 void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void cast(Ray *ray, int object_index, V3 position, LIGHT *light, double max_dist, V3 color) {
        V3 normal;
        V3 object_diff_color;
        V3 object_spec_color;

        // find normal and color
        if (objects[object_index].type == PLAN) {
                v3_copy(objects[object_index].plane.normal, normal);
                v3_copy(objects[object_index].plane.diff_color, object_diff_color);
                v3_copy(objects[object_index].plane.spec_color, object_spec_color);
        } else if (objects[object_index].type == SPH) {
                // find normal of our current intersection on the sphere
                v3_subtract(ray->origin, objects[object_index].sphere.position, normal);
                // copy the colors into temp variables
                v3_copy(objects[object_index].sphere.diff_color, object_diff_color);
                v3_copy(objects[object_index].sphere.spec_color, object_spec_color);
        } else {
                fprintf(stderr, "Error: shade: Trying to shade unsupported type of object\n");
                exit(1);
        }
        normalize(normal);
        // find light, reflection and camera vectors
        V3 L;
        V3 R;
        V3 V;
        v3_copy(ray->direction, L);
        normalize(L);
        v3_reflect(L, normal, R);
        v3_copy(position, V);
        V3 diffuse_color;
        //Vector_scale(diffuse_color, 0, diffuse_color);
        V3 specular_color;
        //Vector_scale(specular_color, 0, specular_color);

        get_diffuse(normal, L, light->color, object_diff_color, diffuse_color);
        get_specular(SHININESS, L, R, normal, V, object_spec_color, light->color, specular_color);

        // calculate the angular and radial attenuation
        double fang;
        double frad;
        // get the vector from the object to the light
        V3 light_to_object_dir;
        v3_copy(L, light_to_object_dir);
        v3_scale(light_to_object_dir, -1, light_to_object_dir);
        if (light->type == TRACING) {
                fang = 1;
                frad = 1;
        }else{
                fang = calculate_angular_att(light, light_to_object_dir);
                frad = calculate_radial_att(light, max_dist);
        }
        color[0] += frad * fang * (specular_color[0] + diffuse_color[0]);
        color[1] += frad * fang * (specular_color[1] + diffuse_color[1]);
        color[2] += frad * fang * (specular_color[2] + diffuse_color[2]);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 Simple raycast function, computes intersection between a ray an an object,
 then computes the final color based on the object's properties

 Parameters
 ==========
 Ray* ray
      The ray used in ray casting
 int object_index
     The index of the object to compute the color for
 double t
     intersection distance to object
 double current_ior
        The current index of refraction
 int rec_level
     How many times for the ray to bounce off
 V3 color
    Where the final color is stored

 Returns
 =======
 void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void trace(Ray *ray, int object_index, double t, double current_ior, int rec_level, V3 color) {
        if (rec_level > MAX_REC_LEVEL) {
                // if the recursive times more than MAX Times, then we return black color
                v3_scale(color, 0, color);
                return;
        }
        // base case, no intersection object, so set color is black
        if (object_index == -1) {
                v3_scale(color, 0, color);
                return;
        }
        if (ray == NULL) {
                fprintf(stderr, "Error: shade: Ray had is Empty\n");
                exit(1);
        }
        V3 new_origin = {0,0,0};
        V3 new_direction ={0,0,0};

        v3_scale(ray->direction, t, new_origin);
        v3_add(new_origin, ray->origin, new_origin);

        Ray new_ray={
                .origin = {new_origin[0], new_origin[1], new_origin[2]},
                .direction = {new_direction[0], new_direction[1], new_direction[2]}
        };

        V3 reflection ={0,0,0};
        V3 refraction ={0,0,0};
        normalize(ray->direction);
        reflection_vector(ray->direction, new_ray.origin, object_index, reflection);
        refraction_vector(ray->direction, new_ray.origin, object_index, current_ior, refraction);


        int best_reflection_object_index;
        double best_reflection_t;
        int best_refraction_object_index;
        double best_refraction_t;

        Ray ray_reflected = {
                .origin = {new_origin[0], new_origin[1], new_origin[2]},
                .direction = {reflection[0], reflection[1], reflection[2]}
        };

        Ray ray_refracted = {
                .origin = {new_origin[0], new_origin[1], new_origin[2]},
                .direction = {refraction[0], refraction[1], refraction[2]}
        };
        // add some tiny things
        V3 reflected_offset = {0,0,0};
        v3_scale(ray_reflected.direction, 0.01, reflected_offset);
        v3_add(ray_reflected.origin, reflected_offset, ray_reflected.origin);
        normalize(ray_reflected.direction);

        V3 refracted_offset ={0,0,0};
        v3_scale(ray_refracted.direction, 0.01, refracted_offset);
        v3_add(ray_refracted.origin,refracted_offset, ray_refracted.origin);
        normalize(ray_refracted.direction);

        get_best_solution(&ray_reflected,object_index, INFINITY, &best_reflection_object_index, &best_reflection_t);
        if (objects[object_index].type == PLAN) {
                get_best_solution(&ray_refracted, object_index, INFINITY, &best_refraction_object_index, &best_refraction_t);
        }
        else{
                get_best_solution(&ray_refracted, -1, INFINITY, &best_refraction_object_index, &best_refraction_t);
        }
        if (best_reflection_object_index == -1 && best_refraction_object_index == -1) {
                v3_scale(color, 0, color);

        }
        else {
                // we have an intersection, so we use recursive
                V3 reflection_color ={0,0,0};
                V3 refraction_color ={0,0,0};
                double reflect_coefficient = get_reflectivity(object_index);
                double refract_coefficient = get_refractivity(object_index);
                double reflect_ior =1;
                double refract_ior =1;


                LIGHT reflection_light;
                reflection_light.type = TRACING;
                //at very begining, forget to malloc for light direction and color
                reflection_light.direction = malloc(sizeof(V3));
                reflection_light.color = malloc(sizeof(V3));

                LIGHT refraction_light;
                refraction_light.type = TRACING;
                //at very begining, forget to malloc for light direction and color
                refraction_light.direction = malloc(sizeof(V3));
                refraction_light.color = malloc(sizeof(V3));

                if (best_reflection_object_index >= 0) {
                        reflect_ior = get_ior(best_reflection_object_index);
                        trace(&ray_reflected, best_reflection_object_index, best_reflection_t, reflect_ior, rec_level+1, reflection_color);
                        v3_scale(reflection_color, reflect_coefficient, reflection_color);
                        v3_scale(reflection, -1, reflection_light.direction);
                        v3_copy(reflection_color, reflection_light.color);

                        v3_scale(ray_reflected.direction, best_reflection_t, ray_reflected.direction);
                        v3_subtract(ray_reflected.direction, new_ray.origin, new_ray.direction);
                        normalize(new_ray.direction);
                        cast(&new_ray, object_index, ray->direction, &reflection_light, INFINITY, color);
                }
                if (best_refraction_object_index >= 0) {
                        refract_ior = get_ior(best_refraction_object_index);
                        trace(&ray_refracted, best_refraction_object_index, best_refraction_t, refract_ior, rec_level+1, refraction_color);
                        v3_scale(refraction_color, refract_coefficient, refraction_color);
                        v3_scale(refraction, -1, refraction_light.direction);
                        v3_copy(refraction_color, refraction_light.color);

                        v3_scale(ray_refracted.direction, best_refraction_t, ray_refracted.direction);
                        v3_subtract(ray_refracted.direction, new_ray.origin, new_ray.direction);
                        normalize(new_ray.direction);

                        cast(&new_ray, object_index, ray->direction, &refraction_light, INFINITY, color);
                        //Vector_add(color, refraction_color, color);
                }
                if (reflect_coefficient == -1) {
                        reflect_coefficient = 0;
                }
                if (refract_coefficient == -1) {
                        refract_coefficient = 0;
                }
                if (fabs(reflect_coefficient) < 0.00001 && fabs(refract_coefficient) < 0.00001) {
                        v3_copy(background_color, color);
                }else{
                        double color_coefficient = 1.0 - reflect_coefficient - refract_coefficient;
                        if (fabs(color_coefficient) < 0.0001) // account for numbers that are really close to 0, but still negative
                                color_coefficient = 0;
                        V3 object_color = {0, 0, 0};
                        v3_copy(objects[object_index].plane.diff_color, object_color);
                        v3_scale(object_color, color_coefficient, object_color);
                        color[0] += object_color[0];
                        color[1] += object_color[1];
                        color[2] += object_color[2];
                }
                //cleanup
                free(reflection_light.direction);
                free(reflection_light.color);
                free(refraction_light.direction);
                free(refraction_light.color);
        }
        for (int i=0; i<num_lights; i++) {
                // find new ray direction
                int best_o;
                double best_t;
                v3_zero(new_ray.direction);
                v3_subtract(lights[i].position, new_ray.origin, new_ray.direction);
                double distance_to_light = v3_len(new_ray.direction);
                normalize(new_ray.direction);

                // new check new ray for intersections with other objects
                get_best_solution(&new_ray, object_index, distance_to_light, &best_o, &best_t);

                if (best_o == -1) { // this means there was no object in the way between the current one and the light
                        cast(&new_ray, object_index, ray->direction, &lights[i], distance_to_light, color);
                }
                // there was an object in the way, so we don't do anything. It's shadow
        }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 Draws a complete scene using a recursive raytracing techinique

 Parameters
 ==========
 Image* img
        The image to draw the scene on
 double cam_width
        The width of the camera
 double cam_height
        The height of the camera
 OBJECT* objects
         The array of objects to draw

 Returns
 =======
 void
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void raytrace_scene(Image *img, double cam_width, double cam_height, OBJECT *objects) {

        V3 vp_pos = {0, 0, 1}; // view plane position

        V3 point = {0, 0, 0}; // point on viewplane where intersection happens

        double pixheight = (double)cam_height / (double)img->height;
        double pixwidth = (double)cam_width / (double)img->width;

        Ray ray = {
                .origin = {0, 0, 0},
                .direction = {0, 0, 0}
        };

        for (int x = 0; x < img->height; x++) {
                for (int y = 0; y < img->width; y++) {
                        v3_zero(ray.origin);
                        v3_zero(ray.direction);
                        point[0] = vp_pos[0] - cam_width/2.0 + pixwidth*(y + 0.5);
                        point[1] = -(vp_pos[1] - cam_height/2.0 + pixheight*(x + 0.5));
                        point[2] = vp_pos[2]; // set intersecting point Z to viewplane Z
                        normalize(point); // normalize the point

                        // store normalized point as our ray direction
                        v3_copy(point, ray.direction);
                        V3 color = {0, 0, 0};

                        int best_o; // index of 'best' or closest object
                        double best_t; // closest distance
                        get_best_solution(&ray, -1, INFINITY, &best_o, &best_t);

                        // set ambient color
                        if (best_t > 0 && best_t != INFINITY && best_o != -1) {// there was an intersection
                                trace(&ray, best_o, best_t,1,0,color);
                                shade_pixel(color, x, y, img);
                        }
                        else {
                                shade_pixel(background_color, x, y, img);
                        }
                }
        }
}

#endif
