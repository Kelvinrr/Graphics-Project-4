#ifndef _JSON_PARSER_H_
#define _JSON_PARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <stdbool.h>

typedef uint8_t boolean;

#include "pplib.h"

const uint16_t MAX_OBJECTS = 256;
const uint16_t MAX_COLOR_VAL = 255;
const uint8_t CAM = 1;
const uint8_t SPH = 2;
const uint8_t PLAN = 3;
const uint8_t LIG = 5;
const uint8_t SPOTLIG = 6;
const uint8_t TRACING = 7;

typedef struct CAMERA {
        double width;
        double height;
} CAMERA;

typedef struct SPHERE {
        double *diff_color;
        double *spec_color;
        double *position;
        double radius;
        double reflect;
        double refract;
        double ior;
} SPHERE;

typedef struct PLANE {
        double *diff_color;
        double *spec_color;
        double *position;
        double *normal;
        double reflect;
        double refract;
        double ior;
} PLANE;

typedef struct LIGHT {
        int type;
        double *color;
        double *position;
        double *direction;
        double theta_deg;
        double rad_att0;
        double rad_att1;
        double rad_att2;
        double ang_att0;
} LIGHT;

typedef struct OBJECT {
        int type;
        union {
                CAMERA camera;
                SPHERE sphere;
                PLANE plane;
        };
} OBJECT;

int line = 1;  // global variable, it will tells us which line is not correct

OBJECT objects[MAX_OBJECTS]; // Allocate an array for All Objects in Json File
LIGHT lights[MAX_OBJECTS];      // allocate space for lights

int num_lights;
int num_objects;

int next_c(FILE* json) {
        int c = fgetc(json);

        if (c == '\n') {
                line++;;
        }
        if (c == EOF) {
                fprintf(stderr, "Error: next_c: Unexpected EOF: %d\n", line);
                exit(1);
        }
        return c;
}

void skip_ws(FILE* json) {
        int c = next_c(json);
        while (isspace(c)) {
                c = next_c(json);
        }
        ungetc(c, json);
}

void expect_c(FILE* json, int d) {
        int c = next_c(json);
        if (c == d) return;
        fprintf(stderr, "Error: Expected '%c': %d\n", d, line);
        exit(1);
}

double next_number(FILE* json) {
        double value;
        int res = fscanf(json, "%lf", &value);
        if (res == EOF) {
                fprintf(stderr, "Error: Expected a number but found EOF: %d\n", line);
                exit(1);
        }
        //printf("next_number: %lf\n", value);
        return value;
}

int check_color_val(double v) {
        if (v < 0.0 || v > MAX_COLOR_VAL)
                return 0;
        return 1;
}


int check_light_color_val(double v) {
        if (v < 0.0)
                return 0;
        return 1;
}

double* next_vector(FILE* json) {
        double* v = malloc(3*sizeof(double));
        expect_c(json, '[');
        skip_ws(json);
        v[0] = next_number(json);
        skip_ws(json);
        expect_c(json, ',');
        skip_ws(json);
        v[1] = next_number(json);
        skip_ws(json);
        expect_c(json, ',');
        skip_ws(json);
        v[2] = next_number(json);
        skip_ws(json);
        expect_c(json, ']');
        return v;
}

double* next_color(FILE* json, boolean is_rgb) {
        double* v = malloc(sizeof(double)*3);
        skip_ws(json);
        expect_c(json, '[');
        skip_ws(json);
        v[0] = next_number(json);
        skip_ws(json);
        expect_c(json, ',');
        skip_ws(json);
        v[1] = next_number(json);
        skip_ws(json);
        expect_c(json, ',');
        skip_ws(json);
        v[2] = next_number(json);
        skip_ws(json);
        expect_c(json, ']');
        // check that all values are valid
        if (is_rgb) {
                if (!check_color_val(v[0]) ||
                    !check_color_val(v[1]) ||
                    !check_color_val(v[2])) {
                        fprintf(stderr, "Error: next_color: rgb value out of range: %d\n", line);
                        exit(1);
                }
        }
        else {
                if (!check_light_color_val(v[0]) ||
                    !check_light_color_val(v[1]) ||
                    !check_light_color_val(v[2])) {
                        fprintf(stderr, "Error: next_color: light value out of range: %d\n", line);
                        exit(1);
                }

        }
        return v;
}

char* next_string(FILE* json) {
        char buffer[129];
        int c = next_c(json);
        if (c != '"') {
                fprintf(stderr, "Error: Expected string on line %d.\n", line);
                exit(1);
        }
        c = next_c(json);
        int i = 0;
        while (c != '"') {
                if (i >= 128) {
                        fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
                        exit(1);
                }
                if (c == '\\') {
                        fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
                        exit(1);
                }
                if (c < 32 || c > 126) {
                        fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
                        exit(1);
                }
                buffer[i] = c;
                i += 1;
                c = next_c(json);
        }
        buffer[i] = 0;
        return strdup(buffer);
}

double* next_coefficient(FILE* json){
        double* v = malloc(10*sizeof(double));

        skip_ws(json);
        expect_c(json, '[');
        skip_ws(json);
        v[0] = next_number(json);
        for(int i=1; i<10; i++) {
                skip_ws(json);
                expect_c(json, ',');
                skip_ws(json);
                v[i] = next_number(json);
        }
        skip_ws(json);
        expect_c(json, ']');
        return v;
}

void read_scene(const char* filename) {
        FILE* json = fopen(filename, "r");

        if (json == NULL) {
                fprintf(stderr, "Error: Could not open file\n");
                exit(1);
        }
        skip_ws(json);

        // find beginning of the list
        int c  = next_c(json);
        if (c != '[') {
                fprintf(stderr, "Error: read_json: JSON file must begin with [\n");
                exit(1);
        }
        skip_ws(json);
        c = next_c(json);

        // check if file empty
        if (c == ']' || c == EOF) {
                fprintf(stderr, "Error: read_json: Empty json file\n");
                exit(1);
        }
        skip_ws(json);

        int object_counter = 0;
        int light_counter = 0;
        int object_type = 0;
        boolean find_object = true;
        // find the objects
        while (find_object) {
                //c  = next_c(json);
                if (object_counter > MAX_OBJECTS) {
                        fprintf(stderr, "Error: read_json: Number of objects is too large: %d\n", line);
                        exit(1);
                }
                if (c == ']') {
                        fprintf(stderr, "Error: read_json: Unexpected ']': %d\n", line);
                        fclose(json);
                        exit(1);
                }
                if (c == '{') { // found an object
                        skip_ws(json);
                        char *key = next_string(json);
                        if (strcmp(key, "type") != 0) {
                                fprintf(stderr, "Error: read_json: First key of an object must be 'type': %d\n", line);
                                exit(1);
                        }
                        skip_ws(json);
                        // get the colon
                        expect_c(json, ':');
                        skip_ws(json);

                        char *type = next_string(json);
                        if (strcmp(type, "camera") == 0) {
                                object_type = CAM;
                                objects[object_counter].type = CAM;
                        }
                        else if (strcmp(type, "sphere") == 0) {
                                object_type = SPH;
                                objects[object_counter].type = SPH;
                        }
                        else if (strcmp(type, "plane") == 0) {
                                object_type = PLAN;
                                objects[object_counter].type = PLAN;
                        }
                        else if (strcmp(type, "light") == 0) {
                                object_type = LIG;
                        }
                        else {
                                exit(1);
                        }

                        skip_ws(json);

                        while (true) {
                                //  , }
                                c = next_c(json);
                                if (c == '}') {
                                        // stop parsing this object
                                        break;
                                }
                                else if (c == ',') {
                                        // read another field
                                        skip_ws(json);
                                        char* key = next_string(json);
                                        skip_ws(json);
                                        expect_c(json, ':');
                                        skip_ws(json);
                                        if (strcmp(key, "width") == 0) {
                                                if (object_type != CAM) {
                                                        fprintf(stderr, "Error: read_json: Width cannot be set on this type: %d\n", line);
                                                        exit(1);
                                                }
                                                double temp = next_number(json);
                                                if (temp <= 0) {
                                                        fprintf(stderr, "Error: read_json: width must be positive: %d\n", line);
                                                        exit(1);
                                                }
                                                objects[object_counter].camera.width = temp;

                                        }
                                        else if (strcmp(key, "height") == 0) {
                                                if (object_type != CAM) {
                                                        fprintf(stderr, "Error: read_json: Width cannot be set on this type: %d\n", line);
                                                        exit(1);
                                                }
                                                double temp = next_number(json);
                                                if (temp <= 0) {
                                                        fprintf(stderr, "Error: read_json: height must be positive: %d\n", line);
                                                        exit(1);
                                                }
                                                objects[object_counter].camera.height = temp;
                                        }
                                        else if (strcmp(key, "radius") == 0) {
                                                if (object_type != SPH) {
                                                        fprintf(stderr, "Error: read_json: Radius cannot be set on this type: %d\n", line);
                                                        exit(1);
                                                }
                                                double temp = next_number(json);
                                                if (temp <= 0) {
                                                        fprintf(stderr, "Error: read_json: radius must be positive: %d\n", line);
                                                        exit(1);
                                                }
                                                objects[object_counter].sphere.radius = temp;
                                        }
                                        else if (strcmp(key, "theta") == 0) {
                                                if (object_type != LIG) {
                                                        fprintf(stderr, "Error: read_json: Theta cannot be set on this type: %d\n", line);
                                                        exit(1);
                                                }
                                                double temp = next_number(json);
                                                if (temp > 0.0) {
                                                        lights[light_counter].type = SPOTLIG;
                                                }
                                                else if (temp < 0.0) {
                                                        fprintf(stderr, "Error: read_json: theta must be >= 0: %d\n", line);
                                                        exit(1);
                                                }
                                                lights[light_counter].theta_deg = temp;
                                        }

                                        else if (strcmp(key, "radial-a0") == 0) {
                                                if (object_type != LIG) {
                                                        fprintf(stderr, "Error: read_json: Radial-a0 cannot be set on this type: %d\n", line);
                                                        exit(1);
                                                }
                                                double temp = next_number(json);
                                                if (temp < 0) { // TODO: find out if this should be <=
                                                        fprintf(stderr, "Error: read_json: radial-a0 must be positive: %d\n", line);
                                                        exit(1);
                                                }
                                                lights[light_counter].rad_att0 = temp;
                                        }
                                        else if (strcmp(key, "radial-a1") == 0) {
                                                if (object_type != LIG) {
                                                        fprintf(stderr, "Error: read_json: Radial-a0 cannot be set on this type: %d\n", line);
                                                        exit(1);
                                                }
                                                double temp = next_number(json);
                                                if (temp < 0) { // TODO: find out if this should be <=
                                                        fprintf(stderr, "Error: read_json: radial-a1 must be positive: %d\n", line);
                                                        exit(1);
                                                }
                                                lights[light_counter].rad_att1 = temp;
                                        }
                                        else if (strcmp(key, "radial-a2") == 0) {
                                                if (object_type != LIG) {
                                                        fprintf(stderr, "Error: read_json: Radial-a0 cannot be set on this type: %d\n", line);
                                                        exit(1);
                                                }
                                                double temp = next_number(json);
                                                if (temp < 0) { // TODO: find out if this should be <=
                                                        fprintf(stderr, "Error: read_json: radial-a2 must be positive: %d\n", line);
                                                        exit(1);
                                                }
                                                lights[light_counter].rad_att2 = temp;
                                        }
                                        else if (strcmp(key, "angular-a0") == 0) {
                                                if (object_type != LIG) {
                                                        fprintf(stderr, "Error: read_json: Radial-a0 cannot be set on this type: %d\n", line);
                                                        exit(1);
                                                }
                                                double temp = next_number(json);
                                                if (temp < 0) { // TODO: find out if this should be <=
                                                        fprintf(stderr, "Error: read_json: angular-a0 must be positive: %d\n", line);
                                                        exit(1);
                                                }
                                                lights[light_counter].ang_att0 = temp;
                                        }
                                        else if (strcmp(key, "color") == 0) {
                                                if (object_type != LIG) {
                                                        fprintf(stderr, "Error: Just plain 'color' vector can only be applied to a light object\n");
                                                        exit(1);
                                                }
                                                lights[light_counter].color = next_color(json, false);
                                        }
                                        else if (strcmp(key, "direction") == 0) {
                                                if (object_type != LIG) {
                                                        fprintf(stderr, "Error: Direction vector can only be applied to a light object\n");
                                                        exit(1);
                                                }
                                                lights[light_counter].type = SPOTLIG;
                                                lights[light_counter].direction = next_vector(json);
                                        }
                                        else if (strcmp(key, "specular_color") == 0) {
                                                if (object_type == SPH) {
                                                        objects[object_counter].sphere.spec_color = next_color(json, true);
                                                }
                                                else if (object_type == PLAN) {
                                                        objects[object_counter].plane.spec_color = next_color(json, true);
                                                }
                                                else {
                                                        fprintf(stderr, "Error: read_json: speculaor_color vector can't be applied here: %d\n", line);
                                                        exit(1);
                                                }
                                        }
                                        else if (strcmp(key, "diffuse_color") == 0) {
                                                if (object_type == SPH) {
                                                        objects[object_counter].sphere.diff_color = next_color(json, true);
                                                }
                                                else if (object_type == PLAN) {
                                                        objects[object_counter].plane.diff_color = next_color(json, true);
                                                }
                                                else {
                                                        fprintf(stderr, "Error: read_json: diffuse_color vector can't be applied here: %d\n", line);
                                                        exit(1);
                                                }
                                        }
                                        else if (strcmp(key, "position") == 0) {
                                                if (object_type == SPH) {
                                                        objects[object_counter].sphere.position = next_vector(json);
                                                }
                                                else if (object_type == PLAN) {
                                                        objects[object_counter].plane.position = next_vector(json);
                                                }
                                                else if (object_type == LIG) {
                                                        lights[light_counter].position = next_vector(json);
                                                }
                                                else {
                                                        fprintf(stderr, "Error: read_json: Position vector can't be applied here: %d\n", line);
                                                        exit(1);
                                                }

                                        }
                                        else if (strcmp(key, "reflectivity") == 0) {
                                                if (object_type == SPH) {
                                                        objects[object_counter].sphere.reflect = next_number(json);
                                                }
                                                else if(object_type == PLAN) {
                                                        objects[object_counter].plane.reflect = next_number(json);
                                                }
                                                else{
                                                        fprintf(stderr, "Error: read_json: Reflectivity can't be applied here: %d\n", line);
                                                        exit(1);
                                                }
                                        }
                                        else if (strcmp(key, "refractivity") == 0) {
                                                if (object_type == SPH) {
                                                        objects[object_counter].sphere.refract = next_number(json);
                                                }
                                                else if(object_type == PLAN) {
                                                        objects[object_counter].plane.refract = next_number(json);
                                                }
                                                else{
                                                        fprintf(stderr, "Error: read_json: Refractivity can't be applied here: %d\n", line);
                                                        exit(1);
                                                }
                                        }
                                        else if (strcmp(key, "ior") == 0) {
                                                if (object_type == PLAN) {
                                                        objects[object_counter].plane.ior = next_number(json);
                                                }
                                                else if(object_type == SPH) {
                                                        objects[object_counter].sphere.ior = next_number(json);
                                                }
                                                else{
                                                        fprintf(stderr, "Error: read_json: ior can't be applied here: %d\n", line);
                                                        exit(1);
                                                }
                                        }
                                        else if (strcmp(key, "normal") == 0) {
                                                if (object_type != PLAN) {
                                                        fprintf(stderr, "Error: read_json: Normal vector can't be applied here: %d\n", line);
                                                        exit(1);
                                                }
                                                else
                                                        objects[object_counter].plane.normal = next_vector(json);
                                        }
                                        else {
                                                fprintf(stderr, "Error: read_json: '%s' not a valid object: %d\n", key, line);
                                                exit(1);
                                        }

                                        skip_ws(json);
                                }
                                else {
                                        fprintf(stderr, "Error: read_json: Unexpected value '%c': %d\n", c, line);
                                        exit(1);
                                }
                        }
                        skip_ws(json);
                        c = next_c(json);
                        if (c == ',') {
                                // noop
                                skip_ws(json);
                        }
                        else if (c == ']') {
                                find_object = false;
                        }
                        else {
                                fprintf(stderr, "Error: read_json: Expecting comma or ]: %d\n", line);
                                exit(1);
                        }
                }
                if (object_type == LIG) {
                        light_counter++;
                }

                else{
                        object_counter++;
                }

                if (find_object) {
                        c = next_c(json);
                }
        }
  
        fclose(json);
        num_lights = light_counter;
        num_objects = object_counter;
}


//Get Objects
void get_objects(OBJECT *object) {
        int i = 0;
        while (i < MAX_OBJECTS && object[i].type > 0) {
                printf("object type: %d\n", object[i].type);
                if (object[i].type == CAM) {
                        printf("height: %lf\n", object[i].camera.height);
                        printf("width: %lf\n", object[i].camera.width);
                }
                else if (object[i].type == SPH) {
                        printf("color: %lf %lf %lf\n", object[i].sphere.spec_color[0],
                               object[i].sphere.spec_color[1],
                               object[i].sphere.spec_color[2]);
                        printf("position: %lf %lf %lf\n", object[i].sphere.position[0],
                               object[i].sphere.position[1],
                               object[i].sphere.position[2]);
                        printf("radius: %lf\n", object[i].sphere.radius);
                }
                else if (object[i].type == PLAN) {
                        printf("color: %lf %lf %lf\n", object[i].plane.spec_color[0],
                               object[i].plane.spec_color[1],
                               object[i].plane.spec_color[2]);
                        printf("position: %lf %lf %lf\n", object[i].plane.position[0],
                               object[i].plane.position[1],
                               object[i].plane.position[2]);
                        printf("normal: %lf %lf %lf\n", object[i].plane.normal[0],
                               object[i].plane.normal[1],
                               object[i].plane.normal[2]);
                }
                else {
                        printf("unsupported value\n");
                }
                i++;
        }
        printf("end at i=%d\n", i);
}

#endif
