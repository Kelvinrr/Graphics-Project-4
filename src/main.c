#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "raytrace.h"
 
#define MAX_SIZE 1024

int main(int argc, const char * argv[]) {
        if (argc < 5) {
                fprintf(stderr, "Error: main: You must have at least 4 arguments\n");
                exit(1);
        }

        if (atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0) {
                fprintf(stderr, "Error: main: width and height parameters must be > 0\n");
                exit(1);
        }

        const char *input_file = argv[3];
        const char *output_file = argv[4];

        read_scene(input_file);
  
        Image *image = (Image *)malloc(sizeof(Image));

        image->width = atoi(argv[1]);
        image->height = atoi(argv[2]);
        image->maxval = 255;
        image->data = malloc(image->width * image->height*4);
        int pos = get_camera(objects);

        printf("Image Width: %d\n", image->width);
        printf("Image Height: %d\n", image->height);
        printf("Image maxval: %d\n", image->maxval);
  
        if (pos == -1) {
                fprintf(stderr, "Error: main: No camera object found in data\n");
                exit(1);
        }

        raytrace_scene(image, objects[pos].camera.width, objects[pos].camera.height, objects);
        ImageWrite(image, output_file, 6);

        return 0;
}
