#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>


#define	type        double
#define	MATRIX      type*


// --MEMORY ALLOCATION--
void* get_block(int size, int elements) {
	return _mm_malloc(elements*size,16);
}
void free_block(void* p) {
	_mm_free(p);
}
MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}
int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}
char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}
void dealloc_matrix(void* mat) {
	free_block(mat);
}

// --ASSEMBLY FUNCTIONS--
extern void prova_ass(double* data, double* r);

// --MAIN--
int main(int argc, char* argv[]) {
    double* data = alloc_matrix(8,1);
    data[0] = 1.0;
    data[1] = 2.0;
    data[2] = 3.0;
    data[3] = 0.0;

    data[4] = 4.0;
    data[5] = 5.0;
    data[6] = 6.0;
    data[7] = 0.0;

    double r = 0.0;
    prova_ass(data, &r);
    double ref = sqrt((data[4]-data[0])*(data[4]-data[0]) + (data[5]-data[1])*(data[5]-data[1]) + (data[6]-data[2])*(data[6]-data[2]));
    printf("---\n");
    printf("%f (%f)\n", r, ref);

    dealloc_matrix(data);
    return 0;
}