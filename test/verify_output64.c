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

// --LOAD FUNCTIONS--
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols;

	fp = fopen(filename, "rb");

	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	fread(&cols, sizeof(int), 1, fp);
	fread(&rows, sizeof(int), 1, fp);

	MATRIX data = alloc_matrix(rows,cols);
	fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);

	*n = rows;
	*k = cols;

	return data;
}

#define	PI	3.14159265358979323846
type rad_to_deg(type rad) {
	return rad * 180.0 / PI;
}

// --MAIN--
int main(int argc, char const *argv[])
{
    printf("---\n");
    printf("[START]\n");
    printf("---\n");

    // -Load output-
    char fname_phi[] = "../src/x64/out64_256_3_phi.ds2";
	char fname_psi[] = "../src/x64/out64_256_3_psi.ds2";

    int r_phi, c_phi, r_psi, c_psi;

    MATRIX phi = load_data(fname_phi, &r_phi, &c_phi);
    MATRIX psi = load_data(fname_psi, &r_psi, &c_psi);

    // -Load reference-
    char fname_phi_ref[] = "./ref_output/phi_double.ds2";
	char fname_psi_ref[] = "./ref_output/psi_double.ds2";

    int r_phi_ref, c_phi_ref, r_psi_ref, c_psi_ref;

    MATRIX phi_ref = load_data(fname_phi_ref, &r_phi_ref, &c_phi_ref);
    MATRIX psi_ref = load_data(fname_psi_ref, &r_psi_ref, &c_psi_ref);

    // -Verify size-
    printf("[SIZE INFO]: %dx%d (%dx%d)", r_phi, c_phi, r_phi_ref, c_phi_ref);
    if(r_phi != r_phi_ref || c_phi != c_phi_ref){
        printf(" [ERROR]\n");
        return -1;
    } else {
        printf(" [OK]\n");
    }
    printf("[SIZE INFO]: %dx%d (%dx%d)", r_psi, c_psi, r_psi_ref, c_psi_ref);
    if(r_psi != r_psi_ref || c_psi != c_psi_ref){
        printf(" [ERROR]\n");
        return -1;
    } else {
        printf(" [OK]\n");
    }
    printf("---\n");

    double qme1 = 0;
    double qme2 = 0;
    for(int i = 0; i < r_phi*c_phi; i++){
        printf("[%d-PHI]: %f (%f)", i, rad_to_deg(phi[i]), rad_to_deg(phi_ref[i]));
        qme1 += (phi[i] - phi_ref[i]) * (phi[i] - phi_ref[i]);
        if(fabs(phi[i] - phi_ref[i]) > 0.00001){
            printf(" [ERROR]\n");
        } else {
            printf(" [OK]\n");
        }
    }
    printf("---\n");
    for(int i = 0; i < r_psi*c_psi; i++){
        printf("[%d-PSI]: %f (%f)", i, rad_to_deg(psi[i]), rad_to_deg(psi_ref[i]));
        qme2 += (psi[i] - psi_ref[i]) * (psi[i] - psi_ref[i]);
        if(fabs(psi[i] - psi_ref[i]) > 0.00001){
            printf(" [ERROR]\n");
        } else {
            printf(" [OK]\n");
        }
    }
    printf("---\n");
    printf("[QUANTITY MEAN ERROR]: %f, %f\n", qme1/(r_phi*c_phi), qme2/(r_psi*c_psi));
    printf("---\n");
    printf("[END]\n");
    printf("---\n");
    return 0;
}
