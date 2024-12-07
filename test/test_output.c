#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include <xmmintrin.h>

char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols; //* i, status (NON USATI)

	fp = fopen(filename, "rb");

	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	fread(&cols, sizeof(int), 1, fp); //* status = ...
	fread(&rows, sizeof(int), 1, fp); //* status = ...


	char* data = alloc_char_matrix(rows,cols);
	fread(data, sizeof(char), rows*cols, fp); //* status = ...
	fclose(fp);

	*n = rows;
	*k = cols;

	return data;
}

void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	//* int i (NON USATO)
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

int main(int argc, char** argv) {
    char fname_phi[256];
    char fname_psi[256];
    int n, sd;

    sprintf(fname_phi, "out32_%d_%d_phi.ds2", n, sd);
	save_out(fname_phi, input->phi, n);
	sprintf(fname_psi, "out32_%d_%d_psi.ds2", n, sd);
	save_out(fname_psi, input->psi, n);
}