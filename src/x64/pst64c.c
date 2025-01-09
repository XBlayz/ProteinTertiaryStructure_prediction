/**************************************************************************************
*
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2024/25
*
* Progetto dell'algoritmo Predizione struttura terziaria proteine 221 231 a
* in linguaggio assembly x86-64 + SSE
*
* F. Angiulli F. Fassetti S. Nisticò, novembre 2024
*
**************************************************************************************/

/*
*
* Software necessario per l'esecuzione:
*
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
*
* entrambi sono disponibili come pacchetti software
* installabili mediante il packaging tool del sistema
* operativo; per esempio, su Ubuntu, mediante i comandi:
*
*    sudo apt-get install nasm
*    sudo apt-get install gcc
*
* potrebbe essere necessario installare le seguenti librerie:
*
*    sudo apt-get install lib64gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
*
* Per generare il file eseguibile:
*
* nasm -f elf64 pst64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o pst64.o pst64c.c -o pst64c -lm && ./pst64c $pars
*
* oppure
*
* ./runpst64
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

//! Definizione della costante M_PI se non definito
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};

typedef struct {
	char* seq;					// sequenza di amminoacidi
	int N;						// lunghezza sequenza
	unsigned int sd; 			// seed per la generazione	casuale
	type to;					// temperatura INIZIALE
	type alpha;					// tasso di raffredamento
	type k;						// costante
	VECTOR hydrophobicity; 		// hydrophobicity
	VECTOR volume;				// volume
	VECTOR charge;				// charge
	VECTOR phi;					// vettore angoli phi
	VECTOR psi;					// vettore angoli psi
	type e;						// energy
	int display;
	int silent;
} params;


/*
*
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere
* 	memorizzate mediante array di array (float**).
*
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
*
*/

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

/*
*
* 	load_data
* 	=========
*
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
*
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
*
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice.
*****************************************************************************
*
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols; //! i, status (NON USATI)

	fp = fopen(filename, "rb");

	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	fread(&cols, sizeof(int), 1, fp); //! status = ...
	fread(&rows, sizeof(int), 1, fp); //! status = ...

	MATRIX data = alloc_matrix(rows,cols);
	fread(data, sizeof(type), rows*cols, fp); //! status = ...
	fclose(fp);

	*n = rows;
	*k = cols;

	return data;
}

/*
*
* 	load_seq
* 	=========
*
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
*
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*1 byte: matrix data in row-major order --> charatteri che compongono la stringa
*
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice.
*****************************************************************************
*
*/
char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols; //! i, status (NON USATI)

	fp = fopen(filename, "rb");

	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	fread(&cols, sizeof(int), 1, fp); //! status = ...
	fread(&rows, sizeof(int), 1, fp); //! status = ...


	char* data = alloc_char_matrix(rows,cols);
	fread(data, sizeof(char), rows*cols, fp); //! status = ...
	fclose(fp);

	*n = rows;
	*k = cols;

	return data;
}

/*
* 	save_data
* 	=========
*
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
*
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
*
*	Salva su file un array lineare composto da k elementi.
*
* 	Codifica del file:
* 	primi 4 byte: contenenti l'intero 1 		--> numero intero a 32 bit
* 	successivi 4 byte: numero di elementi k     --> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> k numero floating-point a precisione singola
*/
void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	//! int i (NON USATO)
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

/*
* 	gen_rnd_mat
* 	=========
*
*	Genera in maniera casuale numeri reali tra -pi e pi
*	per riempire una struttura dati di dimensione Nx1
*
*/
void gen_rnd_mat(VECTOR v, int N){
	int i;

	for(i=0; i<N; i++){
		// Campionamento del valore + scalatura
		v[i] = (random()*2 * M_PI) - M_PI;
	}
}

// PROCEDURE ASSEMBLY
extern void prova(params* input);

// --ROTATION--
// -UTILITY-
extern void sum_quad(VECTOR v, type* r); // r = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]

type modulo(VECTOR v) {
	// Modulo di un vettore 3D

	//* ASSEMBLY
	/* type r = 0;
	sum_quad(v, &r);

	return sqrt(r); */

	//* C
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

VECTOR cos_and_sin(type x) {
    // Calcolo potenze
	type x2 = x*x;
	type x3 = x2*x;
	type x4 = x3*x;
	type x5 = x4*x;
	type x6 = x5*x;
	type x7 = x6*x;

    // Calcolo coseno e seno
	VECTOR r = alloc_matrix(2, 1);
	r[0] = 1 - x2/2 + x4/24 - x6/720;
	r[1] = x - x3/6 + x5/120 - x7/5040;

	return r;
}

// -MAIN-
MATRIX rotation(VECTOR axis, type theta) {
	// Normalizzazione asse
	type mod_axis = modulo(axis);
	type div_axis = mod_axis*mod_axis;
	axis[0] /= div_axis;
	axis[1] /= div_axis;
	axis[2] /= div_axis;

	// Calcolo coseno e seno
	VECTOR cos_and_sin_v = cos_and_sin(theta/2.0);

	// Calcolo parametri matrice risultante
	type a = cos_and_sin_v[0];
	type b = -1 * axis[0] * cos_and_sin_v[1];
	type c = -1 * axis[1] * cos_and_sin_v[1];
	type d = -1 * axis[2] * cos_and_sin_v[1];

	// Deallocazione vettore coseno e seno
	dealloc_matrix(cos_and_sin_v);

	// Costruzione matrice risultante
	MATRIX r = alloc_matrix(3,3);
	r[0] = a*a + b*b - c*c - d*d;
	r[1] = 2*(b*c + a*d);
	r[2] = 2*(b*d - a*c);

	r[3] = 2*(b*c - a*d);
	r[4] = a*a + c*c - b*b - d*d;
	r[5] = 2*(c*d + a*b);

	r[6] = 2*(b*d + a*c);
	r[7] = 2*(c*d - a*b);
	r[8] = a*a + d*d - b*b - c*c;
	return r;
}

// --BACKBONE--
// -UTILITY-
VECTOR prodotto_vet_mat(VECTOR v, MATRIX m) {
    VECTOR r = alloc_matrix(3, 1);
    r[0] = v[0]*m[0] + v[1]*m[3] + v[2]*m[6];
    r[1] = v[0]*m[1] + v[1]*m[4] + v[2]*m[7];
    r[2] = v[0]*m[2] + v[1]*m[5] + v[2]*m[8];
    return r;
}

// -MAIN-
MATRIX backbone(int n, VECTOR phi, VECTOR psi) {
    // Definizione costanti raggi
	type R_CA_N = 1.46;
	type R_CA_C = 1.52;
	type R_C_N  = 1.33;

    // Definizione costanti angoli
	//type THETA_CA_C_N = 2.028; (Non usati)
	type THETA_C_N_CA = 2.124;
	//type THETA_N_CA_C = 1.940; (Non usati)

    // Allocazione matrice (vettore di vettori 3D) risultante
	MATRIX coords = alloc_matrix(3, n*3);

    // Posizionamento primo atomo backbone
	coords[0] = 0.0;
	coords[1] = 0.0;
	coords[2] = 0.0;

    // Posizionamento secondo atomo backbone
	coords[3] = R_CA_N;
	coords[4] = 0.0;
	coords[5] = 0.0;

	// Alloco le variabili da usare
	VECTOR v = alloc_matrix(3, 1);
	type mod_v = -1.0;

	MATRIX rot = NULL;
	VECTOR p = alloc_matrix(3, 1);
	VECTOR newv = NULL;

    // Posizionamento atomi restanti backbone
	for (int i = 0; i < n; i++) {
		int idx = i*3;
		if(i > 0) {
			v[0] = coords[(idx-1)*3] - coords[(idx-2)*3];
			v[1] = coords[(idx-1)*3+1] - coords[(idx-2)*3+1];
			v[2] = coords[(idx-1)*3+2] - coords[(idx-2)*3+2];
			mod_v = modulo(v);
			v[0] /= mod_v;
			v[1] /= mod_v;
			v[2] /= mod_v;

			rot = rotation(v, THETA_C_N_CA);
			p[0] = 0.0;
			p[1] = R_C_N;
			p[2] = 0.0;
			newv = prodotto_vet_mat(p, rot);
			dealloc_matrix(rot);

			coords[(idx)*3] = coords[(idx-1)*3] + newv[0];
			coords[(idx)*3+1] = coords[(idx-1)*3+1] + newv[1];
			coords[(idx)*3+2] = coords[(idx-1)*3+2] + newv[2];
			dealloc_matrix(newv);

			// -----

			v[0] = coords[(idx)*3] - coords[(idx-1)*3];
			v[1] = coords[(idx)*3+1] - coords[(idx-1)*3+1];
			v[2] = coords[(idx)*3+2] - coords[(idx-1)*3+2];
			mod_v = modulo(v);
			v[0] /= mod_v;
			v[1] /= mod_v;
			v[2] /= mod_v;

			rot = rotation(v, phi[i]);
			p[0] = 0.0;
			p[1] = R_CA_N;
			p[2] = 0.0;
			newv = prodotto_vet_mat(p, rot);
			dealloc_matrix(rot);

			coords[(idx+1)*3] = coords[(idx)*3] + newv[0];
			coords[(idx+1)*3+1] = coords[(idx)*3+1] + newv[1];
			coords[(idx+1)*3+2] = coords[(idx)*3+2] + newv[2];
			dealloc_matrix(newv);
		}
		v[0] = coords[(idx+1)*3] - coords[(idx)*3];
		v[1] = coords[(idx+1)*3+1] - coords[(idx)*3+1];
		v[2] = coords[(idx+1)*3+2] - coords[(idx)*3+2];
		mod_v = modulo(v);
		v[0] /= mod_v;
		v[1] /= mod_v;
		v[2] /= mod_v;

		rot = rotation(v, psi[i]);
		p[0] = 0.0;
		p[1] = R_CA_C;
		p[2] = 0.0;
		newv = prodotto_vet_mat(p, rot);
		dealloc_matrix(rot);

		coords[(idx+2)*3] = coords[(idx+1)*3] + newv[0];
		coords[(idx+2)*3+1] = coords[(idx+1)*3+1] + newv[1];
		coords[(idx+2)*3+2] = coords[(idx+1)*3+2] + newv[2];
        dealloc_matrix(newv);
	}

	// Deallocazione
	dealloc_matrix(v);
	dealloc_matrix(p);

	return coords;
}

// --RAMA-ENERGY--
// -MAIN-
type rama_energy(int n, VECTOR phi, VECTOR psi) {
	// Definizione costanti funzione
	type a_phi = -57.8;
	type a_psi = -47.0;
	type b_phi = -119.0;
	type b_psi = 113.0;

	type energy = 0.0;
	for(int i=0; i<n; i++) {
		VECTOR tmp = alloc_matrix(3, 1);
		tmp[0] = phi[i] - a_phi;
		tmp[1] = psi[i] - a_psi;
		tmp[2] = 0;
		type a = modulo(tmp);
		tmp[0] = phi[i] - b_phi;
		tmp[1] = psi[i] - b_psi;
		tmp[2] = 0;
		type b = modulo(tmp);
		dealloc_matrix(tmp);
		// Calcolo energia
		energy += 0.5 * (a<b ? a : b);
	}
	return energy;
}

// --HYDROPHOBIC-ENERGY--
// -UTILITY-
type dist(VECTOR i, VECTOR j) {
	// Calcolo distanza tra due vettori 3D
	VECTOR tmp = alloc_matrix(3, 1);
	tmp[0] = j[0]-i[0];
	tmp[1] = j[1]-i[1];
	tmp[2] = j[2]-i[2];
	type r = modulo(tmp);
	dealloc_matrix(tmp);

	return r;
}

// -MAIN-
extern void h_energy(char* s, int n, MATRIX coords, type* e, type* h);

type hydrophobic_energy(char* s, int n, MATRIX coords) {
	type energy = 0;

	//* ASSEMBLY
	h_energy(s, n, coords, &energy, hydrophobicity);

	//* C
	/* for(int i = 0; i<n; i++) {
		for(int j = i+1; j<n; j++) {
			VECTOR vi = alloc_matrix(3, 1);
			// Indicizzazione atomi Ca [(i*3)+1]
			vi[0] = coords[((i*3)+1)*3];
			vi[1] = coords[((i*3)+1)*3+1];
			vi[2] = coords[((i*3)+1)*3+2];
			VECTOR vj = alloc_matrix(3, 1);
			// Indicizzazione atomi Ca [(i*3)+1]
			vj[0] = coords[((j*3)+1)*3];
			vj[1] = coords[((j*3)+1)*3+1];
			vj[2] = coords[((j*3)+1)*3+2];
			type d = dist(vi, vj);
			dealloc_matrix(vi);
			dealloc_matrix(vj);

			if(d < 10.0) {
				// Calcolo energia
				energy += (hydrophobicity[s[i]-'a'] * hydrophobicity[s[j]-'a']) / d;
			}
		}
	} */
	return energy;
}

// --ELECTROSTATIC-ENERGY--
// -MAIN-
extern void e_energy(char* s, int n, MATRIX coords, type* e, type* c);

type electrostatic_energy(char* s, int n, MATRIX coords) {
	type energy = 0;

	//* ASSEMBLY
	e_energy(s, n, coords, &energy, charge);

	//* C
	/* for(int i = 0; i<n; i++) {
		for(int j = i+1; j<n; j++) {
			VECTOR vi = alloc_matrix(3, 1);
			// Indicizzazione atomi Ca [(i*3)+1]
			vi[0] = coords[((i*3)+1)*3];
			vi[1] = coords[((i*3)+1)*3+1];
			vi[2] = coords[((i*3)+1)*3+2];
			VECTOR vj = alloc_matrix(3, 1);
			// Indicizzazione atomi Ca [(i*3)+1]
			vj[0] = coords[((j*3)+1)*3];
			vj[1] = coords[((j*3)+1)*3+1];
			vj[2] = coords[((j*3)+1)*3+2];
			type d = dist(vi, vj);
			dealloc_matrix(vi);
			dealloc_matrix(vj);

			if(i != j && d < 10.0 && charge[s[i]-'a'] != 0 && charge[s[j]-'a'] != 0) {
				// Calcolo energia
				energy += (charge[s[i]-'a'] * charge[s[j]-'a']) / (d*4.0);
			}
		}
	} */
	return energy;
}

// --PACKING-ENERGY--
// -MAIN-
extern void p_energy(char* s, int n, MATRIX coords, type* e, type* v);

type packing_energy(char* s, int n, MATRIX coords) {
	type energy = 0;

	//* ASSEMBLY
	p_energy(s, n, coords, &energy, volume);

	//* C
	/* for(int i = 0; i<n; i++) {
		type density = 0;
		for(int j = 0; j<n; j++) {
			VECTOR vi = alloc_matrix(3, 1);
			// Indicizzazione atomi Ca [(i*3)+1]
			vi[0] = coords[((i*3)+1)*3];
			vi[1] = coords[((i*3)+1)*3+1];
			vi[2] = coords[((i*3)+1)*3+2];
			VECTOR vj = alloc_matrix(3, 1);
			// Indicizzazione atomi Ca [(j*3)+1]
			vj[0] = coords[((j*3)+1)*3];
			vj[1] = coords[((j*3)+1)*3+1];
			vj[2] = coords[((j*3)+1)*3+2];
			type d = dist(vi, vj);
			dealloc_matrix(vi);
			dealloc_matrix(vj);

			if(i != j && d < 10.0) {
				// Calcolo densità
				density += volume[s[j]-'a'] / (d*d*d);
			}
		}
		// Calcolo energia
		type tmp = volume[s[i]-'a'] - density;
		energy += tmp*tmp;
	} */
	return energy;
}

// --ENERGY--
// -MAIN-
type energy(char* s, int n, VECTOR phi, VECTOR psi) {
	// Calcolo vettore delle coordinate del backbone
	MATRIX coords = backbone(n, phi, psi);

	// Calcolo delle energie
	type rama = rama_energy(n, phi, psi);
	type hydro = hydrophobic_energy(s, n, coords);
	type elec = electrostatic_energy(s, n, coords);
	type pack = packing_energy(s, n, coords);

	// Deallocazione
	dealloc_matrix(coords);

	// Definizione pesi delle energie
	type w_rama = 1.0;
	type w_hydro = 0.5;
	type w_elec = 0.2;
	type w_pack = 0.3;

	// Somma pesata delle energie
	return w_rama*rama + w_hydro*hydro + w_elec*elec + w_pack*pack;
}

void pst(params* input) {
	// Calcolo dell'energia
	type temperatura = input->to;
	input->e = energy(input->seq, input->N, input->phi, input->psi);

	int t = 0;
	do {
		//! int i = random() * input->N;
		int i = rand() % input->N;
		type dev_phi = (random()*2 * M_PI) - M_PI;
		type dev_psi = (random()*2 * M_PI) - M_PI;

		input->phi[i] += dev_phi;
		input->psi[i] += dev_psi;
		type new_e = energy(input->seq, input->N, input->phi, input->psi);
		type delta_e = new_e - input->e;

		if(delta_e <= 0) {
			input->e = new_e;
		} else {
			type p = exp(-1*(delta_e / (input->k * temperatura)));
			type r = random();

			if(r <= p) {
				input->e = new_e;
			}
			input->phi[i] -= dev_phi;
			input->psi[i] -= dev_psi;
		}

		// Aggiorno la temperatura
		t++;
		temperatura = input->to - sqrt(input->alpha * t);
	} while(temperatura > 0);
}

int main(int argc, char** argv) {
	char fname_phi[256];
	char fname_psi[256];
	char* seqfilename = NULL;
	clock_t t;
	float time;
	int d;

	//
	// Imposta i valori di default dei parametri
	//
	params* input = malloc(sizeof(params));
	input->seq = NULL;
	input->N = -1;
	input->to = -1;
	input->alpha = -1;
	input->k = -1;
	input->sd = -1;
	input->phi = NULL;
	input->psi = NULL;
	input->e = -1;
	input->silent = 0;
	input->display = 0;
	input->hydrophobicity = hydrophobicity;
	input->volume = volume;
	input->charge = charge;


	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//
	if(argc <= 1){
		printf("%s -seq <SEQ> -to <to> -alpha <alpha> -k <k> -sd <sd> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tSEQ: il nome del file ds2 contenente la sequenza amminoacidica\n");
		printf("\tto: parametro di temperatura\n");
		printf("\talpha: tasso di raffredamento\n");
		printf("\tk: costante\n");
		printf("\tsd: seed per la generazione casuale\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//
	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-seq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			seqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-to") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing to value!\n");
				exit(1);
			}
			input->to = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-alpha") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing alpha value!\n");
				exit(1);
			}
			input->alpha = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sd") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing seed value!\n");
				exit(1);
			}
			input->sd = atoi(argv[par]);
			par++;
		}else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(seqfilename == NULL || strlen(seqfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	input->seq = load_seq(seqfilename, &input->N, &d);


	if(d != 1){
		printf("Invalid size of sequence file, should be %ix1!\n", input->N);
		exit(1);
	}

	if(input->to <= 0){
		printf("Invalid value of to parameter!\n");
		exit(1);
	}

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	if(input->alpha <= 0){
		printf("Invalid value of alpha parameter!\n");
		exit(1);
	}

	input->phi = alloc_matrix(input->N, 1);
	input->psi = alloc_matrix(input->N, 1);
	// Impostazione seed
	srand(input->sd);
	// Inizializzazione dei valori
	gen_rnd_mat(input->phi, input->N);
	gen_rnd_mat(input->psi, input->N);

	//
	// Visualizza il valore dei parametri
	//
	if(!input->silent){
		printf("Dataset file name: '%s'\n", seqfilename);
		printf("Sequence lenght: %d\n", input->N);
	}

	// COMMENTARE QUESTA RIGA!
	// prova(input);
	//

	//
	// Predizione struttura terziaria
	//
	t = clock();
	pst(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("PST time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname_phi, "out64_%d_%d_phi.ds2", input->N, input->sd);
	save_out(fname_phi, input->phi, input->N);
	sprintf(fname_psi, "out64_%d_%d_psi.ds2", input->N, input->sd);
	save_out(fname_psi, input->psi, input->N);
	if(input->display){
		if(input->phi == NULL || input->psi == NULL)
			printf("out: NULL\n");
		else{
			int i; //! j (NON USATO)
			printf("energy: %f, phi: [", input->e);
			for(i=0; i<input->N; i++){
				printf("%f,", input->phi[i]);
			}
			printf("]\n");
			printf("psi: [");
			for(i=0; i<input->N; i++){
				printf("%f,", input->psi[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->phi);
	dealloc_matrix(input->psi);
	free(input);

	return 0;
}