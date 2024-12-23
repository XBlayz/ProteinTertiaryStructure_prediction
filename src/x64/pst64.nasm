; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli, F. Fassetti, S. Nisticò
; 12/11/2024
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software
; installabili mediante il packaging tool del sistema
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf64 pst64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

section .bss			; Sezione contenente dati non inizializzati

alignb 32
e		resq		1

section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro


; ------------------------------------------------------------
; Funzione prova
; ------------------------------------------------------------
global prova

msg	db 'e:',32,0
nl	db 10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------
		; rdi = indirizzo della struct input

		; esempio: stampa input->e
       	; [RDI] input->seq; 			    // sequenza
		; [RDI + 8]  input->N;			    // lunghezza della sequenza
		; [RDI + 12] input->sd; 		    // tasso raffredamento
		; [RDI + 16] input->to;			    // temperatura
		; [RDI + 24] input->alpha;		    // tasso raffredamento
		; [RDI + 32] input->k; 		        // numero di features da estrarre
		; [RDI + 40] input->hydrophobicity;	// hydrophobicity
		; [RDI + 48] input->volume;		    // volume
		; [RDI + 56] input->charge;		    // charge
		; [RDI + 64] input->phi;		    // vettore angoli phi
		; [RDI + 72] input-> psi;		    // vettore angoli psi
		; [RDI + 80] input->e:			    // energy
		; [RDI + 88] input->display;
		; [RDI + 92] input->silent;

		VMOVSD		XMM0, [RDI+80]
		VMOVSD		[e], XMM0
		prints 		msg
		printsd		e
		prints 		nl
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		    ; ripristina il Base Pointer
		ret				    ; torna alla funzione C chiamante


;; ------------------------------------------------------------
;; Funzione sum_quad
;; ------------------------------------------------------------
;global sum_quad
;
;mask	dq	0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x0
;
;sum_quad:
;	; ------------------------------------------------------------
;	; Sequenza di ingresso nella funzione
;	; ------------------------------------------------------------
;		push		rbp				; salva il Base Pointer
;		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
;		pushaq						; salva i registri generali
;
;	; ------------------------------------------------------------
;	; Funzione
;	; ------------------------------------------------------------
;		VMOVAPD			XMM0, [RDI]					; XMM0 <- v[0], v[1]
;		VMOVSD			XMM1, [RDI+16]				; XMM1 <- v[2]
;		VPERM2F128		YMM0, YMM0, YMM1, 0x20		; YMM0 <- | xxx | v[2] | v[1] | v[0] |
;		VANDPD			YMM0, YMM0, [mask]			; YMM0 <- |  0  | v[2] | v[1] | v[0] |
;
;		VMULPD			YMM0, YMM0, YMM0			; YMM0 <- |    0   | v[2]^2 | v[1]^2 | v[0]^2 |
;		VPERM2F128		YMM1, YMM0, YMM0, 0x01		; YMM1 <- | v[1]^2 | v[0]^2 |    0   | v[2]^2 |
;		VADDPD			XMM0, XMM0, XMM1			; XMM0 <- | v[1]^2 | v[0]^2 + v[2]^2 |
;		VPERMILPD 		XMM1, XMM0, 0x03			; XMM1 <- | v[0]^2 + v[2]^2 | v[1]^2 |
;		VADDPD			XMM0, XMM0, XMM1			; XMM0 <- | v[0]^2 + v[1]^2 + v[2]^2 | v[0]^2 + v[1]^2 + v[2]^2 |
;		VMOVSD			[RSI], XMM0					; r = v[0]^2 + v[1]^2 + v[2]^2
;	; ------------------------------------------------------------
;	; Sequenza di uscita dalla funzione
;	; ------------------------------------------------------------
;		popaq				; ripristina i registri generali
;		mov		rsp, rbp	; ripristina lo Stack Pointer
;		pop		rbp		    ; ripristina il Base Pointer
;		ret				    ; torna alla funzione C chiamante

; ------------------------------------------------------------
; Funzione p_energy
; ------------------------------------------------------------
global p_energy

mask	dq	0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x0

p_energy:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; Funzione
	; ------------------------------------------------------------
		; RDI: *s, RSI: n, RDX: *coords, RCX: *e
		XOR 		R10, R10				; i
		XOR			R11, R11				; j
		VXORPD		XMM0, XMM0				; | energy | density |

		MOV			R8, 3*8					; R8 <- 3*8
loopi:	; i
		MOV			RAX, 3					; RAX <- 3
		IMUL		RAX, R10				; RAX <- i*3
		ADD			RAX, 1					; RAX <- i*3 + 1
		IMUL		RAX, R8					; RAX <- (i*3 + 1)*3*8

		; if (i < n)
		CMP			RAX, RSI
		JGE			endiloop

loopj:	; j
		MOV			R9, 3					; R9 <- 3
		IMUL		R9, R11					; R9 <- j*3
		ADD			R9, 1					; R9 <- j*3 + 1
		IMUL		R9, R8					; R9 <- (j*3 + 1)*3*8

		; if (j < n)
		CMP			R11, RSI
		JGE			endjloop

		; if (i != j)
		CMP			RAX, R9
		JNE			loopj

		; v[i]
		VMOVAPD			XMM1, [RDX+RAX]				; XMM1 <- v[(i*3 + 1)*3], v[(i*3 + 1)*3+1]
		VMOVSD			XMM2, [RDX+RAX+16]			; XMM2 <- v[(i*3 + 1)*3+2]
		VPERM2F128		YMM1, YMM1, YMM2, 0x20		; YMM1 <- | xxx | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		VANDPD			YMM1, YMM1, [mask]			; YMM1 <- |  0  | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		; v[i]
		VMOVAPD			XMM2, [RDX+R9]				; XMM2 <- v[(j*3 + 1)*3], v[(j*3 + 1)*3+1]
		VMOVSD			XMM3, [RDX+R9+16]			; XMM3 <- v[(j*3 + 1)*3+2]
		VPERM2F128		YMM2, YMM2, YMM3, 0x20		; YMM2 <- | xxx | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
		VANDPD			YMM2, YMM2, [mask]			; YMM2 <- |  0  | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |

		VSUBPD			YMM1, YMM2, YMM1			; YMM1 <- | v[(j*3 + 1)*3] - v[(i*3 + 1)*3] | v[(j*3 + 1)*3+1] - v[(i*3 + 1)*3+1] | v[(j*3 + 1)*3+2] - v[(i*3 + 1)*3+2] |

		VMULPD			YMM1, YMM1, YMM1
		VPERM2F128		YMM2, YMM1, YMM1, 0x01
		VADDPD			XMM1, XMM1, XMM2
		VPERMILPD 		XMM2, XMM1, 0x03
		VADDPD			XMM1, XMM1, XMM2
		VSQRTSD			XMM1, XMM1, XMM1

		CMP				XMM1, 10
		JL				loopj

		; CALCOLA DENSITà #TODO

		JMP				loopj

endjloop:
		; CALCOLA ENERGIA #TODO
		JMP				loopi

endiloop:
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		    ; ripristina il Base Pointer
		ret				    ; torna alla funzione C chiamante