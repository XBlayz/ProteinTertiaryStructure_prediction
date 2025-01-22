%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
alignb  4
mask	dd	0xffffffff, 0xffffffff, 0xffffffff, 0x0
ten		dd	10.0
four	dd	4.0

nl	db	10, 0
sep	db	'|', 0

section .bss			; Sezione contenente dati non inizializzati
__print_xmm_temp__	resd	4

section .text			; Sezione contenente il codice macchina
; ------------------------------------------------------------
; Macro debug
; ------------------------------------------------------------
%macro print_int 1
	printi	%1
	prints	nl
%endmacro

%macro print_float 1
	sprint	%1
	prints	nl
%endmacro

%macro print_xmm 1
	MOVUPD		[__print_xmm_temp__], %1
	prints		sep
	sprint		__print_xmm_temp__+12
	prints		sep
	sprint		__print_xmm_temp__+8
	prints		sep
	sprint		__print_xmm_temp__+4
	prints		sep
	sprint		__print_xmm_temp__
	prints		sep
	prints		nl
%endmacro

; ------------------------------------------------------------
; Macro utility
; ------------------------------------------------------------
; -Formula distanza-
%macro dist 2
	; %1: Va(XMMa), %2: Vb(XMMb)
    ; Differenza
    SUBPS			%1, %2          	; %1: |0|Az-Bz|Ay-By|Ax-Bx|
    ; Quadrati
    MULPS			%1, %1          	; %1: |0|(Az-Bz)^2|(Ay-By)^2|(Ax-Bx)^2|

    ; Somma
    HADDPS         %1, %1            	; %1: |(Az-Bz)^2|(Az-Bz)^2|(Ay-By)^2 + (Ax-Bx)^2|(Ay-By)^2 + (Ax-Bx)^2|
    HADDPS         %1, %1            	; %1: |---|---|---|(Az-Bz)^2 + (Ay-By)^2 + (Ax-Bx)^2|

    ; Radice quadrata
    SQRTSS         %1, %1               ; %1: |---|sqrt((Az-Bz)^2 + (Ay-By)^2 + (Ax-Bx)^2)|
%endmacro

; ------------------------------------------------------------
; Funzione p_energy
; ------------------------------------------------------------
global p_energy

s		equ		8		; char*
n  		equ		12		; int
coords	equ		16		; float*
energy	equ		20		; float*
volume	equ		24		; float*

p_energy:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi

	; ------------------------------------------------------------
	; Funzione
	; ------------------------------------------------------------
		MOV			edi, [ebp+n]			; edi: n

		XORPS		XMM7, XMM7				; XMM7: |0|0|0|Energy=0|
		XOR 		eax, eax				; eax: i=0
p_start_loop_i:
		; if (i < n)
		CMP			eax, edi
		JGE			p_end_i_loop

		XORPS		XMM6, XMM6				; XMM6: |0|0|0|Density=0|
		XOR			ebx, ebx				; ebx: j=0
p_start_loop_j:
		; if (j < n)
		CMP			ebx, edi
		JGE			p_end_j_loop
		; if (i != j)
		CMP			eax, ebx
		JE			p_loop_j

		; ix
		mov  		ecx, eax			; ecx = i
		shl  		ecx, 1      		; ecx = i*2
		add  		ecx, eax     		; ecx = i*3 (= i*2 + i)
		add  		ecx, 1      		; ecx = i*3 + 1
		IMUL		ecx, ecx, 12		; ecx = (i*3 + 1)*3*4 //Cooedinata
		; jx
		mov  		edx, ebx			; edx = j
		shl  		edx, 1      		; edx = j*2
		add  		edx, ebx     		; edx = j*3 (= j*2 + j)
		add  		edx, 1      		; edx = j*3 + 1
		IMUL		edx, edx, 12		; edx = (j*3 + 1)*3*4 //Cooedinata

		MOV				esi, [ebp+coords]
		MOVUPS			XMM5, [mask]
		; v[ix]
		MOVUPS			XMM0, [esi+ecx]				; XMM0: |[(i*3 + 1)*3+3]|v[(i*3 + 1)*3+1+2]|[(i*3 + 1)*3+1]|v[(i*3 + 1)*3]|
		ANDPS			XMM0, XMM5					; XMM0: |0|v[(i*3 + 1)*3+1+2]|v[(i*3 + 1)*3+1]|v[(i*3 + 1)*3]|
		; v[jx]
		MOVUPS			XMM1, [esi+edx]				; XMM1: |[(j*3 + 1)*3+3]|v[(j*3 + 1)*3+1+2]|[(j*3 + 1)*3+1]|v[(j*3 + 1)*3]|
		ANDPS			XMM1, XMM5					; XMM1: |0|v[(j*3 + 1)*3+1+2]|v[(j*3 + 1)*3+1]|v[(j*3 + 1)*3]|

		; Distanza
		dist			XMM0, XMM1

		MOVSS			XMM1, [ten]
		; if (d > 10.0)
		COMISS			XMM0, XMM1
		JBE				p_loop_j

		; Calcolo densità
		MOVSS			XMM1, XMM0					; XMM1 = d
		MULSS			XMM0, XMM0					; XMM0 = d^2
		MULSS			XMM0, XMM1					; XMM0 = d^3
		XOR				ecx, ecx					; ecx = 0
		MOV				edx, [ebp+s]				; edx = s
		MOV				cl, [edx+ebx]				; CL = s[j]
		SUB				cl, 65						; CL = s[j]-'A'
		MOV				edx, [ebp+volume]			; edx = volume
		MOVSS			XMM1, [edx+ecx*4]			; XMM1 = volume[s[j]-'a']
		DIVSS			XMM1, XMM0					; XMM1 = volume[s[j]-'a'] / d^3
		ADDSS			XMM6, XMM1					; Density += volume[s[j]-'a'] / d^3

p_loop_j:
		INC 	ebx									; j++
		JMP		p_start_loop_j

p_end_j_loop:
		XOR				ecx, ecx					; ecx = 0
		MOV				edx, [ebp+s]				; edx = s
		MOV				cl, [edx+eax]				; CL = s[i]
		SUB				cl, 65						; CL = s[i]-'A'
		MOV				edx, [ebp+volume]			; edx = volume
		MOVSS			XMM1, [edx+ecx*4]			; XMM1 = volume[s[j]-'A']
		SUBSS			XMM1, XMM6					; XMM1 = volume[s[j]-'A'] - Density
		MULSS			XMM1, XMM1					; XMM1 = (volume[s[j]-'A'] - Density)^2
		ADDSS			XMM7, XMM1					; Energy += (volume[s[j]-'A'] - Density)^2

p_loop_i:
		INC 	eax									; i++
		JMP		p_start_loop_i

p_end_i_loop:
		MOV				eax, [ebp+energy]
		MOVSS			[eax], XMM7			; e = Energy

	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		pop	edi			; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp			; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante



;; ------------------------------------------------------------
;; Funzione e_energy
;; ------------------------------------------------------------
;global e_energy
;
;four	dq	4.0
;
;e_energy:
;	; ------------------------------------------------------------
;	; Sequenza di ingresso nella funzione
;	; ------------------------------------------------------------
;		push		ebp		; salva il Base Pointer
;		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
;		push		ebx		; salva i registri da preservare
;		push		esi
;		push		edi					; salva i registri generali
;
;	; ------------------------------------------------------------
;	; Funzione
;	; ------------------------------------------------------------
;		; edi: *s, esi: n, RDX: *coords, ECX: *e, R8: *charge
;		MOV			[eEPointer], ECX		; ePointer = *e
;		XOR 		R9, R9					; R9: i=0
;		VXORPD		XMM0, XMM0				; XMM0: | 0 | energy=0 |
;
;		MOV			R11, 3*8				; R11 = 3*8
;		VMOVSD		XMM4, [ten]				; XMM4 = 10.0
;e_start_loop_i:
;		; if (i < n)
;		CMP			R9, esi
;		JGE			e_end_i_loop
;
;		MOV			R10, R9				; R10: j=i
;		INC			R10					; R10: j=i+1
;e_start_loop_j:
;		; if (j < n)
;		CMP			R10, esi
;		JGE			e_end_j_loop
;		; if (i != j)
;		CMP			R9, R10
;		JE			e_j_loop
;
;		; charge[s[i]-'a']
;		XOR			eax, eax				; eax = 0
;		MOV			AL, [edi+R9]			; AL = s[i]
;		SUB			AL, 65					; AL = s[i]-'a'
;		VMOVSD		XMM6, [R8+eax*8]		; XMM6 = charge[s[i]-'a']
;		; if (charge[s[i]-'a'] != 0)
;		VXORPD			XMM2, XMM2
;		VCMPSD			XMM5, XMM6, XMM2, 0
;		MOVMSKPD		ECX, XMM5
;		CMP				ECX, 0
;		JG				e_j_loop
;
;		; charge[s[j]-'a']
;		XOR			eax, eax				; eax = 0
;		MOV			AL, [edi+R10]			; AL = s[j]
;		SUB			AL, 65					; AL = s[j]-'a'
;		VMOVSD		XMM7, [R8+eax*8]		; XMM7 = charge[s[i]-'a']
;		; if (charge[s[j]-'a'] != 0)
;		VXORPD			XMM2, XMM2
;		VCMPSD			XMM5, XMM7, XMM2, 0
;		MOVMSKPD		ECX, XMM5
;		CMP				ECX, 0
;		JG				e_j_loop
;
;		; ix
;		MOV			eax, 3					; eax = 3
;		IMUL		eax, R9					; eax = i*3
;		ADD			eax, 1					; eax = i*3 + 1 //Atomo Ca
;		IMUL		eax, R11				; eax = (i*3 + 1)*3*8 //Cooedinata
;		; jx
;		MOV			ECX, 3					; ECX = 3
;		IMUL		ECX, R10				; ECX = j*3
;		ADD			ECX, 1					; ECX = j*3 + 1 //Atomo Ca
;		IMUL		ECX, R11				; ECX = (j*3 + 1)*3*8 //Cooedinata
;
;		; v[ix]
;		VMOVUPD			XMM1, [RDX+eax]				; XMM1 <- v[(i*3 + 1)*3], v[(i*3 + 1)*3+1]
;		VMOVSD			XMM2, [RDX+eax+16]			; XMM2 <- v[(i*3 + 1)*3+2]
;		VPERM2F128		YMM1, YMM1, YMM2, 0x20		; YMM1 <- | xxx | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
;		VANDPD			YMM1, YMM1, [mask]			; YMM1 <- |  0  | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
;		; v[jx]
;		VMOVUPD			XMM2, [RDX+ECX]				; XMM2 <- v[(j*3 + 1)*3], v[(j*3 + 1)*3+1]
;		VMOVSD			XMM3, [RDX+ECX+16]			; XMM3 <- v[(j*3 + 1)*3+2]
;		VPERM2F128		YMM2, YMM2, YMM3, 0x20		; YMM2 <- | xxx | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
;		VANDPD			YMM2, YMM2, [mask]			; YMM2 <- |  0  | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
;
;		; Distance
;		dist			YMM1, YMM2, XMM1, XMM2
;
;		; if (d < 10.0)
;		VCMPSD			XMM5, XMM1, XMM4, 5
;		MOVMSKPD		ECX, XMM5
;		CMP				ECX, 0
;		JG				e_j_loop
;
;		; ENERGIA
;		VMULSD			XMM6, XMM6, XMM7			; XMM6 = charge[s[i]-'a'] * charge[s[j]-'a']
;		VMULSD			XMM7, XMM1, [four]			; XMM7 = d * 4.0
;		VDIVSD			XMM6, XMM6, XMM7			; XMM6 = XMM6 / XMM7
;		VADDSD			XMM0, XMM0, XMM6			; XMM0 += XMM6
;
;e_j_loop:
;		INC			R10						; j++
;		JMP			e_start_loop_j
;
;e_end_j_loop:
;		INC 		R9						; i++
;		JMP			e_start_loop_i
;
;e_end_i_loop:
;		MOV			ECX, [eEPointer]
;		VMOVSD		[ECX], XMM0			; e = XMM7
;	; ------------------------------------------------------------
;	; Sequenza di uscita dalla funzione
;	; ------------------------------------------------------------
;		pop	edi			; ripristina i registri da preservare
;		pop	esi
;		pop	ebx
;		mov	esp, ebp	; ripristina lo Stack Pointer
;		pop	ebp			; ripristina il Base Pointer
;		ret				; torna alla funzione C chiamante
;
;
;; ------------------------------------------------------------
;; Funzione h_energy
;; ------------------------------------------------------------
;global h_energy
;
;h_energy:
;	; ------------------------------------------------------------
;	; Sequenza di ingresso nella funzione
;	; ------------------------------------------------------------
;		push		ebp		; salva il Base Pointer
;		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
;		push		ebx		; salva i registri da preservare
;		push		esi
;		push		edi		; salva i registri generali
;
;	; ------------------------------------------------------------
;	; Funzione
;	; ------------------------------------------------------------
;		; edi: *s, esi: n, RDX: *coords, ECX: *e, R8: *hydrophobicity
;		MOV			[hEPointer], ECX		; ePointer = *e
;		XOR 		R9, R9					; R9: i=0
;		VXORPD		XMM0, XMM0				; XMM0: | 0 | energy=0 |
;
;		MOV			R11, 3*8				; R11 = 3*8
;		VMOVSD		XMM4, [ten]				; XMM4 = 10.0
;h_start_loop_i:
;		; if (i < n)
;		CMP			R9, esi
;		JGE			h_end_i_loop
;
;		MOV			R10, R9				; R10: j=i
;		INC			R10					; R10: j=i+1
;h_start_loop_j:
;		; if (j < n)
;		CMP			R10, esi
;		JGE			h_end_j_loop
;
;		; hydrophobicity[s[i]-'a']
;		XOR			eax, eax				; eax = 0
;		MOV			AL, [edi+R9]			; AL = s[i]
;		SUB			AL, 65					; AL = s[i]-'a'
;		VMOVSD		XMM6, [R8+eax*8]		; XMM6 = hydrophobicity[s[i]-'a']
;
;		; hydrophobicity[s[j]-'a']
;		XOR			eax, eax				; eax = 0
;		MOV			AL, [edi+R10]			; AL = s[j]
;		SUB			AL, 65					; AL = s[j]-'a'
;		VMOVSD		XMM7, [R8+eax*8]		; XMM7 = hydrophobicity[s[i]-'a']
;
;		; ix
;		MOV			eax, 3					; eax = 3
;		IMUL		eax, R9					; eax = i*3
;		ADD			eax, 1					; eax = i*3 + 1 //Atomo Ca
;		IMUL		eax, R11				; eax = (i*3 + 1)*3*8 //Cooedinata
;		; jx
;		MOV			ECX, 3					; ECX = 3
;		IMUL		ECX, R10				; ECX = j*3
;		ADD			ECX, 1					; ECX = j*3 + 1 //Atomo Ca
;		IMUL		ECX, R11				; ECX = (j*3 + 1)*3*8 //Cooedinata
;
;		; v[ix]
;		VMOVUPD			XMM1, [RDX+eax]				; XMM1 <- v[(i*3 + 1)*3], v[(i*3 + 1)*3+1]
;		VMOVSD			XMM2, [RDX+eax+16]			; XMM2 <- v[(i*3 + 1)*3+2]
;		VPERM2F128		YMM1, YMM1, YMM2, 0x20		; YMM1 <- | xxx | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
;		VANDPD			YMM1, YMM1, [mask]			; YMM1 <- |  0  | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
;		; v[jx]
;		VMOVUPD			XMM2, [RDX+ECX]				; XMM2 <- v[(j*3 + 1)*3], v[(j*3 + 1)*3+1]
;		VMOVSD			XMM3, [RDX+ECX+16]			; XMM3 <- v[(j*3 + 1)*3+2]
;		VPERM2F128		YMM2, YMM2, YMM3, 0x20		; YMM2 <- | xxx | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
;		VANDPD			YMM2, YMM2, [mask]			; YMM2 <- |  0  | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
;
;		; Distance
;		dist			YMM1, YMM2, XMM1, XMM2
;
;		; if (d < 10.0)
;		VCMPSD			XMM5, XMM1, XMM4, 5
;		MOVMSKPD		ECX, XMM5
;		CMP				ECX, 0
;		JG				h_j_loop
;
;		; ENERGIA
;		VMULSD			XMM6, XMM6, XMM7			; XMM6 = hydrophobicity[s[i]-'a'] * hydrophobicity[s[j]-'a']
;		VDIVSD			XMM6, XMM6, XMM1			; XMM6 = XMM6 / d
;		VADDSD			XMM0, XMM0, XMM6			; XMM0 += XMM6
;
;h_j_loop:
;		INC			R10						; j++
;		JMP			h_start_loop_j
;
;h_end_j_loop:
;		INC 		R9						; i++
;		JMP			h_start_loop_i
;
;h_end_i_loop:
;		MOV			ECX, [hEPointer]
;		VMOVSD		[ECX], XMM0			; e = XMM7
;	; ------------------------------------------------------------
;	; Sequenza di uscita dalla funzione
;	; ------------------------------------------------------------
;		pop	edi			; ripristina i registri da preservare
;		pop	esi
;		pop	ebx
;		mov	esp, ebp	; ripristina lo Stack Pointer
;		pop	ebp			; ripristina il Base Pointer
;		ret				; torna alla funzione C chiamante