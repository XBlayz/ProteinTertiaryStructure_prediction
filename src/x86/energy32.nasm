%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
alignb  16
mask	dd	0xffffffff, 0xffffffff, 0xffffffff, 0x0
ten		dd	10.0
four	dd	4.0

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina
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
		shl			ecx, 4				; ecx = (i*3 +1)*4*4
		;IMUL		ecx, ecx, 12		; ecx = (i*3 + 1)*3*4 //Cooedinata
		; jx
		mov  		edx, ebx			; edx = j
		shl  		edx, 1      		; edx = j*2
		add  		edx, ebx     		; edx = j*3 (= j*2 + j)
		add  		edx, 1      		; edx = j*3 + 1
		shl			edx, 4				; edx = (j*3 +1)*4*4
		;IMUL		edx, edx, 12		; edx = (j*3 + 1)*3*4 //Cooedinata

		MOV				esi, [ebp+coords]
		;MOVUPS			XMM5, [mask]
		; v[ix]
		MOVAPS			XMM0, [esi+ecx]				; XMM0: |[(i*3 + 1)*3+3]|v[(i*3 + 1)*3+1+2]|[(i*3 + 1)*3+1]|v[(i*3 + 1)*3]|
		;ANDPS			XMM0, XMM5					; XMM0: |0|v[(i*3 + 1)*3+1+2]|v[(i*3 + 1)*3+1]|v[(i*3 + 1)*3]|
		; v[jx]
		MOVAPS			XMM1, [esi+edx]				; XMM1: |[(j*3 + 1)*3+3]|v[(j*3 + 1)*3+1+2]|[(j*3 + 1)*3+1]|v[(j*3 + 1)*3]|
		;ANDPS			XMM1, XMM5					; XMM1: |0|v[(j*3 + 1)*3+1+2]|v[(j*3 + 1)*3+1]|v[(j*3 + 1)*3]|

		; Distanza
		dist			XMM0, XMM1

		MOVSS			XMM1, [ten]
		; if (d < 10.0)
		COMISS			XMM0, XMM1
		JAE				p_loop_j

		; Calcolo densità
		MOVSS			XMM1, XMM0					; XMM1 = d
		MULSS			XMM0, XMM0					; XMM0 = d^2
		MULSS			XMM0, XMM1					; XMM0 = d^3
		XOR				ecx, ecx					; ecx = 0
		MOV				edx, [ebp+s]				; edx = s
		MOV				cl, [edx+ebx]				; CL = s[j]
		SUB				cl, 65						; CL = s[j]-'A'
		MOV				edx, [ebp+volume]			; edx = volume
		MOVSS			XMM1, [edx+ecx*4]			; XMM1 = volume[s[j]-'A']
		DIVSS			XMM1, XMM0					; XMM1 = volume[s[j]-'A'] / d^3
		ADDSS			XMM6, XMM1					; Density += volume[s[j]-'a'] / d^3

p_loop_j:
		INC 	ebx									; j++
		JMP		p_start_loop_j

p_end_j_loop:
		; Calcolo energia
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



; ------------------------------------------------------------
; Funzione e_energy
; ------------------------------------------------------------
global e_energy

charge	equ		24		; float*

e_energy:
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
		MOV				edi, [ebp+n]			; edi: n

		XORPS			XMM7, XMM7				; XMM7: |0|0|0|Energy=0|
		XOR 			eax, eax				; eax: i=0

		MOVSS			XMM6, [ten]				; XMM6: 10.0
e_start_loop_i:
		; if (i < n)
		CMP				eax, edi
		JGE				e_end_i_loop

		MOV				ebx, eax				; ebx: j=i
		INC				ebx						; ebx: j=i+1
e_start_loop_j:
		; if (j < n)
		CMP				ebx, edi
		JGE				e_end_j_loop
		; if (i != j)
		CMP				eax, ebx
		JE				e_j_loop

		MOV				edx, [ebp+s]				; edx = s
		; charge[s[i]-'A']
		XOR				ecx, ecx					; ecx = 0
		MOV				cl, [edx+eax]				; cl = s[i]
		SUB				cl, 65						; cl = s[i]-'A'
		MOV				edx, [ebp+charge]			; edx = charge
		MOVSS			XMM0, [edx+ecx*4]			; XMM0 = charge[s[i]-'A']
		; if (charge[s[i]-'A'] != 0)
		XORPS			XMM1, XMM1
		COMISS			XMM0, XMM1
		JE				e_j_loop

		MOV				edx, [ebp+s]				; edx = s
		; charge[s[j]-'A']
		XOR				ecx, ecx					; ecx = 0
		MOV				cl, [edx+ebx]				; cl = s[j]
		SUB				cl, 65						; cl = s[j]-'A'
		MOV				edx, [ebp+charge]			; edx = charge
		MOVSS			XMM2, [edx+ecx*4]			; XMM2 = charge[s[j]-'A']
		; if (charge[s[j]-'A'] != 0)
		XORPS			XMM1, XMM1
		COMISS			XMM2, XMM1
		JE				e_j_loop

		; ix
		mov  		ecx, eax			; ecx = i
		shl  		ecx, 1      		; ecx = i*2
		add  		ecx, eax     		; ecx = i*3 (= i*2 + i)
		add  		ecx, 1      		; ecx = i*3 + 1
		shl			ecx, 4				; ecx = (i*3 +1)*4*4
		;IMUL		ecx, ecx, 12		; ecx = (i*3 + 1)*3*4 //Cooedinata
		; jx
		mov  		edx, ebx			; edx = j
		shl  		edx, 1      		; edx = j*2
		add  		edx, ebx     		; edx = j*3 (= j*2 + j)
		add  		edx, 1      		; edx = j*3 + 1
		shl			edx, 4				; edx = (j*3 +1)*4*4
		;IMUL		edx, edx, 12		; edx = (j*3 + 1)*3*4 //Cooedinata

		MOV				esi, [ebp+coords]
		;MOVUPS			XMM5, [mask]
		; v[ix]
		MOVAPS			XMM3, [esi+ecx]				; XMM3: |[(i*3 + 1)*3+3]|v[(i*3 + 1)*3+1+2]|[(i*3 + 1)*3+1]|v[(i*3 + 1)*3]|
		;ANDPS			XMM3, XMM5					; XMM3: |0|v[(i*3 + 1)*3+1+2]|v[(i*3 + 1)*3+1]|v[(i*3 + 1)*3]|
		; v[jx]
		MOVAPS			XMM4, [esi+edx]				; XMM4: |[(j*3 + 1)*3+3]|v[(j*3 + 1)*3+1+2]|[(j*3 + 1)*3+1]|v[(j*3 + 1)*3]|
		;ANDPS			XMM4, XMM5					; XMM4: |0|v[(j*3 + 1)*3+1+2]|v[(j*3 + 1)*3+1]|v[(j*3 + 1)*3]|

		; Distance
		dist			XMM3, XMM4

		; if (d < 10.0)
		COMISS			XMM3, XMM6
		JAE				e_j_loop

		; Calcolo energia
		MULSS			XMM0, XMM2					; XMM0 = charge[s[i]-'a'] * charge[s[j]-'a']
		MULSS			XMM3, [four]				; XMM3 = d * 4.0
		DIVSS			XMM0, XMM3					; XMM0 = (charge[s[i]-'a'] * charge[s[j]-'a']) / (d * 4.0)
		ADDSS			XMM7, XMM0					; Energy += (charge[s[i]-'a'] * charge[s[j]-'a']) / (d * 4.0)

e_j_loop:
		INC			ebx								; j++
		JMP			e_start_loop_j

e_end_j_loop:
		INC 		eax								; i++
		JMP			e_start_loop_i

e_end_i_loop:
		MOV			eax, [ebp+energy]
		MOVSS		[eax], XMM7						; e = Energy
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		pop	edi			; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp			; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante


; ------------------------------------------------------------
; Funzione h_energy
; ------------------------------------------------------------
global h_energy

hydrophobicity	equ		24		; float*

h_energy:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi		; salva i registri generali

	; ------------------------------------------------------------
	; Funzione
	; ------------------------------------------------------------
		MOV				edi, [ebp+n]			; edi: n

		XORPS			XMM7, XMM7				; XMM7: |0|0|0|Energy=0|
		XOR 			eax, eax				; eax: i=0

		MOVSS			XMM6, [ten]				; XMM6: 10.0
h_start_loop_i:
		; if (i < n)
		CMP				eax, edi
		JGE				h_end_i_loop

		MOV				ebx, eax				; ebx: j=i
		INC				ebx						; ebx: j=i+1
h_start_loop_j:
		; if (j < n)
		CMP				ebx, edi
		JGE				h_end_j_loop

		MOV				edx, [ebp+s]				; edx = s
		; hydrophobicity[s[i]-'A']
		XOR				ecx, ecx					; ecx = 0
		MOV				cl, [edx+eax]				; cl = s[i]
		SUB				cl, 65						; cl = s[i]-'A'
		MOV				edx, [ebp+hydrophobicity]	; edx = hydrophobicity
		MOVSS			XMM0, [edx+ecx*4]			; XMM0 = hydrophobicity[s[i]-'A']

		MOV				edx, [ebp+s]				; edx = s
		; hydrophobicity[s[j]-'A']
		XOR				ecx, ecx					; ecx = 0
		MOV				cl, [edx+ebx]				; cl = s[j]
		SUB				cl, 65						; cl = s[j]-'A'
		MOV				edx, [ebp+hydrophobicity]	; edx = hydrophobicity
		MOVSS			XMM1, [edx+ecx*4]			; XMM1 = hydrophobicity[s[i]-'A']

		; ix
		mov  		ecx, eax			; ecx = i
		shl  		ecx, 1      		; ecx = i*2
		add  		ecx, eax     		; ecx = i*3 (= i*2 + i)
		add  		ecx, 1      		; ecx = i*3 + 1
		shl			ecx, 4				; ecx = (i*3 +1)*4*4
		;IMUL		ecx, ecx, 12		; ecx = (i*3 + 1)*3*4 //Cooedinata
		; jx
		mov  		edx, ebx			; edx = j
		shl  		edx, 1      		; edx = j*2
		add  		edx, ebx     		; edx = j*3 (= j*2 + j)
		add  		edx, 1      		; edx = j*3 + 1
		shl			edx, 4				; edx = (j*3 +1)*4*4
		;IMUL		edx, edx, 12		; edx = (j*3 + 1)*3*4 //Cooedinata

		MOV				esi, [ebp+coords]
		;MOVUPS			XMM5, [mask]
		; v[ix]
		MOVAPS			XMM2, [esi+ecx]				; XMM2: |[(i*3 + 1)*3+3]|v[(i*3 + 1)*3+1+2]|[(i*3 + 1)*3+1]|v[(i*3 + 1)*3]|
		;ANDPS			XMM2, XMM5					; XMM2: |0|v[(i*3 + 1)*3+1+2]|v[(i*3 + 1)*3+1]|v[(i*3 + 1)*3]|
		; v[jx]
		MOVAPS			XMM3, [esi+edx]				; XMM3: |[(j*3 + 1)*3+3]|v[(j*3 + 1)*3+1+2]|[(j*3 + 1)*3+1]|v[(j*3 + 1)*3]|
		;ANDPS			XMM3, XMM5					; XMM3: |0|v[(j*3 + 1)*3+1+2]|v[(j*3 + 1)*3+1]|v[(j*3 + 1)*3]|

		; Distance
		dist			XMM2, XMM3

		; if (d < 10.0)
		COMISS			XMM2, XMM6
		JAE				h_j_loop

		; Calcolo energia
		MULSS			XMM1, XMM0					; XMM1 = hydrophobicity[s[i]-'a'] * hydrophobicity[s[j]-'a']
		DIVSS			XMM1, XMM2					; XMM1 = (hydrophobicity[s[i]-'a'] * hydrophobicity[s[j]-'a']) / d
		ADDSS			XMM7, XMM1					; Energy += (hydrophobicity[s[i]-'a'] * hydrophobicity[s[j]-'a']) / d

h_j_loop:
		INC			ebx						; j++
		JMP			h_start_loop_j

h_end_j_loop:
		INC 		eax						; i++
		JMP			h_start_loop_i

h_end_i_loop:
		MOV			eax, [ebp+energy]
		MOVSS		[eax], XMM7						; e = Energy
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		pop	edi			; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp			; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante