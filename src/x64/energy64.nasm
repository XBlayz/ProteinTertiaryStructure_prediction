%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
alignb  32
mask	dq	0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x0
ten		dq	10.0
four	dq	4.0

section .bss			; Sezione contenente dati non inizializzati
pEPointer	resq    1
eEPointer	resq    1
hEPointer	resq    1

section .text			; Sezione contenente il codice macchina
; ------------------------------------------------------------
; Macro utility
; ------------------------------------------------------------
; -Formula distanza-
%macro dist 4
	; %1: Va(YMMa), %2: Vb(YMMb), %3: (Vr)XMMa, %4: (-)XMMb
    ; Differenza
    VSUBPD			%1, %2, %1          	; %1: |0|Az-Bz|Ay-By|Ax-Bx|
    ; Quadrati
    VMULPD			%1, %1, %1          	; %1: |0|(Az-Bz)^2|(Ay-By)^2|(Ax-Bx)^2|

    ; Somma
    VHADDPD         %1, %1, %1            	; %1: |(Az-Bz)^2|(Az-Bz)^2|(Ay-By)^2 + (Ax-Bx)^2|(Ay-By)^2 + (Ax-Bx)^2|
    VEXTRACTF128    %4, %1, 0b00000001      ; %4: |(Az-Bz)^2|(Az-Bz)^2|
    VADDSD          %3, %3, %4            	; %3: |---|(Az-Bz)^2 + (Ay-By)^2 + (Ax-Bx)^2|

    ; Radice quadrata
    VSQRTSD         %3, %3                  ; %3: |---|sqrt((Az-Bz)^2 + (Ay-By)^2 + (Ax-Bx)^2)|
%endmacro

; ------------------------------------------------------------
; Funzione p_energy
; ------------------------------------------------------------
global p_energy

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
		; RDI: *s, RSI: n, RDX: *coords, RCX: *e, R8: *volume
		MOV			[pEPointer], RCX		; ePointer = *e
		MOV			R11, R8					; R11: *volume
		XOR 		R8, R8					; R8: i=0
		VXORPD		XMM7, XMM7				; XMM7: | 0 | energy=0 |

		MOV			R10, 4*8				; R10 = 4*8
		VMOVSD		XMM4, [ten]
start_loop_i:
		; if (i < n)
		CMP			R8, RSI
		JGE			end_i_loop

		VXORPD		XMM0, XMM0				; XMM0: | 0 | density=0 |
		XOR			R9, R9					; R9: j=0
start_loop_j:
		; if (j < n)
		CMP			R9, RSI
		JGE			end_j_loop
		; if (i != j)
		CMP			R8, R9
		JE			loop_j

		; ix
		MOV			RAX, 3					; RAX = 3
		IMUL		RAX, R8					; RAX = i*3
		ADD			RAX, 1					; RAX = i*3 + 1 //Atomo Ca
		IMUL		RAX, R10				; RAX = (i*3 + 1)*4*8 //Coordinata
		; jx
		MOV			RCX, 3					; RCX = 3
		IMUL		RCX, R9					; RCX = j*3
		ADD			RCX, 1					; RCX = j*3 + 1 //Atomo Ca
		IMUL		RCX, R10				; RCX = (j*3 + 1)*4*8 //Coordinata

		; v[ix]
		VMOVAPD			YMM1, [RDX+RAX]				; YMM1 <- | 0 | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		;VANDPD			YMM1, YMM1, [mask]			; YMM1 <- | 0 | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		; v[jx]
		VMOVAPD			YMM2, [RDX+RCX]				; YMM2 <- | 0 | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
		;VANDPD			YMM2, YMM2, [mask]			; YMM2 <- | 0 | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |

		; Distanza
		dist			YMM1, YMM2, XMM1, XMM2

		; if (d < 10.0)
		VCMPSD			XMM5, XMM1, XMM4, 5
		MOVMSKPD		RCX, XMM5
		CMP				RCX, 0
		JG				loop_j

		; CALCOLA DENSITà
		VMULSD			XMM2, XMM1, XMM1			; XMM2 = d^2
		VMULSD			XMM2, XMM2, XMM1			; XMM2 = d^3
		XOR				RAX, RAX					; RAX = 0
		MOV				AL, [RDI+R9]				; AL = s[j]
		SUB				AL, 65						; AL = s[j]-'a'
		VMOVSD			XMM3, [R11+RAX*8]			; XMM3 = volume[s[j]-'a']
		VDIVSD			XMM3, XMM3, XMM2			; XMM3 = volume[s[j]-'a'] / d^3
		VADDSD			XMM0, XMM0, XMM3			; XMM0 += volume[s[j]-'a'] / d^3

loop_j:
		INC 	R9				; j++
		JMP		start_loop_j

end_j_loop:
		XOR				RAX, RAX					; RAX = 0
		MOV				AL, [RDI+R8]				; AL = s[i]
		SUB				AL, 65						; AL = s[i]-'a'
		VMOVSD			XMM3, [R11+RAX*8]			; XMM3 = volume[s[i]-'a']
		VSUBSD			XMM3, XMM3, XMM0			; XMM3 = volume[s[i]-'a'] - XMM0
		VMULSD			XMM3, XMM3, XMM3			; XMM3 = (volume[s[i]-'a'] - XMM0)^2
		VADDSD			XMM7, XMM7, XMM3			; XMM7 += (volume[s[i]-'a'] - XMM0)^2
loop_i:
		INC 	R8				; i++
		JMP		start_loop_i

end_i_loop:
		MOV				RCX, [pEPointer]
		VMOVSD			[RCX], XMM7					; e = XMM7
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		    ; ripristina il Base Pointer
		ret				    ; torna alla funzione C chiamante


; ------------------------------------------------------------
; Funzione e_energy
; ------------------------------------------------------------
global e_energy

e_energy:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; Funzione
	; ------------------------------------------------------------
		; RDI: *s, RSI: n, RDX: *coords, RCX: *e, R8: *charge
		MOV			[eEPointer], RCX		; ePointer = *e
		XOR 		R9, R9					; R9: i=0
		VXORPD		XMM0, XMM0				; XMM0: | 0 | energy=0 |

		MOV			R11, 4*8				; R11 = 4*8
		VMOVSD		XMM4, [ten]				; XMM4 = 10.0
e_start_loop_i:
		; if (i < n)
		CMP			R9, RSI
		JGE			e_end_i_loop

		MOV			R10, R9				; R10: j=i
		INC			R10					; R10: j=i+1
e_start_loop_j:
		; if (j < n)
		CMP			R10, RSI
		JGE			e_end_j_loop
		; if (i != j)
		CMP			R9, R10
		JE			e_j_loop

		; charge[s[i]-'a']
		XOR			RAX, RAX				; RAX = 0
		MOV			AL, [RDI+R9]			; AL = s[i]
		SUB			AL, 65					; AL = s[i]-'a'
		VMOVSD		XMM6, [R8+RAX*8]		; XMM6 = charge[s[i]-'a']
		; if (charge[s[i]-'a'] != 0)
		VXORPD			XMM2, XMM2
		VCMPSD			XMM5, XMM6, XMM2, 0
		MOVMSKPD		RCX, XMM5
		CMP				RCX, 0
		JG				e_j_loop

		; charge[s[j]-'a']
		XOR			RAX, RAX				; RAX = 0
		MOV			AL, [RDI+R10]			; AL = s[j]
		SUB			AL, 65					; AL = s[j]-'a'
		VMOVSD		XMM7, [R8+RAX*8]		; XMM7 = charge[s[i]-'a']
		; if (charge[s[j]-'a'] != 0)
		VXORPD			XMM2, XMM2
		VCMPSD			XMM5, XMM7, XMM2, 0
		MOVMSKPD		RCX, XMM5
		CMP				RCX, 0
		JG				e_j_loop

		; ix
		MOV			RAX, 3					; RAX = 3
		IMUL		RAX, R9					; RAX = i*3
		ADD			RAX, 1					; RAX = i*3 + 1 //Atomo Ca
		IMUL		RAX, R11				; RAX = (i*3 + 1)*4*8 //Coordinata
		; jx
		MOV			RCX, 3					; RCX = 3
		IMUL		RCX, R10				; RCX = j*3
		ADD			RCX, 1					; RCX = j*3 + 1 //Atomo Ca
		IMUL		RCX, R11				; RCX = (j*3 + 1)*4*8 //Coordinata

		; v[ix]
		VMOVAPD			YMM1, [RDX+RAX]				; YMM1 <- | 0 |v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		;VANDPD			YMM1, YMM1, [mask]			; YMM1 <- | 0 |v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		; v[jx]
		VMOVAPD			YMM2, [RDX+RCX]				; YMM2 <- | 0 |v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
		;VANDPD			YMM2, YMM2, [mask]			; YMM2 <- | 0 | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |

		; Distance
		dist			YMM1, YMM2, XMM1, XMM2

		; if (d < 10.0)
		VCMPSD			XMM5, XMM1, XMM4, 5
		MOVMSKPD		RCX, XMM5
		CMP				RCX, 0
		JG				e_j_loop

		; ENERGIA
		VMULSD			XMM6, XMM6, XMM7			; XMM6 = charge[s[i]-'a'] * charge[s[j]-'a']
		VMULSD			XMM7, XMM1, [four]			; XMM7 = d * 4.0
		VDIVSD			XMM6, XMM6, XMM7			; XMM6 = XMM6 / XMM7
		VADDSD			XMM0, XMM0, XMM6			; XMM0 += XMM6

e_j_loop:
		INC			R10						; j++
		JMP			e_start_loop_j

e_end_j_loop:
		INC 		R9						; i++
		JMP			e_start_loop_i

e_end_i_loop:
		MOV			RCX, [eEPointer]
		VMOVSD		[RCX], XMM0			; e = XMM7
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		    ; ripristina il Base Pointer
		ret				    ; torna alla funzione C chiamante


; ------------------------------------------------------------
; Funzione h_energy
; ------------------------------------------------------------
global h_energy

h_energy:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; Funzione
	; ------------------------------------------------------------
		; RDI: *s, RSI: n, RDX: *coords, RCX: *e, R8: *hydrophobicity
		MOV			[hEPointer], RCX		; ePointer = *e
		XOR 		R9, R9					; R9: i=0
		VXORPD		XMM0, XMM0				; XMM0: | 0 | energy=0 |

		MOV			R11, 4*8				; R11 = 4*8
		VMOVSD		XMM4, [ten]				; XMM4 = 10.0
h_start_loop_i:
		; if (i < n)
		CMP			R9, RSI
		JGE			h_end_i_loop

		MOV			R10, R9				; R10: j=i
		INC			R10					; R10: j=i+1
h_start_loop_j:
		; if (j < n)
		CMP			R10, RSI
		JGE			h_end_j_loop

		; hydrophobicity[s[i]-'a']
		XOR			RAX, RAX				; RAX = 0
		MOV			AL, [RDI+R9]			; AL = s[i]
		SUB			AL, 65					; AL = s[i]-'a'
		VMOVSD		XMM6, [R8+RAX*8]		; XMM6 = hydrophobicity[s[i]-'a']

		; hydrophobicity[s[j]-'a']
		XOR			RAX, RAX				; RAX = 0
		MOV			AL, [RDI+R10]			; AL = s[j]
		SUB			AL, 65					; AL = s[j]-'a'
		VMOVSD		XMM7, [R8+RAX*8]		; XMM7 = hydrophobicity[s[i]-'a']

		; ix
		MOV			RAX, 3					; RAX = 3
		IMUL		RAX, R9					; RAX = i*3
		ADD			RAX, 1					; RAX = i*3 + 1 //Atomo Ca
		IMUL		RAX, R11				; RAX = (i*3 + 1)*4*8 //Coordinata
		; jx
		MOV			RCX, 3					; RCX = 3
		IMUL		RCX, R10				; RCX = j*3
		ADD			RCX, 1					; RCX = j*3 + 1 //Atomo Ca
		IMUL		RCX, R11				; RCX = (j*3 + 1)*4*8 //Coordinata

		; v[ix]
		VMOVAPD			YMM1, [RDX+RAX]				; YMM1 <- |0 |v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		;VANDPD			YMM1, YMM1, [mask]			; YMM1 <- | 0 | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		; v[jx]
		VMOVAPD			YMM2, [RDX+RCX]				; YMM2 <- | 0 |v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
		;VANDPD			YMM2, YMM2, [mask]			; YMM2 <- | 0 | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |

		; Distance
		dist			YMM1, YMM2, XMM1, XMM2

		; if (d < 10.0)
		VCMPSD			XMM5, XMM1, XMM4, 5
		MOVMSKPD		RCX, XMM5
		CMP				RCX, 0
		JG				h_j_loop

		; ENERGIA
		VMULSD			XMM6, XMM6, XMM7			; XMM6 = hydrophobicity[s[i]-'a'] * hydrophobicity[s[j]-'a']
		VDIVSD			XMM6, XMM6, XMM1			; XMM6 = XMM6 / d
		VADDSD			XMM0, XMM0, XMM6			; XMM0 += XMM6

h_j_loop:
		INC			R10						; j++
		JMP			h_start_loop_j

h_end_j_loop:
		INC 		R9						; i++
		JMP			h_start_loop_i

h_end_i_loop:
		MOV			RCX, [hEPointer]
		VMOVSD		[RCX], XMM0			; e = XMM7
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		    ; ripristina il Base Pointer
		ret				    ; torna alla funzione C chiamante