%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
alignb  32
mask	dq	0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x0

section .bss			; Sezione contenente dati non inizializzati

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

global fn_energy

fn_energy:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali
	; ------------------------------------------------------------
	; Funzione
	; ------------------------------------------------------------
		;RDI: *coords, RSI: ix, RDX: jx, RCX: *d
		VMOVAPD			YMM7, [mask]

		; v[ix]
		VMOVUPD			YMM0, [RDI+RSI*8]			; YMM0 <- | v[(i*3 + 1)*3+3] | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		VANDPD			YMM0, YMM0, YMM7			; YMM0 <- | 0 | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		; v[jx]
		VMOVUPD			YMM1, [RDI+RDX*8]			; YMM1 <- | v[(j*3 + 1)*3+3] | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
		VANDPD			YMM1, YMM1, YMM7			; YMM1 <- | 0 | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |

		; Distanza
		dist			YMM0, YMM1, XMM0, XMM1

		VMOVSD			[RCX], XMM0
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		    ; ripristina il Base Pointer
		ret				    ; torna alla funzione C chiamante