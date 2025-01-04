%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

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
; Funzione sum_quad
; ------------------------------------------------------------
global sum_quad

alignb  32
mask	dq	0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x0

sum_quad:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; Funzione
	; ------------------------------------------------------------
		VMOVAPD			XMM0, [RDI]					; XMM0 <- v[0], v[1]
		VMOVSD			XMM1, [RDI+16]				; XMM1 <- v[2]
		VPERM2F128		YMM0, YMM0, YMM1, 0x20		; YMM0 <- | xxx | v[2] | v[1] | v[0] |
		VANDPD			YMM0, YMM0, [mask]			; YMM0 <- |  0  | v[2] | v[1] | v[0] |

		VMULPD			YMM0, YMM0, YMM0			; YMM0 <- |    0   | v[2]^2 | v[1]^2 | v[0]^2 |
		VPERM2F128		YMM1, YMM0, YMM0, 0x01		; YMM1 <- | v[1]^2 | v[0]^2 |    0   | v[2]^2 |
		VADDPD			XMM0, XMM0, XMM1			; XMM0 <- | v[1]^2 | v[0]^2 + v[2]^2 |
		VPERMILPD 		XMM1, XMM0, 0x03			; XMM1 <- | v[0]^2 + v[2]^2 | v[1]^2 |
		VADDPD			XMM0, XMM0, XMM1			; XMM0 <- | v[0]^2 + v[1]^2 + v[2]^2 | v[0]^2 + v[1]^2 + v[2]^2 |
		VMOVSD			[RSI], XMM0					; r = v[0]^2 + v[1]^2 + v[2]^2
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		    ; ripristina il Base Pointer
		ret				    ; torna alla funzione C chiamante