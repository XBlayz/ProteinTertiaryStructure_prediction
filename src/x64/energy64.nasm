%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
int_fmt db '%d', 0
msgNl	db  0x0a, 0

section .bss			; Sezione contenente dati non inizializzati
alignb  32
pEPointer	resq    1
eEPointer	resq    1

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

%macro iprint 1
    pushaq                      ; Save all registers
    mov rdi, int_fmt            ; Load format string into rdi
    mov rsi, %1                 ; Load the integer value into rsi
    xor rax, rax                ; Clear rax as required by printf
    call printf                 ; Call printf
    popaq                       ; Restore all registers

	prints msgNl
%endmacro

%macro vprint 2
	prints %1
	iprint %2
%endmacro

; ------------------------------------------------------------
; Funzione p_energy
; ------------------------------------------------------------
global p_energy

mask	dq	0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x0
ten		dq	10.0

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
		VXORPD		XMM0, XMM0				; XMM0: | 0 | density=0 |
		VXORPD		XMM7, XMM7				; XMM7: | 0 | energy=0 |

		MOV			R10, 3*8				; R10 = 3*8
start_loop_i:
		; if (i < n)
		CMP			R8, RSI
		JGE			end_i_loop

		XOR			R9, R9					; R9: j=0
start_loop_j:
		; if (j < n)
		CMP			R9, RSI
		JGE			end_j_loop
		; if (i == j)
		CMP			R8, R9
		JE			loop_j

		; ix
		MOV			RAX, 3					; RAX = 3
		IMUL		RAX, R8					; RAX = i*3
		ADD			RAX, 1					; RAX = i*3 + 1 //Atomo Ca
		IMUL		RAX, R10				; RAX = (i*3 + 1)*3*8 //Coordinata
		; jx
		MOV			RCX, 3					; RCX = 3
		IMUL		RCX, R9					; RCX = j*3
		ADD			RCX, 1					; RCX = j*3 + 1 //Atomo Ca
		IMUL		RCX, R10				; RCX = (j*3 + 1)*3*8 //Coordinata

		; v[ix]
		VMOVUPD			XMM1, [RDX+RAX]				; XMM1 <- v[(i*3 + 1)*3], v[(i*3 + 1)*3+1]
		VMOVSD			XMM2, [RDX+RAX+16]			; XMM2 <- v[(i*3 + 1)*3+2]
		VPERM2F128		YMM1, YMM1, YMM2, 0x20		; YMM1 <- | xxx | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		VANDPD			YMM1, YMM1, [mask]			; YMM1 <- |  0  | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		; v[jx]
		VMOVUPD			XMM2, [RDX+RCX]				; XMM2 <- v[(j*3 + 1)*3], v[(j*3 + 1)*3+1]
		VMOVSD			XMM3, [RDX+RCX+16]			; XMM3 <- v[(j*3 + 1)*3+2]
		VPERM2F128		YMM2, YMM2, YMM3, 0x20		; YMM2 <- | xxx | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
		VANDPD			YMM2, YMM2, [mask]			; YMM2 <- |  0  | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |

		VSUBPD			YMM1, YMM2, YMM1			; YMM1 <- | v[(j*3 + 1)*3] - v[(i*3 + 1)*3] | v[(j*3 + 1)*3+1] - v[(i*3 + 1)*3+1] | v[(j*3 + 1)*3+2] - v[(i*3 + 1)*3+2] |

		VMULPD			YMM1, YMM1, YMM1
		VPERM2F128		YMM2, YMM1, YMM1, 0x01
		VADDPD			XMM1, XMM1, XMM2
		VPERMILPD 		XMM2, XMM1, 0x03
		VADDPD			XMM1, XMM1, XMM2
		VSQRTSD			XMM1, XMM1, XMM1

		; if (d < 10.0)
		VMOVSD			XMM4, [ten]
		VCMPSD			XMM5, XMM1, XMM4, 5
		MOVMSKPD		RCX, XMM5
		CMP				RCX, 0
		JG				loop_j

		; CALCOLA DENSITà
		VMULSD			XMM2, XMM1, XMM1			; XMM2 = d^2
		VMULSD			XMM2, XMM2, XMM1			; XMM2 = d^3
		XOR				RAX, RAX					; RAX = 0
		MOV				AL, [RDI+R9]				; AL = s[j]
		SUB				AL, 97						; AL = s[j]-'a'
		VMOVSD			XMM3, [R11+RAX]				; XMM3 = volume[s[j]-'a']
		VDIVSD			XMM3, XMM3, XMM2			; XMM3 = volume[s[j]-'a'] / d^3
		VADDSD			XMM0, XMM0, XMM3			; XMM0 += volume[s[j]-'a'] / d^3

loop_j:
		INC 	R9				; j++
		JMP		start_loop_j

end_j_loop:
		XOR				RAX, RAX					; RAX = 0
		MOV				AL, [RDI+R9]				; AL = s[j]
		SUB				AL, 97						; AL = s[j]-'a'
		VMOVSD			XMM3, [R11+RAX]				; XMM3 = volume[s[j]-'a']
		VSUBSD			XMM3, XMM3, XMM0			; XMM3 = volume[s[j]-'a'] - XMM0
		VMULSD			XMM3, XMM3, XMM3			; XMM3 = (volume[s[j]-'a'] - XMM0)^2
		VADDSD			XMM7, XMM7, XMM3			; XMM7 += (volume[s[j]-'a'] - XMM0)^2
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

four	dq	4.0
msgOk	db	"ok", 10, 0
msgI	db	"i:", 0
msgJ	db	"j:", 0

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

		MOV			R11, 3*8				; R11 = 3*8
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
		SUB			AL, 97					; AL = s[i]-'a'
		VMOVSD		XMM6, [R8+RAX]			; XMM6 = charge[s[i]-'a']
		; if (charge[s[i]-'a'] != 0)
		VXORPD			XMM2, XMM2
		VCMPSD			XMM5, XMM6, XMM2, 0
		MOVMSKPD		RCX, XMM5
		CMP				RCX, 0
		JG				e_j_loop

		; charge[s[j]-'a']
		XOR			RAX, RAX				; RAX = 0
		MOV			AL, [RDI+R10]			; AL = s[j]
		SUB			AL, 97					; AL = s[j]-'a'
		VMOVSD		XMM7, [R8+RAX]			; XMM7 = charge[s[i]-'a']
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
		IMUL		RAX, R11				; RAX = (i*3 + 1)*3*8 //Coordinata
		; jx
		MOV			RCX, 3					; RCX = 3
		IMUL		RCX, R10				; RCX = j*3
		ADD			RCX, 1					; RCX = j*3 + 1 //Atomo Ca
		IMUL		RCX, R11				; RCX = (j*3 + 1)*3*8 //Coordinata

		; v[ix]
		VMOVUPD			XMM1, [RDX+RAX]				; XMM1 <- v[(i*3 + 1)*3], v[(i*3 + 1)*3+1]
		VMOVSD			XMM2, [RDX+RAX+16]			; XMM2 <- v[(i*3 + 1)*3+2]
		VPERM2F128		YMM1, YMM1, YMM2, 0x20		; YMM1 <- | xxx | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		VANDPD			YMM1, YMM1, [mask]			; YMM1 <- |  0  | v[(i*3 + 1)*3+2] | v[(i*3 + 1)*3+1] | v[(i*3 + 1)*3] |
		; v[jx]
		VMOVUPD			XMM2, [RDX+RCX]				; XMM2 <- v[(j*3 + 1)*3], v[(j*3 + 1)*3+1]
		VMOVSD			XMM3, [RDX+RCX+16]			; XMM3 <- v[(j*3 + 1)*3+2]
		VPERM2F128		YMM2, YMM2, YMM3, 0x20		; YMM2 <- | xxx | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |
		VANDPD			YMM2, YMM2, [mask]			; YMM2 <- |  0  | v[(j*3 + 1)*3+2] | v[(j*3 + 1)*3+1] | v[(j*3 + 1)*3] |

		; Distance
		VSUBPD			YMM1, YMM2, YMM1			; YMM1 <- | v[(j*3 + 1)*3] - v[(i*3 + 1)*3] | v[(j*3 + 1)*3+1] - v[(i*3 + 1)*3+1] | v[(j*3 + 1)*3+2] - v[(i*3 + 1)*3+2] |
		VMULPD			YMM1, YMM1, YMM1
		VPERM2F128		YMM2, YMM1, YMM1, 0x01
		VADDPD			XMM1, XMM1, XMM2
		VPERMILPD 		XMM2, XMM1, 0x03
		VADDPD			XMM1, XMM1, XMM2
		VSQRTSD			XMM1, XMM1, XMM1

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