%include "sseutils64.nasm"

section .data
; DEBUG String
msg_d   db  '%f', 10, 0
msg_end db  '-', 10, 0

data_t    dq  1.0, 2.0, 3.0, 0.0, 4.0, 5.0, 6.0, 0.0

section .bss
; Utility memory address
tmp1    resq    4

section .text
global prova_ass

; -Print value macro-
; Print double
extern printf
%macro	print_double	1
		pushaq
		mov		  rax, 1
		mov		  rdi, msg_d
		vmovsd  xmm0, [%1]
		call		printf
		popaq
%endmacro
; Print end
%macro	print_end 0
    pushaq
    mov		  rax, 1
    mov		  rdi, msg_end
    call		printf
    popaq
%endmacro

; -Print register macro-
; Print ymm
%macro	print_ymm	1
    VMOVUPD   [tmp1], %1
    print_double tmp1
    print_double tmp1+8
    print_double tmp1+16
    print_double tmp1+24
    print_end
%endmacro
; Print xmm
%macro	print_xmm	1
    VMOVUPD   [tmp1], %1
    print_double tmp1
    print_double tmp1+8
    print_end
%endmacro

; -TEST Procedure-
prova_ass:
    ; ---CALL---
		push		rbp				  ; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						  ; salva i registri generali


    ; RDI: *data, RSI: *r
    VMOVUPD         YMM4, [data_t]
    VMOVUPD         YMM5, [data_t+32]

    ; -Formula distanza-
    ; Differenza
    VSUBPD			    YMM4, YMM5, YMM4            ; YMM4: |0|Az-Bz|Ay-By|Ax-Bx|
    ; Quadrati
    VMULPD			    YMM4, YMM4, YMM4            ; YMM4: |0|(Az-Bz)^2|(Ay-By)^2|(Ax-Bx)^2|

    ; Somma
    VHADDPD         YMM4, YMM4, YMM4            ; YMM4: |(Az-Bz)^2|(Az-Bz)^2|(Ay-By)^2 + (Ax-Bx)^2|(Ay-By)^2 + (Ax-Bx)^2|
    VEXTRACTF128    XMM5, YMM4, 0b00000001      ; XMM5: |(Az-Bz)^2|(Az-Bz)^2|
    VADDSD          XMM4, XMM4, XMM5            ; XMM4: |---|(Az-Bz)^2 + (Ay-By)^2 + (Ax-Bx)^2|

    ; Radice quadrata
    VSQRTSD         XMM4, XMM4                  ; XMM4: |---|sqrt((Az-Bz)^2 + (Ay-By)^2 + (Ax-Bx)^2)|

    VMOVSD          [RSI], XMM4


    ; ---RETURN---
    popaq				        ; ripristina i registri generali
		mov		rsp, rbp      ; ripristina lo Stack Pointer
		pop		rbp		        ; ripristina il Base Pointer
		ret				          ; torna alla funzione C chiamante