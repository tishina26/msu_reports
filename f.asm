global f1
global f2
global f3
global df1
global df2
global df3

section  .bss

section .data
; константы - коэффициенты в функциях
    c1 dq 0.35
    c2 dq -0.95
    c3 dq 2.7
    c4 dq 3.0
    c5 dq 1.0
    c6 dq 2.0
    c7 dq 0.7

section .text


;__________f______ФУНКЦИИ____________   
   
; 0.35x*x-0.95x+2.7 
f1:
    push ebp
    mov ebp, esp
    
    fld qword[ebp+8]
    fld qword[ebp+8]
    fmulp               ; x*x
    fld qword[c1]     ;0.35
    fmulp               ;0.35x*x
    fld qword[ebp+8]
    fld qword[c2]     ;-0.95
    fmulp               ;-0.95x
    faddp               ;0.35x*x-0.95x
    fld qword[c3]     ;2.7
    faddp               ;0.35x*x-0.95x +2.7
    
    leave
    ret
    
; 3x+1
f2:
    push ebp
    mov ebp, esp
    
    fld qword[ebp+8]
    fld qword[c4]     ;3
    fmulp               ;3x
    fld qword[c5]     ;1
    faddp               ;3x+1
    
    leave
    ret

; 1/(x+2)
f3:
    push ebp
    mov ebp, esp
    
    fld qword[ebp+8]
    fld qword[c6]     ;2
    faddp               ;x+2
    fld qword[c5]     ;1
    fdivrp              ;1/(x+2)
    
    leave
    ret
    
    

;_________dfi______ПРОИЗВОДНЫЕ_________

; 0.7x-0.95
df1:
    push ebp
    mov ebp, esp
    
    fld qword[ebp+8]
    fld qword[c7]     ;0.7
    fmulp               ;0.7x
    fld qword[c2]     ;-0.95
    faddp               ;0.7x-0.95
    
    leave
    ret
    
; 3 always
df2:
    push ebp
    mov ebp, esp
    
    fld qword[c4]     ;3
    
    leave
    ret
    
; -1/(x+2)^2
df3:
    push ebp
    mov ebp, esp
    
    fld qword[ebp+8]
    fld qword[c6]     ;2
    faddp               ;x+2
    fld qword[ebp+8]
    fld qword[c6]     ;2
    faddp               ;x+2
    fmulp               ; (x+2)^2
    fld qword[c5]     ;1
    fdivrp              ;1/(x+2)^2
    fchs                ;-1/(x+2)^2
    
    leave
    ret