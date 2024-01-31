# An-lise-estrutural-por-elementos-finitos
Projeto Técnicas e análise estrutural
Código implmentado no Jupyter notebook


import sympy as sp
sp.init_printing

import sympy as sp
I,E,L,A,theta = sp.symbols(['I','E','L','A','theta'])

#Simbolos dos esforços
symbols_force = ['s1','s2','s3','s4','s5','s5','s6','s7','s8','s9', 's10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36','s37','s38','s39']
#simbolos de deslocamentos
symbols_disp = ['g1','g2','g3','g4','g5','g5','g6','g7','g8','g9','g10','g11','g12','g13','g14','g15','g16','g17','g18','g19','g20','g21','g22','g23','g24','g25','g26','g27','g28','g29','g30','g31','g32','g33','g34','g35','g36','g37','g38','g39']

s1, s2, s3, s4, s5, s5, s6, s7, s8, s9,s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36, s37, s38, s39 = sp.symbols(symbols_force)
g1, g2,g3 , g4 , g5, g5, g6, g7 , g8, g9, g10, g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21, g22, g23, g24, g25, g26, g27, g28, g29, g30, g31, g32, g33, g34, g35, g36, g37, g38, g39 = sp.symbols(symbols_disp)


Ke_ = sp.Matrix([
    [E*A/L,0,0,-E*A/L,0,0],
    [0,12*E*I/(L**3), 6*E*I/(L**2),0,-12*E*I/(L**3),6*E*I/(L**2)],
    [0,6*E*I/(L**2),4*E*I/L,0,-6*E*I/(L**2),2*E*I/L],
    [-E*A/L,0,0,E*A/L,0,0],
    [0,-12*E*I/(L**3), -6*E*I/(L**2), 0, 12*E*I/(L**3), -6*E*I/(L**2)],
    [0,6*E*I/(L**2),2*E*I/L,0, -6*E*I/(L**2),4*E*I/L]])
Ke_

c = sp.cos(theta)
s = sp.sin(theta)

T = sp.Matrix([
    [c,-s,0,0,0,0],
    [s,c,0,0,0,0],
    [0,0,1,0,0,0],
    [0,0,0,c,-s,0],
    [0,0,0,s,c,0],
    [0,0,0,0,0,1]
]).T

ke = T.T*Ke_*T
T
ke


def barra(I_,E_,L_,A_,theta_,q1,q2,q3,q4,q5,q6):
    theta_= theta_ * 3.1415 / 180.
    n_element = ke.subs([
        (I, I_),
        (E, E_),
        (L, L_),
        (A, A_),
        (theta, theta_)
        ])
    
    n_element = n_element.row_insert(0,sp.Matrix([[q1,q2,q3,q4,q5,q6]]))
    n_element = n_element.col_insert(0,sp.Matrix([0,q1,q2,q3,q4,q5,q6]))
    
    return n_element 
    
    #I_,E_[MPa],L[mm]_,A_[mm^2],theta_,q1,q2,q3,q4,q5,q6
#módulo de elasticidade em MPa
    
barra1 = barra(30*10**7, 113.8e3, 950, 6000, 0, 1, 2, 3, 4, 5, 6)
barra2 = barra(30*10**7, 113.8e3, 950, 6000, 0, 4, 5, 6, 7, 8, 9)
barra3 = barra(30*10**7, 113.8e3, 750, 6000, 350, 7, 8, 9, 10, 11, 12)
barra4 = barra(30*10**7, 113.8e3, 750, 6000, 350, 10, 11, 12, 13, 14, 15)
barra5 = barra(30*10**7, 113.8e3, 500, 6000, 0, 13, 14, 15, 16, 17, 18)
barra6 = barra(30*10**7, 113.8e3, 500, 6000, 0, 16, 17, 18, 19, 20, 21)
barra7 = barra(30*10**7, 113.8e3, 500, 6000, 90, 19, 20, 21, 22, 23, 24)
barra8 = barra(30*10**7, 113.8e3, 500, 6000, 90, 22, 23, 24, 25, 26, 27)
barra9 = barra(30*10**7, 113.8e3, 800, 6000, 0, 25, 26, 27, 28, 29, 30)
barra10 = barra(30*10**7,113.8e3, 800, 6000, 0, 28, 29, 30, 31, 32, 33)
barra11 = barra(30*10**7,113.8e3, 500, 6000, 0, 31, 32, 33, 34, 35, 36)
barra12 = barra(30*10**7,113.8e3,1212.719,6000,18, 34, 35, 36, 37, 38, 39)
barra13 = barra(30*10**7,113.8e3,1212.719,6000,18, 37, 38, 39, 1, 2, 3)
barra14 = barra(30*10**7,113.8e3, 740,6000,46, 16, 17, 18, 22, 23, 24)
barra15 = barra(30*10**7,113.8e3,1000, 6000,90, 13, 14, 15, 28, 29, 30)
barra16 = barra(30*10**7,113.8e3,880,6000,88, 10, 11,12,31, 32, 33)
barra17 = barra(30*10**7,113.8e3,820,6000,63, 7, 8, 9,34, 35, 36)
barra18 = barra(30*10**7,113.8e3,380,6000,68, 4, 5,6,37, 38, 39)

barra1


#Armazena matrizes de rigidez de todos os elementos (usando método append)
M_elementos = []

M_elementos.append(barra1)
M_elementos.append(barra2)
M_elementos.append(barra3)
M_elementos.append(barra4)
M_elementos.append(barra5)
M_elementos.append(barra6)
M_elementos.append(barra7)
M_elementos.append(barra8)
M_elementos.append(barra9)
M_elementos.append(barra10)
M_elementos.append(barra11)
M_elementos.append(barra12)
M_elementos.append(barra13)
M_elementos.append(barra14)
M_elementos.append(barra15)
M_elementos.append(barra16)
M_elementos.append(barra17)
M_elementos.append(barra18)


k_de_elementos = len(M_elementos)
n_de_GL = 39
n_Struct = sp.Matrix.zeros(n_de_GL, n_de_GL)


for n in range(k_de_elementos):
    # Extrai os graus de liberdade do elemento n: no caso dessa serão 39 GDL matriz 39x39
    
    GL = []
    GL.append(M_elementos[n][0,1] - 1)
    GL.append(M_elementos[n][0,2] - 1)
    GL.append(M_elementos[n][0,3] - 1)
    GL.append(M_elementos[n][0,4] - 1)
    GL.append(M_elementos[n][0,5] - 1)
    GL.append(M_elementos[n][0,6] - 1)
   
    
   
    for i in range(len(GL)):
        for j in range(len(GL)):
            n_Struct[GL[i],GL[j]] += M_elementos[n][i+1,j+1]
            
           
           
       
#g deslocamentos que a estrutura está submetida

F = sp.Matrix([0,-4.75e3,-2.26e3,0,-14.25e3,-2.26e3,0,-14.25e3,-2.26e3,s10,s11,s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, 0, 0, 0, 0, 0, 0, 0, 0, 0, s34, s35, s36, s37, s38, s39])
#deslocamento
u = sp.Matrix([g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21, g22, g23, g24, 0, 0, 0, 0, 0, 0, 0, 0, 0, g34, g35, g36, g37, g38, g39 ])

resultado
