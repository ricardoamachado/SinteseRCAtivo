clear
clc

%Especificações do filtro.
freq_passagem = 1e6;
freq_rejeicao = 250e3;
A_max = 0.1;
A_min = 85;

%Obtenção de T(s) com base na aproximação de Chebyshev.
omega_passagem = 2*pi*freq_passagem;
omega_rejeicao = 2*pi*freq_rejeicao;
[n,omega_ne] = cheb1ord(omega_passagem,omega_rejeicao,A_max,A_min,'s');
[b,a] = cheby1(n,A_max,omega_ne,'high','s');
T = tf(b,a)

%Cálculo das raízes e polos de T(s).
[zn,pn,kn] = tf2zp(b,a);

%Síntese dos biquads.
n1 = poly([zn(1) zn(2)]);
n2 = poly([zn(3) zn(4)]);
n3 = poly([zn(5) zn(6)]);
p1 = poly([pn(1) pn(2)]);
p2 = poly([pn(3) pn(4)]);
p3 = poly([pn(5) pn(6)]);
t1 = kn*tf(n1,p1)
t2 = tf(n2,p2)
t3 = tf(n3,p3)

%Cálculos relativos ao biquad 1:
omega0_1 = sqrt(p1(3));
Q_1 = omega0_1/p1(2);
%Considerei R1 = R2 = R e C1 = C2 = 10pF
C_1 = 10e-12;
syms R k  
eq1_1 = 1/(R*C_1) == omega0_1;
eq2_1 = R*C_1 / (R*2*C_1+R*C_1*(1-k)) == Q_1;
resultado_1 = vpasolve(eq1_1,eq2_1);
Rnum_1 = resultado_1.R;
knum_1 = resultado_1.k;
% Valor escolhido para R3: 10 kohm.
R3_1 = 10000;
R4_1 = (knum_1 - 1)*R3_1;

%Cálculos relativos ao biquad 2:
omega0_2 = sqrt(p2(3));
Q_2 = omega0_2/p2(2);
%Considerei R1 = R2 = R e C1 = C2 = 10pF
C_2 = 10e-12; 
eq1_2 = 1/(R*C_2) == omega0_2;
eq2_2 = R*C_2 / (R*2*C_2+R*C_2*(1-k)) == Q_2;
resultado_2 = vpasolve(eq1_2,eq2_2);
Rnum_2 = resultado_2.R;
knum_2 = resultado_2.k;
% Valor escolhido para R3: 10 kohm.
R3_2 = 10000;
R4_2 = (knum_2 - 1)*R3_2;

%Cálculos relativos ao biquad 3:
omega0_3 = sqrt(p3(3));
Q_3 = omega0_3/p3(2);
%Considerei R1 = R2 = R e C1 = C2 = 10pF
C_3 = 10e-12;
eq1_3 = 1/(R*C_3) == omega0_3;
eq2_3 = R*C_3 / (R*2*C_3+R*C_3*(1-k)) == Q_3;
resultado_3 = vpasolve(eq1_3,eq2_3);
Rnum_3 = resultado_3.R;
knum_3 = resultado_3.k;
% Valor escolhido para R3: 10 kohm.
R3_3 = 10000;
R4_3 = (knum_3 - 1)*R3_3;

%Print valores dos biquads
fprintf(['Valores obtidos para o biquad 1:\n' ...
    'R1 = R2 = %.2f\nK = %.2f\nR4 = %.2f\n'],Rnum_1,knum_1,R4_1)
fprintf(['Valores obtidos para o biquad 2:\n' ...
    'R1 = R2 = %.2f\nK = %.2f\nR4 = %.2f\n'],Rnum_2,knum_2,R4_2)
fprintf(['Valores obtidos para o biquad 3:\n' ...
    'R1 = R2 = %.2f\nK = %.2f\nR4 = %.2f\n'],Rnum_3,knum_3,R4_3)

%Carremento dados simulação
dados = readmatrix('Circuito.txt');

%Gráficos
figure(1)
bode(T)
figure(2)
f = logspace(3,8,1000);
h = freqs(b,a,2*pi*f);
semilogx(f,20*log10(abs(h)));
xlabel('Frequência (Hz)')
ylabel('|T(s)|')
figure(3)
ganho_ajustado = knum_1*knum_2*knum_3/kn;
semilogx(f,20*log10(abs(ganho_ajustado*h)))
hold on
semilogx(dados(:,1),dados(:,2))
legend('Gráfico Calculado', 'Gráfico Simulado')
hold off
xlabel('Frequência (Hz)')
ylabel('|T(s)|')
