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
ganho = 1.5; %Ganho escolhido para a banda de passagem do filtro.
T = ganho*tf(b,a)

%Cálculo das raízes e polos de T(s).
[zn,pn,kn] = tf2zp(b,a);
kn = ganho*kn;

%Síntese dos biquads.
n1 = poly([zn(1) zn(2)]);
n2 = poly([zn(3) zn(4)]);
n3 = poly([zn(5) zn(6)]);
p1 = poly([pn(1) pn(2)]);
p2 = poly([pn(3) pn(4)]);
p3 = poly([pn(5) pn(6)]);
kn_biquad = nthroot(kn,3); %Ganho de cada um dos biquads - Raiz cúbica.
t1 = kn_biquad*tf(n1,p1)
t2 = kn_biquad*tf(n2,p2)
t3 = kn_biquad*tf(n3,p3)

%Valores dos componentes escolhidos.
r3_1 = 10000; 
r2_1 = 20000;
c2_1 = 220e-12;
%Cálculo dos componentes do biquad 1 (Sallen-&-Key Passa Altas).
syms r1_1 r4_1 c1_1
eqk_1 = r4_1/r3_1 + 1 == kn_biquad;
eq1_1 = 1/(r2_1*c1_1) + 1/(r2_1*c2_1) + (1-kn_biquad)/(r1_1*c1_1) == p1(2);
eq2_1 = 1/(r1_1*r2_1*c1_1*c2_1) == p1(3);
r4_1_num = vpasolve(eqk_1, r4_1)
resultado = vpasolve(eq1_1,eq2_1,r1_1,c1_1);
r1_1_num = resultado.r1_1
c1_1_num = resultado.c1_1
%Escolha de valores númericos para os componentes:


%Gráficos
figure(1)
bode(T)
figure(2)
f = logspace(3,8,1000);
h = freqs(b,a,2*pi*f);
semilogx(f,20*log10(abs(h)));
xlabel('Frequência (Hz)')
ylabel('|T(s)|')



