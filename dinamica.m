clear all
close all
clc

syms lc2 lc3 a2 a3 d1 m2 m3 I2 I3 d1(t) th2(t) th3(t) g real

Jl2 = [1 -lc2*sin(th2) 0;
       0 lc2*cos(th2)  0];

Ja2 = [0 1 0];

Jl3 = [1 -a2*sin(th2)-lc3*sin(th2+th3) -lc3*sin(th2+th3);
       0 a2*cos(th2)+lc2*cos(th2+th3) lc2*cos(th2+th3)];

Ja3 = [0 1 1];

H = m2*Jl2'*Jl2 + I2*Ja2'*Ja2 + m3*Jl3'*Jl3 + I3*Ja3'*Ja3;

dq = [diff(d1(t),t); diff(th2(t),t); diff(th3(t),t)];

T = (dq'*H*dq)/2;

U = m2*g*lc2*sin(th2)+m3*g*(a2*sin(th2)+lc3*sin(th2+th3));

L = T-U;

D1 = diff(L, dq(1));

D2 = diff(L, dq(2));

D3 = diff(L, dq(3));

F1 = diff(D1, t) - diff(L, d1)

T2 = diff(D2, t) - diff(L, th2)

T3 = diff(D3, t) - diff(L, th3)

%% Arrumando as equações em função dos termos de aceleração
clear
clc

syms f1 t2 t3 ddd1 ddtheta_2 ddtheta_3 dtheta_2 dtheta_3 dd1 real
syms lc2 lc3 a2 a3 d1 m2 m3 I2 I3 s2 s3 s23 c2 c3 c23 g real


%Equação F1 sem os conjugados
%(m2 + m3)*diff(d1(t), t, t) - (m3*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))) + lc2*m2*sin(th2(t)))*diff(th2(t), t, t) - diff(th2(t), t)*(m3*(lc3*cos(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t)) + a2*cos(th2(t))*diff(th2(t), t)) + lc2*m2*cos(th2(t))*diff(th2(t), t)) - lc3*m3*sin(th2(t) + th3(t))*diff(th3(t), t, t) - lc3*m3*cos(th2(t) + th3(t))*diff(th3(t), t)*(diff(th2(t), t) + diff(th3(t), t))

%Equação T2 sem os conjugados
%diff(th2(t), t, t)*(I2 + I3 + m3*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t)))^2 + m3*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t)))^2 + lc2^2*m2*cos(th2(t))^2 + lc2^2*m2*sin(th2(t))^2) - (m3*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))) + lc2*m2*sin(th2(t)))*diff(d1(t), t, t) + diff(th3(t), t)*(lc3*m3*sin(th2(t) + th3(t))*(lc3*cos(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t)) + a2*cos(th2(t))*diff(th2(t), t)) - lc2*m3*cos(th2(t) + th3(t))*(lc2*sin(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t)) + a2*sin(th2(t))*diff(th2(t), t)) - lc2*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t)))*(diff(th2(t), t) + diff(th3(t), t)) + lc3*m3*cos(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t)))) + (diff(th3(t), t)*(diff(th3(t), t)*(2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t))*lc2^2 - 2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t))*lc3^2) + diff(th2(t), t)*(lc2*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t))) - lc3*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc3*cos(th2(t) + th3(t))) + lc2*m3*cos(th2(t) + th3(t))*(a2*sin(th2(t)) + lc2*sin(th2(t) + th3(t))) - lc3*m3*cos(th2(t) + th3(t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t)))) + lc3*m3*cos(th2(t) + th3(t))*diff(d1(t), t)))/2 + (((m3*(a2*cos(th2(t)) + lc3*cos(th2(t) + th3(t))) + lc2*m2*cos(th2(t)))*diff(th2(t), t) + lc3*m3*cos(th2(t) + th3(t))*diff(th3(t), t))*diff(d1(t), t))/2 + (I3 + lc2*m3*cos(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t))) + lc3*m3*sin(th2(t) + th3(t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))))*diff(th3(t), t, t) - (2*m3*(lc2*sin(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t)) + a2*sin(th2(t))*diff(th2(t), t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t))) - 2*m3*(lc3*cos(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t)) + a2*cos(th2(t))*diff(th2(t), t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))))*diff(th2(t), t) + (diff(th2(t), t)*(diff(th3(t), t)*(lc2*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t))) - lc3*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc3*cos(th2(t) + th3(t))) + lc2*m3*cos(th2(t) + th3(t))*(a2*sin(th2(t)) + lc2*sin(th2(t) + th3(t))) - lc3*m3*cos(th2(t) + th3(t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t)))) + (2*m3*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t)))*(a2*sin(th2(t)) + lc2*sin(th2(t) + th3(t))) - 2*m3*(a2*cos(th2(t)) + lc3*cos(th2(t) + th3(t)))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))))*diff(th2(t), t) + (m3*(a2*cos(th2(t)) + lc3*cos(th2(t) + th3(t))) + lc2*m2*cos(th2(t)))*diff(d1(t), t)))/2 - diff(d1(t), t)*(m3*(lc3*cos(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t)) + a2*cos(th2(t))*diff(th2(t), t)) + lc2*m2*cos(th2(t))*diff(th2(t), t)) + g*m3*(a2*cos(th2(t)) + lc3*cos(th2(t) + th3(t))) + g*lc2*m2*cos(th2(t))

%Equação T3 sem os conjugados
%diff(th2(t), t)*(lc3*m3*sin(th2(t) + th3(t))*(lc3*cos(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t)) + a2*cos(th2(t))*diff(th2(t), t)) - lc2*m3*cos(th2(t) + th3(t))*(lc2*sin(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t)) + a2*sin(th2(t))*diff(th2(t), t)) - lc2*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t)))*(diff(th2(t), t) + diff(th3(t), t)) + lc3*m3*cos(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t)))) + (diff(th2(t), t)*(diff(th3(t), t)*(lc2*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t))) - lc3*m3*cos(th2(t) + th3(t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))) + lc2^2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t)) - lc3^2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t))) + (2*lc2*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t))) - 2*lc3*m3*cos(th2(t) + th3(t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))))*diff(th2(t), t) + lc3*m3*cos(th2(t) + th3(t))*diff(d1(t), t)))/2 + diff(th3(t), t, t)*(m3*lc2^2*cos(th2(t) + th3(t))^2 + m3*lc3^2*sin(th2(t) + th3(t))^2 + I3) + ((lc3*m3*cos(th2(t) + th3(t))*diff(th2(t), t) + lc3*m3*cos(th2(t) + th3(t))*diff(th3(t), t))*diff(d1(t), t))/2 + (diff(th3(t), t)*(diff(th3(t), t)*(2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t))*lc2^2 - 2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t))*lc3^2) + diff(th2(t), t)*(lc2*m3*sin(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t))) - lc3*m3*cos(th2(t) + th3(t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))) + lc2^2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t)) - lc3^2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t))) + lc3*m3*cos(th2(t) + th3(t))*diff(d1(t), t)))/2 + (I3 + lc2*m3*cos(th2(t) + th3(t))*(a2*cos(th2(t)) + lc2*cos(th2(t) + th3(t))) + lc3*m3*sin(th2(t) + th3(t))*(a2*sin(th2(t)) + lc3*sin(th2(t) + th3(t))))*diff(th2(t), t, t) - (2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t))*lc2^2 - 2*m3*cos(th2(t) + th3(t))*sin(th2(t) + th3(t))*(diff(th2(t), t) + diff(th3(t), t))*lc3^2)*diff(th3(t), t) - lc3*m3*sin(th2(t) + th3(t))*diff(d1(t), t, t) + g*lc3*m3*cos(th2(t) + th3(t)) - lc3*m3*cos(th2(t) + th3(t))*diff(d1(t), t)*(diff(th2(t), t) + diff(th3(t), t))

%Equação F1 arrumada em termos da variáveis simbólicas
f1 = (m2 + m3)*ddd1 - (m3*(a2*s2 + lc3*s23) + lc2*m2*s2)*ddtheta_2 - dtheta_2*(m3*(lc3*c23*(dtheta_2 + dtheta_3) + a2*c2*dtheta_2) + lc2*m2*c2*dtheta_2) - lc3*m3*s23*ddtheta_3 - lc3*m3*c23*dtheta_3*(dtheta_2 + dtheta_3);

%Equação T2 arrumada em termos da variáveis simbólicas
t2 = ddtheta_2*(I2 + I3 + m3*(a2*s2 + lc3*s23)^2 + m3*(a2*c2 + lc2*c23)^2 + lc2^2*m2*c2^2 + lc2^2*m2*s2^2) - (m3*(a2*s2 + lc3*s23) + lc2*m2*s2)*ddd1 + dtheta_3*(lc3*m3*s23*(lc3*c23*(dtheta_2 + dtheta_3) + a2*c2*dtheta_2) - lc2*m3*c23*(lc2*s23*(dtheta_2 + dtheta_3) + a2*s2*dtheta_2) - lc2*m3*s23*(a2*c2 + lc2*c23)*(dtheta_2 + dtheta_3) + lc3*m3*c23*(dtheta_2 + dtheta_3)*(a2*s2 + lc3*s23)) + (dtheta_3*(dtheta_3*(2*m3*c23*s23*lc2^2 - 2*m3*c23*s23*lc3^2) + dtheta_2*(lc2*m3*s23*(a2*c2 + lc2*c23) - lc3*m3*s23*(a2*c2 + lc3*c23) + lc2*m3*c23*(a2*s2 + lc2*s23) - lc3*m3*c23*(a2*s2 + lc3*s23)) + lc3*m3*c23*dd1))/2 + (((m3*(a2*c2 + lc3*c23) + lc2*m2*c2)*dtheta_2 + lc3*m3*c23*dtheta_3)*dd1)/2 + (I3 + lc2*m3*c23*(a2*c2 + lc2*c23) + lc3*m3*s23*(a2*s2 + lc3*s23))*ddtheta_3 - (2*m3*(lc2*s23*(dtheta_2 + dtheta_3) + a2*s2*dtheta_2)*(a2*c2 + lc2*c23) - 2*m3*(lc3*c23*(dtheta_2 + dtheta_3) + a2*c2*dtheta_2)*(a2*s2 + lc3*s23))*dtheta_2 + (dtheta_2*(dtheta_3*(lc2*m3*s23*(a2*c2 + lc2*c23) - lc3*m3*s23*(a2*c2 + lc3*c23) + lc2*m3*c23*(a2*s2 + lc2*s23) - lc3*m3*c23*(a2*s2 + lc3*s23)) + (2*m3*(a2*c2 + lc2*c23)*(a2*s2 + lc2*s23) - 2*m3*(a2*c2 + lc3*c23)*(a2*s2 + lc3*s23))*dtheta_2 + (m3*(a2*c2 + lc3*c23) + lc2*m2*c2)*dd1))/2 - dd1*(m3*(lc3*c23*(dtheta_2 + dtheta_3) + a2*c2*dtheta_2) + lc2*m2*c2*dtheta_2) + g*m3*(a2*c2 + lc3*c23) + g*lc2*m2*c2;

%Equação T3 arrumada em termos da variáveis simbólicas
t3 = dtheta_2*(lc3*m3*s23*(lc3*c23*(dtheta_2 + dtheta_3) + a2*c2*dtheta_2) - lc2*m3*c23*(lc2*s23*(dtheta_2 + dtheta_3) + a2*s2*dtheta_2) - lc2*m3*s23*(a2*c2 + lc2*c23)*(dtheta_2 + dtheta_3) + lc3*m3*c23*(dtheta_2 + dtheta_3)*(a2*s2 + lc3*s23)) + (dtheta_2*(dtheta_3*(lc2*m3*s23*(a2*c2 + lc2*c23) - lc3*m3*c23*(a2*s2 + lc3*s23) + lc2^2*m3*c23*s23 - lc3^2*m3*c23*s23) + (2*lc2*m3*s23*(a2*c2 + lc2*c23) - 2*lc3*m3*c23*(a2*s2 + lc3*s23))*dtheta_2 + lc3*m3*c23*dd1))/2 + ddtheta_3*(m3*lc2^2*c23^2 + m3*lc3^2*s23^2 + I3) + ((lc3*m3*c23*dtheta_2 + lc3*m3*c23*dtheta_3)*dd1)/2 + (dtheta_3*(dtheta_3*(2*m3*c23*s23*lc2^2 - 2*m3*c23*s23*lc3^2) + dtheta_2*(lc2*m3*s23*(a2*c2 + lc2*c23) - lc3*m3*c23*(a2*s2 + lc3*s23) + lc2^2*m3*c23*s23 - lc3^2*m3*c23*s23) + lc3*m3*c23*dd1))/2 + (I3 + lc2*m3*c23*(a2*c2 + lc2*c23) + lc3*m3*s23*(a2*s2 + lc3*s23))*ddtheta_2 - (2*m3*c23*s23*(dtheta_2 + dtheta_3)*lc2^2 - 2*m3*c23*s23*(dtheta_2 + dtheta_3)*lc3^2)*dtheta_3 - lc3*m3*s23*ddd1 + g*lc3*m3*c23 - lc3*m3*c23*dd1*(dtheta_2 + dtheta_3);

%Isolando o termo de aceleração na primeira
% isolate(eqn1, ddd1);
% isolate(eqn2, ddtheta_2);
% isolate(eqn3, ddtheta_3);



%Equações prontas para o arquivo model2.m
% ddd1 == (f1 + dtheta_2*(m3*(a2*dtheta_2*c2 + lc3*c23*(dtheta_2 + dtheta_3)) + dtheta_2*lc2*m2*c2) + ddtheta_2*(m3*(a2*s2 + lc3*s23) + lc2*m2*s2) + ddtheta_3*lc3*m3*s23 + dtheta_3*lc3*m3*c23*(dtheta_2 + dtheta_3))/(m2 + m3)
% 
% ddtheta_2 == -((dtheta_3*(dtheta_3*(2*m3*c23*s23*lc2^2 - 2*m3*c23*s23*lc3^2) + dtheta_2*(lc2*m3*s23*(a2*c2 + lc2*c23) - lc3*m3*s23*(a2*c2 + lc3*c23) + lc2*m3*c23*(a2*s2 + lc2*s23) - lc3*m3*c23*(a2*s2 + lc3*s23)) + lc3*m3*c23*diff(d1(t), t)))/2 - dtheta_2*(2*m3*(a2*dtheta_2*s2 + lc2*s23*(dtheta_2 + dtheta_3))*(a2*c2 + lc2*c23) - 2*m3*(a2*dtheta_2*c2 + lc3*c23*(dtheta_2 + dtheta_3))*(a2*s2 + lc3*s23)) - t2 + dtheta_3*(lc3*m3*s23*(a2*dtheta_2*c2 + lc3*c23*(dtheta_2 + dtheta_3)) - lc2*m3*c23*(a2*dtheta_2*s2 + lc2*s23*(dtheta_2 + dtheta_3)) - lc2*m3*s23*(dtheta_2 + dtheta_3)*(a2*c2 + lc2*c23) + lc3*m3*c23*(dtheta_2 + dtheta_3)*(a2*s2 + lc3*s23)) - ddd1*(m3*(a2*s2 + lc3*s23) + lc2*m2*s2) - (m3*(a2*dtheta_2*c2 + lc3*c23*(dtheta_2 + dtheta_3)) + dtheta_2*lc2*m2*c2)*diff(d1(t), t) + (dtheta_2*(dtheta_2*(2*m3*(a2*c2 + lc2*c23)*(a2*s2 + lc2*s23) - 2*m3*(a2*c2 + lc3*c23)*(a2*s2 + lc3*s23)) + (m3*(a2*c2 + lc3*c23) + lc2*m2*c2)*diff(d1(t), t) + dtheta_3*(lc2*m3*s23*(a2*c2 + lc2*c23) - lc3*m3*s23*(a2*c2 + lc3*c23) + lc2*m3*c23*(a2*s2 + lc2*s23) - lc3*m3*c23*(a2*s2 + lc3*s23))))/2 + ddtheta_3*(I3 + lc2*m3*c23*(a2*c2 + lc2*c23) + lc3*m3*s23*(a2*s2 + lc3*s23)) + (diff(d1(t), t)*(dtheta_2*(m3*(a2*c2 + lc3*c23) + lc2*m2*c2) + dtheta_3*lc3*m3*c23))/2 + g*m3*(a2*c2 + lc3*c23) + g*lc2*m2*c2)/(I2 + I3 + m3*(a2*s2 + lc3*s23)^2 + m3*(a2*c2 + lc2*c23)^2 + lc2^2*m2*c2^2 + lc2^2*m2*s2^2)
% 
% ddtheta_3 == -((dtheta_2*(dtheta_3*(lc2*m3*s23*(a2*c2 + lc2*c23) - lc3*m3*c23*(a2*s2 + lc3*s23) + lc2^2*m3*c23*s23 - lc3^2*m3*c23*s23) + dtheta_2*(2*lc2*m3*s23*(a2*c2 + lc2*c23) - 2*lc3*m3*c23*(a2*s2 + lc3*s23)) + lc3*m3*c23*diff(d1(t), t)))/2 - dtheta_3*(2*m3*c23*s23*(dtheta_2 + dtheta_3)*lc2^2 - 2*m3*c23*s23*(dtheta_2 + dtheta_3)*lc3^2) - t3 + dtheta_2*(lc3*m3*s23*(a2*dtheta_2*c2 + lc3*c23*(dtheta_2 + dtheta_3)) - lc2*m3*c23*(a2*dtheta_2*s2 + lc2*s23*(dtheta_2 + dtheta_3)) - lc2*m3*s23*(dtheta_2 + dtheta_3)*(a2*c2 + lc2*c23) + lc3*m3*c23*(dtheta_2 + dtheta_3)*(a2*s2 + lc3*s23)) + (dtheta_3*(dtheta_2*(lc2*m3*s23*(a2*c2 + lc2*c23) - lc3*m3*c23*(a2*s2 + lc3*s23) + lc2^2*m3*c23*s23 - lc3^2*m3*c23*s23) + dtheta_3*(2*m3*c23*s23*lc2^2 - 2*m3*c23*s23*lc3^2) + lc3*m3*c23*diff(d1(t), t)))/2 + ddtheta_2*(I3 + lc2*m3*c23*(a2*c2 + lc2*c23) + lc3*m3*s23*(a2*s2 + lc3*s23)) + (diff(d1(t), t)*(dtheta_2*lc3*m3*c23 + dtheta_3*lc3*m3*c23))/2 + g*lc3*m3*c23 - ddd1*lc3*m3*s23 - lc3*m3*c23*(dtheta_2 + dtheta_3)*diff(d1(t), t))/(m3*lc2^2*c23^2 + m3*lc3^2*s23^2 + I3)
