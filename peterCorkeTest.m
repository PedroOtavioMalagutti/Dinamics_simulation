
clc; clear; close all;


% Ordem dos parâmetros: theta d a alpha
DH = [0  0    0    -pi/2 ;
      0  0   -1    0    ;
      0  0   -1   0 ];


% Links
% theta d a alpha tipo offset
% tipo: 0 - rotativa, 1 - primatica
L(1) = Link([DH(1,:) 1 0],'standard');
L(2) = Link([DH(2,:) 0 0],'standard');
L(3) = Link([DH(3,:) 0 0],'standard');

L(1).m=0.8;
L(2).m=1.5;
L(3).m=0.6;

L(1).r = [ 0.0 0.0 0.0];
L(2).r = [ -1/2 0 0];
L(3).r = [ -1/2 0 0];


L(1).I = [ 0 0 0 0 0 0];
L(2).I = [ 0.03125  0.03125  0.03125  0 0 0];
L(3).I = [ 0.125  0.125  0.125  0 0 0];

L(1).Jm=0;
L(2).Jm=0;
L(3).Jm=0;

L(1).qlim = [0 2];
L(2).qlim= [-pi/2 pi/2];
L(3).qlim= [-pi/2 pi/2];

L(1).B=0;
L(2).B=0;
L(3).B=0;

% robot
q0= [1 0 0];
rpy = [0, -pi/2, 0];
T_des = rt2tr(rpy2r(rpy(1),rpy(2),rpy(3)),[0;0;0]);

rob = SerialLink(L,'name','rob');
rob.base =T_des;

[t,q,qd] = rob.fdyn(6, @my_torque_function, q0);
figure
plot(t,q)

% clf
figure
rob.plot(q)

function tau = my_torque_function(rob, t, q, qd)
    qref=[1.5 pi/2 pi/2]; %posicao de junta desejada
    %ganhos do controle pd
    Kp=17;
    Kd=18;
    e=qref-q;
    %parametros do robo
    m1 = 1;
    m2 = 1.5;
    m3 = 0.6;
    i2 = 1.5/12;
    i3 = 0.6/12;
    a2=1;
    l2 = 1;
    l3 = 1;
    q2=q(2);
    q3=q(3);
    %controle de gravidade
    G(1)= 0;
    G(2)=(2943*l2*m2*sin(conj(q2)))/200 + (981*l2*m3*sin(conj(q2)))/100 + (2943*l3*m3*sin(conj(q2) + conj(q3)))/200;
    G(3)=(2943*l3*m3*sin(conj(q2) + conj(q3)))/200;
    %controle pd
    tau = G + Kp.*e-Kd*qd; 

end
