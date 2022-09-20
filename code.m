%%振型叠加法

clear
% clc

M = [2 0 0;
    0 1.5 0
    0 0 1];
K = 600* [5 -2 0;
    -2 3 -1;
    0 -1 1];

[V,D]=eig(inv(M)*K);
freq=diag(D).^0.5;
[Bc,ord] = sort(freq);                  %ord为记录顺序的向量
wsc=freq(ord);                          %角（圆）频率 rad/s
fsc=wsc/2/pi;                           %频率 Hz
V=V(:,ord);                             %振型按频率阶数排序  一阶振型是第一列

VV = V(:,1:3);                                                      %调整阶数
VV(:,1) = VV(:,1)/VV(3,1);
VV(:,2) = VV(:,2)/VV(3,2);
VV(:,3) = VV(:,3)/VV(3,3);

Mn = VV' * M * VV;
M1 = VV(:,1)' * M * VV(:,1);
M2 = VV(:,2)' * M * VV(:,2);
M3 = VV(:,3)' * M * VV(:,3);

Kn  = VV' * K * VV;
K1 = VV(:,1)' * K * VV(:,1);
K2 = VV(:,2)' * K * VV(:,2);
K3 = VV(:,3)' * K * VV(:,3);

wn = sqrt(Kn/Mn);
w1 = sqrt(K1/M1);
w2 = sqrt(K2/M2);
w3 = sqrt(K3/M3);

wn = real(wn);
Phi_ = Mn^(-1) * VV' * M;
PhiT_ = M * VV * Mn^(-1);
Cn = 2 * 0.005 * wn * Mn;
C = PhiT_ * Cn * Phi_;
C = zeros(3,3);
load("ACC_el.mat");


% Pn1 = (VV(:,1)' * (ACC_el(i, 2) * [1;1;1]));
% Pn2 = (VV(:,2)' * (ACC_el(i, 2) * [1;1;1]));
% Pn3 = (VV(:,3)' * (ACC_el(i, 2) * [1;1;1]));


qn1 = zeros(length(ACC_el),1);
qn2 = zeros(length(ACC_el),1);
qn3 = zeros(length(ACC_el),1);
uu = zeros(3,length(ACC_el));

syms t

for i = 2:55000
%     qn1(i) = qn1(i-1) + (1 / (M1*w1)) * (0.001 * ((VV(:,1)' * (ACC_el(i, 2) * [1;1;1])) * sin(w1 * (55 - i * 0.001))));
%     qn2(i) = qn2(i-1) + (1 / (M2*w2)) * (0.001 * ((VV(:,2)' * (ACC_el(i, 2) * [1;1;1])) * sin(w2 * (55 - i * 0.001))));
%     qn3(i) = qn3(i-1) + (1 / (M3*w3)) * (0.001 * ((VV(:,3)' * (ACC_el(i, 2) * [1;1;1])) * sin(w3 * (55 - i * 0.001))));
%     qn1(i) = qn1(i) + (1 / (M1*w1)) * (0.001 * ((VV(:,1)' * (ACC_el(i, 2) * [1;1;1])) * sin(w1 * (i * 0.001))));
%     qn2(i) = qn2(i) + (1 / (M2*w2)) * (0.001 * ((VV(:,2)' * (ACC_el(i, 2) * [1;1;1])) * sin(w2 * (i * 0.001))));
%     qn3(i) = qn3(i) + (1 / (M3*w3)) * (0.001 * ((VV(:,3)' * (ACC_el(i, 2) * [1;1;1])) * sin(w3 * (i * 0.001))));
%     qn1(i) = qn1(i-1) + ((1 / (M1*w1)) * (0.001 * ((VV(:,1)' * (ACC_el(i, 2) * [1;1;1])) * sin(w1 * (55 - i * 0.001)))));
%     qn2(i) = qn2(i-1) + ((1 / (M2*w2)) * (0.001 * ((VV(:,2)' * (ACC_el(i, 2) * [1;1;1])) * sin(w2 * (55 - i * 0.001)))));
%     qn3(i) = qn3(i-1) + ((1 / (M3*w3)) * (0.001 * ((VV(:,3)' * (ACC_el(i, 2) * [1;1;1])) * sin(w3 * (55 - i * 0.001)))));
    qn1(i) = (1 / (M1*w1)) * int( sin(2 * t) * sin(w1*(55-t)), 0, i);
    qn2(i) = (1 / (M2*w2)) * int( sin(2 * t) * sin(w2*(55-t)), 0, i);
    qn3(i) = (1 / (M3*w3)) * int( sin(2 * t) * sin(w3*(55-t)), 0, i);
    
    uu(:,i) = V(:,1) * qn1(i) + V(:,2) * qn2(i) + V(:,3) * qn3(i);
end

plot(uu(1,:));
hold on;

%***********************************************************************************%

dofs = length(M);
diagM = diag(M);
alpha = 0.25;
beta = 0.5;
dt = 0.001;
a0 = (1 / alpha / dt / dt);
a1 = (beta /alpha / dt);
a2 = (1 / alpha / dt);
a3 = (1 / 2 / alpha -1);
a4 = (beta / alpha - 1);
a5 = (dt / 2 * (beta / alpha -2));
a6 = (dt * (1 - beta));
a7 = (dt * beta);
Ke = K + a0 * M + a1 * C;

u = zeros(dofs , length(ACC_el));
v = zeros(dofs , length(ACC_el));
ac = zeros(dofs , length(ACC_el));
uuu = u;
for i = 2:55000
    PP = -ACC_el(i,2) * [1;1;1] + M * (a0 * u(:,i-1) + a2 * v(:,i-1) + a3 * ac(:,i-1)) + C * (a1 * u(:,i-1) + a4 * v(:,i-1) + a5 * ac(:,i-1));
    u(:,i) = Ke \ PP;
    uuu(:,i) = u(:,i) + uuu(:,i-1);
    ac(: , i) = a0 * (u(: , i) - u(: , i-1)) - a3 * ac(: , i-1) - a2 * v(: , i-1);
    v(: , i)= v(: , i-1) + a6 * ac(: , i-1) + a7 * ac(: , i);
end

% plot(u(1,:));
% hold on; 
% plot(uuu(1,:));
