%%振型叠加法
% 使用num2str(i)严重拖慢了运行速度，此脚本中手动编程，取消自动更改变量名

clear
% clc
load("MCK1215.mat","C","K","M");
load("ACC_el.mat");
% M = [2 0 0;
%     0 1.5 0
%     0 0 1];
% K = 600* [5 -2 0;
%     -2 3 -1;
%     0 -1 1];

ACC_el = ACC_el(1:1000,:);

diagM = diag(M);
[V,D]=eig(inv(M)*K);
freq=diag(D).^0.5;
[Bc,ord] = sort(freq);                  %ord为记录顺序的向量
wsc=freq(ord);                          %角（圆）频率 rad/s
fsc=wsc/2/pi;                           %频率 Hz
V=V(:,ord);                             %振型按频率阶数排序  一阶振型是第一列
V = real(V);
VV = V(:,1:length(M));                                                      %调整阶数
for i = 1:length(M)
    VV(:,i) = VV(:,i)/VV(length(M),i);
end

order = 10;

Mn = VV' * M * VV;
for i = 1:order
    Mn_ = VV(:,i)' * M * VV(:,i);
    eval(['M',num2str(i),'=Mn_']);
end

Kn  = VV' * K * VV;
for i = 1:order
    Kn_ = VV(:,i)' * K * VV(:,i);
    eval(['K',num2str(i),'=Kn_;']);
end

wn = sqrt(Kn/Mn);
for i = 1:order
    wn_ = sqrt(K1/M1);
    eval(['w',num2str(i),'=sqrt(K',num2str(i),'/M',num2str(i),');']);
end

wn = real(wn);
Phi_ = Mn^(-1) * VV' * M;
PhiT_ = M * VV * Mn^(-1);
Cn = 2 * 0.005 * wn * Mn;
C = PhiT_ * Cn * Phi_;
C = zeros(length(M),length(M));

% Pn1 = (VV(:,1)' * (ACC_el(i, 2) * [1;1;1]));
% Pn2 = (VV(:,2)' * (ACC_el(i, 2) * [1;1;1]));
% Pn3 = (VV(:,3)' * (ACC_el(i, 2) * [1;1;1]));


for i = 1:order
    eval(['qn',num2str(i),'=zeros(length(ACC_el),1);']);
end

uu = zeros(length(M),length(ACC_el));

for i = 1:length(ACC_el)
    for ii = 1:order
        eval(['P',num2str(ii),'(:,i) = VV(:,',num2str(ii),')''* (ACC_el(i, 2) * diagM);']);
    end
end

for i = 0.001:0.001:length(ACC_el)*0.001
i
    for ii = 0.001:0.001:i
        qn1(round(ii*1e3)+1) = qn1(round(ii*1e3)) + 0.001 * (1 / (M1*w1)) * P1(round(ii*1e3)) * sin(w1*(i-ii));
        qn2(round(ii*1e3)+1) = qn2(round(ii*1e3)) + 0.001 * (1 / (M2*w2)) * P2(round(ii*1e3)) * sin(w2*(i-ii));
        qn3(round(ii*1e3)+1) = qn3(round(ii*1e3)) + 0.001 * (1 / (M3*w3)) * P3(round(ii*1e3)) * sin(w3*(i-ii));
        qn4(round(ii*1e3)+1) = qn4(round(ii*1e3)) + 0.001 * (1 / (M4*w4)) * P4(round(ii*1e3)) * sin(w4*(i-ii));
        qn5(round(ii*1e3)+1) = qn5(round(ii*1e3)) + 0.001 * (1 / (M5*w5)) * P5(round(ii*1e3)) * sin(w5*(i-ii));
        qn6(round(ii*1e3)+1) = qn6(round(ii*1e3)) + 0.001 * (1 / (M6*w6)) * P6(round(ii*1e3)) * sin(w6*(i-ii));
        qn7(round(ii*1e3)+1) = qn7(round(ii*1e3)) + 0.001 * (1 / (M7*w7)) * P7(round(ii*1e3)) * sin(w7*(i-ii));
        qn8(round(ii*1e3)+1) = qn8(round(ii*1e3)) + 0.001 * (1 / (M8*w8)) * P8(round(ii*1e3)) * sin(w8*(i-ii));
        qn9(round(ii*1e3)+1) = qn9(round(ii*1e3)) + 0.001 * (1 / (M9*w9)) * P9(round(ii*1e3)) * sin(w9*(i-ii));
        qn10(round(ii*1e3)+1) = qn10(round(ii*1e3)) + 0.001 * (1 / (M10*w10)) * P10(round(ii*1e3)) * sin(w10*(i-ii));
    end
    uu(:,round(ii*1e3)) = VV(:,1) * qn1(round(ii*1e3)+1) + VV(:,2) * qn2(round(ii*1e3)+1) + VV(:,3) * qn3(round(ii*1e3)+1)...
        + VV(:,4) * qn4(round(ii*1e3)+1)+ VV(:,5) * qn5(round(ii*1e3)+1)+ VV(:,6) * qn6(round(ii*1e3)+1)...
        + VV(:,7) * qn7(round(ii*1e3)+1)+ VV(:,8) * qn8(round(ii*1e3)+1)+ VV(:,9) * qn9(round(ii*1e3)+1)...
        + VV(:,10) * qn10(round(ii*1e3)+1);
end

%***********************************************************************************%
dofs = length(M);
diagM = diag(M);
dt = 0.001;
Ke=M/(dt^2)+((C)/(2*dt));                       
a = K - (2 * M) / (dt)^2;
b=M/dt^2 - C/(2*dt);
u = zeros(dofs , 2000);
v = zeros(dofs , 2000);
ac = zeros(dofs , 2000);

for i = 2 : length(ACC_el)
    PP = ACC_el(i,2)* diagM  - a * u(: , i) - b * u(: , i-1);
    u(:,i+1)=Ke \ PP;                            
    v(: , i) = (u(: , i+1) - u(: , i-1)) / (dt*2);
    ac(: , i) = (u(: , i+1) - 2 * u(: , i) + u(: , i-1)) / (dt^2);
end
ucdm = u;
%******************************************************************************%
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
for i = 2:length(ACC_el)
    PP = ACC_el(i,2)* diagM  + M * (a0 * u(:,i-1) + a2 * v(:,i-1) + a3 * ac(:,i-1)) + C * (a1 * u(:,i-1) + a4 * v(:,i-1) + a5 * ac(:,i-1));
    u(:,i) = Ke \ PP;
    uuu(:,i) = u(:,i) + uuu(:,i-1);
    ac(: , i) = a0 * (u(: , i) - u(: , i-1)) - a3 * ac(: , i-1) - a2 * v(: , i-1);
    v(: , i)= v(: , i-1) + a6 * ac(: , i-1) + a7 * ac(: , i);
end

unmb = u;
%**********************************************************************%

plot(uu(1,:));
hold on;
plot(ucdm(1,:));
hold on;
plot(unmb(1,:));

% hold on; 
% plot(uuu(1,:));
