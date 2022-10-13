%%振型叠加法-含阻尼

clear
% clc

load("ACC_el.mat");
ACC_el = ACC_el(1:1000,:);

M = [2 0 0;
    0 1.5 0
    0 0 1];
K = 600* [5 -2 0;
    -2 3 -1;
    0 -1 1];
ksi = 0.005;
dt = 0.001;
diagM = diag(M);
[V,D]=eig(inv(M)*K);
freq=diag(D).^0.5;
[Bc,ord] = sort(freq);                  %ord为记录顺序的向量
wsc=freq(ord);                          %角（圆）频率 rad/s
fsc=wsc/2/pi;                           %频率 Hz
V=V(:,ord);                             %振型按频率阶数排序  一阶振型是第一列

VV = V(:,1:3);                                                      %调整阶数
% VV(:,1) = VV(:,1)/VV(3,1);
% VV(:,2) = VV(:,2)/VV(3,2);
% VV(:,3) = VV(:,3)/VV(3,3);

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
Cn = 2 * ksi * wn * Mn;

% C = PhiT_ * Cn * Phi_;
Rayleigh_A0 = ((2 * ksi) * (w1 * w2)) / (w1 + w2);
Rayleigh_A1 = ((2 * ksi) * 1) / (w1 + w2);
C = Rayleigh_A0 * M +  Rayleigh_A1 * K;                            %用瑞丽阻尼
wDn = wn * sqrt(1-ksi^2);
Cn = VV' * C * VV;

wDn1 = w1 * sqrt(1-ksi^2);
wDn2 = w2 * sqrt(1-ksi^2);
wDn3 = w3 * sqrt(1-ksi^2);

qn= zeros(3,length(ACC_el));
qn1 = zeros(length(ACC_el),1);
qn2 = zeros(length(ACC_el),1);
qn3 = zeros(length(ACC_el),1);
uu = zeros(3,length(ACC_el));

for i = 1:length(ACC_el)
    P1(:,i) = VV(:,1)' * (ACC_el(i, 2) * diagM);
    P2(:,i) = VV(:,2)' * (ACC_el(i, 2) * diagM);
    P3(:,i) = VV(:,3)' * (ACC_el(i, 2) * diagM);
end

Ken = Mn/(dt^2)+((Cn)/(2*dt));                       
an = Kn - (2 * Mn) / (dt)^2;
bn = Mn / dt^2 - Cn / (2*dt);
diagMn = diag(Mn);
for i = 2 : length(ACC_el)
    PPn = [P1(i);P2(i);P3(i)] - an * qn(: , i) - bn * qn(: , i-1);
    qn(:,i+1)=Ken \ PPn;                            
    uu(:,i) = VV(:,1) * qn(1,i) + VV(:,2) * qn(2,i) + VV(:,3) * qn(3,i);
end

% for i = 0.001:0.001:length(ACC_el)*0.001
%     for ii = 0.001:0.001:i
%         qn1(round(ii*1e3)+1) = qn1(round(ii*1e3)) + 0.001 * (1 / (M1*wDn1)) * exp(-ksi * w1 * (i-ii)) * P1(round(ii*1e3)) * sin(wDn1*(i-ii));
%         qn2(round(ii*1e3)+1) = qn2(round(ii*1e3)) + 0.001 * (1 / (M2*wDn2)) * exp(-ksi * w2 * (i-ii)) * P2(round(ii*1e3)) * sin(wDn2*(i-ii));
%         qn3(round(ii*1e3)+1) = qn3(round(ii*1e3)) + 0.001 * (1 / (M3*wDn3)) * exp(-ksi * w3 * (i-ii)) * P3(round(ii*1e3)) * sin(wDn3*(i-ii));
%     end
%     uu(:,round(ii*1e3)) = VV(:,1) * qn1(round(ii*1e3)+1) + VV(:,2) * qn2(round(ii*1e3)+1) + VV(:,3) * qn3(round(ii*1e3)+1);
% end

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
