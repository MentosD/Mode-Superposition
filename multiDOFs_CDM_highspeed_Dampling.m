%%振型叠加法
% 使用num2str(i)严重拖慢了运行速度，此脚本中手动编程，取消自动更改变量名

clear
% clc
load("MCK1215.mat","K","M");
load("ACC_el.mat");
% M = [2 0 0;
%     0 1.5 0
%     0 0 1];
% K = 600* [5 -2 0;
%     -2 3 -1;
%     0 -1 1];
ACC_el = ACC_el(1:10000,:);
dt = 0.05;
ksi = 0.00;
order = 10;

diagM = diag(M);
% [V,D]=eig(inv(M)*K);
[V,D]=eig(M\K);
freq=diag(D).^0.5;
[Bc,ord] = sort(freq);                                                       %ord为记录顺序的向量
wsc=freq(ord);                                                               %角（圆）频率 rad/s
fsc=wsc/2/pi;                                                                %频率 Hz
V=V(:,ord);                                                                  %振型按频率阶数排序  一阶振型是第一列
V = real(V);
VV = V(:,1:order);                                                           %调整阶数
% for i = 1:order
%     VV(:,i) = VV(:,i)/VV(length(M),i);
% end

for i = 1:order
    i = i;
    Mn_ = VV(:,i)' * M * VV(:,i);
    i = i;
    eval(['M',num2str(i),'=Mn_;']);
end

Kn  = VV' * K * VV;
for i = 1:order
    i = i;
    Kn_ = VV(:,i)' * K * VV(:,i);
    i = i;
    eval(['K',num2str(i),'=Kn_;']);
end
Mn = VV(:,1:order)' * M * VV(:,1:order);
Kn  = VV(:,1:order)' * K * VV(:,1:order);
wn = sqrt(Kn/Mn);
for i = 1:order
    eval(['w',num2str(i),'=sqrt(K',num2str(i),'/M',num2str(i),');']);
end

for i = 1:order
    eval(['wDn',num2str(i),'=w',num2str(i) ,' * sqrt(1-ksi^2);;']);
end

wn = real(wn);
Phi_ = Mn^(-1) * VV' * M;
PhiT_ = M * VV * Mn^(-1);
Cn = 2 * ksi * wn * Mn;

Rayleigh_A0 = ((2 * ksi) * (w1 * w2)) / (w1 + w2);
% Rayleigh_A0 = ((2 * ksi) * (17.5711 * 17.7796)) / (17.5711 + 17.7796);
Rayleigh_A1 = ((2 * ksi) * 1) / (w1 + w2);
% Rayleigh_A1 = ((2 * ksi) * 1) / (17.5711 + 17.7796);

% C = PhiT_ * Cn * Phi_;              
C = (Rayleigh_A0 * M +  Rayleigh_A1 * K);                            %用瑞丽阻尼

load("MCK1215.mat","C");
Cn = VV(:,1:order)' * C * VV(:,1:order);

for i = 1:order
    eval(['qn',num2str(i),'=zeros(length(ACC_el),1);']);
end

uu = zeros(length(M),length(ACC_el));

for i = 1:length(ACC_el)
    for ii = 1:order
        eval(['P',num2str(ii),'(:,i) = VV(:,',num2str(ii),')''* (ACC_el(i, 2) * diagM);']);
        PACC(ii,i) =  VV(:,ii)' * (ACC_el(i, 2) * diagM);
    end
end

Ken = Mn/(dt^2)+((Cn)/(2*dt));                       
an = Kn - (2 * Mn) / (dt)^2;
bn = Mn / dt^2 - Cn / (2*dt);
diagMn = diag(Mn);
qn= zeros(order,length(ACC_el));

for i = 1:order
    eval(['PPP(i,:) = P',num2str(i),';']);
end

for i = 2 : length(ACC_el)
    PPn = PPP(:, i) - an * qn(: , i) - bn * qn(: , i-1);
    qn(:,i+1)=Ken \ PPn;                            
    for iiii = 1:order
        eval(['uu(:,i) = uu(:,i) + VV(:,',num2str(iiii),') * qn(',num2str(iiii),',i);']);
    end
end
plot(real(uu(1,:)),'linewidth',2);


%***********************************************************************************%
tic;
dofs = length(M);
diagM = diag(M);
dt = 0.001;
Ke=M/(dt^2)+((C)/(2*dt));
Kev = inv(Ke);
a = K - (2 * M) / (dt)^2;
b=M/dt^2 - C/(2*dt);
u = zeros(dofs , length(ACC_el));
v = zeros(dofs , length(ACC_el));
ac = zeros(dofs , length(ACC_el));
% ACC_el = gpuArray(ACC_el);
% diagM = gpuArray(diagM);
% a = gpuArray(a);
% b = gpuArray(b);
% u = gpuArray(u);
% v = gpuArray(v);
% ac = gpuArray(ac);
% Ke = gpuArray(Ke);
% PP = gpuArray(zeros(length(Mn),1));

for i = 2 : length(ACC_el)
    i
    PP = ACC_el(i,2)* diagM  - a * u(: , i) - b * u(: , i-1);
    u(:,i+1)=Ke \ PP;                            
    v(: , i) = (u(: , i+1) - u(: , i-1)) / (dt*2);
    ac(: , i) = (u(: , i+1) - 2 * u(: , i) + u(: , i-1)) / (dt^2);
end
ucdm = u;
% ucdm = gather(ucdm);
tCDM = toc;
%******************************************************************************%
tic;
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
Kev = inv(Ke);
% alpha = gpuArray(alpha);
% beta = gpuArray(beta);
% a1 = gpuArray(a1);
% a2 = gpuArray(a2);
% a3 = gpuArray(a3);
% a4 = gpuArray(a4);
% a5 = gpuArray(a5);
% a6 = gpuArray(a6);
% a7 = gpuArray(a7);
% Ke = gpuArray(Ke);

for i = 2:length(ACC_el)
    i
    PP = ACC_el(i,2)* diagM  + M * (a0 * u(:,i-1) + a2 * v(:,i-1) + a3 * ac(:,i-1)) + C * (a1 * u(:,i-1) + a4 * v(:,i-1) + a5 * ac(:,i-1));
    u(:,i) = Ke \ PP;
    ac(: , i) = a0 * (u(: , i) - u(: , i-1)) - a3 * ac(: , i-1) - a2 * v(: , i-1);
    v(: , i)= v(: , i-1) + a6 * ac(: , i-1) + a7 * ac(: , i);
end

unmb = u;
% unmb = gather(u);
tNMB = toc;
%**********************************************************************%

plot(real(uu(1,:)/207),'linewidth',2);
hold on;
plot(real(ucdm(1,:)),'linewidth',2);
hold on;
plot(real(unmb(1,:)),'linewidth',2);

% hold on; 
% plot(uuu(1,:));
