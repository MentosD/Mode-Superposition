clear
clear
close all
load("ACC_el.mat");
%%  结构参数
m1 = 1.01E6; m2 = 9.89E5; m3 = 9.89E5; m4 = 9.89E5; m5 = 9.89E5;
m6 = 9.89E5; m7 = 9.89E5; m8 = 9.89E5; m9 = 1.07E6;                        % 单位：kg
k1 = 3.1836E8; k2 = 3.5761E8; k3 = 3.3031E8; k4 = 2.9785E8; k5 = 2.801E8;
k6 = 2.5567E8; k7 = 2.003E8; k8 = 1.4945E8; k9 = 1.2114E8;                 % 单位：kN/m
M_ns = [[m1, 0, 0, 0, 0, 0, 0, 0, 0];
    [0, m2, 0, 0, 0, 0, 0, 0, 0];
    [0, 0, m3, 0, 0, 0, 0, 0, 0];
    [0, 0, 0, m4, 0, 0, 0, 0, 0];
    [0, 0, 0, 0, m5, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, m6, 0, 0, 0];
    [0, 0, 0, 0, 0, 0, m7, 0, 0];
    [0, 0, 0, 0, 0, 0, 0, m8, 0];
    [0, 0, 0, 0, 0, 0, 0, 0, m9]];
K_ns = [[k1+k2, -k2, 0, 0, 0, 0, 0, 0, 0];
    [-k2, k2+k3, -k3, 0, 0, 0, 0, 0, 0];
    [0, -k3, k3+k4, -k4, 0, 0, 0, 0, 0];
    [0, 0, -k4, k4+k5, -k5, 0, 0, 0, 0];
    [0, 0, 0, -k5, k5+k6, -k6, 0, 0, 0];
    [0, 0, 0, 0, -k6, k6+k7, -k7, 0, 0];
    [0, 0, 0, 0, 0, -k7, k7+k8, -k8, 0];
    [0, 0, 0, 0, 0, 0, -k8, k8+k9, -k9];
    [0, 0, 0, 0, 0, 0, 0, -k9, k9]];
zeta = 0.02; w1 = 0.4479; w2 = 1.16;                                       % 定义结构阻尼比，1、2阶频率 9层
a_coe = 2*w1*w2*zeta/(w1+w2); b_coe = 2*zeta/(w1+w2);
C_ns = a_coe*M_ns + b_coe*K_ns;                                            % 计算瑞利阻尼

M_ps = m9 * 1e-5;
K_ps = k9 * 1e-5;
C_ps = 3.0136e+06 * 1e-5;

M = [[m1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    [0, m2, 0, 0, 0, 0, 0, 0, 0, 0];
    [0, 0, m3, 0, 0, 0, 0, 0, 0, 0];
    [0, 0, 0, m4, 0, 0, 0, 0, 0, 0];
    [0, 0, 0, 0, m5, 0, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, m6, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, 0, m7, 0, 0, 0];
    [0, 0, 0, 0, 0, 0, 0, m8, 0, 0];
    [0, 0, 0, 0, 0, 0, 0, 0, m9, 0];
    [0, 0, 0, 0, 0, 0, 0, 0, 0, M_ps]];
K = [[k1+k2, -k2, 0, 0, 0, 0, 0, 0, 0, 0];
    [-k2, k2+k3, -k3, 0, 0, 0, 0, 0, 0, 0];
    [0, -k3, k3+k4, -k4, 0, 0, 0, 0, 0, 0];
    [0, 0, -k4, k4+k5, -k5, 0, 0, 0, 0, 0];
    [0, 0, 0, -k5, k5+k6, -k6, 0, 0, 0, 0];
    [0, 0, 0, 0, -k6, k6+k7, -k7, 0, 0, 0];
    [0, 0, 0, 0, 0, -k7, k7+k8, -k8, 0, 0];
    [0, 0, 0, 0, 0, 0, -k8, k8+k9, -k9, 0];
    [0, 0, 0, 0, 0, 0, 0, -k9, k9+K_ps, -K_ps];
    [0, 0, 0, 0, 0, 0, 0, 0, -K_ps, K_ps]];
C = C_ns;
C(9, 9)   = C(9,9)+C_ps;
C(10, 9)  = -C_ps;
C(9, 10)  = -C_ps;
C(10, 10) = C_ps;
dt = 0.001;
dofs = length(M);
diagM = diag(M);
%% 整体结构计算
Ke=M/(dt^2)+((C)/(2*dt));
a = K - (2 * M) / (dt)^2;
b=M/dt^2 - C/(2*dt);
u = zeros(dofs , length(ACC_el));
v = zeros(dofs , length(ACC_el));
ac = zeros(dofs , length(ACC_el));

for i = 2 : length(ACC_el)
    PP = -ACC_el(i,2)* diagM  - a * u(: , i) - b * u(: , i-1);
    u(:,i+1)=Ke \ PP;
    v(: , i) = (u(: , i+1) - u(: , i-1)) / (dt*2);
    ac(: , i) = (u(: , i+1) - 2 * u(: , i) + u(: , i-1)) / (dt^2);
end
ucdm = u;


%% 子结构计算
dofs_ns = length(M_ns);
dofs_ps = length(M_ps);
Ke_ns=M_ns / (dt^2) + ((C_ns) / (2*dt));
a_ns = K_ns - (2 * M_ns) / (dt)^2;
b_ns=M_ns / dt^2 - C_ns / (2*dt);
Ke_ps=M_ps / (dt^2) + ((C_ps) / (2*dt));
a_ps = K_ps - (2 * M_ps) / (dt)^2;
b_ps=M_ps / dt^2 - C_ps / (2*dt);
diagM_ns = diag(M_ns);
u_ns = zeros(dofs_ns , length(ACC_el));
v_ns = zeros(dofs_ns , length(ACC_el));
ac_ns = zeros(dofs_ns , length(ACC_el));
u_ps = zeros(dofs_ps , length(ACC_el));
v_ps = zeros(dofs_ps , length(ACC_el));
ac_ps = zeros(dofs_ps , length(ACC_el));
F_ps = zeros(1, length(ACC_el));
F_weizhi = zeros(1, length(M_ns));
F_weizhi(9) = 1;
for i = 2 : length(ACC_el)
    PP_ns = -ACC_el(i,2)* diagM_ns - F_ps(i-1) * F_weizhi' - a_ns * u_ns(: , i) - b_ns * u_ns(: , i-1);
    u_ns(:,i+1) = Ke_ns \ PP_ns;
    v_ns(: , i) = (u_ns(: , i+1) - u_ns(: , i-1)) / (dt*2);
    ac_ns(: , i) = (u_ns(: , i+1) - 2 * u_ns(: , i) + u_ns(: , i-1)) / (dt^2);

    PP_ps = -(ac_ns(length(M_ns),i) + ACC_el(i,2)) * M_ps - a_ps * u_ps(i) - b_ps * u_ps(i-1);
    u_ps(i+1) = Ke_ps \ PP_ps;
    v_ps(i) = (u_ps(i+1) - u_ps(i-1)) / (dt*2);
    ac_ps(i) = (u_ps(i+1) - 2 * u_ps(i) + u_ps(i-1)) / (dt^2);
    F_ps(i) = (ac_ps(i) + ac_ns(length(M_ns),i) + ACC_el(i,2)) * M_ps ;
end
ucdm_ns = u_ns;
ucdm_ps = u_ps;


%% 振型叠加法—整体结构
[V,D]=eig(inv(M)*K);
freq=diag(D).^0.5;
[Bc,ord] = sort(freq);                                                     %ord为记录顺序的向量
wsc=freq(ord);                                                             %角（圆）频率 rad/s
fsc=wsc/2/pi;                                                              %频率 Hz
V=V(:,ord);                                                                %振型按频率阶数排序  一阶振型是第一列
VV = V(:,1:length(M));                                                     %调整阶数
for i = 1:length(M)
    VV(:,i) = VV(:,i)/VV(length(M),i);
end

order = 10;
ksi = zeta;
Mn = VV' * M * VV;
for i = 1:order
    Mn_ = VV(:,i)' * M * VV(:,i);
    eval(['M',num2str(i),'=Mn_;']);
end

Kn  = VV' * K * VV;
for i = 1:order
    Kn_ = VV(:,i)' * K * VV(:,i);
    eval(['K',num2str(i),'=Kn_;']);
end

Cn  = VV' * C * VV;
for i = 1:order
    Cn_ = VV(:,i)' * C * VV(:,i);
    eval(['C',num2str(i),'=Cn_;']);
end

wn = sqrt(Kn/Mn);
for i = 1:order
    eval(['w',num2str(i),'=sqrt(K',num2str(i),'/M',num2str(i),');']);
end

for i = 1:order
    eval(['wDn',num2str(i),'=w',num2str(i) ,' * sqrt(1-ksi^2);;']);
end

for i = 1:order
    eval(['qn',num2str(i),'=zeros(length(ACC_el),1);']);
end

uu = zeros(length(M),length(ACC_el));

for i = 1:length(ACC_el)
    for ii = 1:order
        eval(['P',num2str(ii),'(:,i) = VV(:,',num2str(ii),')''* (ACC_el(i, 2) * diagM);']);
    end
end

Ken = Mn/(dt^2)+((Cn)/(2*dt));
an = Kn - (2 * Mn) / (dt)^2;
bn = Mn / dt^2 - Cn / (2*dt);
qn= zeros(order,length(ACC_el));

for i = 1:order
    eval(['PPP(i,:) = P',num2str(i),';']);
end

for i = 2 : length(ACC_el)-1
    PPn = -PPP(:, i) - an * qn(: , i) - bn * qn(: , i-1);
    qn(:,i+1)=Ken \ PPn;
    for iiii = 1:order
        eval(['uu(:,i+1) = uu(:,i+1) + VV(:,',num2str(iiii),') * qn(',num2str(iiii),',i+1);']);
    end
end


%% 振型叠加法—子结构
% close all
M = M_ns;
K = K_ns;
C = C_ns;
diagM = diag(M);
[V,D]=eig(inv(M)*K);
freq=diag(D).^0.5;
[Bc,ord] = sort(freq);                                                     %ord为记录顺序的向量
wsc=freq(ord);                                                             %角（圆）频率 rad/s
fsc=wsc/2/pi;                                                              %频率 Hz
V=V(:,ord);                                                                %振型按频率阶数排序  一阶振型是第一列
VV = V(:,1:length(M));                                                     %调整阶数
for i = 1:length(M)
    VV(:,i) = VV(:,i)/VV(length(M),i);
end

order = 9;
ksi = zeta;
Mn = VV' * M * VV;
for i = 1:order
    Mn_ = VV(:,i)' * M * VV(:,i);
    eval(['M',num2str(i),'=Mn_;']);
end

Kn  = VV' * K * VV;
for i = 1:order
    Kn_ = VV(:,i)' * K * VV(:,i);
    eval(['K',num2str(i),'=Kn_;']);
end

Cn  = VV' * C * VV;
for i = 1:order
    Cn_ = VV(:,i)' * C * VV(:,i);
    eval(['C',num2str(i),'=Cn_;']);
end

wn = sqrt(Kn/Mn);
for i = 1:order
    eval(['w',num2str(i),'=sqrt(K',num2str(i),'/M',num2str(i),');']);
end

for i = 1:order
    eval(['wDn',num2str(i),'=w',num2str(i) ,' * sqrt(1-ksi^2);;']);
end

for i = 1:order
    eval(['qn',num2str(i),'=zeros(length(ACC_el),1);']);
end

u_ModeS_ns = zeros(length(M),length(ACC_el));
v_ModeS_ns = zeros(length(M),length(ACC_el));
ac_ModeS_ns = zeros(length(M),length(ACC_el));
u_ModeS_ps = zeros(1, length(ACC_el));
v_ModeS_ps = zeros(1, length(ACC_el));
ac_ModeS_ps = zeros(1, length(ACC_el));
F_ModeS_ps = zeros(1, length(ACC_el));


Ken = Mn/(dt^2)+((Cn)/(2*dt));
an = Kn - (2 * Mn) / (dt)^2;
bn = Mn / dt^2 - Cn / (2*dt);
Ke_ps=M_ps / (dt^2) + ((C_ps) / (2*dt));
a_ps = K_ps - (2 * M_ps) / (dt)^2;
b_ps=M_ps / dt^2 - C_ps / (2*dt);
diagMn = diag(Mn);
qn= zeros(order,length(ACC_el));


ps_weizhi = zeros(length(M_ns),1);
ps_weizhi(length(M_ns),1) = 1;
PACC_el = zeros(length(M_ns),1);

for i = 2 : length(ACC_el)-1

    for ii = 1:order
        PACC_el(ii,i) = VV(:,ii)' * (ACC_el(i, 2) * diagM + F_ModeS_ps(i-1) * ps_weizhi);
    end

    PPn = -PACC_el(:, i) - an * qn(: , i) - bn * qn(: , i-1);
    qn(:,i+1)=Ken \ PPn;

    for iiii = 1:order
        eval(['u_ModeS_ns(:,i+1) = u_ModeS_ns(:,i+1) + VV(:,',num2str(iiii),') * qn(',num2str(iiii),',i+1);']);
    end

    v_ModeS_ns(:,i) = (u_ModeS_ns(: , i+1) - u_ModeS_ns(: , i-1)) / (dt*2);
    ac_ModeS_ns(: , i) = (u_ModeS_ns(: , i+1) - 2 * u_ModeS_ns(: , i) + u_ModeS_ns(: , i-1)) / (dt^2);

    PP_ModeS_ps = -(ac_ModeS_ns(length(M_ns),i) + ACC_el(i,2)) * M_ps...
                  - a_ps * u_ModeS_ps(i) - b_ps * u_ModeS_ps(i-1);

    u_ModeS_ps(i+1) = Ke_ps \ PP_ModeS_ps;
    v_ModeS_ps(i) = (u_ModeS_ps(i+1) - u_ModeS_ps(i-1)) / (dt*2);
    ac_ModeS_ps(i) = (u_ModeS_ps(i+1) - 2 * u_ModeS_ps(i) + u_ModeS_ps(i-1)) / (dt^2);

    F_ModeS_ps(i) = (ac_ModeS_ps(i) + ac_ModeS_ns(2,i) + ACC_el(i,2)) * M_ps ;
end

figure(1)
plot(u(10,:));
title('整体结构与子结构10层响应对比')
hold on;
plot(ucdm_ps + ucdm_ns(9,:));
hold on;
plot(real(uu(10,:)));                                                      %振型叠加法与CDM整体
hold on;
plot(u_ModeS_ps + u_ModeS_ns(9,:));
