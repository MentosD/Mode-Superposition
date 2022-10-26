% CDM-2dofs 子结构

clear
close all
load("ACC_el.mat");

M_ns = 2;
K_ns = 1800;
C_ns = 300;

M_ps = 0.5*1e-5;
K_ps = 1200*1e-5;
C_ps = 200*1e-5;

M = [M_ns,0;0,M_ps];
K = [K_ns+K_ps, -K_ps; -K_ps, K_ps];
C = [C_ns+C_ps, -C_ps; -C_ps, C_ps];

dofs = length(M);
diagM = diag(M);
dt = 0.001;
Ke=M/(dt^2)+((C)/(2*dt));
a = K - (2 * M) / (dt)^2;
b=M/dt^2 - C/(2*dt);
u = zeros(dofs , length(ACC_el));
v = zeros(dofs , length(ACC_el));
ac = zeros(dofs , length(ACC_el));
Ke_inv = inv(Ke);
for i = 2 : length(ACC_el)
    PP = -ACC_el(i,2)* diagM  - a * u(: , i) - b * u(: , i-1);
    u(:,i+1)=Ke_inv * PP;
    v(: , i) = (u(: , i+1) - u(: , i-1)) / (dt*2);
    ac(: , i) = (u(: , i+1) - 2 * u(: , i) + u(: , i-1)) / (dt^2);
end

ucdm = u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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

for i = 2 : length(ACC_el)
    PP_ns = (-ACC_el(i,2) * M_ns) - F_ps(i-1) - a_ns * u_ns(i) - b_ns * u_ns(i-1);  %数值子结构部分
    u_ns(i+1) = Ke_ns \ PP_ns;
    v_ns(i) = (u_ns(i+1) - u_ns(i-1)) / (dt*2);
    ac_ns(i) = (u_ns(i+1) - 2 * u_ns(i) + u_ns(i-1)) / (dt^2);
    
    PP_ps = -(ac_ns(i) + ACC_el(i,2)) * M_ps - a_ps * u_ps(i) - b_ps * u_ps(i-1);   %物理子结构部分
    u_ps(i+1) = Ke_ps \ PP_ps;
    v_ps(i) = (u_ps(i+1) - u_ps(i-1)) / (dt*2);
    ac_ps(i) = (u_ps(i+1) - 2 * u_ps(i) + u_ps(i-1)) / (dt^2);
    
    F_ps(i) = (ac_ps(i) + ac_ns(i) + ACC_el(i,2)) * M_ps ;                          %物理子结构对数值子结构的反力
end

ucdm_ns = u_ns;
ucdm_ps = u_ps;
% plot(ucdm_ns(2, :));
%%
% G = tf([-M_ps,0,0],[M_ps,C_ps,K_ps]); 
% zzz = ac_ns(2,:)' + ACC_el(:,2);
% uG = lsim(G, ac_ns(2,:)'+ACC_el(:,2), [0.001:0.001:0.001*55001]);

% plot(uG);
% hold on;
% plot(ac_ps);

figure(1)
plot(u(2,:));
hold on;
plot(u_ns+u_ps);

% figure(2)
% plot(ac(1,:));
% hold on;
% plot(ac_ns)
