% CDM-3dofs 子结构

clear
close all
load("ACC_el.mat");

M_ns = [2 0;
        0 1.5];
K_ns = 600* [5 -2;
            -2 2];
C_ns = 0.25 * M_ns + 0.25 * K_ns;

M_ps = 1;
K_ps = 600;
C_ps = 0.25 * M_ps + 0.25 * K_ps;

M = [2 0 0;
    0 1.5 0
    0 0 M_ps];
K = 600* [5 -2 0;
    -2 2+K_ps/600 -K_ps/600;
    0 -K_ps/600 K_ps/600];
C = [750.5000 -300.0000 0;
     -300.0000 300.3750+150.2500 -150.2500;
     0 -150.2500 150.2500;]

dofs = length(M);
diagM = diag(M);
dt = 0.001;
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
Ke_ns_inv = inv(Ke_ns);
Ke_ps_inv = inv(Ke_ps);
ps_weizhi = [0 ; 1];
%%
for i = 2 : length(ACC_el)
    PP_ns = -ACC_el(i,2)* diagM_ns - F_ps(i-1) * ps_weizhi - a_ns * u_ns(: , i) - b_ns * u_ns(: , i-1);
    u_ns(:,i+1) = Ke_ns_inv * PP_ns;
    v_ns(: , i) = (u_ns(: , i+1) - u_ns(: , i-1)) / (dt*2);
    ac_ns(: , i) = (u_ns(: , i+1) - 2 * u_ns(: , i) + u_ns(: , i-1)) / (dt^2);
    PP_ps = -(ac_ns(2,i) + ACC_el(i,2)) * M_ps - a_ps * u_ps(i) - b_ps * u_ps(i-1);
    u_ps(i+1) = Ke_ps_inv * PP_ps;
    v_ps(i) = (u_ps(i+1) - u_ps(i-1)) / (dt*2);
    ac_ps(i) = (u_ps(i+1) - 2 * u_ps(i) + u_ps(i-1)) / (dt^2);
    F_ps(i) = (ac_ps(i) + ac_ns(2,i) + ACC_el(i,2)) * M_ps ;
end
ucdm_ns = u_ns;
ucdm_ps = u_ps;

%%

plot(u(3,:));
hold on;
plot(u_ps+u_ns(2,:));
