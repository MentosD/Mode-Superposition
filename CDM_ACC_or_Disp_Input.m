clear

load("ACC_el.mat");

M = [2 0 0;
    0 1.5 0
    0 0 1];
K = 600* [5 -2 0;
    -2 3 -1;
    0 -1 1];
C = 0.05 * M + 0.025 * K;
dt = 0.001;

dofs = length(M);
diagM = diag(M);
dt = 0.001;
Ke=M/(dt^2)+((C)/(2*dt));                       
a = K - (2 * M) / (dt)^2;
b=M/dt^2 - C/(2*dt);
u = zeros(dofs , length(ACC_el));
v = zeros(dofs , length(ACC_el));
ac = zeros(dofs , length(ACC_el));
jiazaiweizhi = ones(length(M), 1);

for i = 2 : length(ACC_el)
    PP = -ACC_el(i,2) * 9.8 * diagM  - a * u(: , i) - b * u(: , i-1);
    u(:,i+1)=Ke \ PP;                            
    v(: , i) = (u(: , i+1) - u(: , i-1)) / (dt*2);
    ac(: , i) = (u(: , i+1) - 2 * u(: , i) + u(: , i-1)) / (dt^2);
end
plot(u(3,1:54950)+ACC_el(:,4)');hold on;

for i = 2 : length(ACC_el)
    PP = (K * (ACC_el(i,4) * jiazaiweizhi) + C * (ACC_el(i,3) * jiazaiweizhi)) - a * u(: , i) - b * u(: , i-1);
    u(:,i+1)=Ke \ PP;                            
    v(: , i) = (u(: , i+1) - u(: , i-1)) / (dt*2);
    ac(: , i) = (u(: , i+1) - 2 * u(: , i) + u(: , i-1)) / (dt^2);
end
plot(u(3,:));

