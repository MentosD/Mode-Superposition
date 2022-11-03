clear
close all

format long

load("data.mat")
x = Input;
y = Output;
x = x-x(1);
y = y-y(1);
t=0.001;
[Hs,www]=tfestimate(x,y,[] , [] , length(x) ,1000);
mag=abs(Hs);
pha=angle(Hs)/pi*180;

data=iddata(y,x,0.001);
G=tfest(data,3,0);

Gcontrol = tf([50*pi],[1 , 50*pi]);
Gcompensate = (1/G) ;
t=[1:length(x)]*0.001-0.001;
y1=lsim(G,x,t);
y2=lsim(G*Gcompensate* (Gcontrol),x,t);
y2=lsim(G * (Gcontrol),x,t);

figure(2)
plot(t,[x,y,y1,y2],'linewidth',1)
legend('com','resp','3^t^hmodel','location','best')
xlabel('t(s)');
ylabel('Dis.(mm)')


[mag1,pha1]=bode(G,www*2*pi);
[mag2,pha2]=bode(G*Gcompensate* (Gcontrol^3),www*2*pi);
[mag3,pha3]=bode(Gcontrol^3,www*2*pi);

%
figure(1)
subplot(121)
mag = mag(:); mag1 = mag1(:); mag2 = mag2(:);
mag = mag(1:4000); mag1 = mag1(1:4000); mag2 = mag2(1:4000);
pha = pha(:); pha1 = pha1(:); pha2 = pha2(:);
pha = pha(1:4000); pha1 = pha1(1:4000); pha2 = pha2(1:4000);
www = www(1:4000);
plot(www,[mag],'linewidth',1);hold on;
plot(www,[mag1],'linewidth',1);hold on;
plot(www,[mag2],'linewidth',1);hold on;
xlim([0 20]);
legend('试验','3阶传递函数','location','best')
xlabel('w(Hz)');
ylabel('|Gst|')
subplot(122)
plot(www,[pha],'linewidth',1);hold on;
plot(www,[pha1(:)],'linewidth',1);hold on;
plot(www,[pha2(:)],'linewidth',1);hold off;
xlabel('w(Hz)'); xlim([0 40]);
ylabel('<Gst(^o)')
magg =  mag(1:4000);
magg(:,2) = mag1(1:4000);
magg(:,3) = mag2(1:4000);
phaa = pha(1:4000);
phaa(:,2) = pha1(1:4000);
phaa(:,3) = pha2(1:4000);


%
% figure(1)
% subplot(121)
% mag = mag(:); mag1 = mag1(:);
% mag = mag(1:4000); mag1 = mag1(1:4000);
% pha = pha(:); pha1 = pha1(:);
% pha = pha(1:4000); pha1 = pha1(1:4000);
% www = www(1:4000);
% plot(www,[mag],'linewidth',1);hold on;
% plot(www,[mag1],'linewidth',1);hold on;
% xlim([0 20]);
% legend('试验','3阶传递函数','location','best')
% xlabel('w(Hz)');
% ylabel('|Gst|')
% subplot(122)
% plot(www,[pha],'linewidth',1);hold on;
% plot(www,[pha1(:)],'linewidth',1);hold on;
% xlabel('w(Hz)'); xlim([0 40]);
% ylabel('<Gst(^o)')
% magg =  mag(1:4000);
% phaa = pha(1:4000);

