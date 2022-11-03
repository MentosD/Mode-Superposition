clear;
filename = 'whitenoise08.tsi';
[cloumn1,cloumn2,cloumn3,cloumn4,cloumn5,cloumn6]=textread(filename,'%f%f%f%f%f%f','headerlines', 4);
data=[cloumn1 cloumn2 cloumn3 cloumn4 cloumn5 cloumn6];

%
AAccelEXC = cloumn1;
AAccelEXA = cloumn2;
AAccelEXB = cloumn3;
APlayedDriveX = cloumn4;
ADispX = cloumn5;
APosition_XDispCmdInput = cloumn6;

plot(APosition_XDispCmdInput*100,'DisplayName','APosition_XDispCmdInput');hold on;plot(ADispX,'DisplayName','ADispX');hold off;


Input = APosition_XDispCmdInput*100;
Output = ADispX;


%
% AAccelEXC = cloumn5;
% AAccelEXA = cloumn3;
% AAccelEXB = cloumn4;
% Acctable = cloumn4;
% ADispX = cloumn1;
% UDPcmd = cloumn6;