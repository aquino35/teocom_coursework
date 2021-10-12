%s096g01dc14dtmf 
%clear all 
%close all 
Fs=8192; 
Ts=1/Fs; 
%V=10;           %Time Transmission window (TW) 
t=0:Ts:(81920*Ts)-Ts;    %Time axis for TW 
tt=0:Ts:(2048*Ts)-Ts;    %Pulse time tone “ON” width 
tz=0:Ts:(6144*Ts)-Ts;    %Pulse time zero “OFF” width 
s1=cos(2*pi*697*tt)+cos(2*pi*1209*tt); 
s3=cos(2*pi*697*tt)+cos(2*pi*1477*tt); 
s8=cos(2*pi*852*tt)+cos(2*pi*1336*tt); 
s0=cos(2*pi*941*tt)+cos(2*pi*1336*tt); 
sz=0*tz;          %Pulse signal segment “OFF” 
s=[s8 sz s3 sz s1 sz s8 sz s3 sz s1 sz s0 sz s3 sz s1 sz s8 sz]; 
figure
plot(t,s)
axis([0 max(t) -max(abs(s)) max(abs(s))]) 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('DTFM Message or Word of Order or Length L=10 and Duty Cycle=1/2 to 1/4') 
sound(s,Fs)


%s096g01dc38dtmf 
%clear all 
%close all 
Fs=8192; 
Ts=1/Fs; 
%V=10;                  %Time Transmission window (TW) 
t=0:Ts:(81920*Ts)-Ts;    %Time axis for TW 
tt=0:Ts:2048*Ts-Ts;    %Pulse time tone “ON” width 
tz=0:Ts:6144*Ts-Ts;    %Pulse time zero “OFF” width 
s1=cos(2*pi*697*tt)+cos(2*pi*1209*tt); 
s3=cos(2*pi*697*tt)+cos(2*pi*1477*tt); 
s8=cos(2*pi*852*tt)+cos(2*pi*1336*tt); 
s0=cos(2*pi*941*tt)+cos(2*pi*1336*tt); 
sz=0*tz;          %Pulse signal segment “OFF” 
s=[s8 sz s3 sz s1 sz s8 sz s3 sz s1 sz s0 sz s3 sz s1 sz s8 sz]; 
figure
plot(t,s) 
axis([0 max(t) -max(abs(s)) max(abs(s))]) 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('DTFM Message or Word of Order or Length L=10 and Duty Cycle=3/8 to 1/4') 
sound(s,Fs) 
