

Fs=40000; 
Ts=1/Fs; 
%V=10;           %Time Transmission window (TW) 
tt=0:Ts:(2048*Ts)-Ts;    %Time axis for TW 
tz=0:Ts:(1024*Ts)-Ts;    %Pulse time tone “ON” width 
s1=cos(2*pi*697*tt)+cos(2*pi*1209*tt); 
s3=cos(2*pi*697*tt)+cos(2*pi*1477*tt); 
s8=cos(2*pi*852*tt)+cos(2*pi*1336*tt); 
s0=cos(2*pi*941*tt)+cos(2*pi*1336*tt); 
sz=0*tz;          %Pulse signal segment “OFF” 
s=[s1 sz s3 sz s8 sz s8 sz s3 sz s3 sz s1 sz s0 sz s0 sz s1 sz];

sound(s, Fs);

audiowrite('dtmf01.wav', s, Fs);