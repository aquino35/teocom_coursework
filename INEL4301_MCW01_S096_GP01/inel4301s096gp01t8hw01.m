%inel4301s096gp01t8hw01 
%Inel 4301 - Communication Theory I 
%MATLAB CLASS ASSIGNAMENT 01(MCA01) 
%Double Sideband Suppressed Carrier Systems 
%Prof. Domingo Antonio Rodríguez 
% 
%% 
%%PARAMETER SETTINGS****************************** 
%% 
clear all 
close all 
Fs=16800;           %Sampling frequency 
Ts=1/Fs;            %Sampling time or temporal resolution 
f1=131;             %First frequency component (C3) of s(t) 
f2=175;             %Second frequency component (F3) of s(t) 
f3=247;             %Third frequency component of (B3) s(t) 
fc=1000;            %Frequency of carrier signal > 2*fm 
fm=300;             %Bandwidth of desired signal s(t) 
Ns=1200;             %Number of samples of s(t) 
Tw=Ns*Ts;           %Time observation window for s(t) 
t=0:Ts:Tw-Ts;       %Time axis for desired signal s(t)  
%% 
%%TRANSMITTER SYSTEM****************************** 
%% 
s=(-7*cos(2*pi*f1*t)+2*cos(2*pi*f2*t)+5*cos(2*pi*f3*t))/8; 
g=1*cos(2*pi*800*t)+1*cos(2*pi*500*t);    %Interference signal g(t) 
max_sg=max(s+g);        %Absolute value of signal s(t)+g(t) 
xm=(s+g)/max_sg;        %Normalized Modulating signal xm(t)         
ct=cos(2*pi*fc*t);      %Transmitter carrier signal ct(t) 
xc=xm.*ct;              %Modulated signal xc(t) 
%% 
%%Transmitter's Band-Pass Filter Design 
%% 
fu=fc+fm;           %Upper cut-off frequency 
fl=fc-fm;           %Lower cut-off frequency 
wu=(2*fu)/Fs;       %Normalized upper frequency       
wl=(2*fl)/Fs;       %Normalized lower frequency 
wb=[wl wu];         %Bandwidth of band-pass filter 
Mb=501;             %Order of band-pass filter 
thB=0:Ts:Mb*Ts-Ts;  %Time axis of impulse response signal hB(t) 
hB=fir1(Mb-1,wb,'bandpass');    %Filter design function 
yci=conv(xc,hB);        %Discrete-time convolution operation 
lyci=length(yci);       %Length of yci=Ns+Mb-1; 
tyci=0:Ts:(lyci-1)*Ts;  %Time axis of signal yci(t) 
%% 
%% COMMUNICATION CHANNEL************************** 
%%
nn=randn(lyci,1);           %AWGN with mean=0 and variance=1.

n=(1/1000)*transpose(nn);   %Attenuated noise signal n(t) 
yco=yci+0*n; 
%% 
%% RECEIVER SYSTEM******************************** 
%% 
cr=cos(2*pi*fc*tyci); 
yd=2*yco.*cr;       %Demodulated signal 
% 
%Low-Pass Filter Design*************************** 
% 
Wn=2*fm/Fs;         %Normalized cut-off frequency 
M=091;              %Order of the impulse response signal hL(t) 
thL=0:Ts:M*Ts-Ts;   %Time axis to plot impulse response signal 
hL=fir1(M-1,Wn);    %Impulse response signal 
% 
%Signal Filtering at the Receiver***************** 
% 
tr=0:Ts:(lyci+M-1)*Ts-Ts;   %Time axis for received signal 
xr=conv(yd,hL);             %Convolution operation for filtering 
Af=1;                       %Signal amplification factor 
xr=Af*xr; 
%% 
%%SIGNALS SPECTRA********************************** 
%% 
fres_s=1/(length(s)*Ts); 
faxis_s=-(Fs/2):fres_s:(Fs/2)-fres_s; 
Fss=fft(s); 
sFss=fftshift(Fss); 
asFss=abs(sFss); 
% 
fres_g=1/(length(g)*Ts); 
faxis_g=-(Fs/2):fres_g:(Fs/2)-fres_g; 
Fg=fft(g); 
sFg=fftshift(Fg); 
asFg=abs(sFg); 
% 
fres_xm=1/(length(xm)*Ts); 
faxis_xm=-(Fs/2):fres_xm:(Fs/2)-fres_xm; 
Fxm=fft(xm); 
sFxm=fftshift(Fxm); 
asFxm=abs(sFxm); 
% 
fres_ct=1/(length(ct)*Ts); 
faxis_ct=-(Fs/2):fres_ct:(Fs/2)-fres_ct; 
Fct=fft(ct); 
sFct=fftshift(Fct); 
asFct=abs(sFct); 
% 
fres_xc=1/(length(xc)*Ts); 

faxis_xc=-(Fs/2):fres_xc:(Fs/2)-fres_xc; 
Fxc=fft(xc); 
sFxc=fftshift(Fxc); 
asFxc=abs(sFxc); 
% 
fres_hB=1/(length(hB)*Ts); 
faxis_hB=-(Fs/2):fres_hB:(Fs/2)-fres_hB; 
FhB=fft(hB); 
sFhB=fftshift(FhB); 
asFhB=abs(sFhB); 
% 
fres_yci=1/(length(yci)*Ts); 
faxis_yci=-(Fs/2):fres_yci:(Fs/2)-fres_yci; 
Fyci=fft(yci); 
sFyci=fftshift(Fyci); 
asFyci=abs(sFyci); 
% 
fres_n=1/(length(n)*Ts); 
faxis_n=-(Fs/2):fres_n:(Fs/2)-fres_n; 
Fn=fft(n); 
sFn=fftshift(Fn); 
asFn=abs(sFn); 
% 
fres_yco=1/(length(yco)*Ts); 
faxis_yco=-(Fs/2):fres_yco:(Fs/2)-fres_yco; 
Fyco=fft(yco); 
sFyco=fftshift(Fyco); 
asFyco=abs(sFyco); 
% 
fres_cr=1/(length(cr)*Ts); 
faxis_cr=-(Fs/2):fres_cr:(Fs/2)-fres_cr; 
Fcr=fft(cr); 
sFcr=fftshift(Fcr); 
asFcr=abs(sFcr); 
% 
fres_yd=1/(length(yd)*Ts); 
faxis_yd=-(Fs/2):fres_yd:(Fs/2)-fres_yd; 
Fyd=fft(yd); 
sFyd=fftshift(Fyd); 
asFyd=abs(sFyd); 
% 
fres_hL=1/(length(hL)*Ts); 
faxis_hL=-(Fs/2):fres_hL:(Fs/2)-fres_hL; 
FhL=fft(hL); 
sFhL=fftshift(FhL); 
asFhL=abs(sFhL); 
% 
fres_xr=1/(length(xr)*Ts); 
faxis_xr=-(Fs/2):fres_xr:(Fs/2)-fres_xr;

Fxr=fft(xr); 
sFxr=fftshift(Fxr); 
asFxr=abs(sFxr); 
% 
%% 
%%PLOTS******************************************* 
%% 
% 
figure 
plot(t,s)               
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Desired Signal s(t)') 
% 
figure 
plot(t,g)               
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Interference Signal g(t)') 
% 
figure 
plot(t,xm)               
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Modulating Signal xm(t)') 
% 
figure 
plot(t,ct) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Transmitter Carrier Signal ct(t)') 
% 
figure 
plot(t,xc) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Modulated Signal xc(t)') 
% 
figure 
plot(thB,hB,thB,hB,'o') 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Impulse Response Signal hB(t)')

% 
figure 
plot(tyci,yci) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Channel Input Signal yci(t)') 
% 
figure 
plot(tyci,n) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Channel Noise Signal n(t)') 
% 
figure 
plot(tyci,yco) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Channel Output Signal yco(t)') 
% 
figure 
plot(tyci,cr) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Receiver Carrier Signal cr(t)') 
% 
figure 
plot(tyci,yd) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Demodulated Signal yd(t)') 
% 
figure 
plot(thL,hL,thL,hL,'o') 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Impulse Response Signal hL(t)') 
% 
figure 
plot(tr,xr) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Receiver Output Signal Signal xr')

% 
figure 
plot(tr,xr) 
grid 
xlabel('Time in Seconds') 
ylabel('Amplitude') 
title('Signals s(t) vs. xr(t)') 
hold on 
plot(t,s) 
hold off 
% 
figure 
plot(faxis_s,asFss) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Desired Signal s(t)') 
% 
figure 
plot(faxis_g,asFg) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Interference Signa g(t)') 
% 
figure 
plot(faxis_xm,asFxm) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Modulating Signal xm(t)') 
% 
figure 
plot(faxis_ct,asFct) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Carrier Signal ct(t)') 
% 
figure 
plot(faxis_xc,asFxc) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Modulated Signal xc(t)') 
% 
figure 
plot(faxis_hB,asFhB,faxis_hB,asFhB,'o') 
grid 

xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Band-Pass Filter Signal hB(t)') 
% 
figure 
plot(faxis_yci,asFyci) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Channel Input Signal yci(t)') 
% 
figure 
plot(faxis_n,asFn) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of AWGN Signal n(t)') 
% 
figure 
plot(faxis_yco,asFyco) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Channel Output Signal yco(t)') 
% 
figure 
plot(faxis_cr,asFcr) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Receiver Carrier cr(t)') 
% 
figure 
plot(faxis_yd,asFyd) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Demodulated Signal yd(t)') 
% 
figure 
plot(faxis_hL,asFhL,faxis_hL,asFhL,'o') 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Low-Pass Filter Signal hL(t)') 
% 
figure 
plot(faxis_xr,asFxr) 
grid 

xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectrum of Receiver Output Signal xr(t)') 
% 
figure 
plot(faxis_xr,asFxr) 
grid 
xlabel('Frequency in Hertz') 
ylabel('Magnitude') 
title('Spectra of Signals s(t) vs. xr(t)') 
hold on 
plot(faxis_s,asFss) 
hold off 
% 
sound(s,Fs) 
pause(2) 
sound(xm,Fs) 
pause(2) 
sound(xr,Fs) 
% 
clc 
disp('******************************') 
disp('Simulation ended succesfully') 
disp('******************************') 