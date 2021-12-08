%**************************Prof. D. Rodriguez***********************
%*********DSB-SC COMMUNICATIONS SYSTEMS AND GAUSSIAN NOISE**********
%*******************************************************************
clear all                            %Clear working space
close all
%*******************************************************************
%********************Wanted Signal Input/Output*********************
%*******************************************************************
[sig,Fs]=audioread('hw03gp01t05.wav');  %Get overall sampling frequency
save dsig.txt sig -ascii               %Save signal as dsig.txt file
load -ascii dsig.txt;                %Load ascii file dsig.txt
wsiz=length(dsig);                   %Get the length of the "dsig" signal
%*******************************************************************
%**************************Parameter Setting************************
%*******************************************************************
Ts=1/Fs;                          %Sampling time or temporal resolution
Nq=wsiz;                          %Number of samples of wanted signal
Vq=Nq*Ts;                         %Observation window or signal duration
fincq=1/Vq;                       %Frequency resolution of wanted signal
fkq=-(Fs/2):fincq:+(Fs/2)-fincq;  %Frequency axis of wanted signal
tq=0:Ts:Vq-Ts;                    %Time axis of wanted signal
%
fV=1800;                          %Cut-off frequency of wanted signal
Wn=(2*fV/Fs);                     %Normalized cut-off frequency
Lp=1/500;                         %Length proportionality factor
Nh=floor(Lp*Nq);                  %Low-pass filter length
hL=fir1(Nh,Wn);                   %Impulse response of low-pass filter
dsigL=conv(dsig,hL);              %Filtered wanted signal up to fV=4000 Hz
%
f1=4400;                             %Input interference frequency
f2=4900;                             %Input interference frequency
f3=5400;                             %Input interference frequency
g1=cos(2*pi*f1*tq);                  %Unwanted input signal 
g2=cos(2*pi*f2*tq);                  %Unwanted input signal
g3=cos(2*pi*f3*tq);                  %Unwanted input signal
g=(1/80)*g1+(1/120)*g2+(1/240)*g3;	 %Sum of unwanted signals
gmax=max(g);                         %Maximum value of sum
g=transpose((1/gmax)*g);             %Normalized interference signal
%********************************************************************
%*********************SIGNALS AND SYSTEMS MODELING*******************
fc=8640;                        %Carrier frequency
c=transpose((cos(2*pi*fc*tq))); %Carrier signal
wsig=dsigL(1:Nq,1);             %Wanted input signal
xm=wsig+g;                      %Sum of wanted and unwanted signals
xc=xm.*c;                       %Modulated wanted signal plus interference
nsig=randn(size(xm));           %Generation of AWGN signal: noise signal
nATT=sqrt(1/2000);              %Noise attenuation about 33dB SNR
nsig=nATT*nsig;                 %Addttenuated noise present in the channel
%*********************************************************************
%*************************BAND-PASS FILTER****************************
%Nb=;                           %Order of band-pass filter
%Wb=;                           %Vector of normilized cut-off frequencies
%hB=fir1(Nc,Wc,'bandpass');     %Impulse response of band-pass filter
%yci=conv(xc,hB);               %Removal of unwanted interference signals
%...                            %Other Additional instructions needed...
%*********************************************************************
yci=wsig.*c;            %Assuming bandpass filter worked correctly
yco=yci+nsig;           %Channel output signal with Gaussian noise
%********************************SPECTRA*****************************
%********************************************************************
fwsig=fft(wsig);                  %Spectrum of input signal (Fund. Reg.)
sfwsig=fftshift(fwsig);           %Spectrum of input signal (Princ. Reg.)
asfwsig=abs(sfwsig);              %Magnitude of spectrum of input signal
%
fxc=fft(xc);                  %Spectrum of input signal (Fund. Reg.)
sfxc=fftshift(fxc);           %Spectrum of input signal (Princ. Reg.)
asfxc=abs(sfxc);              %Magnitude of spectrum of input signal
%
fyci=fft(yci);                  %Spectrum of input signal (Fund. Reg.)
sfyci=fftshift(fyci);           %Spectrum of input signal (Princ. Reg.)
asfyci=abs(sfyci);              %Magnitude of spectrum of input signal
%
fyco=fft(yco);                  %Spectrum of output signal (Fund. Reg.)
sfyco=fftshift(fyco);           %Spectrum of output signal (Princ. Reg.)
asfyco=abs(sfyco);              %Magnitude of spectrum of output signal
%
fnsig=fft(nsig);                %Spectrum of noise signal (Fund. Reg.)
sfnsig=fftshift(fnsig);         %Spectrum of noise signal (Princ. Reg.)
asfnsig=abs(sfnsig);            %Magnitude of spectrum of noise signal
%******************************PLOTS*********************************
plot(tq,wsig);
axis([0 (length(wsig)*Ts) min(wsig)  max(wsig)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Wanted signal in the Time Domain')
grid
%
figure
plot(tq,xm)
axis([0 (length(xm)*Ts) min(xm)  max(xm)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Modulating Signal')
grid
%
figure
plot(tq,xc)
axis([0 (length(xc)*Ts) min(xc)  max(xc)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Modulated Signal')
grid
%
figure
plot(tq,yci)
axis([0 (length(yci)*Ts) min(yci) max(yci)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Channel Input Signal')
grid
%
figure
plot(tq,yco)
axis([0 (length(yco)*Ts) min(yco)  max(yco)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Channel Output Signal')
grid
%
figure
plot(fkq,asfwsig)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Wanted Signal')
grid
%
figure
plot(fkq,asfxc)
xlabel('Freq in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Modulated Signal')
grid
%
figure
plot(fkq,asfyci)
xlabel('Freq in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Channel Input Signal')
grid
%
figure
plot(fkq,asfyco)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Channel Output Signal')
grid
%
%*******************************OUTPUTS******************************
sound(wsig,Fs)  %Sound of Wanted Signal
pause(20)
sound(xc,Fs)    %Sound of Modulated Signal    
pause(20)
sound(yci,Fs)   %Sound of Cnannel's Input Signal 
pause(20)
sound(yco,Fs)   %Sound of Channel's Output Signal   

%********************************************************************
clc
disp('******************************')
disp('Simulation ended succesfully.')
disp('******************************')