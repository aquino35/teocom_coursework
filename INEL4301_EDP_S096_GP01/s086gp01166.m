%*******************************************************************
%********************Wanted Signal Input/Output*********************
%*******************************************************************
[st, Fss] = audioread('S086GP01_EDP_INP1s/dtmf01.wav');       % s(t) signal
[gt, Fsg] = audioread('S086GP01_EDP_INP1g/uwtd01.wav');       % g(t) signal
[nt,Fsn]=audioread('S086GP01_EDP_INP1n/nois01.wav');         % n(t) signal


%
%%****Processing Input Signals*****
%

save st01.txt st -ascii                   	%Save signal as s.txt file
load -ascii st01.txt;                     	%Load ascii file s.txt
ssiz=length(st);                        	%Get the length of the "s" signal

save gt01.txt gt -ascii                   	%Save signal as g.txt file
load -ascii gt01.txt;                     	%Load ascii file g.txt
gsiz=length(gt);                       	%Get the length of the "g" signal

save nt01.txt nt -ascii                   	%Save signal as n.txt file
load -ascii nt01.txt;                      %Load ascii file n.txt
nsiz=length(nt);                         %Get the length of the "n" signal

%*******************************************************************
%**************************Parameter Settings***********************
%*******************************************************************
Fs=40000;                               %Sampling rate
Ts=1/Fs;
fc=8000;                                %Carrier frequency
N=ssiz;
V=N*Ts;
fincq=1/V;                       %Frequency resolution of wanted signal
fkq=-((Fs+1)/2):fincq:+((Fs+1)/2)-fincq;  %Frequency axis of wanted signal

t=0:Ts:V-Ts;
M=100;                                   %Order of the filters
Fv=2000;
 

%% 
%% TRANSMITTER SYSTEM
%% 
c=cos(2*pi*fc*t);                       %Carrier signal
xm=st+gt;
c1=transpose(c);
xc=xm.*c1;                               %Modulated wanted signal plus interference

%% Part #1 task #1
%% 
%% BAND_PASS FILTER DESIGN (hB)
%%
fm=Fv;                                    %Cuttoff frequency
fUpper = fc + fm;
fLower = fc - fm;
Wnl = 2*(fLower)/Fs;
Wnl2 = 2*(fc)/Fs;
Wnu = 2*(fUpper)/Fs;
Fb = [Wnl2 Wnu];                             %band limits
hB = fir1(M-1,Fb,'bandpass');               %filter design
yci=conv(xc,hB,'same');                     %Removal of unwanted interference signals

ycisiz=length(yci);
Nyci=ycisiz;
Vyci=Nyci*Ts;
fincqyci=1/Vyci;                       %Frequency resolution of channel input signal
fkqyci=-((Fs+1)/2):fincqyci:+((Fs+1)/2)-fincqyci;  %Frequency axis of channel input signal


%% Part #1 task #2
%% 
%% LOW-PASS FILTER DESIGN (Tc) 
%%
B = 12000;
Wn=2*B/Fs;                              %Normilized frequency
th=0:Ts:M*Ts-Ts;                        %Time axis to plot impulse response signal
hc=fir1(M-1,Wn);                        %Impulse response signal
yc=conv(yci,hc,'same');

ycsiz=length(yc);
Nyc=ycsiz;
Vyc=Nyc*Ts;
fincqyc=1/Vyc;                       %Frequency resolution of channel filter output signal
fkqyc=-((Fs+1)/2):fincqyc:+((Fs+1)/2)-fincqyc;  %Frequency axis of channel filter output signal

%% Part #1 task #3
%% 
%% ADDITION OF NOISE SIGNAL (Tn)
yco = yc+nt;

ycosiz=length(yco);
Nyco=ycosiz;
Vyco=Nyco*Ts;
fincqyco=1/Vyco;                       %Frequency resolution of channel output signal
fkqyco=-((Fs+1)/2):fincqyco:+((Fs+1)/2)-fincqyco;  %Frequency axis of channel output signal

xmSNR = snr(st,gt);
ycPdB = pow2db(bandpower(yc));
nPdB = pow2db(bandpower(nt));
SNR = snr(yc,nt);

nsig_avg = mean(nt)                    %Mean of noise 
nsig_var = var(nt)                      %Variance of noise 

%% Part #1 task #4
%% 
%% BANDPASS FILTER (Tr)
%%
fR = [Wnl Wnu];                         %band limits
hR = fir1(15,fR,'bandpass');           %filter design
ydi = conv(yco,hR,'same');                       %Removal of unwanted interference signals

ydisiz = length(ydi);
Nydi = ydisiz;
Vydi = Nydi*Ts;
fincqydi = 1/Vydi;                       %Frequency resolution of demodulator input signal
fkqydi = -((Fs+1)/2):fincqydi:+((Fs+1)/2)-fincqydi;  %Frequency axis of demodulator input signal

ydiPdB = pow2db(bandpower(ydi));
SNR2 = snr(ydi,nt)

%***** Demodulator Td*******
ydo = ydi.*(2*c1);

ydosiz=length(ydo);
Nydo=ydosiz;
Vydo=Nydo*Ts;
fincqydo=1/Vydo;                       %Frequency resolution of demodulator output signal
fkqydo=-((Fs+1)/2):fincqydo:+((Fs+1)/2)-fincqydo;  %Frequency axis of demodulator output signal

%% Part #1 Task #5
%% 
%% LOW_PASS FILTER DESIGN (Tr)
%%
Wn=2*fm/Fs;         %Normalized cut-off frequency 
M=50;              %Order of the impulse response signal hL(t) 
thL=0:Ts:M*Ts-Ts;   %Time axis to plot impulse response signal 
hL=fir1(M-1,Wn);    %Impulse response signal 
xr=conv(ydo,hL);             %Convolution operation for filtering 

audiowrite("S086GP01_EDP_OUP166/dtmf01.wav",xr,Fs);

Dr=((M-1)/2)*Ts
xrsiz=length(xr);
Nxr=xrsiz;
Vxr=Nxr*Ts;
fincqxr=1/Vxr;                       %Frequency resolution of recieved signal
fkqxr=-((Fs+1)/2):fincqxr:+((Fs+1)/2)-fincqxr;  %Frequency axis of recieved signal

%%
%% ********PLOTS*********
%%
tyci=0:Ts:(length(yci)-1)*Ts;
tyc=0:Ts:(length(yc)-1)*Ts;
txr=0:Ts:(length(xr)-1)*Ts;


% % Time Domain

plot(t,st);
axis([0 (length(st)*Ts) min(st)  max(st)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Wanted signal s(t)')
grid

% %

figure
plot(t,gt);
axis([0 (length(gt)*Ts) min(gt)  max(gt)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Interference Signal g(t)')
grid

% %

figure
plot(t,xm)
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Modulating Signal (xm(t) = s(t) + g(t))')
grid

% %

figure
plot(t,xc)
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Modulated Signal xc(t)')
grid

% % 

figure
plot(tyci,yci)
axis([0, (length(yci)*Ts), min(yci), max(yci)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Channel Input Signal yci(t)')
grid

% %

figure
plot(tyci,yc)
axis([0 (length(yc)*Ts) min(yc)  max(yc)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Channel Filter Output Signal yc(t)')
grid

% %

figure
plot(tyc,nt)
axis([0 (length(nt)*Ts) min(nt)  max(nt)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Channel Noise Signal n(t)') 
grid

% %

figure
plot(tyc,yco)
axis([0 (length(yco)*Ts) min(yco)  max(yco)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Channel Output Signal yco(t)')
grid

% %

figure
plot(tyc,ydi)
axis([0 (length(ydi)*Ts) min(ydi)  max(ydi)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Demodulator Input Signal ydi(t)')
grid

% %

figure
plot(tyc,ydo)
%axis([0, (length(ydo)*Ts), min(ydo),  max(ydo)]);	****
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Demodulator Output Signal ydo(t)')
grid

% %

figure
plot(txr,xr)
%axis([0 (length(xr)*Ts) min(xr)  max(xr)]);    
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Received Signal xr(t)')
grid

% % Spectrums

fw=fft(st);                  %Spectrum of wanted signal (Fund. Reg.)
sfs=fftshift(fw);           %Spectrum of wanted signal (Princ. Reg.)
asfs=abs(sfs);              %Magnitude of spectrum of wanted signal

%

fg=fft(gt);                  %Spectrum of interference signal (Fund. Reg.)
sfg=fftshift(fg);           %Spectrum of interference signal (Princ. Reg.)
asfg=abs(sfg);              %Magnitude of spectrum of wanted signal

%

fxm=fft(xm);                  %Spectrum of input signal (Fund. Reg.)
sfxm=fftshift(fxm);           %Spectrum of input signal (Princ. Reg.)
asfxm=abs(sfxm);              %Magnitude of spectrum of input signal

%

fxc=fft(xc);                  %Spectrum of modulated signal (Fund. Reg.)
sfxc=fftshift(fxc);           %Spectrum of modulated signal (Princ. Reg.)
asfxc=abs(sfxc);              %Magnitude of spectrum of modulated signal

%

fyci=fft(yci);                  %Spectrum of channel input signal (Fund. Reg.)
sfyci=fftshift(fyci);           %Spectrum of channel input signal (Princ. Reg.)
asfyci=abs(sfyci);              %Magnitude of spectrum of channel input signal

%

fyc=fft(yc);                  %Spectrum of channel filter output signal (Fund. Reg.)
sfyc=fftshift(fyc);           %Spectrum of channel filter output signal (Princ. Reg.)
asfyc=abs(sfyc);              %Magnitude of spectrum of channel filter output signal

%

fn=fft(nt);                  %Spectrum of AWGN signal (Fund. Reg.)
sfn=fftshift(fn);           %Spectrum of AWGN signal (Princ. Reg.)
asfn=abs(sfn);              %Magnitude of spectrum of wanted signal

%

fyco=fft(yco);                  %Spectrum of channel output signal (Fund. Reg.)
sfyco=fftshift(fyco);           %Spectrum of channel output signal (Princ. Reg.)
asfyco=abs(sfyco);              %Magnitude of spectrum of channel output signal

%

fydi=fft(ydi);                  %Spectrum of demodulator input signal (Fund. Reg.)
sfydi=fftshift(fydi);           %Spectrum of demodulator input signal (Princ. Reg.)
asfydi=abs(sfydi);              %Magnitude of spectrum of demodulator input signal

%

fydo=fft(ydo);                  %Spectrum of demodulator output signal (Fund. Reg.)
sfydo=fftshift(fydo);           %Spectrum of demodulator output signal (Princ. Reg.)
asfydo=abs(sfydo);              %Magnitude of spectrum of demodulator output signal

%

fxr=fft(xr);                  %Spectrum of recieved signal (Fund. Reg.)
sfxr=fftshift(fxr);           %Spectrum of recieved signal (Princ. Reg.)
asfxr=abs(sfxr);              %Magnitude of spectrum of recieved signal

%

figure
plot(fkq,asfs)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Wanted Signal S(f)')
grid

% %

figure
plot(fkq,asfg)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Interference Signal G(f)')
grid

% %

figure
plot(fkq,asfxm)
xlabel('Freq in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Modulating Signal Xm(f)')
grid

% %

figure
plot(fkq,asfxc)
xlabel('Freq in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Modulated Signal Xc(f)')
grid

% %

asfyci=transpose(asfyci);
figure
plot(fkqyci,asfyci)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Channel Input Signal Yci(f)')
grid

% %

asfyc=transpose(asfyc);
figure
plot(fkqyc,asfyc)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Channel Filter Output Signal Yc(f)')
grid

% %

figure
plot(fkq,asfn)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of AWGN Signal N(f)')
grid

% %

asfyco=transpose(asfyco);
figure
plot(fkqyco,asfyco)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Channel Output Signal Yco(f)')
grid

% %

asfydi=transpose(asfydi);
figure
plot(fkqydi,asfydi)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Demodulator Input Signal Ydi(f)')
grid

% %

figure
plot(fkqydo,asfydo)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Demodulator Output Signal Ydo(f)')
grid

% %

figure
plot(fkqxr,asfxr)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Recieved Signal Xr(f)')
grid

% %

sound(st, Fs);
pause(3);
sound(xr, Fs);