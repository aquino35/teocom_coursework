clear all                            %Clear working space
close all
%%******Input Signals********
%
Fs=40000
n = randn(30720,1);
variance = sqrt(1/5000);
n= variance* n;
audiowrite('nois01.wav', n, Fs);
