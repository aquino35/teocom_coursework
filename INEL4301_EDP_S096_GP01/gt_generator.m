Fs = 40000;
Ts = 1/Fs;
f1 = 4800;
f2 = 5200;
f3 = 5800;
N = 30720;
V = N*Ts;
t = 0:Ts:V-Ts;

%   Signals 

g1 = cos(2*pi*f1*t);
g2 = cos(2*pi*f2*t);
g3 = cos(2*pi*f3*t);
g=(1/80)*g1+(1/120)*g2+(1/240)*g3;

audiowrite('uwtd01.wav', g, Fs);
