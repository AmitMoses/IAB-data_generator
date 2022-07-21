close all
clear all
clc

d = 0:0.1:5; % km
f = 3.5e9/1e6; %Mhz
Ptx = 30; % dBm
L = 32.44 + 20*log10(d) + 20*log10(f);
Psens = -L + 3 + 3 +  Ptx;
N0 = -174 + 10*log10(100e6);
SNR = Psens - N0;

figure()
plot(d,SNR)
grid on
grid minor
