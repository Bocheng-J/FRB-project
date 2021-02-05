clc; clear all
close all;
clear all;
clc;

fs = 256e6;
t = 0:1/fs:0.000004-1/fs;
A = 3;
Ar = 2;
Vs1 = [zeros(1,1000), A*chirp(t,50e6,0.000004,40e6)];

figure;
t1 = 0:1/fs:0.000004*2-1/fs;
plot(Vs1);
N = length(Vs1);
N_2 = ceil(N/2);
fax_Hz = (0:N-1)*(fs/N);
fftVs1 = fft(Vs1);
figure;
plot(fax_Hz(1:N_2)/1e6,abs(fftVs1(1:N_2)));


refVs1 = Ar*chirp(t,50e6,0.000004,40e6);
energy_refVs1 = sum(refVs1.^2);
refVs1 = refVs1/sqrt(energy_refVs1);
figure;
plot(refVs1);


Vs1 = Vs1+0.1*randn(1,length(Vs1));
corr = conv(Vs1,fliplr(refVs1));
%corr = conv(fliplr(refVs1),Vs1);
figure; plot(corr);title('corr');

[tmp, Tmax] = max(corr);
Tx_hat = Tmax-length(refVs1);


N = length(corr);
N_2 = ceil(N/2);
fax_Hz = (0:N-1)*(fs/N);
fftcorr = fft(corr);
figure;
plot(fax_Hz(1:N_2)/1e6,abs(fftcorr(1:N_2)));




