fs = 256e6;
N = 1024;
noise1 = rand(1,N)*0.1;
noise2 = randn(1,N)*0.1;

figure;
plot(noise1);hold on;
plot(noise2);hold on;
hold off;

fax_Hz = (0:N-1)*(fs/N);
N_2 = ceil(N/2);
fftnoise1 = fft(noise1);
fftnoise2 = fft(noise2);

figure;
plot(fax_Hz(1:N_2),abs(fftnoise1(1:N_2))); hold on;
plot(fax_Hz(1:N_2),abs(fftnoise2(1:N_2))); hold on;
hold off;
