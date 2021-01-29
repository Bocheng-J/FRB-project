close all;
clear;
clc; 

%% Generate signal Vs1
% parameter
N=1000;                         % number of samples
tscale= 0.05;                   % the length of time displayed on the x-axis, unit in second
dt=tscale/N;                    % time interval of each sampling point
t=0:dt:(N-1)*dt;        
f=100;                          % frequency of Vs1 in Hz
A=10;                            % amplitude of Vs1
Vs1 = A*sin(2*pi*f*t);
% x1 = max(Vs1)./sqrt(2);         % root mean square of Vs1
% Vs1 = Vs1./x1;                  % normalize Vs0

% plot
figure(1);
plot(t.*1000,Vs1);              % times 1000 to change the unit from s to ms
str_title1 = sprintf('Input Vs1, f = %d Hz, A = %d V',f, A);
title(str_title1);
%axis([-inf,+inf,-1,+1]);
xlabel('t/ms');
ylabel('amplitude/V');
hold on;


%% Amplifier A1
% parameter
G1 = 10;                        % voltage gain of A1 in dB
Ph1 = 45;                      % phase introduced by A1
Vn1 = 1;                        % noise introduced by A1 in dB

Ph1 = Ph1*pi/180; 
Vs1 = A*sin(2*pi*f*t+Ph1);

sigPower_Lin = sum(Vs1.^2)/N;       % signal power in watt
sigPower = 10*log10(sigPower_Lin);  % signal power in dB
noisePower_Lin = 10^(Vn1/10);           % noise power in watt
noise = sqrt(noisePower_Lin)*randn(1,N);    % generate noise signal

A1_Out = Vs1+noise;                 % add noise to signal

G1_Lin = 10^(G1/20); % convert voltage gain from dB to linear

A1_Out = A1_Out.* G1_Lin;           % Output of amplifier A1

% plot
%figure(2);
plot(t.*1000,A1_Out);            % times 1000 to change the unit from s to ms
% str_title2 = sprintf('Vs1 after amplifier A1, with G = %d dB, Ph1 = %d and Vn1 = %d dB', G1, Ph1, Vn1);
% title(str_title2);
% xlabel('t/ms');
% ylabel('amplitude/V');



%% FFT
y = fft(A1_Out);
fs = 1000;                               % sample frequency 100Mhz
F_range = (0:N-1)*(fs/N);
power = abs(y).^2/N;


figure(3)
plot(F_range,power);
xlabel('Frequency');
ylabel('Power');







