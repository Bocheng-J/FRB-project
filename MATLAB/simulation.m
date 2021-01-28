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
phase=0;                        % phase of Vs1
A=1;                            % amplitude of Vs1
Vs1 = A*sin(2*pi*f*t+phase);
x1 = max(Vs1)./sqrt(2);         % root mean square of Vs1
Vs1 = Vs1./x1;                  % normalize Vs1

% plot
figure(1);
plot(t.*1000,Vs1);              % times 1000 to change the unit from s to ms
str_title2 = sprintf('Original signal Vs1 with frequency %dHz',f);
title(str_title2);
%axis([-inf,+inf,-1,+1]);
xlabel('t/ms');
ylabel('amplitude/V');


%% Amplifier A1
% parameter
G1 = 10;                        % voltage gain of A1 in dB
Ph1 = 0;                      % phase introduced by A1
Vn1 = 40;                        % noise introduced by A1 in dB
Vs1_Power = sum(abs(Vs1).^2);
Vs1_Power = 20*log10(Vs1_Power);
SNR = Vs1_Power-Vn1;

% convert dB to linear
G1_Lin = 10^(G1/20);

% apply to Vs1
% A1_Out = awgn(Vs1,Vn1);
% A1_Out = A1_Out.* G1_Lin;

A1_Out = Vs1.*G1_Lin;
A1_Out = awgn(A1_Out,SNR);

% plot
figure(2);
plot(t.*1000,A1_Out);            % times 1000 to change the unit from s to ms
str_title2 = sprintf('Vs1 after amplifier A1, with G = %d dB, Ph1 = %d dB and Vn1 = %d dB', G1, Ph1, Vn1);
title(str_title2);
%axis([-inf,+inf,-1,+1]);
xlabel('t/ms');
ylabel('amplitude/V');





