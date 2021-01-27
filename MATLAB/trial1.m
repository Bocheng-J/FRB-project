clear;
clc;

N = 1000;                      % sample points
tscale= 0.05;                  % x axis range 
dt = tscale/N;                 % time slot

A = 20;                         % amplitude
f = 100;                       % frequency
phase = 0;                     % phase
t = 0:dt:(N-1)*dt;             % time range
V_s1 = A*sin(2*pi*f*t+phase);    % define a sine wave

%V_s1 = V_s1./max(abs(V_s1));                     % set amplitude rms to always 1
x1 = max(V_s1)./sqrt(2);
V_s1 = V_s1./x1;

% plot
figure(1);
plot(t.*1000,V_s1); %t times 1000 change unit from s to ms
title('100Hz sine wave');
%axis([-inf,+inf,-1,+1]);          % set x&y axis range
xlabel('t/ms');
ylabel('amplitude/V');


A1 = amplifier;
A1 = amplifier('Gain', 1);

c = circuit([V_s1 A1]);


