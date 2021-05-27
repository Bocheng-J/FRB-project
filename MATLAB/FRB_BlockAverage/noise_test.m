clc;
clear;
close all;


fs = 5*1e9;            % 5G Hz Sampling rate
clk_cycle_number = 32000; % how many clock cycles does the target signal have
T_tot = clk_cycle_number * 1/fs;       %[s], signal duration: 5 us
t_samp = 1/fs;         % resulting time sample
t_chirp = 0:t_samp:T_tot-1/fs; %time vector

% generate noise
noise = normrnd(0,0.01257,[1,length(t_chirp)]);

figure;
plot(noise);

% generate chirp
f_start = 1000*1e6;   % start frequency of chirp
f_stop = 500*1e6;     % end frequency of chirp
chirp_sig = chirp(t_chirp,f_start,T_tot,f_stop);     % generate chirp
figure;
plot(t_chirp,chirp_sig);
max_chirp = max(chirp_sig)
min_chirp = min(chirp_sig)


% normalization
[PN]=tramnmx(chirp_sig,min_chirp,max_chirp);
figure;
plot(t_chirp,PN);