clc;
clear;
close all;


%% Create chirped signal
clf, clearvars
disp('--Beginning signal generation (step 1/7)--')
f_dev = 2000*10^6; %[Hz], frequency deviation -f/2,...,f/2
DM = 10; %[pc cm^-3], dispersion measure
fc = 2000*10^6; %[Hz], centre frequency
 deltaT_max = 4.148808*(((fc-f_dev/2)*10^-9)^-2 - ((fc+f_dev/2)*10^-9)^-2)*DM; %[ms]
% deltaT_max = 4.148808*(((fc)*10^-9)^-2 - ((fc+f_dev)*10^-9)^-2)*DM; %[ms]

f_samp = 5000*10^6; %[Hz], sampling frequency
%T_tot = 10*10^-6; %[s], signal duration manually set
% T_tot = deltaT_max/10*10^-3; %[s], signal duration from DM
T_tot = 2 *10^-3; %[s], signal duration from DM

disp('--Generating sampling points (step 2/7)--')
t_samp = 1/f_samp; %resulting time sample
t = 0:t_samp:T_tot; %time vector

disp('--Generating amplitude modulation (step 3/7)--')
am = ones(1, length(t)); %amplitude modulation (constant)

disp('--Generating frequency modulations (step 4/7)--')
disp('-----Quad-----')
fm_quad = linspace(f_dev/2, f_dev/2, length(t));
for i=1:length(t)
    f_shift = (((t(i)*10^3)/(4.148808*DM) + 1/((fc+f_dev/2)*10^-9)^2)^(-1/2))*10^9 - (fc+f_dev/2);
%    f_shift = (((t(i)*10^3)/(4.148808*DM) + 1/((fc+f_dev/2)*10^-9)^2)^(-1/2))*10^9 - (fc-f_dev/2);
    fm_quad(i) = fm_quad(i) + f_shift;
end


disp('--Converting to phase (step 5/7)--')
fm = fm_quad; %Choose modulation
phase = 2.0*pi/f_samp*cumsum(fm); %convert to phase


disp('--Generating I_data (step 6/7)--')
I_data = sqrt(2) * am .* cos(phase);
disp('--Generating Q_data (step 7/7)--')
Q_data = sqrt(2) * am .* sin(phase);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%This is for generating waveform data%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQInfo.I_data         = I_data;
IQInfo.Q_data         = Q_data;
IQInfo.comment        = 'Quadratic Chirp';     % optional comment
IQInfo.clock          = f_samp;            % sample rate
IQInfo.filename       = strcat('quadratic_', num2str(f_dev/10^6), 'M_dm', num2str(DM), '_fc', num2str(fc*10^-6), 'M.wv');
IQInfo.no_scaling     = 0;


%status_gen = rs_generate_wave(0, IQInfo, 0, 1);
%status_vis = rs_visualize(f_samp, I_data, Q_data)

% figure;
% spectrogram(I_data+1i*Q_data, 256*4, 250*4, 256*4, 'centered', f_samp, 'yaxis');
sig = I_data+1i*Q_data;

 
 plotPlus(sig,t,f_samp,'input chirp');
 
 % plots
figure;
spectrogram(sig,256,128,256,f_samp,'yaxis');         % Time vs frequency     
title('Frequency vs Time: raw signal');