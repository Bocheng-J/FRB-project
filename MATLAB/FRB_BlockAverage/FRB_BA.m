%% Script description
% Script: FRB capture system simulation

% Feature: Do digital signal processing for the input signal,find target
% FRB/chirp, save snapshot of the candidate

% Script outline: 1.generate input signal; 2. Downsampling; 3.Use PFB to
% generate spectrum; 4.Moving average; 5.Incoherent-dedispersion; 6.Sum
% channels to find trigger; 7.Triggering and save snapshots.


%% Script initialization
clc;
clear;
close all;
load('I_data_DM10_1000to500_fs5G');
load('Q_data_DM10_1000to500_fs5G');


%% General settings
noise_level = 1450;           % SNR -2: 1.585; SNR 0: 1.261; SNR 3 : 0.8603; SNR 5: 0.715
trigger_level = 10;             % in which SNR level to trigger (in dB) 
trigger_time = 0;               % create a variable to save the time of triggering, default: 0 
trigger_cnt = 0;               % indicate trigger or not, when trigger_cnt>= 3, we have the trigger
snapshot_flag = 0;              % snapshot flag
no_trigger_cnt = 0;            % flag for no trigger 
trigger_number = 0;             % record how many times it was triggered
shift_number = 1;               % check how many times of shift is performed
average_factor = 500;            % moving average span. (only allow odd number), better to be close to power of 2


%% generate input signal
fs = 5*1e9;            % 5G Hz Sampling rate
clk_cycle_number = 10000000; % how many clock cycles does the target signal have
T_tot = clk_cycle_number * 1/fs;       %[s], signal duration: 5 us
t_samp = 1/fs;         % resulting time sample
t_FRB = 0:t_samp:T_tot-1/fs; %time vector

DM = 10; %[pc cm^-3], dispersion measure
fc = 2000*10^6; %[Hz], centre frequency
sig_frb = I_data+1i*Q_data;
sig_original = real(sig_frb);                          % get real signal
f_start = 1000*1e6;   % start frequency of FRB
f_stop = 500*1e6;     % end frequency of FRB

length_extra =  length(t_FRB);                                
sig_extra = zeros(1,length_extra);                % add an extra length of zeros after chirp signal
sig = [sig_extra sig_original sig_extra];
t = 0:t_samp:(2*length_extra+length(sig_original)-1)*1/fs;                        % time 

% add noise
 input_SNR = 3;                                    % Set input SNR
sig = awgn(sig,input_SNR,'measured');             % Add noise
%sig = sig;

% convert floating point data to binary integer data
actual_analog_range = 1;                            % Vpp = 1
analog_level = (actual_analog_range/2)*sig;
sig = round(analog_level/(actual_analog_range/2)*(2^14)); % should be 2^15, but there is a bug when use 2^15 to generate test signal

% plots
figure;
spectrogram(sig,256,128,256,fs,'yaxis');         % Time vs frequency     
title('Frequency vs Time: raw signal');
plotPlus(sig,t,fs,'chirp');                       % plot in time domain & frequency domain



%% downsampling 
Sample_skip = 2;
sig_downsp =downsample(sig,Sample_skip);                   % downsampling input signal
t_1 = downsample(t,Sample_skip);                           % downsampling time
fs2 = fs/Sample_skip;                                      % New fs 
sig_length = length(sig)/fs2;                              % total time length of the input signal 
plotPlus(sig_downsp,t_1,fs2,'input chirp (After downsampling, fs = 2.5 GHz)');    


%% 32 channel PFB setting
M = 32;                                                   % channel number
D = 16;                                                   % decimation number
load('lp_125G_32ch_2');                                   % load prototype filter coefficients
h = reshape(lp_125G_32ch_2,M,[]);
h = round(h * 2^15 -1);

x = [sig_downsp ];                                        % get the downsampled signal as PFB input

% PFB initialization
inputDat = zeros(D,1);
inputDatBuf = zeros(M,size(h,2));
filtOutBuf = zeros(M,1);
chanOut = zeros(D,1);
chanOutBuf = zeros(D,round(length(x)/D));
flag = 0;

fs3 = fs2/D;                                              % New fampling frequency after PFB
t_2 = 0:1/fs3:(length(chanOutBuf)-1)*1/fs3;               


%% average time

% average time
t_sum = zeros(1,1);
t_average = zeros(1,floor(size(chanOutBuf,2)/average_factor));
for i = 1: size(chanOutBuf,2)
    if mod(i,average_factor)~=0 
      t_sum = t_2(i) + t_sum;
    else
        t_sum = t_2(i) + t_sum;
        t_average(i/average_factor) = t_sum/average_factor;
        t_sum = zeros(1,1);
    end
end

f_average = 1/(t_average(2)-t_average(1));



%% Preparation for incoherent dedispersion
% Dedispersion: delay higher frequencies to f_stop

frequencies = zeros(1,M/2);
for i = 1: M/2
    frequencies(i) = ((i-1)*(fs3/2/1e6));                    %calculate centre frequencies for each channel
end

% delay calculation: 
delays = zeros(1,M);
for i = 1:M/2
   delays(i) = 4.15*DM*((((fc/1e6+f_stop/1e6)/1000)^-2)-(((fc/1e6 + frequencies(i))/1000)^-2))/1000;   %FRB delay
end
delays_in_us = delays/1e-6;
delay_units = zeros(1,M);
for i = 1:M/2                                                              % convert time of delay to clock cycles of delay 
     delay_units(i) = round(delays(i)/(t_average(2)-t_average(1))); 
     if delay_units(i) <0
         delay_units(i)=0;
     end
end
delay_units = delay_units +1;


%% Pipeline 2.0
nof_samples = length(sig_downsp);
nof_clk = floor(nof_samples/D);
de_dispersed_data = zeros(D,1);
de_dispersed_data_plot = zeros(D,length(t_average));

% average filterbank
sum_average = zeros(D,1);
fifo_average = cell(D,1);
for k=1:D
    fifo_average{k,1} = zeros(1,average_factor);
end
chanOutBuf_average_plot = zeros(D,length(t_average));
fifo_delay = cell(D,1);
for i=1:D
    fifo_delay{i,1} = zeros(1,delay_units(i));
end
chanOutBuf_sum = zeros(D,1);
flux_plot = zeros(1,length(t_average));
trigger_SNR = zeros(1,length(t_average));                         % vector to store SNR of matched filter result
cnt = 1;


%% PFB loop implementation
for clk=1:nof_clk
    inputDat = fliplr(x((clk-1)*D+1:clk*D)).';
    inputDatBuf = [inputDatBuf(D+1:M,:);inputDatBuf(1:D,:)];
    inputDatBuf(1:D,:) = [inputDat inputDatBuf(1:D,1:size(h,2)-1)];
    for k=1:M
        filtOutBuf(k) = inputDatBuf(k,:)*h(k,:)';
    end
    if(flag == 0)
        flag = 1;
    else
        flag = 0;
        filtOutBuf = [filtOutBuf(D+1:M);filtOutBuf(1:D)];
    end
    filtOutBuf = filtOutBuf/(2^10);                 % we have 10bits truncation here in FPGA
    chanOut = fft(filtOutBuf,M)/(2^8);           % we have 8bits truncation here in FPGA    
    chanOut = abs(round(real(chanOut(1:D))));    % get the real part
    chanOutBuf(:,clk) = (chanOut);

        % average PFB out
    if mod(clk,average_factor)~=0 
      chanOutBuf_sum = chanOut + chanOutBuf_sum;
    else
        chanOutBuf_sum = chanOut + chanOutBuf_sum;
        chanOutBuf_average_data = floor(chanOutBuf_sum/average_factor);
        chanOutBuf_average_plot(:,clk/average_factor) = chanOutBuf_average_data;
        chanOutBuf_sum = zeros(D,1);
        
    % De-dispersion
    for k=1:D
        fifo_delay{k,1}(length(fifo_delay{k,1})+1) = chanOutBuf_average_data(k);
        fifo_delay{k,1}(1) = [];
        de_dispersed_data(k) = fifo_delay{k,1}(1);
    end    
    cnt = cnt+1;
    de_dispersed_data_plot(:,clk/average_factor) = de_dispersed_data;
    
    
        % Channel Sum
    flux = sum(de_dispersed_data);
    flux_plot(:,clk/average_factor) = flux;
  
    %% triggering

% triggering
trigger_SNR(clk/average_factor) = 10*log10((flux^2)/(noise_level)^2);
if trigger_SNR(clk/average_factor) >= trigger_level
    trigger_cnt = trigger_cnt +1;
    %trigger_number = trigger_number +1;              % check how many times of triggers we have 
    no_trigger_cnt = 0;
else 
    trigger_cnt = 0;
    no_trigger_cnt = no_trigger_cnt +1;
end

if trigger_cnt == 3                             % only when we have 3 consecutive triggers, we consider there is a valid trigger
        trigger_time = clk*(1/fs3);               % save the trigger time, unit: s
        candidate1 = snapshot(sig_downsp,T_tot,sig_length,0.2,trigger_time,fs2,snapshot_flag);  % get candidate
        str = sprintf('-----Possible candidate detected at %d s-----',trigger_time);
        disp(str);
        %figure;
        %spectrogram(candidate1,256,128,256,fs2,'yaxis');         % Time vs frequency     
        %title('Frequency vs Time: candidate snapshot');
end

if no_trigger_cnt == 3
    trigger_time = 0;
    disp('-----No candidate detected-----');
end  

    end

% show percentage of completion
if mod(clk,5/100*nof_clk)==0
percentage = clk/nof_clk * 100;
str_p = sprintf('-----Loop completion: %d%%-----',round(percentage));
disp(str_p);
end

end


%% plot results
figure;
for i = 1:D
   plot(t_average*1e6,de_dispersed_data_plot(i,:));
   hold on;
end 
title('Arrival time of 32 channel signals: averaged de-dispersed signal');
xlabel('Time/us'); 
ylabel('Magnitude');

figure;
for k=1:D
    plot(t_2,chanOutBuf(k,:)); hold on;
end
xlabel('Time/us');ylabel('digital level');
title('polyphase filterbank output');

figure;
plot(t_average*1e6,(flux_plot));
str = sprintf('Dedispersed flux when input SNR = %d dB', input_SNR);
title(str);
xlabel('Time/us'); 
ylabel('Magnitude');

figure;
plot(t_average,(trigger_SNR));
str2 = sprintf('trigger SNR when input SNR = %d dB', input_SNR);
title(str2);
xlabel('Time/us'); 
ylabel('SNR/dB');

% plot original signal (dispersed)
figure;
for i = 1:D
   plot(t_average*1e6,(chanOutBuf_average_plot(i,:)));
   hold on;
end
title('Arrival time of 32 channels: averaged dispersed signal');
xlabel('Time/us'); 
ylabel('Magnitude');


% Final Plot: dispersion vs de-dispersion
figure;
subplot(2,2,1);
% plot dispersed sum

% plot(t_2*1e6,(chanOutBuf_sum));
% title('Dispersed sum');
% xlabel('Time/us'); 
% ylabel('Magnitude');
% xlim([0 15]);
% ylim([0 1.8*1e9]);

% subplot(2,2,2);
% plot(t_3*1e6,ddp_chanOutBuf_sum);
% str = sprintf('De-dispersed sum');
% title(str);
% xlabel('Time/us'); 
% ylabel('Magnitude');
% xlim([0 22]);
% ylim([0 1.8*1e9]);

subplot(2,2,3);
imagesc(t_average*1e6,frequencies,chanOutBuf_average_plot(1:16,:));
set(gca,'YDir','normal') ;
title('Dispersed signal');
xlabel('Time/us'); 
ylabel('Frequency/MHz');
%xlim([0 10]);

subplot(2,2,4);
imagesc(t_average*1e6,frequencies,de_dispersed_data_plot);
set(gca,'YDir','normal') ;
title('De-dispersed signal');
xlabel('Time/us'); 
ylabel('Frequency/MHz');
%xlim([0 10]);

%% plot dispersion vs dedispersion separately
% figure;
% imagesc(t_3*1e6,frequencies,abs(dedisp_sig));
% set(gca,'YDir','normal') ;
% title('De-dispersed signal');
% xlabel('Time/us'); 
% ylabel('Frequency/MHz');
% xlim([0 22]);
% 
% 
% figure;
% plot(t_3*1e6,ddp_chanOutBuf_sum);
% str = sprintf('De-dispersed sum');
% title(str);
% xlabel('Time/us'); 
% ylabel('Magnitude');
% xlim([0 22]);
% ylim([0 1.8*1e9]);
% 
% figure;
% imagesc(t_2*1e6,frequencies,abs(chanOutBuf_average(1:16,:)));
% set(gca,'YDir','normal') ;
% title('Dispersed signal');
% xlabel('Time/us'); 
% ylabel('Frequency/MHz');
% xlim([0 15]);
% 
% figure;
% plot(t_2*1e6,(chanOutBuf_sum));
% title('Dispersed sum');
% xlabel('Time/us'); 
% ylabel('Magnitude');
% xlim([0 15]);
% ylim([0 1.8*1e9]);


%% check particulate channel result (if necessary)
% 
%  figure;
%  subplot(2,3,1);
%  plot(t_average*1e6,(de_dispersed_data_plot(14,:)));       % check channel 14
%  title('Check channel 14');
%  
%   subplot(2,3,2);
%  plot(t_average*1e6,(de_dispersed_data_plot(13,:)));       % check channel 13
%  title('Check channel 13');
%   
%    subplot(2,3,3);
%  plot(t_average*1e6,abs(de_dispersed_data_plot(12,:)));       % check channel 12
%  title('Check channel 12');
%   
%    subplot(2,3,4);
%  plot(t_average*1e6,(de_dispersed_data_plot(11,:)));       % check channel 11
%  title('Check channel 11');
%   
%    subplot(2,3,5);
%  plot(t_average*1e6,(de_dispersed_data_plot(10,:)));       % check channel 10
%  title('Check channel 10');
%   
%    subplot(2,3,6);
%  plot(t_average*1e6,(de_dispersed_data_plot(9,:)));       % check channel 9
%  title('Check channel 9');
%  
