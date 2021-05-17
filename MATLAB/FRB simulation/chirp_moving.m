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
% load('I_data_DM20_fs5G');
% load('Q_data_DM20_fs5G');


%% General settings
noise_level = 0.8603;           % SNR -2: 1.585; SNR 0: 1.261; SNR 3 : 0.8603; SNR 5: 0.715
trigger_level = 10;             % in which SNR level to trigger (in dB) 
trigger_time = 0;               % create a variable to save the time of triggering, default: 0 
trigger_flag = 0;               % indicate trigger or not, when trigger_flag~= 0, we have the trigger
trigger_number = 0;             % record how many times it was triggered
shift_number = 1;               % check how many times of shift is performed
average_factor = 55;            % moving average span. (only allow odd number)


%% generate input signal
fs = 5*1e9;            % 5G Hz Sampling rate
T_tot = 5*10^-6;       %[s], signal duration: 5 us
t_samp = 1/fs;         % resulting time sample
t_chirp = 0:t_samp:T_tot-1/fs; %time vector

f_start = 1000*1e6;   % start frequency of chirp
f_stop = 600*1e6;     % end frequency of chirp
chirp_sig = chirp(t_chirp,f_start,T_tot,f_stop);     % generate chirp

length_extra =  length(t_chirp);                                
sig_extra = zeros(1,length_extra);                % add an extra length of zeros after chirp signal
sig = [chirp_sig sig_extra];
t = 0:t_samp:(length_extra+length(chirp_sig)-1)*1/fs;                        % time 

% add noise
input_SNR = 3;                                    % Set input SNR
sig = awgn(sig,input_SNR,'measured');             % Add noise

% plots
%myspectro(sig,128,f_start,f_stop);                % Time vs frequency
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
x = [sig_downsp ];                                        % get the downsampled signal as PFB input

% PFB initialization
inputDat = zeros(D,1);
inputDatBuf = zeros(M,size(h,2));
filtOutBuf = zeros(M,1);
chanOut = zeros(M,1);
chanOutBuf = zeros(M,round(length(x)/D));
flag = 0;
mm = 1;

fs3 = fs2/D;                                              % New fampling frequency after PFB
t_2 = 0:1/fs3:(length(chanOutBuf)-1)*1/fs3;               


%% Preparation for incoherent dedispersion
% Dedispersion: delay higher frequencies to f_stop

frequencies = zeros(1,M/2);
for i = 1: M/2
    frequencies(i) = ((i-1)*(fs3/2/1e6));                    %calculate centre frequencies for each channel
end

% delay calculation: 
delays = zeros(1,M);
chirp_coefficient = (f_stop/1e6 -f_start/1e6)/T_tot; 
for i = 1:M/2
 %   delays(i) =
 %   4.15*DM*(((2000/1000)^-2)-((frequencies(i)/1000)^-2))/1000;   %FRB delay
 delays(i)= (T_tot)-((frequencies(i)-1000)/chirp_coefficient);             % calculate the delay for linear chirp
end
delays_in_us = delays/1e-6;
delay_units = zeros(1,M);
for i = 1:M/2                                                              % convert time of delay to clock cycles of delay 
     delay_units(i) = round(delays(i)/(t_2(2)-t_2(1))); 
     if delay_units(i) <0
         delay_units(i)=0;
     end
end
delay_units = delay_units +1;

% create shift registers to implement delay
scale = 1;
dedisp_sig_16 = zeros(1,round(scale*(delay_units(16)))+size(chanOutBuf,2));
dedisp_sig_15 = zeros(1,round(scale*(delay_units(15)))+size(chanOutBuf,2));
dedisp_sig_14 = zeros(1,round(scale*(delay_units(14)))+size(chanOutBuf,2));
dedisp_sig_13 = zeros(1,round(scale*(delay_units(13)))+size(chanOutBuf,2));
dedisp_sig_12 = zeros(1,round(scale*(delay_units(12)))+size(chanOutBuf,2));
dedisp_sig_11 = zeros(1,round(scale*(delay_units(11)))+size(chanOutBuf,2));
dedisp_sig_10 = zeros(1,round(scale*(delay_units(10)))+size(chanOutBuf,2));
dedisp_sig_9 = zeros(1,round(scale*(delay_units(9)))+size(chanOutBuf,2));
dedisp_sig_8 = zeros(1,round(scale*(delay_units(8)))+size(chanOutBuf,2));
dedisp_sig_7 = zeros(1,round(scale*(delay_units(7)))+size(chanOutBuf,2));
dedisp_sig_6 = zeros(1,round(scale*(delay_units(6)))+size(chanOutBuf,2));
dedisp_sig_5 = zeros(1,round(scale*(delay_units(5)))+size(chanOutBuf,2));
dedisp_sig_4 = zeros(1,round(scale*(delay_units(4)))+size(chanOutBuf,2));
dedisp_sig_3 = zeros(1,round(scale*(delay_units(3)))+size(chanOutBuf,2));
dedisp_sig_2 = zeros(1,round(scale*(delay_units(2)))+size(chanOutBuf,2));
dedisp_sig_1 = zeros(1,round(scale*(delay_units(1)))+size(chanOutBuf,2));
dedisp_sig = zeros(D,size(dedisp_sig_16,2));
dedisp_sum = zeros(1,length(dedisp_sig_16));                            % vector to store the de-dispersed sum of channels
t_3 = 0:1/fs3:(size(dedisp_sig,2)-1)*1/fs3;                             % time axis for the dedispersed signal


%% pfb implementation
MF = zeros(1,length(dedisp_sig_16));                                  % vector to store matched filter result
trigger_SNR = zeros(1,length(dedisp_sig_16));                         % vector to store SNR of matched filter result

% initialization for moving average
p = (average_factor-1)/2;
q = p+1; 
chanOutBuf_average = zeros(M,length(chanOutBuf));
count = 1;

% dynamic plots initialization
% %warning: enable dynamic plots would make the simulation run in a
% %extremely slow speed, so do not enable if really necessary
% pfb_out = 1;                                                         % dynamic plot 1
% fHz = 0;
% figure;
% subplot(4,1,1);
% p1 = plot(fHz,pfb_out,'EraseMode','background','MarkerSize',5);
% 
% pfb_out2 = 1;                                                       % dynamic plot 2
% subplot(4,1,2);
% p2 = plot(fHz,pfb_out2,'EraseMode','background','MarkerSize',5);
% 
% pfb_dedisp_out = 1;                                                   % dynamic plot 3
% subplot(4,1,3);
% p3 = plot(t_2*1e6,pfb_dedisp_out,'EraseMode','background','MarkerSize',5);


% PFB implementation
for n=1:D:length(x)-M
    inputDat = fliplr(x(n:n+D-1)).';
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
    chanOut = M*ifft(filtOutBuf);
    chanOutBuf(:,mm) = (chanOut);
    
    % moving average PFB output
input = abs(chanOutBuf);
if (mm<= average_factor) 
    if (count ==1 )
        temp_vec = zeros(M,mm);
            for j = 1:mm
                temp_vec(:,j) = input(:,j);
            end   
        temp_addition = cumsum(temp_vec,2);
        chanOutBuf_average(:,(mm+1)/2) = temp_addition(:,mm)/(mm);
        count = count + 1;
        
     %% shift
     dedisp_sig_16 = ArrayDisplacement(dedisp_sig_16,2,1,2);
     dedisp_sig_15 = ArrayDisplacement(dedisp_sig_15,2,1,2);
     dedisp_sig_14 = ArrayDisplacement(dedisp_sig_14,2,1,2);
     dedisp_sig_13 = ArrayDisplacement(dedisp_sig_13,2,1,2);
     dedisp_sig_12 = ArrayDisplacement(dedisp_sig_12,2,1,2);
     dedisp_sig_11 = ArrayDisplacement(dedisp_sig_11,2,1,2);
     dedisp_sig_10 = ArrayDisplacement(dedisp_sig_10,2,1,2);
     dedisp_sig_9 = ArrayDisplacement(dedisp_sig_9,2,1,2);
     dedisp_sig_8 = ArrayDisplacement(dedisp_sig_8,2,1,2);
     dedisp_sig_7 = ArrayDisplacement(dedisp_sig_7,2,1,2);
     dedisp_sig_6 = ArrayDisplacement(dedisp_sig_6,2,1,2);
     dedisp_sig_5 = ArrayDisplacement(dedisp_sig_5,2,1,2);
     dedisp_sig_4 = ArrayDisplacement(dedisp_sig_4,2,1,2);
     dedisp_sig_3 = ArrayDisplacement(dedisp_sig_3,2,1,2);
     dedisp_sig_2 = ArrayDisplacement(dedisp_sig_2,2,1,2);
     dedisp_sig_1 = ArrayDisplacement(dedisp_sig_1,2,1,2);
     shift_number = shift_number + 1;
     %% get new data 
     dedisp_sig_16(1) = chanOutBuf_average(16,(mm+1)/2);
     dedisp_sig_15(1) = chanOutBuf_average(15,(mm+1)/2);
     dedisp_sig_14(1) = chanOutBuf_average(14,(mm+1)/2);
     dedisp_sig_13(1) = chanOutBuf_average(13,(mm+1)/2);
     dedisp_sig_12(1) = chanOutBuf_average(12,(mm+1)/2);
     dedisp_sig_11(1) = chanOutBuf_average(11,(mm+1)/2);
     dedisp_sig_10(1) = chanOutBuf_average(10,(mm+1)/2);
     dedisp_sig_9(1) = chanOutBuf_average(9,(mm+1)/2);
     dedisp_sig_8(1) = chanOutBuf_average(8,(mm+1)/2);
     dedisp_sig_7(1) = chanOutBuf_average(7,(mm+1)/2);
     dedisp_sig_6(1) = chanOutBuf_average(6,(mm+1)/2);
     dedisp_sig_5(1) = chanOutBuf_average(5,(mm+1)/2);
     dedisp_sig_4(1) = chanOutBuf_average(4,(mm+1)/2);
     dedisp_sig_3(1) = chanOutBuf_average(3,(mm+1)/2);
     dedisp_sig_2(1) = chanOutBuf_average(2,(mm+1)/2);
     dedisp_sig_1(1) = chanOutBuf_average(1,(mm+1)/2);
    elseif( count ~=1 )
    count = 1;
    end
else                                                                       
    current_index = mm-(average_factor-1)/2;
    chanOutBuf_average(:,current_index) = (average_factor*chanOutBuf_average(:,current_index-1)+input(:,current_index+p)-input(:,current_index-q))/average_factor;
     %% shift
     dedisp_sig_16 = ArrayDisplacement(dedisp_sig_16,2,1,2);
     dedisp_sig_15 = ArrayDisplacement(dedisp_sig_15,2,1,2);
     dedisp_sig_14 = ArrayDisplacement(dedisp_sig_14,2,1,2);
     dedisp_sig_13 = ArrayDisplacement(dedisp_sig_13,2,1,2);
     dedisp_sig_12 = ArrayDisplacement(dedisp_sig_12,2,1,2);
     dedisp_sig_11 = ArrayDisplacement(dedisp_sig_11,2,1,2);
     dedisp_sig_10 = ArrayDisplacement(dedisp_sig_10,2,1,2);
     dedisp_sig_9 = ArrayDisplacement(dedisp_sig_9,2,1,2);
     dedisp_sig_8 = ArrayDisplacement(dedisp_sig_8,2,1,2);
     dedisp_sig_7 = ArrayDisplacement(dedisp_sig_7,2,1,2);
     dedisp_sig_6 = ArrayDisplacement(dedisp_sig_6,2,1,2);
     dedisp_sig_5 = ArrayDisplacement(dedisp_sig_5,2,1,2);
     dedisp_sig_4 = ArrayDisplacement(dedisp_sig_4,2,1,2);
     dedisp_sig_3 = ArrayDisplacement(dedisp_sig_3,2,1,2);
     dedisp_sig_2 = ArrayDisplacement(dedisp_sig_2,2,1,2);
     dedisp_sig_1 = ArrayDisplacement(dedisp_sig_1,2,1,2);
     shift_number = shift_number + 1;
     %% get new data 
     dedisp_sig_16(1) = chanOutBuf_average(16,current_index);
     dedisp_sig_15(1) = chanOutBuf_average(15,current_index);
     dedisp_sig_14(1) = chanOutBuf_average(14,current_index);
     dedisp_sig_13(1) = chanOutBuf_average(13,current_index);
     dedisp_sig_12(1) = chanOutBuf_average(12,current_index);
     dedisp_sig_11(1) = chanOutBuf_average(11,current_index);
     dedisp_sig_10(1) = chanOutBuf_average(10,current_index);
     dedisp_sig_9(1) = chanOutBuf_average(9,current_index);
     dedisp_sig_8(1) = chanOutBuf_average(8,current_index);
     dedisp_sig_7(1) = chanOutBuf_average(7,current_index);
     dedisp_sig_6(1) = chanOutBuf_average(6,current_index);
     dedisp_sig_5(1) = chanOutBuf_average(5,current_index);
     dedisp_sig_4(1) = chanOutBuf_average(4,current_index);
     dedisp_sig_3(1) = chanOutBuf_average(3,current_index);
     dedisp_sig_2(1) = chanOutBuf_average(2,current_index);
     dedisp_sig_1(1) = chanOutBuf_average(1,current_index);
end
    
  %% add channels together
dedisp_sig(16,:) = dedisp_sig_16;
dedisp_sig(15,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_15)) dedisp_sig_15];
dedisp_sig(14,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_14)) dedisp_sig_14];
dedisp_sig(13,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_13)) dedisp_sig_13];
dedisp_sig(12,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_12)) dedisp_sig_12];
dedisp_sig(11,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_11)) dedisp_sig_11];
dedisp_sig(10,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_10)) dedisp_sig_10];
dedisp_sig(9,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_9)) dedisp_sig_9];
dedisp_sig(8,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_8)) dedisp_sig_8];
dedisp_sig(7,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_7)) dedisp_sig_7];
dedisp_sig(6,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_6)) dedisp_sig_6];
dedisp_sig(5,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_5)) dedisp_sig_5];
dedisp_sig(4,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_4)) dedisp_sig_4];
dedisp_sig(3,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_3)) dedisp_sig_3];
dedisp_sig(2,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_2)) dedisp_sig_2];
dedisp_sig(1,:) = [ zeros(1,length(dedisp_sig_16)-length(dedisp_sig_1)) dedisp_sig_1];
dedisp_sum = zeros(1,length(dedisp_sig_16));
for i = 1:D
    dedisp_sum = dedisp_sum + dedisp_sig(i,:);
end

%% triggering
% generate MF for triggering
MF(mm)=dedisp_sig_16(delay_units(16))+dedisp_sig_15(delay_units(15))+dedisp_sig_14(delay_units(14))+dedisp_sig_13(delay_units(13))+dedisp_sig_12(delay_units(12)) ...
    +dedisp_sig_11(delay_units(11))+dedisp_sig_10(delay_units(10))+dedisp_sig_9(delay_units(9))+dedisp_sig_8(delay_units(8))+dedisp_sig_7(delay_units(7)) ...
    +dedisp_sig_6(delay_units(6))+dedisp_sig_5(delay_units(5))+dedisp_sig_4(delay_units(4))+dedisp_sig_3(delay_units(3))+dedisp_sig_2(delay_units(2))+dedisp_sig_1(delay_units(1));

% triggering
trigger_SNR(mm) = 10*log10((MF(mm)^2)/(noise_level)^2);
if trigger_SNR(mm) > trigger_level
    trigger_flag = trigger_flag +1;
    trigger_number = trigger_number +1;
else 
    trigger_flag = 0;
end

if trigger_flag == 1
        trigger_time = mm*(1/fs3)/1e-6;               % save the first trigger time, unit: us
end

%% dynamic plot
% %warning: enable dynamic plots would make the simulation run in a
% %extremely slow speed, so do not enable if really necessary
    
   % %dynamic_plot
%     N_pfb_plot = size(chanOutBuf,1);
%     N_pfb_plot_2 = ceil(N_pfb_plot/2);
%     fHz = (0:N_pfb_plot-1)*(fs2/N_pfb_plot);
%     spectrum_pfb =[abs(chanOutBuf(:,mm))]';
%     subplot(4,1,1);
%     title('real-time PFB spectrum');
% xlabel('Frequency / MHz');
% ylabel('Magnitude');
%     set(p1,'XData',fHz(1:N_pfb_plot_2)/1e6,'YData',spectrum_pfb(1:N_pfb_plot_2))
%     drawnow
%     hold on;
     
    % %dynamic_plot: average
%     if (flag_average)
%     N_pfb_plot = size(chanOutBuf,1);
%     N_pfb_plot_2 = ceil(N_pfb_plot/2);
%     fHz = (0:N_pfb_plot-1)*(fs2/N_pfb_plot);
%     spectrum_pfb2 =[abs(chanOutBuf_average(:,mm/average_factor))]';
%     subplot(4,1,2);
%     title('real-time PFB spectrum (averaged)');
% xlabel('Frequency / MHz');
% ylabel('Magnitude');
%     set(p2,'XData',fHz(1:N_pfb_plot_2)/1e6,'YData',spectrum_pfb2(1:N_pfb_plot_2))
%     drawnow
%     hold on;
%         end

    
   % %dynamic_plot dedispersed output
%     subplot(4,1,3);
%     title('Dedispersed output: whole signal length');
% xlabel('Time/us');
% ylabel('Magnitude');
%     set(p3,'XData',t_3*1e6,'YData',abs(dedisp_sum))
%     drawnow
%     hold on;
   
    mm = mm+1;   
end


%% capture candidate snapshot
% get start time of snapshot
snapshot_start = trigger_time*1e-6 - T_tot;
start_index = ceil(snapshot_start/(1/fs2));               % transfer start time to vector index
if snapshot_start <= 0                                   
    snapshot_start = 0;
    start_index = 1;
end

% get end time of snapshot
extra_factor= 0.2;                                         % to save a snapshot longer than target signal length (to add some buffer so as not to miss any target signal)
snapshot_end = trigger_time*1e-6 + T_tot*extra_factor;
end_index = ceil(snapshot_end/(1/fs2));                    % transfer end time to vector index
if snapshot_end > sig_length                               % if snapshot_end exceed the length of the raw signal
    snapshot_end = sig_length;
    end_index = snapshot_end/(1/fs2);
end

% get snapshot
snapshot = sig_downsp(start_index:end_index);
snapshot_t = 0:1/fs2:snapshot_end-1/fs2;
plotPlus(snapshot,snapshot_t,fs2,'Candidate snapshot');


%% plot results
figure;
for i = 1:D
   plot(t_3*1e6,dedisp_sig(i,:));
   hold on;
end 
title('Arrival time of 32 channel signals: averaged de-dispersed signal');
xlabel('Time/us'); 
ylabel('Magnitude');

%plot dedispersed sum
% figure;
% plot(t_3*1e6,(dedisp_sum));
% title('dedispersion sum');
% xlabel('Time/us'); 
% ylabel('Magnitude');

figure;
plot(t_3*1e6,(MF));
str = sprintf('MF trigger when input SNR = %d dB', input_SNR)
title(str)
xlabel('Time/us'); 
ylabel('Magnitude');

figure;
plot(t_3*1e6,(trigger_SNR));
str2 = sprintf('trigger SNR when input SNR = %d dB', input_SNR)
title(str2)
xlabel('Time/us'); 
ylabel('SNR/dB');

% plot original signal (dispersed)
figure;
for i = 1:D
   plot(t_2*1e6,(chanOutBuf_average(i,:)));
   hold on;
end
title('Arrival time of 32 channels: averaged dispersed signal');
xlabel('Time/us'); 
ylabel('Magnitude');

% % check particulate channel result (if necessary)
%
%  figure;
%  subplot(2,3,1);
%  plot(t_2*1e6,(chanOutBuf_average(14,:)));       % check channel 14
%  title('Check channel 14');
%  
%   subplot(2,3,2);
%  plot(t_2*1e6,(chanOutBuf_average(13,:)));       % check channel 13
%  title('Check channel 13');
%   
%    subplot(2,3,3);
%  plot(t_2*1e6,abs(chanOutBuf_average(12,:)));       % check channel 12
%  title('Check channel 12');
%   
%    subplot(2,3,4);
%  plot(t_2*1e6,(chanOutBuf_average(11,:)));       % check channel 11
%  title('Check channel 11');
%   
%    subplot(2,3,5);
%  plot(t_2*1e6,(chanOutBuf_average(10,:)));       % check channel 10
%  title('Check channel 10');
%   
%    subplot(2,3,6);
%  plot(t_2*1e6,(chanOutBuf_average(9,:)));       % check channel 9
%  title('Check channel 9');
 
