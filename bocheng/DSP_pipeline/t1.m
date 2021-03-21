% Simulation of FRB detection pipelines
% Created: Bocheng Jia, 2021/3/18
close all;
clc;
clear;
trigger_flag = 0;                             % set the trigger flag to 0

%% simulate input sine/chirp signals
fs = 24000;                                   % sample frequency 
tMax = 240e-3;                               % time duration: 0.24s (240 ms)
t_1 = 0:1/fs:tMax-1/fs;
 t_length = length(t_1);

% input chirp
fchirp0 = 2200;
fchirp1 = 3600;
tchirp1 = 0.2;
sig = chirp(t_1,fchirp0,tchirp1,fchirp1);  
sig_complex = hilbert(sig);
sig_chirp = sig_complex;

% input sine
f1 = 1500;
f2 = 2100;
f3 = 2500;
f4 = 5100;

s1_com = 0.4*exp(1j*2*pi*t_1*f1);     % complex sine input waves
s2_com = 0.6*exp(1j*2*pi*t_1*f2);
s3_com = 0.8*exp(1j*2*pi*t_1*f3);
s4_com = 1*exp(1j*2*pi*t_1*f4);

s4_com((t_length/4+1):t_length)=0;
s3_com(1:(t_length/4))=0;
s3_com((t_length*2/4)+1:t_length)=0;
s2_com(1:(t_length*2/4))=0;
s2_com((t_length*3/4)+1:t_length)=0;
s1_com(1:(t_length*3/4))=0;

%sig = s1_com + s2_com +s3_com+s4_com;           % sum
signal = s1_com + sig_chirp;           % sum
% sig = signal;
% t =t2;

tMax2 = 440e-3;                               % time duration: 0.44s (440 ms)
t2 = 0:1/fs:tMax2-1/fs;
sig = zeros(1,length(t2));

sig(1,(length(t2)-length(t_1))/2:((length(t2)-length(t_1))/2)+length(t_1)-1)=signal;
t = t2;
sig_title5 = sprintf('Input signal ');
plotPlus(sig,t,fs,sig_title5);



%% quantization 4-bit
bit_Width = 4;

% quantize input signal
[quantized_sig] = quantization(sig,bit_Width);
plotPlus(quantized_sig,t,fs,'4-bit quantized signal');

% quantize tamplate
[quantized_template] = quantization(sig_chirp,bit_Width);


%% pepeline 1: pfb channelizer to 6 channels
M = 6;                                %% channel number
D = 3;                                % decimation number
load('lp_12000_6ch');                 % load prototype filter coefficients

% channelization for quantized input signal
channelizer_out = PFB_channelizer(quantized_sig,M,D,lp_12000_6ch);
intermediate_data = channelizer_out(1:M/2,:);                                % use intermediate data to swap first half channels and last half channels, to make channels placed in order
channelizer_out(1:M/2,:) = channelizer_out(M/2+1: M,:);                      
channelizer_out(M/2+1: M,:) = intermediate_data;

% calculate the max value of each channel
for k = 1:6                                                              
fft_value(k,:) = fftshift(20*log10(abs(fft(channelizer_out(k,:)))));
fft_max = 1.2*max(max(fft_value)); 
end

% plot channelization result in frequency domain
N = length(channelizer_out);
N_2 = ceil(N/2);
fax_Hz = ((0:N-1)-N_2)*(fs/D/N);
figure;
sgtitle('pipeline 1:  channelization result (spectrum)');
for k=1:M
    subplot(2,3,k);
    plot((fax_Hz),fftshift(20*log10(abs(fft(channelizer_out(k,:))))),'linewidth',2);grid on;
    title(['Spectrum, Channel(',num2str(k),')']);
    xlabel('Frequency');
    ylabel('Log Mag (dB)');
    axis([-(fs/D)/2 (fs/D)/2 -0.2*fft_max 1.2*fft_max])
end

% plot channelization result in time domain
figure;
sgtitle('pipeline 1: channelization result (time domain)');
for k=1:M
    subplot(2,3,k);
    plot(real(channelizer_out(k,:)));grid on;
    title(['Channel(',num2str(k),')']);
    xlabel('samples');
    ylabel('amplitude');
    ylim([-1 1])
end


% channelization for template signal
ch_out_template = PFB_channelizer(quantized_template,M,D,lp_12000_6ch);
intermediate_data_temp = ch_out_template(1:M/2,:);
ch_out_template(1:M/2,:) = ch_out_template(M/2+1: M,:);
ch_out_template(M/2+1: M,:) = intermediate_data_temp;

% find the template signal's channel
[x y]=find(ch_out_template==max(max(ch_out_template)));
template_ch = x;


%% pepeline 1: matched filter

MF_Out = zeros(M,length(channelizer_out)+length(ch_out_template)-1);
MF_Index = zeros(M,1);
MF_Max = zeros(M,1);

% do matched filter for each channel 
for k = 1:M
[MF_Out(k,:),MF_Index(k),MF_Max(k)] = matchedFiltering(real(channelizer_out(k,:)),real(ch_out_template(template_ch,:)));
end

% calculate the max of matched filter out
for k = 1:6                                                              
mf_max_data = 1.2*max(max(MF_Out)); 
end

% plot matched filter result
figure;
for k = 1:M
    sgtitle('pipeline 1: Matched filter output');
    subplot(2,M/2,k);
    plot(MF_Out(k,:));
    title(['channel',num2str(k)]);
    xlabel('sample');ylabel('amplitude');
    ylim([-0.2*mf_max_data 1.1*mf_max_data])
end



%% pepeline 1: find the target signal's channel and arrival time
[x1 y1]=find(MF_Out==max(max(MF_Out)));
target_ch = x1;
arrival_index = MF_Index(target_ch);
arrival_time = arrival_index/length(channelizer_out) * tMax2 * 1e3;       % arrival time in ms
disp(['Target signal is in channel ', num2str(target_ch)]);
disp(['Arrival time :', num2str(arrival_time),' ms']);


%% pepeline 1: generate trigger
threshold = 5;

for k = 1:M 
if MF_Max(k)>threshold                                         % when the target signal is found, set trigger to 1
    trigger_flag = 1;
end
end


%% pipeline 2: get the snapshot of target signal
Tmax_snapshot = 240e-3;                                          % snapshot length: 240 ms
t_snapshot = 0:1/(fs/D):Tmax_snapshot-1/(fs/D);

% decimation
decimated_data = reshape(quantized_sig,D,[]);
t_decimated = reshape(t,D,[]);

% create a 240ms long snapshot
 snapshot = zeros(1,length(t_snapshot));                                        
if trigger_flag == 1
   arrival_index_decimated = round(arrival_index*length(t_decimated)/length(channelizer_out));
   snapshot = decimated_data(1,arrival_index_decimated:arrival_index_decimated+length(t_snapshot)-1);
   trigger_flag = 0;
end

sig_title1 = sprintf('pipeline 2: target snapshot ');
plotPlus(snapshot,t_snapshot,fs/D,sig_title1);


%% pipeline 3: record orginal signal

recording_data = decimated_data(1,:);
% always recording
plotPlus(recording_data(1,:),t_decimated,fs/D,'pipeline 3: recording mode');

