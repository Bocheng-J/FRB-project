close all;
clear;
clc;
load coef_LP_2M.mat
load coef_LP_32M.mat

%% generate test signal
fs = 256e6;
tMax = 4e-5;
t = 0:1/fs:tMax-1/fs;

fchirp0 = 90e6;
fchirp1 = 80e6;
tchirp1 = 4e-5;
sig = chirp(t,fchirp0,tchirp1,fchirp1);

MF_template = downsample(sig,4);
plotPlus(sig,t,fs,'orginal signal');

%% add delay and rfi

% add delay
delay = randi([1,1024], 1);       % Generate a random delay between 1 and 1024
sig = [zeros(1,delay), sig, zeros(1,1024-delay)];

% add noise
noise = 0.1*randn(1,length(sig));

% add RFI
fRFI = 40e6;
t = 0:1/fs:(1/fs)*(length(sig)-1);
rfi = 0.1*sin(2*pi*fRFI*t);
sig = sig+noise+rfi;

plotPlus(sig,t,fs,'delayed signal with noise and RFI');


%% quantization

sig = quantization(sig,12);
plotPlus(sig,t,fs,'12-bit quantized signal');


%% channelization
D = 4;
h = coef_LP_32M;
sig = polyPhaseChannelizer(sig,D,h);

plotChannelizer(sig,t,D,fs,0);
plotChannelizer(sig,t,D,fs,1);                  % 我感觉可能是因为我们画频谱图的方法有问题


%% matched filtering
figure;
sgtitle('Matched filter output');
MF_Out = zeros(D,length(sig)+length(MF_template)-1);
MF_Index = zeros(D,1);
MF_Max = zeros(D,1);
for i = 1:D
    [MF_Out(i,:),MF_Index(i),MF_Max(i)] = matchedFiltering(real(sig(i,:)),real(MF_template));
    subplot(4,D/4,i);
    plot(MF_Out(i,:));
    title(['channel',num2str(i)]);
    xlabel('sample');ylabel('amplitude');
end

[~,numChannel] = max(MF_Max(1:D/2));            % find the channel number of the target signal
sig = sig(numChannel,MF_Index(numChannel):end); % throw away the delay before the target signal
fs = fs/D;
t = 0:1/fs:(1/fs)*(length(sig)-1);
plotPlus(sig,t,fs,'target signal');

%% channelization for de-dispersion
D = 16;
h = coef_LP_2M;
pad = D*ceil(length(sig)/D)-length(sig);
sig = [sig,zeros(1,pad)];
t = 0:1/fs:(1/fs)*(length(sig)-1);
sig = polyPhaseChannelizer(sig,D,h);

plotChannelizer(sig,t,D,fs,0);
plotChannelizer(sig,t,D,fs,1);
