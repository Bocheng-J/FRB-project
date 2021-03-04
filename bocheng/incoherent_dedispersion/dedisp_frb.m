close all;
clear;
clc;
load coef_LP_1250_4.mat;
format long g;

%% generate FRB
f_start = 1000;  % start frequency           
f_stop = 5000;  % stop frequency
DM = 349;      % dispersion measure
SNR = 29;       % SNR in dB
nFiltSize = 2^10;   % size of FRB filter
alpha = 1;    % FRB curvature index
frbwidth = 35;    % frb instantaneous frequency bandwidth index


BW = f_stop - f_start;  % total signal bandwidth
t_start = 4140*DM*f_start^(-2); % FRB delay at start frequency
t_stop = 4140*DM*f_stop^(-2);   % FRB delay at stop frequency
t_d = t_stop - t_start;         % delay between start and stop frequencies
t_d_samples = round(abs(t_d*BW));   % number of samples corresponding to t_d
%nTotSam = round(t_d_samples/0.6);   % total sample point
nTotSam = 9648;

if t_d_samples + round(nTotSam/10) >= nTotSam        % signal starts from 10%
    disp('increase number of samples');
    disp(['need at least ',num2str(t_d_samples + round(nTotSam/10)),' samples in total']);
    %break;
end

t = 0:t_d_samples-1;    % time....

fs1 = 2.34/24000;
% FRB filter
filt = ifft(exp(-(5*linspace(-(nFiltSize-1)/nFiltSize,(nFiltSize-1)/nFiltSize,nFiltSize)).^2/2*frbwidth)).*exp(-1i*2*pi*0.5*(0:nFiltSize-1));
puls = 10^(SNR/20)*(randn(1,t_d_samples)+1i*randn(1,t_d_samples))/sqrt(2);  % pulse before filtering
puls = conv(puls,filt); puls = puls(1:t_d_samples); % narrow pulse (i.e. after filtering)
% apply fade-in/-out and dispersion
puls  = puls.*exp(-(1*linspace(-(t_d_samples-1)/t_d_samples,(t_d_samples-1)/t_d_samples,t_d_samples)).^2/2*5).*exp(1i*2*pi*( alpha/(3*t_d_samples^2)*t.^3 - (alpha+1)*t.^2/(2*t_d_samples) - t));
signal = (randn(1,nTotSam)+1i*randn(1,nTotSam))/sqrt(2);    % background noise
% background noise + FRB starting at 10%
signal(round(nTotSam/10):round(nTotSam/10)+t_d_samples-1) = puls + signal(round(nTotSam/10):round(nTotSam/10)+t_d_samples-1);


%% plot signal in frequency domain
figure(1);
subplot(1,2,1);
myspectro(signal,128,f_start,f_stop);
title('initial signal');

%% plot signal in time domain 
nSam = length(signal);
nResol = 128;
nTime = floor(nSam/nResol);
sig = reshape(signal(1:nTime*nResol),nResol,nTime);
sig = abs(fft(sig,[],1)).^2;
BW = f_stop - f_start;
t2 = linspace(0,nSam/BW,length(randn(1,nTotSam)));

   subplot(1,2,2);
plot(t2,abs(signal));
xlabel('time [s]');
ylabel('magnitude');
title('initial signal in time domain');


%% generate the original input signal
fs = 10000;                                   % sample frequency 1 GHz
tMax = 2400e-3;                               % time duration: 2.4 s
t = 0:1/fs:tMax-1/fs;
t_length = length(t);
f1 = 1000;
f2 = 1700;
f3 = 2500;
f4 = 5000;

% s1 = 0.4*sin(2*pi*f1*t);                        % real sine input waves
% s2 = 0.6*sin(2*pi*f2*t);
% s3 = 0.8*sin(2*pi*f3*t);
% s4 = sin(2*pi*f4*t);

s1_com = 0.4*exp(1j*2*pi*t*f1);     % complex sine input waves
s2_com = 0.6*exp(1j*2*pi*t*f2);
s3_com = 0.8*exp(1j*2*pi*t*f3);
s4_com = 1*exp(1j*2*pi*t*f4);

%make the high frequency singals come earlier, low frequency signals come
%later
s4_com((t_length/4+1):t_length)=0;
s3_com(1:(t_length/4))=0;
s3_com((t_length*2/4)+1:t_length)=0;
s2_com(1:(t_length*2/4))=0;
s2_com((t_length*3/4)+1:t_length)=0;
s1_com(1:(t_length*3/4))=0;


% fchirp0 = 1700;
% fchirp1 = 3800;
% tchirp1 = 2;
% sig = chirp(t,fchirp0,tchirp1,fchirp1);


%sig = s1_com+s2_com+s3_com+s4_com;           % sum
sig = signal;
sig_title5 = sprintf('Input signal ');
plotPlus(sig,t2,fs,sig_title5);



%% plot time vs frequency
figure;
window = hamming(128);  			% haming window
noverlap = 120; 					%overlap length
nfft = 2^nextpow2(length(window)); 	% dft number
spectrogram(sig,window,noverlap,nfft,fs1,'yaxis');
% [S,F,T,P]= spectrogram(X,window,noverlap,nfft,fs,'yaxis');
title('Input signal: Frequency vs Time');
xlabel('t/s'); ylabel('Frequency/Hz'); 


%% PFB channelization 
D = 4;
h = coef_LP_1250_4;

channelizer_output = polyPhaseChannelizer(sig,D,h);

plotChannelizer(channelizer_output,t2,D,fs,0);
plotChannelizer(channelizer_output,t2,D,fs,1);  


%% Incoherent de-dispersion
DM = 200;                                   %  unit:(pc cm-3) 
k_dm = 4150;                                 % unit: MHz2 pc?1 cm3 s 

% calculate center frequencies
channels_Fc = zeros(1,D+1);
 for i = 1:D+1                          
     if(i == 1)
        channels_Fc(i) = fs/(4*D);
     elseif(i ~= D+1)
         channels_Fc(i) = (i-1)*(fs/D);
     else
         channels_Fc(i) = fs-(fs/(4*D));
     end
 end

 % calculate center frequency square difference
shift_t = zeros(1,D);
for i = 1:D+1
     if i ~= D+1
   shift_t(i)= (1/(channels_Fc(1)^2))-(1/(channels_Fc(i+1)^2));
     else 
   shift_t(i) = 1/(channels_Fc(1)^2);
     end
end

% calcluate delay time
delay_t = zeros(1,D+1);
for i = 1:D+1
delay_t(i) = k_dm * DM * shift_t(i);
end


%% plot signal spectrum
N = length(sig);
N_2 = ceil(N/2);
fax_Hz = (0:N-1)*(fs/N);
fftsig = fft(sig);
% figure;
% plot(fax_Hz(1:N_2),abs(fftsig(1:N_2))); grid on;
% title('magnitude spectrum');
% xlabel('frequency/Hz'); ylabel('magnitude');


%% add delay 
delayed_sig = zeros(D,N/D);
delayLength_bin = zeros(1,D+1);                                   
  for i=1:D+1
      delayLength_bin(i) = round(delay_t(i)/(1/(fs/D)));
  end
 
dedispersed_sig = zeros(1,10*N/D);
channel_bin_length = length(sig)/D;
for i = 1:D
dedispersed_sig(1+delayLength_bin(i):channel_bin_length+delayLength_bin(i)) = dedispersed_sig(1+delayLength_bin(i):channel_bin_length+delayLength_bin(i)) + channelizer_output(i,:);
end 


%% View result
N_new = length(dedispersed_sig);
fs_2 = 10000/D;            
tmax_2 = N_new / fs_2;
t_2 = 0:1/fs_2:tmax_2-1/fs_2;
plotPlus_fullband(dedispersed_sig,t_2,fs_2,'Dedispersed signal');

% figure;
% window = hamming(128);  			% haming window
% noverlap = 120; 					% overlap length
% nfft = 2^nextpow2(length(window)); 	% dft number
% spectrogram(dedispersed_sig,window,noverlap,nfft,fs_2,'yaxis');
% %[S,F,T,P]= spectrogram(X,window,noverlap,nfft,fs,'yaxis');
% title('dedispersed signal: Frequency vs Time');
% xlabel('t/s'); ylabel('Frequency/Hz'); 

figure;
myspectro(dedispersed_sig,128,f_start,f_stop);
title('de-dispersed signal');

%% coherent dedispersion

% CD=coherent_dedispersion(sig,200, 1000, 5000, 1);
% 
% N_new_CD = length(CD);         
% 
% plotPlus_fullband(CD,t2,fs,'coherent Dedispersed signal');
% 
% figure;
% window = hamming(128);  			% haming window
% noverlap = 120; 					% overlap length
% nfft = 2^nextpow2(length(window)); 	% dft number
% spectrogram(CD,window,noverlap,nfft,fs,'yaxis');
% %[S,F,T,P]= spectrogram(X,window,noverlap,nfft,fs,'yaxis');
% title('coherent dedispersed signal: Frequency vs Time');
% xlabel('t/s'); ylabel('Frequency/Hz'); 
