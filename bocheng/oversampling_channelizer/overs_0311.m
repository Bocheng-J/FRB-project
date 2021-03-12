clc;
clear;
close all;

%% prototype filter
h0=remez(70,[0 0.5 1.5 6]/6,{'myfrf',[1 1 0 0]},[1 1]);  % prototype filter
scl=max(h0);
hh2=reshape([h0/scl 0],6,12);
f1=0.5;
f2=1.5;

figure(1)
% subplot(2,1,1);
fh0=fftshift(20*log10(abs(fft(h0,2048))));
plot((-0.5:1/2048:0.5-1/2048)*12,fh0,'linewidth',2)
hold on
plot( [-f1 -f1 +f1 +f1],[-90 -0.1 -0.1 -90],'r--','linewidth',2)
plot( [f2 f2 20],[-20 -80 -80],'r--','linewidth',2)
plot(-[f2 f2 20],[-20 -80 -80],'r--','linewidth',2)
hold off
grid on
axis([-6 6 -150 10])
title(['Spectrum; Passband 0-to-',num2str(f1),' kHz, Stopband ',num2str(f2),'-to-6 kHz, 80 db Attenuation'],'fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)


f_h0=fftshift(abs(fft(h0,128)));

% subplot(2,1,2);
% for zz=-1:1
% plot((-0.5:1/128:0.5-1/128)*12+zz*6,f_h0,'--','linewidth',2)
% end
% axis([-6 6 0 1.1]);

%% simulate FRB

f_start = 1500;  % start frequency           
f_stop = 9000;  % stop frequency
DM = 349;       % dispersion measure
SNR = 39;       % SNR in dB
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
figure;
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



%% simulate input sine/chirp signals
fs = 12000;                                   % sample frequency 
tMax = 240e-3;                               % time duration: 36 us
t = 0:1/fs:tMax-1/fs;
t_length = length(t);

% input chirp
fchirp0 = 1700;
fchirp1 = 3800;
tchirp1 = 0.2;
sig = chirp(t,fchirp0,tchirp1,fchirp1);
 sig_complex = hilbert(sig);
 sig = sig_complex;

% input sine
f1 = 2000;
f2 = 2100;
f3 = 2500;
f4 = 5100;

% s1 = 0.4*sin(2*pi*f1*t);                        % real sine input waves
% s2 = 0.6*sin(2*pi*f2*t);
% s3 = 0.8*sin(2*pi*f3*t);
% s4 = sin(2*pi*f4*t);

s1_com = 0.4*exp(1j*2*pi*t*f1);     % complex sine input waves
s2_com = 0.6*exp(1j*2*pi*t*f2);
s3_com = 0.8*exp(1j*2*pi*t*f3);
s4_com = 1*exp(1j*2*pi*t*f4);

s4_com((t_length/4+1):t_length)=0;
s3_com(1:(t_length/4))=0;
s3_com((t_length*2/4)+1:t_length)=0;
s2_com(1:(t_length*2/4))=0;
s2_com((t_length*3/4)+1:t_length)=0;
s1_com(1:(t_length*3/4))=0;

%sig = s1_com + s2_com +s3_com+s4_com;           % sum
%sig = s1_com ;           % sum
sig = signal;
t =t2;
sig_title5 = sprintf('Input signal ');
plotPlus(sig,t,fs,sig_title5);




%% pfb 
reg=zeros(3,24);
v0=zeros(1,3)';
v1=zeros(1,6)';
v2=zeros(6,27);


x0 = sig; 
figure;
plot(t,real(x0));

m=1;
flg=0;
for n=1:3:2880-3
    v0=fliplr(x0(n:n+2)).';
    reg=[v0 reg(:,1:23)];
        for k=1:3
          v1(k)=reg(k,1:2:24)*hh2(k,:)';
          v1(k+3)=reg(k,2:2:24)*hh2(k+3,:)';
        end
        if flg==0
            flg=1;
        else
            flg=0;
            v1=[v1(4:6);v1(1:3)];
        end
        v2(:,m)=ifft(v1);
        m=m+1;
end






%% filter response plot
figure;
for k=1:6
  subplot(3,2,k)
plot((-0.5:1/128:0.5-1/128)*4,fftshift(20*log10(abs(3*fft(v2(k,:),128)))),'linewidth',2)
hold on
plot([-1 -1 +1 +1],[-100 0 0 -100],'--r','linewidth',2)
hold off
grid on
%axis([-6 6 -100 10])
title(['Frequency Response, Chan(',num2str(k-1),')'],'fontsize',14)
xlabel('Frequency/KHz','fontsize',14)
ylabel('Amplitude','fontsize',14)
end


%% plot channelizer result



for k = 1:6                                                              
fft_value(k,:) = fftshift(20*log10(abs(fft(v2(k,:)))));
fft_max = 1.2*max(max(fft_value)); 
end

figure;
sgtitle('channelizer output (magnitude spectrum)');

for k=1:6
subplot(2,3,k)
plot((-0.5:1/959:0.5-1/959)*4,fftshift(20*log10(abs(fft(v2(k,:))))))
hold on
%plot((-0.5:1/128:0.5-1/128)*4,20*log10(f_h0),'r','linewidth',2.5)
hold off
axis([-2 2 -0.2*fft_max 1.2*fft_max])
grid on
title(['Channel (',num2str(k),')'],'fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)
xlabel('Frequency','fontsize',14)
end

% result in time domain
figure;
sgtitle('channelizer output (time domain)');
D=6;
       for i=1:D
           subplot(2,D/2,i);
           t2 = downsample(t,D/2);
           t2 = t2(1:959);
           plot(t2,real(v2(i,:)));grid on;
           ylim([-1 1])
           title(['channel',num2str(i)]);
           xlabel('t/us');ylabel('voltage/v');
       end

       
 %% time vs frequency
 
 figure;
window = hamming(128);  			% haming window
noverlap = 120; 					% overlap length
nfft = 2^nextpow2(length(window)); 	% dft number
spectrogram(v2(1,:),window,noverlap,nfft,fs/3,'yaxis');
xlabel('t/s'); ylabel('Frequency/KHz'); 
title('channel 1 de-dispersed signal');

figure;
spectrogram(v2(2,:),window,noverlap,nfft,fs/3,'yaxis');
xlabel('t/s'); ylabel('Frequency/KHz'); 
title('channel 2 de-dispersed signal');

figure;
spectrogram(v2(3,:),window,noverlap,nfft,fs/3,'yaxis');
xlabel('t/s'); ylabel('Frequency/KHz'); 
title('channel 3 de-dispersed signal');

figure;
spectrogram(v2(4,:),window,noverlap,nfft,fs/3,'yaxis');
%myspectro(v2(4,:),128,f_start,f_stop);
xlabel('t/s'); ylabel('Frequency/KHz'); 
title('channel 4 de-dispersed signal');

figure;
spectrogram(v2(5,:),window,noverlap,nfft,fs/3,'yaxis');
xlabel('t/s'); ylabel('Frequency/KHz'); 
title('channel 5 de-dispersed signal');

figure;
spectrogram(v2(6,:),window,noverlap,nfft,fs/3,'yaxis');
xlabel('t/s'); ylabel('Frequency/KHz'); 
title('channel 6 de-dispersed signal');
