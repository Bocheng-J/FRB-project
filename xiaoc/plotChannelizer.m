function [] = plotChannelizer(sig,t,D,fs,flag)
% plot channelizer output in time domain figure and magnitude spectrum
% sig: channelizer output
% t: time array before channelizer
% fs: sampling frequency before channelizer
% flag: 0 for time domain figure; 1 for magnitude spectrum
figure;
switch flag
    case 0
       sgtitle('channelizer output (time domain)'); 
       for i=1:D
           subplot(4,D/4,i);
           plot(downsample(t,D)*1e6,real(sig(i,:)));grid on;
           title(['channel',num2str(i)]);
           xlabel('t/us');ylabel('voltage/v');
       end
    case 1
        sgtitle('channelizer output (magnitude spectrum)');
        for i=1:D
            N = length(sig(i,:));
            %N_2 = ceil(N/2);
            fax_Hz = (0:N-1)*(fs/D/N);
            fftsig = fft(sig(i,:));
            subplot(4,D/4,i);
            plot(fax_Hz/1e6,abs(fftsig));
            title(['channel',num2str(i)]);
            xlabel('frequency/MHz');ylabel('magnitude');
        end
end
       
end

