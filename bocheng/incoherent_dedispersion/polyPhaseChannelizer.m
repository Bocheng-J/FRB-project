function [sig_out] = polyPhaseChannelizer(sig_in,D,h)
% sig_in: input signal
% D: number of channel
% h: coef of prototype lowpass filter
% sig_out: output signal

% decimation
sig_pf = flipud(reshape(sig_in,D,length(sig_in)/D));
h_pf = reshape(h,D,length(h)/D);


% filtering
sig_channel = zeros(D,length(sig_in)/D);
for i=1:D
    sig_channel(i,:) = filter(h_pf(i,:),1,sig_pf(i,:));
end


% ifft
sig_fft = zeros(D,length(sig_channel));
for m=1:length(sig_fft)
    sig_fft(:,m) = ifft(sig_channel(:,m));
end

sig_out = D*sig_fft;         % output of channelizer

end
