function channelizer_out = PFB_channelizer(sig_in,M,D,filter_coef)
% function: oversampling polyphase channelizer filter bank
% sig_in : input signal
% M : number of channels
% D : decimation ratio
% filter_coef : prototype filter coefficient
% channelizer_out : channelizer output of M channels

h = reshape(filter_coef,M,[]);
x = sig_in;

%% OPFB
inputDat = zeros(D,1);
inputDatBuf = zeros(M,size(h,2));
filtOutBuf = zeros(M,1);
chanOut = zeros(M,1);
chanOutBuf = zeros(M,length(x)/D);
flag = 0;

mm = 1;
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
    chanOut = ifft(filtOutBuf);
    chanOutBuf(:,mm) = fftshift(chanOut);
    mm = mm+1;
end

channelizer_out = chanOutBuf;

end