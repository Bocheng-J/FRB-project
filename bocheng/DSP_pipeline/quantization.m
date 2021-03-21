function [sig_out] = quantization(sig_in,bitWidth)
% quantitized input signal
% sig_in: input signal
% bitWidth: bit width

%sig_in = sig_in+abs(min(sig_in));
sig_out = round(sig_in*(bitWidth-1));

end