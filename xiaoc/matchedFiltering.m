function [MF_Out,MF_Index,MF_Max] = matchedFiltering(MF_In,MF_template)
% MF_In: input signal
% MF_template: matched filter template
% MF_Out: matched filter output
% MF_Index: position of the most matched point
% MF_Max: value of the most matched point


energyTemplate = sum(MF_template.^2);
MF_template = MF_template/sqrt(energyTemplate);
MF_Out = conv(MF_In,fliplr(MF_template));

[MF_Max, Tmax] = max(MF_Out);                                             
MF_Index = Tmax-length(MF_template);

end