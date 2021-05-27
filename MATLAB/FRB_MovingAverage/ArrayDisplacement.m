%**********************************************
% Function: rignt shift
% BeforeArray: the input array
% len : right shift length
% AfterArray : output array
%********************************************
function [AfterArray] = ArrayDisplacement(BeforeArray)
% ArrayEnd = size(BeforeArray,2);
%         BeforeArrayL = BeforeArray(1:ArrayEnd-len);
%         BeforeArrayR = zeros(1,len);
%         AfterArray = fliplr([fliplr(BeforeArrayL) fliplr(BeforeArrayR)]);

        AfterArray=[0 BeforeArray(1:(end-1))];
end