%**********************************************
%对一组数组进行左、右位移
%TowardSpring
%2019年7月3日20:28:55
%输入：BeforeArray:需要反转的数组；direction：反转的方向，1向左，2向右;len:位移长度；type:数据空缺补齐方式，1端点移位，2补零；
%输出：AfterArray：移位后的数据
%********************************************
function [AfterArray] = ArrayDisplacement(BeforeArray,direction,len,type)
ArrayEnd = length(BeforeArray);
if type == 1
    if direction == 1
        BeforeArrayL = BeforeArray(1:len);
        BeforeArrayR = BeforeArray(len+1:ArrayEnd);
        BeforeArrayL = fliplr(BeforeArrayL);
        BeforeArrayR = fliplr(BeforeArrayR);
        BeforeArray = [BeforeArrayL BeforeArrayR];
        AfterArray = fliplr(BeforeArray);
    else
        BeforeArrayL = BeforeArray(1:ArrayEnd-len-1);
        BeforeArrayR = BeforeArray(ArrayEnd-len:ArrayEnd);
        BeforeArrayL = fliplr(BeforeArrayL);
        BeforeArrayR = fliplr(BeforeArrayR);
        BeforeArray = [BeforeArrayL BeforeArrayR];
        AfterArray = fliplr(BeforeArray); 
    end
else
    if direction == 1
        BeforeArrayL = zeros(1,len);
        BeforeArrayR = BeforeArray(len+1:ArrayEnd);
        BeforeArrayL = fliplr(BeforeArrayL);
        BeforeArrayR = fliplr(BeforeArrayR);
        BeforeArray = [BeforeArrayL BeforeArrayR];
        AfterArray = fliplr(BeforeArray);
    else
        BeforeArrayL = BeforeArray(1:ArrayEnd-len);
        BeforeArrayR = zeros(1,len);
        BeforeArrayL = fliplr(BeforeArrayL);
        BeforeArrayR = fliplr(BeforeArrayR);
        BeforeArray = [BeforeArrayL BeforeArrayR];
        AfterArray = fliplr(BeforeArray);
    end
end