%**********************************************
%��һ�������������λ��
%TowardSpring
%2019��7��3��20:28:55
%���룺BeforeArray:��Ҫ��ת�����飻direction����ת�ķ���1����2����;len:λ�Ƴ��ȣ�type:���ݿ�ȱ���뷽ʽ��1�˵���λ��2���㣻
%�����AfterArray����λ�������
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