% function: moving averaging
% input: input signal
% M: averaging length
%outpout: averaging output 

function output = moving_average(input,M)
% input = x;
% M = 3;             % averaging length 
p = (M-1)/2;
q = p+1;
output = zeros(1,length(input));
for i = 1:length(input)-p
    if i <= q
        if i == 1
            output(i) = input(1);
        else
            temp_vec = zeros(1,2*i-1);
            for j = 1:2*i-1
                temp_vec(j) = input(j);
            end
            temp_addition = cumsum(temp_vec);
            output(i) = temp_addition(2*i-1)/(2*i-1);
        end
    else
        output(i) = (M*output(i-1)+input(i+p)-input(i-q))/M;
    end
end
end