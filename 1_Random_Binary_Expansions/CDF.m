function [Cumulative_DF] = CDF(p,n)
x_interval=linspace(0,1,2^n);
Cumulative_DF=zeros(1,2^n);
for i=1:2^n-1
    binary_expansion=binary(x_interval(i),n);
    Cumulative_DF(i)=(1-p)*binary_expansion(1);
    for j=1:n-1
        total=0;
        for k=1:j
            total=total+binary_expansion(k);
        end
        Cumulative_DF(i)=Cumulative_DF(i)+(p)^(total)*(1-p)^(j+1-total)*(binary_expansion(j+1));
    end
end
Cumulative_DF(2^n)=1;
plot(x_interval,Cumulative_DF,'b-')
xlabel('x')
ylabel('F(x)')
title('Cumulative Distribution function for x')