function [Cumulative_DF] = CDF_2(p,x)
binary_expansion=binary(x,11);
Cumulative_DF=(1-p)*binary_expansion(1);
for j=1:10
    total=0;
    for k=1:j
        total=total+binary_expansion(k);
    end
    Cumulative_DF=Cumulative_DF+(p)^(total)*(1-p)^(j+1-total)*(binary_expansion(j+1));
end
