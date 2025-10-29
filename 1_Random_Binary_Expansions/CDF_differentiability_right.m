function [Gradient_CDF] = CDF_differentiability_right(p,c)
delta_interval=zeros(1,8);
for j=4:11
    delta_interval(j-3)=1/2^j;
end
Gradient_CDF=zeros(1,8);
for i=1:8
    Gradient_CDF(i)=(CDF_2(p,c+delta_interval(i))-CDF_2(p,c))/delta_interval(i);
end
plot(delta_interval,Gradient_CDF,'b-')
xlabel('delta')
ylabel('(F(c+delta)-F(c))/delta')
title('A graph showing gradient as delta varies (from the right)')