function [Gradient_CDF_left] = CDF_differentiability_left(p,c)
delta_interval=zeros(1,8);
for j=4:11
    delta_interval(j-3)=-1/2^j;
end
Gradient_CDF_left=zeros(1,8);
for i=1:8
    Gradient_CDF_left(i)=(CDF_2(p,c+delta_interval(i))-CDF_2(p,c))/delta_interval(i);
end
plot(delta_interval,Gradient_CDF_left,'b-')
xlabel('delta')
ylabel('(F(c+delta)-F(c))/delta')
title('A graph showing gradient as delta varies (from the left)')