function [vector_f] = Monte_Carlo(p,n)
vector_U=zeros(1,n);
vector_x=zeros(1,n);
vector_X=zeros(1,100000);
vector_f=zeros(1,10000);
for k=1:100000
    for i=1:n
        random_number=rand;
        if random_number<=p
            result=1;
        else
            result=0;
        end
    vector_U(i)=result;
    vector_x(i)=result/(2^i);
    end
    X_n=sum(vector_x);
    vector_X(k)=X_n;
end
x_interval=linspace(0,1,10000);
for j=1:10000
    total_result=0;
    for a=1:100000
        if vector_X(a)<=x_interval(j)
            result_2=1;
        else
            result_2=0;
        end
        total_result=total_result + result_2;
    end
    vector_f(j)=(total_result)/100000;
end
plot(x_interval,vector_f,'b-')
xlabel('x')
ylabel('empirical F')
title('A plot of the empirical distribution function for p = 2/3 and n = 30')
