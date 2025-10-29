function [EULER] = euler_2(x_max)
vector_h=zeros(1,16);
vector_E=zeros(1,16);
EULER=zeros(16,2);
for k=0:15
    n=2^k;
    h=x_max/n;
    Y_n=0;
    x_n=0;
    for i=1:n
        Y_n=Y_n + h*f(x_n,Y_n);
        x_n=x_n+h;
    end
    Z_n=exp(-0.1)-exp(-16*0.1);
    E_n=Y_n-Z_n;
    EULER(k+1,1)=h;
    EULER(k+1,2)=E_n;
    vector_h(k+1)=log(h);
    vector_E(k+1)=log(abs(E_n));
    fprintf('h=%.16f, E_n=%.16f\n',log(h), log(abs(E_n)))
end 
plot(vector_h,vector_E,'b-')
xlabel('log(h)')
ylabel('log(E_n)')
title('a graph to show how the error changes with step size - Euler')