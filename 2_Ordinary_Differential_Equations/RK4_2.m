function [RUNGE_KUTTA] = RK4_2(x_max)
vector_h=zeros(1,16);
vector_E=zeros(1,16);
RUNGE_KUTTA=zeros(16,2);
for k=0:15
    n=2^k;
    h=x_max/n;
    Y_n=0;
    x_n=0;
    for i=1:n
        k1=h*f(x_n,Y_n);
        k2=h*f((x_n+(h/2)), (Y_n+((k1)/2)));
        k3=h*f((x_n+(h/2)), (Y_n+((k2)/2)));
        k4=h*f((x_n+h), (Y_n+k3));
        Y_n=Y_n + (k1+2*k2+2*k3+k4)/6;
        x_n=x_n+h;
    end
    Z_n=exp(-0.1)-exp(-16*0.1);
    E_n=Y_n-Z_n;
    RUNGE_KUTTA(k+1,1)=h;
    RUNGE_KUTTA(k+1,2)=E_n;
    vector_h(k+1)=log(h);
    vector_E(k+1)=log(abs(E_n));
    fprintf('h=%.16f, E_n=%.16f\n',h, E_n)
end 
plot(vector_h,vector_E,'b-')
xlabel('log(h)')
ylabel('log(E_n)')
title('a graph to show how the error changes with step size - RK4')