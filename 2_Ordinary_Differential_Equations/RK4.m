function [RK4_results] = RK4(h) 
Y_n=0;
x_n=0;
RK4_results=zeros(4/h,3);
vector_X=zeros(1,1+(4/h));
vector_Y=zeros(1,1+(4/h));
vector_Z=zeros(1,1+(4/h));
for n=1:4/h
    k1=h*f(x_n,Y_n);
    k2=h*f((x_n+(h/2)), (Y_n+((k1)/2)));
    k3=h*f((x_n+(h/2)), (Y_n+((k2)/2)));
    k4=h*f((x_n+h), (Y_n+k3));
    Y_n=Y_n + (k1+2*k2+2*k3+k4)/6;
    x_n=x_n+h;
    Z_n=exp(-x_n)-exp(-16*x_n);
    vector_X(n+1)=x_n;
    vector_Y(n+1)=Y_n;
    vector_Z(n+1)=Z_n;
    RK4_results(n,1)=x_n;
    RK4_results(n,2)=Y_n;
    RK4_results(n,3)=Z_n;
    fprintf('x_n=%.2f, Y_n=%17.10f, y(x_n)=%.7f\n',x_n, Y_n, Z_n)
end
plot(vector_X,vector_Y,'b-')
hold on;
plot(vector_X,vector_Z,'r-')
xlabel('x')
ylabel('y')
title('estimation using the Runge-Kutta method')
legend('Y_n','y_n','Location','northeast')