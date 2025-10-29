function [RUNGE_KUTTA_2] = RK4_vector(h,t_max,delta)
RUNGE_KUTTA_2 = zeros(t_max/h,2);
Y_n=zeros(1,2);
t_n=0;
vector_T=zeros(1,1+(t_max/h));
vector_Y=zeros(1,1+(t_max/h));
for n=1:t_max/h
    k1=h*F_vector(t_n,Y_n,delta);
    k2=h*F_vector(t_n+(h/2), Y_n+((k1)/2),delta);
    k3=h*F_vector(t_n+(h/2), Y_n+((k2)/2),delta);
    k4=h*F_vector(t_n+h, Y_n+k3,delta);
    Y_n=Y_n + (k1+2*k2+2*k3+k4)/6;
    t_n=t_n + h; 
    vector_T(n+1)=t_n;
    vector_Y(n+1)=Y_n(1);
    RUNGE_KUTTA_2(n,1)=t_n;
    RUNGE_KUTTA_2(n,2)=Y_n(1);
    fprintf('t_n=%.1f, Y_n=%17.10f\n',t_n, Y_n(1))
end
plot(vector_T,vector_Y,'b-')
xlabel('t')
ylabel('y')
title('estimation using RK4 - with coupled first order DEs')