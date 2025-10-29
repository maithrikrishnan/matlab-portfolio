function [RK4_output] = vector_RK4_plot(h,t_max,v,w)
RK4_output = zeros(t_max/h,4);
Y_n=zeros(1,2);
t_n=0;
k=(sqrt(4-v^2))/2;
d=(v*w)^2 +(1-w^2)^2;
B=(v*w)/((v*w)^2 + (1-w^2)^2);
A=((v*B)/2 - (w*(1-w^2))/d)/k;
vector_T=zeros(1,1+(t_max/h));
vector_Y=zeros(1,1+(t_max/h));
vector_Z=zeros(1,1+(t_max/h));
for n=1:t_max/h
    k1=h*vector_F(t_n,Y_n,v,w);
    k2=h*vector_F(t_n+(h/2), Y_n+((k1)/2),v,w);
    k3=h*vector_F(t_n+(h/2), Y_n+((k2)/2),v,w);
    k4=h*vector_F(t_n+h, Y_n+k3,v,w);
    Y_n=Y_n + (k1+2*k2+2*k3+k4)/6;
    t_n=t_n + h;
    Z_n=exp(-(v*t_n)/2)*(A*sin(k*t_n)+B*cos(k*t_n))+(1-w^2)*sin(w*t_n)/d - (v*w)*cos(w*t_n)/d;
    E_n=Y_n(1)-Z_n;
    vector_T(n+1)=t_n;
    vector_Y(n+1)=Y_n(1);
    vector_Z(n+1)=Z_n;
    RK4_output(n,1)=t_n;
    RK4_output(n,2)=Y_n(1);
    RK4_output(n,3)=Z_n;
    RK4_output(n,4)=E_n;
    fprintf('t_n=%.1f, Y_n=%17.10f, Z_n=%17.10f, E_n=%17.10f\n',t_n, Y_n(1), Z_n, E_n)
end
plot(vector_T,vector_Y,'b-')
hold on;
plot(vector_T,vector_Z,'r-')
xlabel('t')
ylabel('y')
title('estimation using RK4 - with coupled first order DEs')
legend('Y_n','y_n','Location','northeast')