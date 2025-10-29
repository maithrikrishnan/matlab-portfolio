function [TABLEAU] = estimation(h)
Y_n=0;
x_n=0;
TABLEAU=zeros(6/h,5);
for n=1:6/h
    Y_n=Y_n + h*f(x_n,Y_n);
    x_n=x_n+h;
    Z_n=exp(-x_n)-exp(-16*x_n);
    E_n=Y_n-Z_n;
    G_n=log(abs(E_n))/x_n;
    TABLEAU(n,1)=x_n;
    TABLEAU(n,2)=Y_n;
    TABLEAU(n,3)=Z_n;
    TABLEAU(n,4)=E_n;
    TABLEAU(n,5)=G_n;    
    fprintf('x_n=%.3f, Y_n=%17.10f, y(x_n)=%.7f, E_n=%17.10f, G_n=%.6f\n',x_n, Y_n, Z_n, E_n, G_n)
end 

