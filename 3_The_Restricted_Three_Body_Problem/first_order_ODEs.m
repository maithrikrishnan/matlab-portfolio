function [derivative] = first_order_ODEs(t,y_val,m)
x=y_val(1);
y=y_val(2); 
u=y_val(3); 
v=y_val(4);
r1=sqrt((x+1-m)^2 +y^2);
r2=sqrt((x-m)^2 +y^2);
dwdx=-(m*(x+1-m))-((1-m)*(x-m))+(m*(x+1-m)/(r1)^3)+((1-m)*(x-m)/(r2)^3);
dwdy=-(m*y)-((1-m)*y)+(m*y/(r1)^3)+((1-m)*y/(r2)^3);
dxdt=u;
dydt=v;
dudt=2*v-dwdx;
dvdt=-2*u-dwdy; 
derivative=[dxdt; dydt; dudt; dvdt];
end