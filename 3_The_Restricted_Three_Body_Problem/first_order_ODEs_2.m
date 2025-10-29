function [derivative] = first_order_ODEs_2(t,y_val,m)
x=y_val(1);
y=y_val(2); 
u=y_val(3); 
v=y_val(4);
r2=sqrt((x-m)^2 +y^2);
dwdx=(x-m)/2*((r2)^3);
dwdy=y/2*((r2)^3);
dxdt=u;
dydt=v;
dudt=2*v-dwdx;
dvdt=-2*u-dwdy; 
derivative=[dxdt; dydt; dudt; dvdt];
end