function[equilibrium] = stability(x0,y0,u0,v0,m)
y0=[x0;y0;u0;v0];
tspan = [0,500]; 
[t, y_val] = ode45(@(t,y_val) first_order_ODEs(t,y_val, m), tspan, y0);
equilibrium=table(t,y_val(:,1),y_val(:,2),'VariableNames',{'time','x','y'});
for i=1:20 
    R1=sqrt(((y_val(i,1)+1-m)^2)+(y_val(i,2))^2);
    R2=sqrt(((y_val(i,1)-m)^2)+(y_val(i,2))^2);
    Omega_val =-(0.5*m*(R1)^2)-(0.5*(1-m)*(R2)^2)-(m/R1)-((1-m)/R2); 
    J_val=0.5*((y_val(i,3))^2)+0.5*((y_val(i,4))^2)+Omega_val;
    sprintf(num2str(J_val))
end
plot(y_val(:,1), y_val(:,2),'-r') 
title('Trajectory of Omega')
xlabel('x')
ylabel('y');
end