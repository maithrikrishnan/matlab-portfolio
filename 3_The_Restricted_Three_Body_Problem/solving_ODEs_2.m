function solving_ODEs_2(m)
y0=[1;0;0;-1];
tspan = [0, 15]; 
[~, y_val] = ode45(@(t,y_val) first_order_ODEs_2(t,y_val, m), tspan, y0);
for i=1:20 
    R1=sqrt(((y_val(i,1)+1-m)^2)+(y_val(i,2))^2);
    R2=sqrt(((y_val(i,1)-m)^2)+(y_val(i,2))^2);
    Omega_val =-(0.5*m*(R1)^2)-(0.5*(1-m)*(R2)^2)-(m/R1)-((1-m)/R2); 
    J_val=0.5*((y_val(i,3))^2)+0.5*((y_val(i,4))^2)+Omega_val;
    sprintf(num2str(J_val))
end
plot(y_val(:,1), y_val(:,2))
title('Trajectory of Omega');
xlabel('x');
ylabel('y');
end