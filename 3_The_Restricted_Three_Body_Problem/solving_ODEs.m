function solving_ODEs(m,v_val)
y0=[0.2;0;0;v_val];
r_1=sqrt((1.2-m)^2);
r_2=sqrt((0.2-m)^2);
omega=-(0.5*m*(r_1)^2)-(0.5*(1-m)*(r_2)^2)-(m/(r_1))-((1-m)/(r_2)); 
J=0.5*(v_val)^2+omega; 
[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2); 
r1=sqrt((X+1-m).^2 +Y.^2);
r2=sqrt((X-m).^2 +Y.^2);
Omega =-(0.5*m*(r1).^2)-(0.5*(1-m)*(r2).^2)-(m./r1)-((1-m)./r2); 
levels=[J,J];
levels_1=[J-1,J-1];
levels_2=[J-2,J-2];
tspan = [0, 15]; 
[~, y_val] = ode45(@(t,y_val) first_order_ODEs(t,y_val, m), tspan, y0);
for i=1:20 
    R1=sqrt(((y_val(i,1)+1-m)^2)+(y_val(i,2))^2);
    R2=sqrt(((y_val(i,1)-m)^2)+(y_val(i,2))^2);
    Omega_val =-(0.5*m*(R1)^2)-(0.5*(1-m)*(R2)^2)-(m/R1)-((1-m)/R2); 
    J_val=0.5*((y_val(i,3))^2)+0.5*((y_val(i,4))^2)+Omega_val;
    sprintf(num2str(J_val))
end
plot(y_val(:,1), y_val(:,2),'-r') 
hold on; 
contour(X, Y, Omega, levels,'-y');
hold on; 
contour(X, Y, Omega, levels_1,'-g');
hold on; 
contour(X, Y, Omega, levels_2,'-b');
legend('trajectory','Omega=J','Omega=J-1','Omega=J-2')
title('Trajectory of Omega')
xlabel('x')
ylabel('y');
end