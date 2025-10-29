function[phi] = SOR_3_efield_plate_q9(L, DX, DY, h, omega, tol)
Nx = DX / h; 
Ny = DY / h; 
N=(Nx+1)*(Ny+1);
x = linspace(0, DX, Nx+1); 
y = linspace(0, DY, Ny+1); 

phi = ones(Nx+1, Ny+1);
phi(Nx+1,:) = 0; 
phi(:,1) = 0; 
phi(:,Ny+1) = 0; 

max_iter = 10000; 
    
for k=1:max_iter
    residual=0;
    phi_old = phi;
    for i = 1:Nx
        for j = 2:Ny
            if i==1
                phi(i,j) = (1 - omega) * phi(i,j) + (omega/4) * (phi(i+1,j) + phi(i+1,j) + phi(i,j-1) + phi(i,j+1));
            else
                phi(i,j) = (1 - omega) * phi(i,j) + (omega/4) * (phi(i-1,j) + phi(i+1,j) + phi(i,j-1) + phi(i,j+1));
            end
            if j==(1/h)+1 && i<=(L/h)+1
                phi(i,j)=1/2;
            end
            residual = residual+abs(phi(i,j) - phi_old(i,j));
        end
    end
    residual=residual/N;
    if residual < tol
        break;
    end
end

psi_values = linspace(-1/2, 0, 11); 
X_vector=zeros(1,11);
e_field_vector=zeros(1,11);
for j = 1:11
    psi = psi_values(j);
    X=(1/pi)*(1-exp(2*pi*psi))+2*psi+L;
    X_vector(j)=X;
    e_field=-0.5*(1-exp(2*pi*psi))^-1;
    e_field_vector(j)=e_field;
end

psi_values_2 = linspace(0, 1/2, 11); 
X_vector_2=zeros(1,11);
e_field_vector_2=zeros(1,11);
for m = 1:11
    psi_2 = psi_values_2(m);
    X_2=(1/pi)*(1-exp(2*pi*psi_2))+2*psi_2+L;
    X_vector_2(m)=X_2;
    e_field_2=-0.5*(1-exp(2*pi*psi_2))^-1;
    e_field_vector_2(m)=e_field_2;
end

e_field_Y1_lower=(1/h)*(phi);
e_field_Y1_lower_2=e_field_Y1_lower(:,((1-h)/h)+1)-e_field_Y1_lower(:,(1/h)+1);
e_field_Y1_higher=-(1/h)*(phi);
e_field_Y1_higher_2=e_field_Y1_higher(:,((1+h)/h)+1)-e_field_Y1_higher(:,(1/h)+1);

subplot(1, 2, 1);
plot(y, e_field_Y1_lower_2, 'b-');
hold on; 
plot(X_vector,e_field_vector)
legend('finite (numerical)','semi-finite (analytic)')
title('Plot of \epsilon_Y(X,1) lower surface');
xlabel('X');
ylabel('\epsilon_Y(X, 1)');

subplot(1, 2, 2);
plot(y, e_field_Y1_higher_2, 'b-');
hold on;
plot(X_vector_2,e_field_vector_2)
legend('finite (numerical)','semi-finite (analytic)')
title('Plot of \epsilon_Y(X,1) upper surface');
xlabel('X');
ylabel('\epsilon_Y(X, 1)');