function[phi] = SOR_3_efield_plate(L, DX, DY, h, omega, tol)
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

e_field_Y1_lower=(1/h)*(phi);
e_field_Y1_lower_2=e_field_Y1_lower(:,((1-h)/h)+1)-e_field_Y1_lower(:,(1/h)+1);
e_field_Y1_higher=-(1/h)*(phi);
e_field_Y1_higher_2=e_field_Y1_higher(:,((1+h)/h)+1)-e_field_Y1_higher(:,(1/h)+1);

subplot(1, 2, 1);
plot(y, e_field_Y1_lower_2, 'b-');
title('Plot of \epsilon_Y(X,1) lower surface');
xlabel('X');
ylabel('\epsilon_Y(X, 1)');
grid on;

subplot(1, 2, 2);
plot(y, e_field_Y1_higher_2, 'r-');
title('Plot of \epsilon_Y(X,1) upper surface');
xlabel('X');
ylabel('\epsilon_Y(X, 1)');
grid on;