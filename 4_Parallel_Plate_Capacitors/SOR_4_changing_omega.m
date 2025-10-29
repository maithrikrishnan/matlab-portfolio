function[phi] = SOR_4_changing_omega(L, DX, DY, h, tol)
Nx = DX / h; 
Ny = DY / h; 
N=(Nx+1)*(Ny+1);
x = linspace(0, DX, Nx+1); 
y = linspace(0, DY, Ny+1); 
k_vector=zeros(1,100);
omega_vector=linspace(1,1.99,100);

for omega_index=1:length(omega_vector)
    omega = omega_vector(omega_index);
    phi = zeros(Nx+1, Ny+1);
    phi(Nx+1,:) = 0; 
    phi(:,1) = 0; 
    phi(:,Ny+1) = 0; 
    phi(1,(L/h)+1)=0.5;
    
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
    k_vector(omega_index)=k;
end
disp(k_vector)    
disp(omega_vector)
plot(omega_vector,k_vector,'-b')
title('How the no. of iterations needed for convergence depends on \omega' )
xlabel('\omega');
ylabel('Number of iterations needed for convergence');
end