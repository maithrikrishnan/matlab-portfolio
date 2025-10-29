function[phi] = SOR(L, DX, DY, h, omega, tol)
Nx = DX / h; 
Ny = DY / h; 
N=(Nx+1)*(Ny+1);
x = linspace(0, DX, Nx+1); 
y = linspace(0, DY, Ny+1); 
   
phi = zeros(Nx+1, Ny+1);
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

fprintf(num2str(k))
    
[Y, X] = meshgrid(y, x);
surf(Y, X,phi); 
title('solution \Phi(x, y)');
xlabel('y');
ylabel('x');
zlabel('\Phi');
colorbar;
end