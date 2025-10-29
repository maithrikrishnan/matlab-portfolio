psi_values = linspace(-1/2, 1/2, 11); 
phi_values = linspace(-1/2,1/2,11);
for i = 1:11
    X_vector=zeros(1,11);
    Y_vector=zeros(1,11);
    psi = psi_values(i);
    for j = 1:11
        phi = phi_values(j);
        X=(1/pi)*(1+(exp(2*pi*psi)*cos(2*pi*phi)))+2*psi;
        Y=(1/pi)*(exp(2*pi*psi)*sin(2*pi*phi))+2*phi;
        X_vector(j)=X;
        Y_vector(j)=Y;
    end
    plot(X_vector,Y_vector)
    hold on;
end
xlabel('X');
ylabel('Y');
title('Field Lines');