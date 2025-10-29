psi_values_1 = linspace(-0.5, -0.05, 10); 
delta_values = linspace(0.0001,0.001,10);
phi_original=1/2;

for m = 1:10
    X_vector=zeros(1,10);
    Y_vector=zeros(1,10);
    phi = phi_original-delta_values(m);
    for n = 1:10
        psi_2 = psi_values_1(n);
        X=(1/pi)*(1+(exp(2*pi*psi_2)*cos(2*pi*phi)))+2*psi_2;
        Y=(1/pi)*(exp(2*pi*psi_2)*sin(2*pi*phi))+2*phi;
        X_vector(n)=X;
        Y_vector(n)=Y;
    end
    plot(X_vector,Y_vector)
    hold on;
end

xlabel('X');
ylabel('Y');
title('Equipotential Lines Below The Plate');