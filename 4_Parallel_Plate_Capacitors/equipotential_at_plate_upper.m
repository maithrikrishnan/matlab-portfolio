psi_values_2 = linspace(0.05, 0.5, 10);
delta_values = linspace(0.0001,0.001,10);
phi_original=1/2;

for p = 1:10
    X_vector=zeros(1,10);
    Y_vector=zeros(1,10);
    phi = phi_original-delta_values(p);
    for q = 1:10
        psi_2 = psi_values_2(q);
        X=(1/pi)*(1+(exp(2*pi*psi_2)*cos(2*pi*phi)))+2*psi_2;
        Y=(1/pi)*(exp(2*pi*psi_2)*sin(2*pi*phi))+2*phi;
        X_vector(q)=X;
        Y_vector(q)=Y;
    end
    plot(X_vector,Y_vector)
    hold on;
end

xlabel('X');
ylabel('Y');
title('Equipotential Lines Above The Plate');