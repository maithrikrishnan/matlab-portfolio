function vdp_phase_plot(h,t_max, a,b, ICs)
num_conditions = size(ICs, 1);

colors = lines(num_conditions);
legend_entries = cell(1, num_conditions);
figure; 
hold on; 

for ic = 1:num_conditions
    Y_n(1)=ICs(ic,1);
    Y_n(2)=ICs(ic,2)+a*((Y_n(1)^3)/3 - Y_n(1));
    t_n=0;
    vector_Y=zeros(1,1+(t_max/h));
    vector_Z=zeros(1,1+(t_max/h));
    vector_Y(1)=Y_n(1);
    vector_Z(1)=ICs(ic,2);
    for n=1:t_max/h
        k1=h*forced_vdp_vector(t_n,Y_n,a,b);
        k2=h*forced_vdp_vector(t_n+(h/2), Y_n+((k1)/2),a,b);
        k3=h*forced_vdp_vector(t_n+(h/2), Y_n+((k2)/2),a,b);
        k4=h*forced_vdp_vector(t_n+h, Y_n+k3,a,b);
        Y_n=Y_n + (k1+2*k2+2*k3+k4)/6;
        t_n=t_n + h; 
        derivative = Y_n(2) - a*((Y_n(1)^3)/3 - Y_n(1));
        vector_Y(n+1)=Y_n(1);
        vector_Z(n+1)=derivative;
    end
    plot(vector_Y,vector_Z,'LineWidth',1.5,'Color',colors(ic,:))
    legend_entries{ic} = sprintf('x₀ = %.1f, ẋ₀ = %.1f', ICs(ic,1), ICs(ic,2));
end
xlabel('x')
ylabel('dx/dt')
title(sprintf('Forced van der Pol oscillator: Phase portrait, a=%d, b=%.3f', a, b))
legend(legend_entries);
end