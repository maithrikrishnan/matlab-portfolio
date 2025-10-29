function duffing_phase_plot_2(h,t_max, a,b, ICs)
num_conditions = size(ICs, 1);
colors = lines(num_conditions);
legend_entries = cell(1, num_conditions);
figure; 
hold on; 

for ic = 1:num_conditions
    Y_n=ICs(ic,:);
    t_n=0;
    vector_Y=[Y_n(1)];
    vector_Z=[Y_n(2)];
    for n=1:t_max/h
        k1=h*fn_vector(t_n,Y_n,a,b);
        k2=h*fn_vector(t_n+(h/2), Y_n+((k1)/2),a,b);
        k3=h*fn_vector(t_n+(h/2), Y_n+((k2)/2),a,b);
        k4=h*fn_vector(t_n+h, Y_n+k3,a,b);
        Y_n=Y_n + (k1+2*k2+2*k3+k4)/6;
        t_n=t_n + h; 
        if mod(t_n, 2*pi) < h/2
            vector_Y=[vector_Y, Y_n(1)];
            vector_Z=[vector_Z, Y_n(2)];
        end
    end
    plot(vector_Y,vector_Z,'x','LineWidth',1.5,'Color',colors(ic,:))
    legend_entries{ic} = sprintf('x₀ = %.1f, ẋ₀ = %.1f', ICs(ic,1), ICs(ic,2));
end
xlabel('x')
ylabel('dx/dt')
title('Phase Portrait of the Duffing Oscillator')
legend(legend_entries);
end