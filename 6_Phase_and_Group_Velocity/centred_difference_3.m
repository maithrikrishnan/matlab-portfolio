function centred_difference_3(dx,lambda,t_final,q,w0)
%lambda is the ratio of the time step and the spatial step
x_final = t_final; 
dt = dx*lambda; 
N = x_final/dx;
x = linspace(0,x_final,N);
u_0 = zeros(size(x));

colors = lines(length(w0));
legend_entries = cell(1, length(w0));

figure;
hold on;
grid on;

u_new = zeros(size(x));

for k = 1:length(w0)
    w_0 = w0(k);
    if t_final == 0 
        u_new = u_0;
    else 
        %first time step using the given symmetric condition 
        u = zeros(size(x));
        for i = 1:N-1
            if i == 1
                u(i) = 0;
            else
                u(i) = u_0(i) + 0.5*((lambda)^2)*(u_0(i+1)-2*u_0(i)+u_0(i-1)) - 0.25*((q*dt)^2)*(u_0(i+1)+u_0(i-1));
            end
        end

        u_prev = u_0;

        for j = 2:t_final/dt
            for i = 1:N-1 
                if i == 1
                %this time if i==1, it's u(0,t) so use the initial condition
                    u_new(i) = sin(w_0*(j-1)*dt);
                else 
                    u_new(i) = 2*u(i) - u_prev(i) + ((lambda)^2)*(u(i+1)-2*u(i)+u(i-1)) - 0.5*((q*dt)^2)*(u(i+1)+u(i-1));
                end
            end
            u_prev = u;
            u = u_new;
        end 
    end

    plot(x, u_new,'LineWidth', 1.5,'Color', colors(k,:));
    legend_entries{k} = sprintf(' Finite Difference, w = %.1f', w_0);

    if q == 0 
        u_exact = sin(w_0*(t_final-x));
        plot(x, u_new,'LineWidth', 1.5)
        plot(x, u_exact, 'r-', 'LineWidth', 1.5)
        ylim([-2 2]);
        legend_entries{end+1} = 'Exact';
    end 
end

xlabel('x'), ylabel('u(x)');
title('Finite Difference Solution of the Wave Equation');
legend(legend_entries);

end 