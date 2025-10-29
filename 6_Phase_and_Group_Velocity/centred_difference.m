function centred_difference(dx,lambda,t_finals,q)
%lambda is the ratio of the time step and the spatial step
dt = dx*lambda;
x_final = 60; 
N = x_final/dx;
x = linspace(0,x_final,N);

u_0 = zeros(size(x));

%set initial condition u(x,0) = f(x)
for n = 1:length(x)
    if x(n) < 1
        u_0(n) = (1 - x(n)^2).^2;
    else 
        u_0(n) = 0;
    end 
end 

colors = lines(length(t_finals));
legend_entries = cell(1, length(t_finals));

figure;
hold on;

u_new = zeros(size(x));

for k = 1:length(t_finals)
    t_final = t_finals(k);
    if t_final == 0 
        u_new = u_0;
    else 
        %first time step using the given symmetric condition 
        u = zeros(size(x));
        for i = 1:N-1
            if i == 1
                u(i) = u_0(i) + ((lambda)^2)*(u_0(i+1)-u_0(i)) - 0.5*((q*dt)^2)*(u_0(i+1));
            else
                u(i) = u_0(i) + 0.5*((lambda)^2)*(u_0(i+1)-2*u_0(i)+u_0(i-1)) - 0.25*((q*dt)^2)*(u_0(i+1)+u_0(i-1));
            end
        end

        u_prev = u_0;

        for j = 2:t_final/dt
            for i = 1:N-1 
                if i == 1
                    u_new(i) = 2*u(i) - u_prev(i) + 2*((lambda)^2)*(u(i+1)-u(i)) - ((q*dt)^2)*(u(i+1));
                else 
                    u_new(i) = 2*u(i) - u_prev(i) + ((lambda)^2)*(u(i+1)-2*u(i)+u(i-1)) - 0.5*((q*dt)^2)*(u(i+1)+u(i-1));
                end
            end
            u_prev = u;
            u = u_new;
        end 
    end
    plot(x, u_new,'LineWidth', 1.5, 'Color', colors(k,:)); 
    legend_entries{k} = sprintf('t = %d', t_final);
end 

xlabel('x'), ylabel('u(x)')
title('Finite Difference Solution of the Wave Equation');
legend(legend_entries);

end