function centred_difference_2(dx,lambda,t_final,q)
%lambda is the ratio of the time step and the spatial step
x_final = t_final + 10;
dt = dx*lambda; 
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

u_new = zeros(size(x));

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

u_approx = zeros(size(x));

for n = 1:length(x)
    x_n = x(n);
    if x_n < t_final
        u_approx_1 = (q^0.5 * t_final)/((2*pi)^0.5 * (t_final^2 - x_n^2).^0.75);
        u_approx_2 = fourier_f((x_n*q)/sqrt(t_final^2 - x_n^2)); 
        u_approx_3 = cos(q*sqrt(t_final^2 - x_n^2) + pi/4);
        u_approx(n) = u_approx_1*u_approx_2*u_approx_3;
    else 
        u_approx(n) = 0;
    end
end 

plot(x, u_new, 'b-','LineWidth', 1.5); hold on; 
plot(x, u_approx, 'r-','LineWidth', 1.5)
legend('Finite Difference', 'Approximation')
xlabel('x'), ylabel('u(x)')
title('Finite Difference Solution of the Wave Equation');

end