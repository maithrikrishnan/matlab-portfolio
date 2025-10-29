function leapfrog_3(delta,t_final,k)
h = 0.01; 
x = 0:h:1;
u_0 = sin(2*pi*x);

%for the first step, use the forward Euler method
u_x0 = 2*pi*cos(2*pi*x);
u_xxx0 = -8*((pi)^3)*cos(2*pi*x); 

u = u_0 - k*(u_0.*u_x0 + ((delta)^2)*u_xxx0);

u_prev = u_0;
u_new = zeros(size(x));

umax = max(abs(u_0));
k_max = h^3 / (4*(delta)^2 + h^2 * umax);
if k > k_max
    fprintf("Your time step k is too big — instability likely!");
end

N = length(x);

for n = 2:round(t_final/k)
    for m = 1:N
        mp2 = mod(m+1, N) + 1;         
        mp1 = mod(m, N) + 1;           
        mm1 = mod(m-2, N) + 1;        
        mm2 = mod(m-3, N) + 1;         

        u_new(m) = u_prev(m) - (k / (3*h)) * (u(mp1) + u(m) + u(mm1)) * (u(mp1) - u(mm1)) - ((k * delta^2) / h^3) * (u(mp2) - 2*u(mp1) + 2*u(mm1) - u(mm2));
    end 
    u_prev = u;
    u = u_new;
end 

plot(x, u_new, 'b-','LineWidth', 1.5)
xlabel('x'), ylabel('u(x)')
title(sprintf('KdV Solution using Leap-Frog, t = %.2f', t_final))

end