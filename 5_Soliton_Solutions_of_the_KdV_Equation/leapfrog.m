function leapfrog(A,delta,x_0,t_final)
h = 0.01; 
k = 0.0001; 
Delta = sqrt(12*(delta)^2/A);
x = 0:h:1;
sech = @(z) 1 ./ cosh(z);
u_0 = A*(sech((x-x_0)/Delta)).^2;

%for the first step, use the forward Euler method
u_x0 = -((2*A)/Delta)*((sech((x-x_0)/Delta)).^2).*tanh((x-x_0)/Delta);
u_xxx0_1 = ((16*A)/(Delta)^3)*((sech((x-x_0)/Delta)).^4).*tanh((x-x_0)/Delta);
u_xxx0_2 = ((8*A)/(Delta)^3)*((sech((x-x_0)/Delta)).^2).*(tanh((x-x_0)/Delta)).^3;
u_xxx0 = u_xxx0_1 - u_xxx0_2; 

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

c = A/3;
u_exact = A*(sech((x-x_0-c*t_final)/Delta).^2);

plot(x, u_new, 'b-','LineWidth', 1.5); hold on;
plot( x, u_exact, 'r-', 'LineWidth', 1.5)
legend('Numerical Solution', 'Exact Solution')
xlabel('x'), ylabel('u(x,0.5)')
title('KdV Solution at t = 0.5 using Leap-Frog')

%finding the root mean square error
error_rmse = sqrt(sum((u_new - u_exact).^2) * h);
fprintf("the estimated (root mean square) error is %.5e\n",error_rmse)

end