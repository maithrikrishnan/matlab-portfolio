function leapfrog_2(t_final)
h = 0.01; 
k = 0.0001; 
delta = 0.03;
x = 0:h:1;

A1 = 2;
x0_1 = 0.25; 
Delta_1 = sqrt(12*(delta)^2/A1);

A2 = 1; 
x0_2 = 0.75; 
Delta_2 = sqrt(12*(delta)^2/A2);

sech = @(z) 1 ./ cosh(z);
u0_1 = A1*(sech((x-x0_1)/Delta_1)).^2;
u0_2 = A2*(sech((x-x0_2)/Delta_2)).^2;

%combined initial condition 
u_0 = u0_1 + u0_2;

%initialise a vector for mass and energy for the masses and energies to be
%calculated at each time step (indexed based on the value of n)
masses = zeros(1, round(t_final/k));
energies = zeros(1, round(t_final/k)); 
timesteps = k:k:t_final;

%for the first step, use the forward Euler method
%have to do separately for first soliton and second soliton and then add 

u_x0_1 = -((2*A1)/Delta_1)*((sech((x-x0_1)/Delta_1)).^2).*tanh((x-x0_1)/Delta_1);
u_x0_2 = -((2*A2)/Delta_2)*((sech((x-x0_2)/Delta_2)).^2).*tanh((x-x0_2)/Delta_2);
u_x0 = u_x0_1 + u_x0_2;

u_xxx0_1_1 = ((16*A1)/(Delta_1)^3)*((sech((x-x0_1)/Delta_1)).^4).*tanh((x-x0_1)/Delta_1);
u_xxx0_2_1 = ((8*A1)/(Delta_1)^3)*((sech((x-x0_1)/Delta_1)).^2).*(tanh((x-x0_1)/Delta_1)).^3;
u_xxx0_1 = u_xxx0_1_1 - u_xxx0_2_1; 

u_xxx0_1_2 = ((16*A2)/(Delta_2)^3)*((sech((x-x0_2)/Delta_2)).^4).*tanh((x-x0_2)/Delta_2);
u_xxx0_2_2 = ((8*A2)/(Delta_2)^3)*((sech((x-x0_2)/Delta_2)).^2).*(tanh((x-x0_2)/Delta_2)).^3;
u_xxx0_2 = u_xxx0_1_2 - u_xxx0_2_2; 

u_xxx0 = u_xxx0_1 + u_xxx0_2; 

u = u_0 - k*(u_0.*u_x0 + ((delta)^2)*u_xxx0);

%need to find the mass and energy for this first time step: index 1 
mass_sum_1 = 0;
for i = 2:length(x)-1
    mass_sum_1 = mass_sum_1 + u(i);
end
masses(1) = (h/2)*((u(1)+u(length(x))) + 2*mass_sum_1);

v = 0.5*(u).^2; 

energy_sum_1 = 0; 
for j = 2:length(x)-1
    energy_sum_1 = energy_sum_1 + v(j); 
end
energies(1) = (h/2)*((v(1)+v(length(x))) + 2*energy_sum_1);

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

    mass_sum = 0;
    for i = 2:length(x)-1
        mass_sum = mass_sum + u(i);
    end

    masses(n) = (h/2)*((u(1)+u(length(x))) + 2*mass_sum); 

    v_new = 0.5*(u).^2; 

    energy_sum = 0; 
    for j = 2:length(x)-1
        energy_sum = energy_sum + v_new(j); 
    end

    energies(n) = (h/2)*((v_new(1)+v_new(length(x))) + 2*energy_sum); 
end 

subplot(1,2,1)
plot(x, u_new, 'b-','LineWidth', 1.5); 
xlabel('x'), ylabel('u(x)')
title(sprintf('KdV Solution using Leap-Frog, t=%.2f', t_final))

subplot(1,2,2)
plot(timesteps, masses, 'b-', 'LineWidth',1.5); hold on;
plot(timesteps, energies, 'r-', 'LineWidth',1.5);
legend('Mass', 'Energy')
xlabel('t'); ylabel('Mass/Energy');
xlim([0 t_final])
ylim([0 1])
title('Mass & Energy vs. Time');

end