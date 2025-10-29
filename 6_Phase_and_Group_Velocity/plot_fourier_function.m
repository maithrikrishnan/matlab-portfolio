k_vector = linspace(0,50,500); 
fourier_f_vector = zeros(size(k_vector));
for k = 1:length(k_vector)
    fourier_f_vector(k) = fourier_f(k_vector(k));
end 

plot(k_vector, fourier_f_vector, 'b-');
xlabel('k'), ylabel('$\tilde{f}$(k)', 'Interpreter', 'latex');
title('Fourier transform of f(x)')
grid on; 