function explore_duffing_a(initial_condition)
a_values = [0.1, 0.2, 0.3, 0.4, 0.5];
h = pi/100;
t_max = 3200*pi;
b = 0.3;

for i = 1:length(a_values)
    figure;
    duffing_phase_plot_2(h, t_max, a_values(i), b, initial_condition);
    title(sprintf('Duffing Oscillator Phase Portrait (a = %.2f, b = 0.3)', a_values(i)));
end