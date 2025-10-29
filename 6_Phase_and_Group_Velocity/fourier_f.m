function[f] = fourier_f(k)
f = -((16/(k^3))*sin(k)) - ((48/(k^4))*cos(k)) + ((48/(k^5))*sin(k));
end