function [vector_derivative] = forced_vdp_vector(t,Y,a,b)
vector_derivative=[Y(2) - a*((Y(1)^3)/3 - Y(1)), -Y(1) + 1 + b];
end