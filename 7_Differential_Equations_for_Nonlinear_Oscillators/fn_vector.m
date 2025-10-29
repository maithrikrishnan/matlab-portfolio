function [vector_derivative] = fn_vector(t,Y,a,b)
vector_derivative=[Y(2), (-a*Y(2) + Y(1) - (Y(1))^3 + b*cos(t))];
end