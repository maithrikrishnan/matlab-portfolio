function [vector_derivative] = fn_vector_2(t,Y,a,b)
vector_derivative=[Y(2), ((b - (Y(1))^2)*Y(2) + a*Y(1) - (Y(1))^3)];
end