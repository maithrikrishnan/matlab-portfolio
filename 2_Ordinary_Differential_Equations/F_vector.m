function [vector_derivative] = F_vector(t,Y,delta)
vector_derivative=[Y(2), (-((delta)^3)*((Y(1))^2)*Y(2)-Y(1)+sin(t))];
end