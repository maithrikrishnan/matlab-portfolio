function [derivative_vector] = vector_F(t,Y,v,w)
derivative_vector=[Y(2), (-v*Y(2)-Y(1)+sin(w*t))];
end