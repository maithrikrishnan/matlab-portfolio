function [contour_plot] = contours_2(m)
[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2); 
r1=sqrt((X+1-m).^2 +Y.^2);
r2=sqrt((X-m).^2 +Y.^2);
Omega =-(0.5*m*(r1).^2)-(0.5*(1-m)*(r2).^2)-(m./r1)-((1-m)./r2); 
J=-5;
levels=zeros(1,200);
for i=1:100
    levels(2*i-1)=J+i/30;
    levels(2*i)=J-i/30;
end
contour_plot=contour(X, Y, Omega,levels);
xlabel('x-axis');
ylabel('y-axis');
title('Contour Plot of Omega');
end