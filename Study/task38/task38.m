[X, Y] = meshgrid(0:0.01:1,0:0.02:2);

T = 0.01;

U = exp(-65/16*pi*pi*T)*sin(2*pi*X).*sin(pi/4*Y);

surf(X,Y,U); axis([0 1 0 2 -1 1]);