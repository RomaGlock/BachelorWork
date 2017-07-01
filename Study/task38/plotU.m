clear all;

fileName = 'mesh.config';
mode = 'r';
formatSpec = ' ';
file = fopen(fileName, mode);
h = fscanf(file, 'Hx = %f N = %f\nHy = %f M = %f');
fclose(file);

hx = h(1); Nx = h(2);
hy = h(3); Ny = h(4);

fileName = 'U0.dat';
mode = 'r';
formatSpec = '%f ';
sizeU = [Nx Ny];
file = fopen(fileName, mode);
U0 = fscanf(file, formatSpec, sizeU);
fclose(file);

fileName = 'U32.dat';
mode = 'r';
formatSpec = '%f ';
sizeU = [Nx Ny];
file = fopen(fileName, mode);
U32 = fscanf(file, formatSpec, sizeU);
fclose(file);

fileName = 'U64.dat';
mode = 'r';
formatSpec = '%f ';
sizeU = [Nx Ny];
file = fopen(fileName, mode);
U64 = fscanf(file, formatSpec, sizeU);
fclose(file);

fileName = 'U128.dat';
mode = 'r';
formatSpec = '%f ';
sizeU = [Nx Ny];
file = fopen(fileName, mode);
U128 = fscanf(file, formatSpec, sizeU);
fclose(file);


[X Y] = meshgrid(0:hx:hx * (Nx - 1),0:hy:hy * (Ny - 1));

figure;
axis([0 1 0 2 -1 1]);
surf(X, Y, U0);
title('U(X,Y), T = 0');
xlabel('X');
ylabel('Y');
zlabel('U(X,Y)');
 

figure;
surf(X, Y, U32); axis([0 1 0 2 -1 1]);
title('U(X,Y), T = 0.025');
xlabel('X');
ylabel('Y');
zlabel('U(X,Y)');

figure
surf(X, Y, U64); axis([0 1 0 2 -1 1]);
title('U(X,Y), T = 0.05');
xlabel('X');
ylabel('Y');
zlabel('U(X,Y)');

figure
surf(X, Y, U128); axis([0 1 0 2 -1 1]);
title('U(X,Y), T = 0.1');
xlabel('X');
ylabel('Y');
zlabel('U(X,Y)');





