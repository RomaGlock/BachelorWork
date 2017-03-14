bottom_point = [40 0 -32];
bottom_point = bottom_point * roty(7);
bottom = bottom_point;
%%

for i = 1 : 17
   bottom = [bottom; bottom_point*rotz(i*20)];
end