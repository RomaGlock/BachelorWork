function line = rotateLine(data1, data2, angle)

    x = vertcat(0, data1, 80);
    y = zeros(81,1);
    z = vertcat(0, (data2 - (max(data2) - 16)), 0);
    
    z = interp1(x,z,((1:81)' - 1));
    x = ((1:81)' - 41);
    line = horzcat(x,y,z);
    
    Rot = rotz(angle);
    line = line * Rot;
    
end