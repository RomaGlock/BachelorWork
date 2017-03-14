function bspline(start_point,end_point)
    fileID = fopen('bspline','w');
    for i = 1 : (length(start_point)-1)
        fprintf(fileID, '//+\n');
        fprintf(fileID, 'BSpline(%d) = {', i*3-2);
        fprintf(fileID, '%d, ', start_point(i):(end_point(i)-1));
        fprintf(fileID, '%d};\n', end_point(i));
        fprintf(fileID, '//+\n');
        fprintf(fileID, 'Circle(%d) = {%d, 365, %d};\n', i*3-1, start_point(i), start_point(i + 1));
        fprintf(fileID, '//+\n');
        fprintf(fileID, 'Circle(%d) = {%d, 365, %d};\n', i*3, end_point(i), end_point(i + 1));
    end
    fclose(fileID);
end