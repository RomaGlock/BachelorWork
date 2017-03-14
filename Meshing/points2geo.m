function result = points2geo(fileName, points)
    fileID = fopen(fileName, 'w');
    fprintf(fileID,'cl__1 = 1;\n');
    for i = 1 : length(points)
        fprintf(fileID,'Point(%d) = {%f, %f, %f, 1};\n', i, points(i,1), points(i, 2), points(i, 3));
    end;
    result = fclose(fileID);
end
