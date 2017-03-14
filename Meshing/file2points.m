function points = file2points(fileName)

    fileID = fopen(fileName, 'r');
    pointsSize = [3 inf];
    points = fscanf(fileID, '%f %f %f', pointsSize);
    fclose(fileID);
    points = points';
end