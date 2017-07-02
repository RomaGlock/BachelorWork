function result = file2points()

%		������ ��� ��������� ����������� ��������� �����������
% �������� ����� ����������� ���������� �� ������ �1-�

		filesList = dir('dataFiles/line*.dat');
		filesCount = length(filesList);
		stepAngle = 180 / filesCount;
		modelRadius = 80;
		
		abletiveSurface = [];

		for i = 1 : filesCount
				
				% ��������� ���� ��� ������ �������
				fileName = strcat('dataFiles/',filesList(i).name);
				fileID = fopen(fileName, 'r');
				fileSize = [2 inf];
				points = fscanf(fileID, '%f %f', fileSize)';
				fclose(fileID);
				% ��������� ������� �� ����������� �����
				z = vertcat(0, points(:,1), 80);
				y = zeros(81,1);
				x = vertcat(0, points(:,2), 0);
				x = -x + max(points(:,2));
				% ������������� ��� ������������ �����
				x = interp1(z,x,(0:80)');
				z = (-40:40)';
				% ��������� ������ �����
				bottom_point1 = [60 0 -40] * roty(-7);
				bottom_point2 = [60 0 40] * roty(7);
				profile = [bottom_point1; horzcat(x,y,z); bottom_point2];			
				% ������������ �� ���� �������������
				rotMatrix = rotx((i - 1) * stepAngle);
				profile = profile * rotMatrix;
				% ��������� � ��������� �������
				abletiveSurface = vertcat(abletiveSurface, profile);
				
		end
		
		points2geo('abletiveSurface.geo', abletiveSurface);
end
		