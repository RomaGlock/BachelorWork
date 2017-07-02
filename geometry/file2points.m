function result = file2points()

%		Скрипт для обработки результатов оцифровки профилограм
% обгарной формы поверхности полученных на стенде В1-А

		filesList = dir('dataFiles/line*.dat');
		filesCount = length(filesList);
		stepAngle = 180 / filesCount;
		modelRadius = 80;
		
		abletiveSurface = [];

		for i = 1 : filesCount
				
				% Открываем файл для чтения профиля
				fileName = strcat('dataFiles/',filesList(i).name);
				fileID = fopen(fileName, 'r');
				fileSize = [2 inf];
				points = fscanf(fileID, '%f %f', fileSize)';
				fclose(fileID);
				% Формируем профиль из прочитанных точек
				z = vertcat(0, points(:,1), 80);
				y = zeros(81,1);
				x = vertcat(0, points(:,2), 0);
				x = -x + max(points(:,2));
				% Интерполируем для выравнивания точек
				x = interp1(z,x,(0:80)');
				z = (-40:40)';
				% Добавляем нижние точки
				bottom_point1 = [60 0 -40] * roty(-7);
				bottom_point2 = [60 0 40] * roty(7);
				profile = [bottom_point1; horzcat(x,y,z); bottom_point2];			
				% Поворачиваем на угол профилограммы
				rotMatrix = rotx((i - 1) * stepAngle);
				profile = profile * rotMatrix;
				% Добовляем к основному массиву
				abletiveSurface = vertcat(abletiveSurface, profile);
				
		end
		
		points2geo('abletiveSurface.geo', abletiveSurface);
end
		