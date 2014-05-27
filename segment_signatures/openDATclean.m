function [xLength, yLength, cubeData_alt, cubeData_orig, linearData_alt, linearData_orig, pathName, fileName] = openDATclean(orderNum)
% Kimberly Chan
% Last edited 4/5/13
% Open DAT file
tic
cd 'C:\FAMILY\Kimberly' % Directory where files are located
[fileName, pathName] = uigetfile('*.dat'); % GUI to select desired file
fID = fopen([pathName, fileName]);

% Read file and isolate relevant data
fprintf('\rReading file and isolating relevant data.\r')
toc
tot_file = fread(fID,inf,'uint8=>char')';
xLength = str2double(tot_file((findstr('* Abscissa points :   ',tot_file) + 22):(findstr('* Ordinate points :   ',tot_file)-3)));
yLength = str2double(tot_file((findstr('* Ordinate points :   ',tot_file) + 22):(findstr('* BLANK LINE', tot_file)-3)));
data = str2num(tot_file((findstr('* DATA',tot_file) + 6):length(tot_file)));
eleData = data(:,9:18);

% Order: Fe, Cu, Zn, Ca, K, S, P, Cl, Si, Mn
cubeData_alt = zeros(yLength,xLength,10);
cubeData_orig = zeros(yLength,xLength,10);
linearData_alt = zeros(10,xLength*yLength);
linearData_orig = zeros(10,xLength*yLength);
fprintf('\rRearranging data into cube form or linear form.')
toc
for i = 1:10
    fprintf('Rearranging element %d of %d.',[i 10]), toc
    cubeData_orig(:,:,i) = rot90(reshape(eleData(:,i),xLength,yLength));
    linearData_orig(i,:) = reshape(cubeData_orig(:,:,i)',1,xLength*yLength);
end

% Stretching the data
maxVals = zeros(1,10);
j = 1:10;
maxVals(j) = max(max(cubeData_orig(:,:,j)));
minVals = min(maxVals);
for k = 1:10
    cubeData_alt(:,:,k) = cubeData_orig(:,:,k)*(minVals/maxVals(k));
end

% Moving average filter and arranging into linear form
H = fspecial('average',11);
for k = 1:10
    cubeData_alt(:,:,k) = imfilter(cubeData_alt(:,:,k),H);
    linearData_alt(k,:) = reshape(cubeData_alt(:,:,k)',1,xLength*yLength);
end

if not(isempty(orderNum))
    figure()
    imagesc(cubeData_alt(:,:,orderNum))
    set(gca,'XTick',[],'YTick',[])
    set(gcf,'Position',[182 208 293 449])
end
end