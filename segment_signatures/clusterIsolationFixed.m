function [allSig] = clusterIsolationFixed()
% Kimberly Chan
% Last edited 4/13/13
% Isolates clusters within synchrotron imaging data from a fixed image and
% produces a list of different clusters and their locations.
% 
% This is a minimalist version of clustIsolation.m meant to be used on
% fixed images of cyanobacteria.
%
% Initialization takes about 5 minutes on my laptop.

% Change me:
minMag = 3.25;  % Minimum magnitude of signature
distBtwnPix = 1.5;  % Maximum distance between contiguous pixels
correlationThresh = 0.975; % Correlation coefficient between pixels in a signature
minPix = 15; % Minimum pixel amount to be considered a signature

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain data from DAT files
% xLength - Number of pixels in the x direction
% yLength - Number of pixels in the y direction
% cubeData - 10 elemental data sets in cube form.
% linearData - 10 elemental data sets in a 10 by xLength*yLength matrix.
fprintf('\rOpening file\r')
tic
[xLength, yLength, cubeData_alt, cubeData_orig, linearData_alt, linearData_orig, pathName, fileName]= openDATclean([]);

% LINEAR TO XY - Change from linear to x,y coordinates
function [picCluster] = linear2xy(tempBlank,xLength,yLength)
picCluster = rot90(rot90(fliplr(rot90(reshape(tempBlank,xLength,yLength)))));
end

% ISOLATE DATA - Isolating the data from the background
function [linTrackingMatrix] = dataIsolation(xLength,yLength,linearData)
fprintf('\rIsolating the data from the background.\r')
toc
linTrackingMatrix = zeros(1,yLength*xLength);
for j = xLength*yLength:-1:1
    curSig = linearData(:,j);
    sigMag = sqrt(sum(curSig.^2));
    if sigMag < minMag
        linTrackingMatrix(j) = 1;
    end
end
end

linTrackingMatrix = dataIsolation(xLength,yLength,linearData_alt);

% SIGNATURE CLASSIFICATION - Bin data into different signature lists.
function [sigList] = sigClassify(xLength,yLength,linTrackingMatrix)
fprintf('\rBinning data into different signature lists.\r')
toc
track = 0;   % Starting tracking cluster number
sigArray = [];   % Initialize sig indices array

for n = 1:xLength*yLength
    fprintf('\rDetermining signature of pixel %d of %d.\r',[n xLength*yLength]), toc
    if linTrackingMatrix(n) == 0     % Continue if not background
        curSig = linearData_alt(:,n);
        sigArray = [sigArray, n];     % Add the first pixel to array
        track = track + 1;       % Signature number
        for next = (n+1):length(linTrackingMatrix) % Go through remaining pixels
            laterSig = linearData_alt(:,next);
            p = dot(curSig, laterSig)/(sqrt(sum(curSig.^2))*sqrt(sum(laterSig.^2)));
            if linTrackingMatrix(next) == 0 && p >correlationThresh  % Bin pixel into list if correlated
                sigArray = [sigArray next];
                linTrackingMatrix(next) = 1; % Mark pixel as "used"
            end
        end
        sigList{track} = sigArray; % Assign an array to a list
        sigArray = []; % Reset the signature array for next list.
    end
end
% Getting rid of signatures with only small amount of pixels in them.
% These are messing up the color schematics.
fprintf('\rGetting rid of signatures with only a small amount of pixels in them.\r')
toc
sigList3 = cell(1,1);
track3 = 0;
for one = 1:length(sigList)
    if length(sigList{one}) > minPix
        track3 = track3 + 1;
        sigList3{track3} = sigList{one};
    end
end
sigList = sigList3;
end

sigList = sigClassify(xLength,yLength,linTrackingMatrix);

% ASSIGNING CLUSTERS TO A SIGNATURE - Assigning clusters to different signatures
function [tempBlank, allSig] = assignClust(xLength,yLength,sigList)
fprintf('\rAssigning clusters to different signatures.\r')
toc
tempBlank = zeros(1,xLength*yLength);
tempBlank2 = zeros(1,xLength*yLength);
allSig = cell(1,1);
for sigNumLoop = 1:length(sigList)
    tempBlank(sigList{sigNumLoop}) = sigNumLoop;
    tempBlank2(sigList{sigNumLoop}) = 1;
    allSig{sigNumLoop} = tempBlank2;
    tempBlank2 = zeros(1,xLength*yLength);
end
end

% Display all signatures.
fprintf('\rDisplaying all signatures.\r')
toc
[tempBlank, allSig] = assignClust(xLength,yLength,sigList);
sigIm = linear2xy(tempBlank,xLength,yLength);
figHandle1 = figure();
subplot(6,19,[1:6 20:25 39:44 58:62 77:82 96:101]);
imagesc(linear2xy(tempBlank,xLength,yLength))
set(gca,'XTick',[],'YTick',[])
title('Cluster Signatures','fontsize',14,'fontweight','b')
set(gcf,'Position',[182 208 1343 449])

% CLUSTER SEPARATION - Separate clusters based on distance between
% contiguous pixels.
function [aClustList2] = clustSeparation(sigList,sigNum,sigIndices)
dist = [];
disTrack = 0;
distTrackingArray = zeros(1,length(sigList{sigNum}));

% Make an array of all the distances from one pixel to the rest.
fprintf('\rMaking an array of all the distances from one pixel to the rest.\r')
toc
for f = 1:length(sigList{sigNum})
    if distTrackingArray(f) == 0
        disTrack = disTrack + 1;    % Track which pixel the program is at.
        for g = (f+1):length(sigList{sigNum})
            dist = [dist sqrt((sigIndices(2,f)- sigIndices(2,g))^2 + (sigIndices(3,f)- sigIndices(3,g))^2)];
            distTrackingArray(f) = 1;
        end
        clusterListDist{disTrack} = dist;  % Assign distance array to a cell
        dist = [];
    end
end

% Isolate contiguous pixels and assign it to a cell
fprintf('\rIsolating contiguous pixels and assign it to a cell\r')
toc
for uqe = 1:length(clusterListDist)
    decIndices{uqe} = uqe + find(clusterListDist{uqe} < distBtwnPix);
end

% Pools all the contiguous pixels together and assign each cluster to a
% cell
fprintf('\rPooling all the contiguous pixels together and assigning each cluster to a cell.\r')
toc
nextDisTrack = 0;
nextDistTrackingArray = zeros(1,length(decIndices));
clustF = [];
for check = 1:length(decIndices)
    if nextDistTrackingArray(check) == 0
        nextDistTrackingArray(check) = 1;
        nextDisTrack = nextDisTrack + 1;
        % Put pixel & contiguous pixels into the array
        clustF = [clustF, check, decIndices{check}];
        for nextCheck = (check+1):(length(decIndices)-1) 
            % Check if other pixels are adjacent to original pixel
            if length(find(ismember(clustF,[nextCheck decIndices{nextCheck}]))) >= 1
                nextDistTrackingArray(nextCheck) = 1;
                % Add that pixel and its contiguous pixels to the cluster
                clustF = [clustF, nextCheck, decIndices{nextCheck}];
            end
        end
        aClustList{nextDisTrack} = unique(clustF);
        clustF = [];
    end
end

% Patching bug.  Code in previous section divides some clusters into 
% separate clusters for some reason.  The following code compares all those
% separate clusters and see if they have any pixels in common and then
% groups those again.
%
% Will fix this more legitimately later.
fprintf('\rPutting separated clusters back together again.\r')
toc
checkArray3 = zeros(1,length(aClustList));
clustF2 = [];
trackNum = 0;
for bettaCheck = 1:length(aClustList)
    if checkArray3(bettaCheck) == 0;
        checkArray3(bettaCheck) = 1;
        trackNum = trackNum + 1;
        for nextBettaCheck = (bettaCheck + 1):length(aClustList)
            if ~isempty(find(ismember(aClustList{bettaCheck}, aClustList{nextBettaCheck}), 1))
                checkArray3(nextBettaCheck) = 1;
                clustF2 = [clustF2, aClustList{bettaCheck}, aClustList{nextBettaCheck}];
            else
                clustF2 = [clustF2, aClustList{bettaCheck}];
            end
        end
        aClustList2{trackNum} = unique(clustF2);
        clustF2 = [];
    end
end
end

% DISPLAY CLUSTERS - Display all clusters in different colors corresponding
% to cluster number.
function [finPicClust] = dispClust(yLength,xLength,aClustList2,sigIndices)
finPicClust = zeros(yLength,xLength);
for binAll = 1:length(aClustList2)
    for w = aClustList2{binAll}
        % Assign parts of the image to different cluster numbers.
        finPicClust(sigIndices(2,w), sigIndices(3,w)) = binAll;
    end
end
end

% Code for starting interactive mode with graph.  Select signature to bring
% up average spectra and graph of all clusters.

% DISPLAY SPECTRA AND CLUSTERS OF SIGNATURE - Displays the spectra and 
% individual clusters of selected signature.
function [] = spectraClust(cursorInfo, cursTrack, xLength, yLength, sigIm, sigList,linearData_orig, allSig)
sigPos = cursorInfo(cursTrack).Position;  % Column, row
eleNum_all = size(linearData_orig);
selectedSig = sigIm(sigPos(2),sigPos(1));

% Clusters of selected signature
[row, col] = find(linear2xy(allSig{selectedSig},xLength,yLength));
sigIndices = [sigList{selectedSig}; row'; col'];
curClustList = clustSeparation(sigList,selectedSig,sigIndices);
finPicClustDisp = dispClust(yLength,xLength,curClustList,sigIndices); % Clusters within sig

% Plot clusters of selected signature
subplot(6,19,[7:12 26:31 45:50 64:69 83:88 102:107]);
imagesc(finPicClustDisp)
set(gca,'XTick',[],'YTick',[])
title(sprintf('Clusters of Signature %s',num2str(selectedSig)),'fontsize',14,'fontweight','b')

% Average spectra of selected signature
spectra = linearData_orig(:,find(allSig{selectedSig}));
spectraSz = size(spectra);
%spectraSzTransposed = size(spectra');
subplot(6,19,[14:19 33:38 52:57 71:76 90:95 109:114]);
lower_range25 = median(spectra') - prctile(spectra', 25);
upper_range75 = prctile(spectra', 75) - median(spectra');

% Write spectra into excel file.
% [fileNamexl, pathNamexl] = uigetfile('*.xlsx','Select the excel spreadsheet you would like to write into.');
% excelFile = [pathNamexl,fileNamexl];
% xlswrite(excelFile,spectra')

% Display spectra.
errorbar(1:eleNum_all(1), median(spectra'),lower_range25, upper_range75, '--k','LineWidth',2)
title('Average Spectra','fontsize',14,'fontweight','b')
ylabel('Pixel Intensity','fontsize',12,'fontweight','b')
set(gca, 'XTick',1:spectraSz(1), 'XTickLabel',{'Fe' 'Cu' 'Zn' 'Ca' 'K' 'S' 'P' 'Cl' 'Si' 'Mn'})
xlabel('Element','fontsize',12,'fontweight','b')
hold on

% Average spectra of selected cluster
selectedClust = finPicClustDisp(sigPos(2),sigPos(1)); % Selected cluster
[clustRow, clustCol] = find(finPicClustDisp == selectedClust);
tempArray = [];
spectraClust = [];
for bin1 = 1:length(clustRow)
    for eleNum = 1:10;
        tempArray = [tempArray; cubeData_orig(clustRow(bin1),clustCol(bin1),eleNum)];
    end
    spectraClust = [spectraClust tempArray];
    tempArray = [];
end
spectraSzClust = size(spectraClust);
lower_range25_clust = median(spectraClust') - prctile(spectraClust', 25);
upper_range75_clust = prctile(spectraClust', 75) - median(spectraClust');
if spectraSzClust > 1
    errorbar(1:eleNum_all(1), median(spectraClust'), lower_range25_clust, upper_range75_clust, 'r','LineWidth',2)
else
    plot(spectraClust,'r','LineWidth',2)
end
title(sprintf('Average Spectra of Cluster %d of Signature %d',[(selectedClust) (selectedSig)]),'fontsize',14,'fontweight','b')
ylabel('Pixel Intensity','fontsize',12,'fontweight','b')
set(gca, 'XTick',1:spectraSzClust(1), 'XTickLabel',{'Fe' 'Cu' 'Zn' 'Ca' 'K' 'S' 'P' 'Cl' 'Si' 'Mn'})
ylim([0 max([prctile(spectraClust', 75), prctile(spectra', 75)])+1])
xlabel('Element','fontsize',12,'fontweight','b')
legend('Signature Average','Cluster Average','Location','NorthEast');
hold off
end

% Initialize interactive mode.
fprintf('\rInitializing interactive mode.\r')
toc
h = datacursormode(figHandle1);
datacursormode on
cursorInfo = [];
cursTrack = 0;
fprintf('\rClick on the image titled "Cluster Signatures" to select a signature\ror a cluster then press Return.  Dark blue is the background.\r')
while(1)
    pause   % Allow user time to select a point.  Press enter to execute
    cursorInfo = [cursorInfo, getCursorInfo(h)];
    cursTrack = cursTrack + 1;
    % Reassigning value to allSig because it keeps disappearing.
    [tempBlank, allSig] = assignClust(xLength,yLength,sigList);
    spectraClust(cursorInfo, cursTrack, xLength, yLength, sigIm, sigList, linearData_orig, allSig)
end
end