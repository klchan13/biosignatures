function varargout = clustIsolation(varargin)
% CLUSTISOLATION M-file for clustIsolation.fig
%      CLUSTISOLATION, by itself, creates a new CLUSTISOLATION or raises the existing
%      singleton*.
%
%      H = CLUSTISOLATION returns the handle to a new CLUSTISOLATION or the handle to
%      the existing singleton*.
%
%      CLUSTISOLATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTISOLATION.M with the given input arguments.
%
%      CLUSTISOLATION('Property','Value',...) creates a new CLUSTISOLATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before clustIsolation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to clustIsolation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help clustIsolation

% Last Modified by GUIDE v2.5 13-Apr-2013 14:20:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @clustIsolation_OpeningFcn, ...
                   'gui_OutputFcn',  @clustIsolation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%
% Kimberly Chan
% Last edited: 4/13/13
% Isolates clusters within synchrotron imaging data and produces a list of
% different clusters and their locations.
%
% This is a more extensive version of clusterIsolationFixed.m meant to
% analyze the live series of cyanobacteria under fluctuating humidity.
%
% Instructions:
% 1. Click on "Select
%
% Warning: Before clicking on button "Selected" on the top right corner of
% Make sure 
%
% To go to a certain line in the code, press ctrl + g and type in the line.
%
% Button Functions and Directory:
% pushbutton1 - Select signature - line
% pushbutton2 - Select cluster - line 
% pushbutton3 - Display spectra for singular cluster. - 
% pushbutton4 - Display spectra for all clusters of the signature. - 
% radiobutton1 - Push button 1 indicator - 
% radiobutton2 - Push button 2 indicator - line
% 
% Inputs:
% correlationThresh - Correlation threshold.
% minMag - Minimum magnitude to isolate data from background


% --- Executes just before clustIsolation is made visible.
function clustIsolation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to clustIsolation (see VARARGIN)

% Choose default command line output for clustIsolation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes clustIsolation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


minMag = 300;  % Minimum magnitude of signature
distBtwnPix = 1.5;  % Maximum distance between contiguous pixels
correlationThresh = 0.975; % Correlation coefficient between pixels in a signature

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain data from DAT files
% xLength - Number of pixels in the x direction
% yLength - Number of pixels in the y direction
% cubeData - 10 elemental data sets in cube form.
% linearData - 10 elemental data sets in a 10 by xLength*yLength matrix.
[xLength, yLength, ~, linearData, pathName, fileName] = openDAT([]);

% LINEAR TO XY - Change from linear to x,y coordinates
function [picCluster] = linear2xy(tempBlank,xLength,yLength)
picCluster = rot90(rot90(fliplr(rot90(reshape(tempBlank,xLength,yLength)))));
end

% ISOLATE DATA - Isolating the data from the background
function [linTrackingMatrix] = dataIsolation(xLength,yLength,linearData)
linTrackingMatrix = zeros(1,yLength*xLength);
for j = xLength*yLength:-1:1
    curSig = linearData(:,j);
    sigMag = sqrt(sum(curSig.^2));
    if sigMag < minMag
        linTrackingMatrix(j) = 1;
    end
end
end

linTrackingMatrix = dataIsolation(xLength,yLength,linearData);
linearData = log(linearData+1);

% SIGNATURE CLASSIFICATION - Bin data into different signature lists.
function [sigList] = sigClassify(xLength,yLength,linTrackingMatrix)
track = 0;   % Starting tracking cluster number
sigArray = [];   % Initialize sig indices array

for n = 1:xLength*yLength
    if linTrackingMatrix(n) == 0     % Continue if not background
        curSig = linearData(:,n);
        sigArray = [sigArray, n];     % Add the first pixel to array
        track = track + 1;       % Signature number
        for next = (n+1):length(linTrackingMatrix) % Go through remaining pixels
            laterSig = linearData(:,next);
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
sigList3 = cell(1,1);
track3 = 0;
for one = 1:length(sigList)
    if length(sigList{one}) > 20
        track3 = track3 + 1;
        sigList3{track3} = sigList{one};
    end
end
sigList = sigList3;
end

sigList = sigClassify(xLength,yLength,linTrackingMatrix);

% ASSIGNING CLUSTERS TO A SIGNATURE - Assigning clusters to different signatures
function [tempBlank, allSig] = assignClust(xLength,yLength,sigList)
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
[tempBlank, allSig] = assignClust(xLength,yLength,sigList);
sigIm = linear2xy(tempBlank,xLength,yLength);
figHandle1 = imagesc(linear2xy(tempBlank,xLength,yLength),'Parent',handles.axes2);
axes(handles.axes2)
axis off
set(gca,'XTick',[],'YTick',[])
title('Cluster Signatures','fontsize',14,'fontweight','b')
axes(handles.axes1)
axis off

% Storing data in invisible GUI boxes
set(handles.edit1,'String',xLength)
set(handles.edit3,'String',yLength)
set(handles.edit4,'String',[pathName, fileName])
set(handles.uitable1,'data',linearData)
set(handles.uitable8,'data',linTrackingMatrix)
end
end


% --- Outputs from this function are returned to the command line.
function varargout = clustIsolation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Variables:
minMag = 300;  % Minimum magnitude of signature
distBtwnPix = 1.5;  % Maximum distance between contiguous pixels
correlationThresh = 0.975; % Correlation coefficient between pixels in a signature

% Getting data stored in invisible GUI boxes.
set(handles.radiobutton2,'Value',0)
set(handles.radiobutton1,'Value',1)
linearData = get(handles.uitable1,'data');
linTrackingMatrix = get(handles.uitable8,'data');
xLength = str2num(get(handles.edit1,'String'));
yLength = str2num(get(handles.edit3,'String'));


% LINEAR TO XY - Change from linear to x,y coordinates
function [picCluster] = linear2xy(tempBlank,xLength,yLength)
picCluster = rot90(rot90(fliplr(rot90(reshape(tempBlank,xLength,yLength)))));
end

% SIGNATURE CLASSIFICATION - Bin data into different signature lists.
function [sigList] = sigClassify(xLength,yLength,linTrackingMatrix)
track = 0;   % Starting tracking cluster number
sigArray = [];   % Initialize sig indices array

for n = 1:xLength*yLength
    if linTrackingMatrix(n) == 0     % Continue if not background
        curSig = linearData(:,n);
        sigArray = [sigArray, n];     % Add the first pixel to array
        track = track + 1;       % Signature number
        for next = (n+1):length(linTrackingMatrix) % Go through remaining pixels
            laterSig = linearData(:,next);
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
sigList3 = cell(1,1);
track3 = 0;
for one = 1:length(sigList)
    if length(sigList{one}) > 20
        track3 = track3 + 1;
        sigList3{track3} = sigList{one};
    end
end
sigList = sigList3;
end

sigList = sigClassify(xLength,yLength,linTrackingMatrix);

% ASSIGNING CLUSTERS TO A SIGNATURE - Assigning clusters to different signatures
function [tempBlank, allSig] = assignClust(xLength,yLength,sigList)
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

[tempBlank, allSig] = assignClust(xLength,yLength,sigList);
sigIm = linear2xy(tempBlank,xLength,yLength);

% CLUSTER SEPARATION - Separate clusters based on distance between
% contiguous pixels.
function [aClustList2] = clustSeparation(sigList,sigNum,sigIndices)
dist = [];
disTrack = 0;
distTrackingArray = zeros(1,length(sigList{sigNum}));

% Make an array of all the distances from one pixel to the rest.
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
for uqe = 1:length(clusterListDist)
    decIndices{uqe} = uqe + find(clusterListDist{uqe} < distBtwnPix);
end

% Pools all the contiguous pixels together and assign each cluster to a
% cell
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
aClustList3 = cell(1,1);
track2 = 0;
for one = 1:length(aClustList2)
    if length(aClustList2{one}) > 1
        track2 = track2 + 1;
        aClustList3{track2} = aClustList2{one};
    end
end
aClustList2 = aClustList3;
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


function [] = spectraClustDisp(cursorInfo, cursTrack, xLength, yLength, sigIm, sigList,allSig,linearData)
    axes(handles.axes3)
    cla %clear axes
    sigPos = cursorInfo(cursTrack).Position;  % Column, row
    selectedSig = sigIm(sigPos(2),sigPos(1));
    set(handles.edit8,'String',selectedSig)
    
    % Clusters of selected signature
    [row, col] = find(linear2xy(allSig{selectedSig},xLength,yLength));
    sigIndices = [sigList{selectedSig}; row'; col'];
    curClustList = clustSeparation(sigList,selectedSig,sigIndices);
    finPicClustDisp = dispClust(yLength,xLength,curClustList,sigIndices); % Clusters within sig
    imagesc(finPicClustDisp,'Parent',handles.axes1)
    axes(handles.axes1)
    set(gca,'XTick',[],'YTick',[])
    title(sprintf('Clusters of Signature %s',num2str(selectedSig)),'fontsize',14,'fontweight','b')
    
    % Average spectra of selected signature
    spectra = linearData(:,find(allSig{selectedSig}));
    spectraSz = size(spectra);
    errorbar(mean(spectra'),std(spectra'),'Parent',handles.axes3,'LineWidth',2)
    axes(handles.axes3)
    ylabel('log(intensity)','fontsize',12,'fontweight','b')
    set(gca, 'XTick',1:spectraSz(1), 'XTickLabel',{'Fe' 'Cu' 'Zn' 'Ca' 'K' 'S' 'P' 'Cl' 'Si' 'Mn'})
    xlabel('Element','fontsize',12,'fontweight','b')
end

h = datacursormode(handles.figure1);
datacursormode on
cursorInfo = [];
cursTrack = 0;
fprintf('\r\rClick on the image titled "Cluster Signatures" to select a signature then press Return.\rDark blue is the background.')
while(1)
    pause   % Allow user time to select a point.  Press enter to execute
    cursorInfo = [cursorInfo, getCursorInfo(h)];
    cursTrack = cursTrack + 1;
    % Reassigning value to allSig because it keeps disappearing.
    [tempBlank, allSig] = assignClust(xLength,yLength,sigList);
    spectraClustDisp(cursorInfo, cursTrack, xLength, yLength, sigIm, sigList,allSig, linearData)
end
end



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% SELECT CLUSTER

% Variables:
minMag = 300;  % Minimum magnitude of signature
distBtwnPix = 1.5;  % Maximum distance between contiguous pixels
correlationThresh = 0.975; % Correlation coefficient between pixels in a signature

% Retrieving hidden data.
set(handles.radiobutton2,'Value',1)
set(handles.radiobutton1,'Value',0)
linearData = get(handles.uitable1,'data');
linTrackingMatrix = get(handles.uitable8,'data');
xLength = str2num(get(handles.edit1,'String'));
yLength = str2num(get(handles.edit3,'String'));

linearData2 = linearData;

cubeData = zeros(yLength,xLength,10);
for i = 1:10
    cubeData(:,:,i) =  rot90(rot90(fliplr(rot90(reshape(linearData2(i,:),xLength,yLength)))));
end

function [picCluster] = linear2xy(tempBlank,xLength,yLength)
picCluster = rot90(rot90(fliplr(rot90(reshape(tempBlank,xLength,yLength)))));
end

% SIGNATURE CLASSIFICATION - Bin data into different signature lists.
function [sigList] = sigClassify(xLength,yLength,linTrackingMatrix)
track = 0;   % Starting tracking cluster number
sigArray = [];   % Initialize sig indices array

for n = 1:xLength*yLength
    if linTrackingMatrix(n) == 0     % Continue if not background
        curSig = linearData(:,n);
        sigArray = [sigArray, n];     % Add the first pixel to array
        track = track + 1;       % Signature number
        for next = (n+1):length(linTrackingMatrix) % Go through remaining pixels
            laterSig = linearData(:,next);
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
sigList3 = cell(1,1);
track3 = 0;
for one = 1:length(sigList)
    if length(sigList{one}) > 20
        track3 = track3 + 1;
        sigList3{track3} = sigList{one};
    end
end
sigList = sigList3;
end

sigList = sigClassify(xLength,yLength,linTrackingMatrix);

% ASSIGNING CLUSTERS TO A SIGNATURE - Assigning clusters to different signatures
function [tempBlank, allSig] = assignClust(xLength,yLength,sigList)
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
[tempBlank, allSig] = assignClust(xLength,yLength,sigList);
sigIm = linear2xy(tempBlank,xLength,yLength);
imagesc(linear2xy(tempBlank,xLength,yLength),'Parent',handles.axes2);
axes(handles.axes2)
axis off
set(gca,'XTick',[],'YTick',[])
title('Cluster Signatures','fontsize',14,'fontweight','b')
% set(gcf,'Position',[182 208 1343 449])

% CLUSTER SEPARATION - Separate clusters based on distance between
% contiguous pixels.
function [aClustList2] = clustSeparation(sigList,sigNum,sigIndices)
dist = [];
disTrack = 0;
distTrackingArray = zeros(1,length(sigList{sigNum}));

% Make an array of all the distances from one pixel to the rest.
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
for uqe = 1:length(clusterListDist)
    decIndices{uqe} = uqe + find(clusterListDist{uqe} < distBtwnPix);
end

% Pools all the contiguous pixels together and assign each cluster to a
% cell
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
aClustList3 = cell(1,1);
track2 = 0;
for one = 1:length(aClustList2)
    if length(aClustList2{one}) > 1
        track2 = track2 + 1;
        aClustList3{track2} = aClustList2{one};
    end
end
aClustList2 = aClustList3;
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

function [oneFinPicClust] = spectOneClust(clustNum,aClustList2,xLength,yLength,sigIndices)
oneClust = aClustList2{clustNum};
oneFinPicClust = zeros(yLength,xLength);
for w2 = oneClust
    oneFinPicClust(sigIndices(2,w2), sigIndices(3,w2)) = 1;
end
end

function [] = spectraClustDisp(cursorInfo, cursTrack, xLength, yLength, sigIm, sigList,aClustList2,linearData,cubeData)
axes(handles.axes3)
cla %clear axes
sigPos = cursorInfo(cursTrack).Position;  % Column, row
selectedSig = sigIm(sigPos(2),sigPos(1));

%Find out which cluster is used.
[row, col] = find(linear2xy(allSig{selectedSig},xLength,yLength));
sigIndices = [sigList{selectedSig}; row'; col'];
curClustList = clustSeparation(sigList,selectedSig,sigIndices);
finPicClustDisp = dispClust(yLength,xLength,curClustList,sigIndices); % Clusters within sig
selectedClust = finPicClustDisp(sigPos(2),sigPos(1));
set(handles.edit9,'String',selectedClust)

% Average spectra of selected signature in dashes
spectraSig = linearData(:,find(allSig{selectedSig}));
spectraSzSig = size(spectraSig);
hold on

% Average spectra of selected cluster
[clustRow, clustCol] = find(finPicClustDisp == selectedClust);
tempArray = [];
spectraClust = [];
for bin1 = 1:length(clustRow)
    for eleNum = 1:10;
        tempArray = [tempArray; cubeData(clustRow(bin1),clustCol(bin1),eleNum)];
    end
    spectraClust = [spectraClust tempArray];
    tempArray = [];
end

spectraSzClust = size(spectraClust);
hold on
errorbar(mean(spectraClust'),std(spectraClust'),'Parent',handles.axes3, 'LineWidth',2)
if get(handles.checkbox1,'Value') == 1
    errorbar(mean(spectraSig'),std(spectraSig'),'Parent',handles.axes3,'--k') % Avg of sig
    legend('Cluster Average','Signature Average','Location','NorthEast');
else
    legend('Cluster Average','Location','NorthEast');
end
ylabel('log(intensity)','fontsize',12,'fontweight','b')
set(gca, 'XTick',1:spectraSzClust(1), 'XTickLabel',{'Fe' 'Cu' 'Zn' 'Ca' 'K' 'S' 'P' 'Cl' 'Si' 'Mn'})
xlabel('Element','fontsize',12,'fontweight','b')
end

h = datacursormode(handles.figure1);
datacursormode on
cursorInfo = [];
cursTrack = 0;
fprintf('\r\rClick on the image titled "Clusters of Signature (number)" to select a cluster then press Return.\rDark blue is the background.')
while(1)
    pause   % Allow user time to select a point.  Press enter to execute
    cursorInfo = [cursorInfo, getCursorInfo(h)];
    cursTrack = cursTrack + 1;
    % Reassigning value to allSig because it keeps disappearing.
    [~, allSig] = assignClust(xLength,yLength,sigList);
    spectraClustDisp(cursorInfo, cursTrack, xLength, yLength, sigIm, sigList,allSig, linearData,cubeData)
end
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
end

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
minMag = 300;  % Minimum magnitude of signature
distBtwnPix = 1.5;  % Maximum distance between contiguous pixels
correlationThresh = 0.975; % Correlation coefficient between pixels in a signature
linearData = get(handles.uitable1,'data');
linTrackingMatrix = get(handles.uitable8,'data');
xLength = str2num(get(handles.edit1,'String'));
yLength = str2num(get(handles.edit3,'String'));

linearData2 = linearData;

cubeData = zeros(yLength,xLength,10);
for i = 1:10
    cubeData(:,:,i) =  rot90(rot90(fliplr(rot90(reshape(linearData2(i,:),xLength,yLength)))));
end

% LINEAR TO XY - Change from linear to x,y coordinates
function [picCluster] = linear2xy(tempBlank,xLength,yLength)
picCluster = rot90(rot90(fliplr(rot90(reshape(tempBlank,xLength,yLength)))));
end

% SIGNATURE CLASSIFICATION - Bin data into different signature lists.
function [sigList] = sigClassify(xLength,yLength,linTrackingMatrix)
track = 0;   % Starting tracking cluster number
sigArray = [];   % Initialize sig indices array

for n = 1:xLength*yLength
    if linTrackingMatrix(n) == 0     % Continue if not background
        curSig = linearData(:,n);
        sigArray = [sigArray, n];     % Add the first pixel to array
        track = track + 1;       % Signature number
        for next = (n+1):length(linTrackingMatrix) % Go through remaining pixels
            laterSig = linearData(:,next);
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
sigList3 = cell(1,1);
track3 = 0;
for one = 1:length(sigList)
    if length(sigList{one}) > 20
        track3 = track3 + 1;
        sigList3{track3} = sigList{one};
    end
end
sigList = sigList3;
end

sigList = sigClassify(xLength,yLength,linTrackingMatrix);

% ASSIGNING CLUSTERS TO A SIGNATURE - Assigning clusters to different signatures
function [tempBlank, allSig] = assignClust(xLength,yLength,sigList)
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
[tempBlank, allSig] = assignClust(xLength,yLength,sigList);
sigIm = linear2xy(tempBlank,xLength,yLength);
figHandle1 = imagesc(linear2xy(tempBlank,xLength,yLength),'Parent',handles.axes2);
axes(handles.axes2)
axis off
set(gca,'XTick',[],'YTick',[])
title('Cluster Signatures','fontsize',14,'fontweight','b')
% set(gcf,'Position',[182 208 1343 449])

% CLUSTER SEPARATION - Separate clusters based on distance between
% contiguous pixels.
function [aClustList2] = clustSeparation(sigList,sigNum,sigIndices)
dist = [];
disTrack = 0;
distTrackingArray = zeros(1,length(sigList{sigNum}));

% Make an array of all the distances from one pixel to the rest.
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
for uqe = 1:length(clusterListDist)
    decIndices{uqe} = uqe + find(clusterListDist{uqe} < distBtwnPix);
end

% Pools all the contiguous pixels together and assign each cluster to a
% cell
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
aClustList3 = cell(1,1);
track2 = 0;
for one = 1:length(aClustList2)
    if length(aClustList2{one}) > 1
        track2 = track2 + 1;
        aClustList3{track2} = aClustList2{one};
    end
end
aClustList2 = aClustList3;
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

axes(handles.axes3)
cla %clear axes
selectedSig = str2num(get(handles.edit8,'String'));

%Find out which cluster is used.
[row, col] = find(linear2xy(allSig{selectedSig},xLength,yLength));
sigIndices = [sigList{selectedSig}; row'; col'];
curClustList = clustSeparation(sigList,selectedSig,sigIndices);
finPicClustDisp = dispClust(yLength,xLength,curClustList,sigIndices); % Clusters within sig
selectedClust = str2num(get(handles.edit9,'String'));

% Average spectra of selected signature in dashes
spectraSig = linearData(:,find(allSig{selectedSig}));
spectraSzSig = size(spectraSig);
hold on

% Average spectra of selected cluster
[clustRow, clustCol] = find(finPicClustDisp == selectedClust);
tempArray = [];
spectraClust = [];
for bin1 = 1:length(clustRow)
    for eleNum = 1:10;
        tempArray = [tempArray; cubeData(clustRow(bin1),clustCol(bin1),eleNum)];
    end
    spectraClust = [spectraClust tempArray];
    tempArray = [];
end

spectraSzClust = size(spectraClust);
errorbar(mean(spectraClust'),std(spectraClust'),'Parent',handles.axes3, 'LineWidth',2)
if get(handles.checkbox1,'Value') == 1
    errorbar(mean(spectraSig'),std(spectraSig'),'Parent',handles.axes3,'--k') % Avg of sig
    legend('Cluster Average','Signature Average','Location','NorthEast');
else
    legend('Cluster Average','Location','NorthEast');
end
ylabel('log(intensity)','fontsize',12,'fontweight','b')
set(gca, 'XTick',1:spectraSzClust(1), 'XTickLabel',{'Fe' 'Cu' 'Zn' 'Ca' 'K' 'S' 'P' 'Cl' 'Si' 'Mn'})
xlabel('Element','fontsize',12,'fontweight','b')
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
minMag = 300;  % Minimum magnitude of signature
distBtwnPix = 1.5;  % Maximum distance between contiguous pixels
correlationThresh = 0.975; % Correlation coefficient between pixels in a signature
linearData = get(handles.uitable1,'data');
linTrackingMatrix = get(handles.uitable8,'data');
xLength = str2num(get(handles.edit1,'String'));
yLength = str2num(get(handles.edit3,'String'));

linearData2 = linearData;

cubeData = zeros(yLength,xLength,10);
for i = 1:10
    cubeData(:,:,i) =  rot90(rot90(fliplr(rot90(reshape(linearData2(i,:),xLength,yLength)))));
end

function [picCluster] = linear2xy(tempBlank,xLength,yLength)
picCluster = rot90(rot90(fliplr(rot90(reshape(tempBlank,xLength,yLength)))));
end

% SIGNATURE CLASSIFICATION - Bin data into different signature lists.
function [sigList] = sigClassify(xLength,yLength,linTrackingMatrix)
track = 0;   % Starting tracking cluster number
sigArray = [];   % Initialize sig indices array

for n = 1:xLength*yLength
    if linTrackingMatrix(n) == 0     % Continue if not background
        curSig = linearData(:,n);
        sigArray = [sigArray, n];     % Add the first pixel to array
        track = track + 1;       % Signature number
        for next = (n+1):length(linTrackingMatrix) % Go through remaining pixels
            laterSig = linearData(:,next);
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
sigList3 = cell(1,1);
track3 = 0;
for one = 1:length(sigList)
    if length(sigList{one}) > 20
        track3 = track3 + 1;
        sigList3{track3} = sigList{one};
    end
end
sigList = sigList3;
end

sigList = sigClassify(xLength,yLength,linTrackingMatrix);

% ASSIGNING CLUSTERS TO A SIGNATURE - Assigning clusters to different signatures
function [tempBlank, allSig] = assignClust(xLength,yLength,sigList)
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
[tempBlank, allSig] = assignClust(xLength,yLength,sigList);
sigIm = linear2xy(tempBlank,xLength,yLength);
figHandle1 = imagesc(linear2xy(tempBlank,xLength,yLength),'Parent',handles.axes2);
axes(handles.axes2)
axis off
set(gca,'XTick',[],'YTick',[])
title('Cluster Signatures','fontsize',14,'fontweight','b')
% set(gcf,'Position',[182 208 1343 449])

% CLUSTER SEPARATION - Separate clusters based on distance between
% contiguous pixels.
function [aClustList2] = clustSeparation(sigList,sigNum,sigIndices)
dist = [];
disTrack = 0;
distTrackingArray = zeros(1,length(sigList{sigNum}));

% Make an array of all the distances from one pixel to the rest.
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
for uqe = 1:length(clusterListDist)
    decIndices{uqe} = uqe + find(clusterListDist{uqe} < distBtwnPix);
end

% Pools all the contiguous pixels together and assign each cluster to a
% cell
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
aClustList3 = cell(1,1);
track2 = 0;
for one = 1:length(aClustList2)
    if length(aClustList2{one}) > 1
        track2 = track2 + 1;
        aClustList3{track2} = aClustList2{one};
    end
end
aClustList2 = aClustList3;
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

axes(handles.axes3)
cla %clear axes
selectedSig = str2num(get(handles.edit8,'String'));

%Find out which cluster is used.
[row, col] = find(linear2xy(allSig{selectedSig},xLength,yLength));
sigIndices = [sigList{selectedSig}; row'; col'];
curClustList = clustSeparation(sigList,selectedSig,sigIndices);
finPicClustDisp = dispClust(yLength,xLength,curClustList,sigIndices); % Clusters within sig

% Average spectra of selected signature in dashes
spectraSig = linearData(:,find(allSig{selectedSig}));
spectraSzClust = size(spectraSig);
hold on
if get(handles.checkbox1,'Value') == 1
    errorbar(mean(spectraSig'),std(spectraSig'),'Parent',handles.axes3,'--k') % Avg of sig
    legend('Signature Average','Location','NorthEast');
end
ylabel('log(intensity)','fontsize',12,'fontweight','b')
set(gca, 'XTick',1:spectraSzClust(1), 'XTickLabel',{'Fe' 'Cu' 'Zn' 'Ca' 'K' 'S' 'P' 'Cl' 'Si' 'Mn'})
xlabel('Element','fontsize',12,'fontweight','b')

% Average spectra of selected cluster
cmap = hsv(length(curClustList));
for allPlot = 1:length(curClustList)
    % Average spectra of selected cluster
    [clustRow, clustCol] = find(finPicClustDisp == allPlot);
    tempArray = [];
    spectraClust = [];
    for bin1 = 1:length(clustRow)
        for eleNum = 1:10;
            tempArray = [tempArray; cubeData(clustRow(bin1),clustCol(bin1),eleNum)];
        end
    spectraClust = [spectraClust tempArray];
    tempArray = [];
    end
    errorbar(mean(spectraClust'),std(spectraClust'),'Parent',handles.axes3, 'LineWidth',2,'Color',cmap(allPlot,:))
    clustLegend{allPlot} = sprintf('Cluster %s',num2str(allPlot));
end
end

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
end

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
end

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
end
