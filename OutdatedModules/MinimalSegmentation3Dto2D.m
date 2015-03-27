function handles = MinimalProjectSegmentation3Dto2D(handles)


% Help for Project3DTo2D
% Category: Image Processing
% This module projects 3D Segmenation to 2D segmentation
% 
% TO DO: NICE VISUALIZATIO,
% add additional projection methods (from script)
% test speed with various different input objects
% make direction optional?
% may call functions inside main function to prevent dublication of memory


drawnow


handles.Settings.VariableValues



[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);
%textVAR01 = What should the projection show?
%choiceVAR01 = OccupiedPlanes
%choiceVAR01 = Unsupported
inProjectionType = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = Which Segementation should be projected?
%infotypeVAR02 = objectgroup
inObjectOfInterest = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = How do you want to call the projected Image?
%defaultVAR03 = Nuclei
%infotypeVAR03 = imagegroup indep
inOutputProjection = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  INITIALIZE SETTINGS  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input format and prepare for following analysis
if isfield(handles.Pipeline.(['Segmented' inObjectOfInterest]),'Format')
    switch handles.Pipeline.(['Segmented' inObjectOfInterest]).Format
        case 'SegmentationCC'
            CC = handles.Pipeline.(['Segmented' inObjectOfInterest]).Label;
        otherwise
            error(['Image processing was canceled in the ',ModuleName,' module because input has unknown format'])
    end
end



%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%

switch inProjectionType
    
    
    case 'OccupiedPlanes'
    [XY, ~, ~] = OccupiedPlanes(CC);
        
    otherwise
        error(['Image processing was canceled in the ',ModuleName,' module because it does not support the desired output format'])
end
        
%%%%%%%%%%%%%%%%%%%%%%
%%% VISUALIZATION %%%%
%%%%%%%%%%%%%%%%%%%%%%

figure; imagesc(XY);

%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE OUTPUT   %%%%
%%%%%%%%%%%%%%%%%%%%%%

%%% Saves the adjusted image to the handles structure so it can be used by
%%% subsequent modules.

handles.Pipeline.(inOutputProjection) = XY;


%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

function [XY XZ YZ] = OccupiedPlanes(CC)
%% make projection which indicates how many z-layers have an object


%initialize
numPixels = cellfun(@numel,CC.PixelIdxList);
countPixels = sum(numPixels,2);
RCZO=NaN(countPixels,4);

%loop through each object and retreive positional information
%(RowColumnZplaneObjectid)
numLastIX=0;
for k=1:CC.NumObjects
    
    numFirstIX=numLastIX+1;
    numLastIX=numLastIX+numPixels(1,k);
    
    [RCZO(numFirstIX:numLastIX,1) RCZO(numFirstIX:numLastIX,2) RCZO(numFirstIX:numLastIX,3)] = ind2sub(CC.ImageSize,CC.PixelIdxList{1,k});
    RCZO(numFirstIX:numLastIX,4) = k;
    
end

% Clear initial pixel mearuremnts
CC.PixelIdxList={};


% for XY Plane
XY = uint16(zeros(CC.ImageSize(1),CC.ImageSize(2)));
SubRC = sub2ind2(size(XY), RCZO(:,[1 2]));
unSubRC = unique(SubRC);
CountPerPx = histc(SubRC,unSubRC);
XY(unSubRC)=CountPerPx;

% for XZ Plane
XZ= uint16(zeros(CC.ImageSize(2),CC.ImageSize(3)));
SubRC = sub2ind2(size(XZ), RCZO(:,[2 3]));
unSubRC = unique(SubRC);
CountPerPx = histc(SubRC,unSubRC);
XZ(unSubRC)=CountPerPx;

% for YZ Plane
YZ = uint16(zeros(CC.ImageSize(1),CC.ImageSize(3)));
SubRC = sub2ind2(size(YZ), RCZO(:,[1 3]));
unSubRC = unique(SubRC);
CountPerPx = histc(SubRC,unSubRC);
YZ(unSubRC)=CountPerPx;


















