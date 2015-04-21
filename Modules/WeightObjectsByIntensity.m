function handles = WeightObjectsByIntensity(handles)

% Help for the Measure Nuclear Spots module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Estimates the number of spots weighting spots by intensity if the area of 
% the given spots in larger than 1 pixel.%
% *************************************************************************
%
%
% How it works:
% Measurement:              Feature Number:
% Estimated number of spots |       1
% 
% All spots which area is one pixel have the weight of one, while spots
% with a larger area have the weight of the total intensity divided by the 
% provided estimate for spot intensity, the result is then rounded up to 
% the closest integer. 
% 
% Authors:
% Nico Battich
%
% $Revision: 4538 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the greyscale images you want to measure?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the parent (Nuclei) objects that you want to measure?
%choiceVAR02 = Image
%choiceVAR02 = Do not use
%infotypeVAR02 = objectgroup
ParentObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What did you call the children (Spots) objects that you want to measure?
%choiceVAR03 = Image
%choiceVAR03 = Do not use
%infotypeVAR03 = objectgroup
ChildrenObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = What is the object (spot) intensity?
%defaultVAR04 = 0.05
SpotIntensity = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,4}));

%textVAR05 = What is the background levels?
%defaultVAR05 = 0.0016
Background = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = What is the run identifier?
%defaultVAR06 = 005
identifier = char(handles.Settings.VariableValues{CurrentModuleNum,6});



%%%VariableRevisionNumber = 2


%%% Set up the window for displaying the results
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    CPfigure(handles,'Text',ThisModuleFigureNumber);
    columns = 1;
end


% Get current Cycle ID
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;


%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
ParentObject = CPretrieveimage(handles,['Segmented', ParentObjectName],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
ChildrenObject = CPretrieveimage(handles,['Segmented', ChildrenObjectName],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieve intensity image
IntensityImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','DontCheckScale');
IntensityImageSub = IntensityImage-Background;
IntensityImageSub(IntensityImageSub<0)=0;

% Estimate the number of children from intensities
ParentsId = unique(ParentObject(:)); ParentsId(ParentsId==0) = [];
spots_area = regionprops(ChildrenObject,'Area');
spots_area = cat(1,spots_area(:).Area);

ix = ismember(ChildrenObject(:),find(spots_area==1));
single_pixel_spots = zeros(size(ChildrenObject));
single_pixel_spots(ix) = 1;

ix = ismember(ChildrenObject(:),find(spots_area>1));
many_pixel_spots = ChildrenObject;
many_pixel_spots(~ix)=0;

spot_num = zeros(size(ParentsId));

for i = 1:length(ParentsId)
    f = ParentObject==ParentsId(i) & many_pixel_spots>0;
    count_large_spots = sum(IntensityImageSub(f(:)))/SpotIntensity;
    f = ParentObject==ParentsId(i) & single_pixel_spots>0;
    count_single_spots = sum(f(:));
    spot_num(i) = count_large_spots+count_single_spots;
end

handles.Measurements.(ParentObjectName).(['EstimatedSpotNumber_',identifier,'_',ParentObjectName,'Features']) = {'Estimated_Number_Of_Spots'};
handles.Measurements.(ParentObjectName).(['EstimatedSpotNumber_',identifier,'_',ParentObjectName]){SetBeingAnalyzed} = round(spot_num);


