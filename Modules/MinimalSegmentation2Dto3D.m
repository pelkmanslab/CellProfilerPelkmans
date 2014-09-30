function handles = MinimalSegmentation2Dto3D(handles)


% Help for Segmentation2Dto3D
% Category: Object Processing

% This module loads individual images of one stack and checks for presence
% of non 0 values. Non 0 values are then connected to individual 3D
% objects. The output of this file corresponds to a CC structure obtained
% by bwconncomp, which uses linear indices per object. 



drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = What did you call the images you want to process?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = FusedSegmentation
%infotypeVAR02 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What is the minimal amount of voxels of one object?
%defaultVAR03 = 1000
numMinimalSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,3}));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  INITIALIZE SETTINGS  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Select List of Files
strFilenames = handles.Pipeline.(ImageName).FileNames;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  PROCESS IMAGES       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sequentially load individual images and create binary mask where the value
%is not 0. In conventional Labelmatrices, this is where a object is

BWImages = cell(1,size(strFilenames,2));
for k = 1:size(strFilenames,2)
    BWImages{k} = imread(strFilenames{1,k})>0;
end
BWImages=cat(3,BWImages{:});

% Create linear indices of objects within 3D stack
CC = bwconncomp(BWImages);

% discard objects with a too small volume
numPixels = cellfun(@numel,CC.PixelIdxList);
TooSmallObject = numPixels<numMinimalSize;
CC.PixelIdxList(TooSmallObject)=[];
CC.NumObjects=CC.NumObjects-sum(TooSmallObject);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  SAVE SEGMENTATION    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldname = ['Segmented', ObjectName];
handles.Pipeline.(fieldname).Label = CC;
handles.Pipeline.(fieldname).Format = 'SegmentationCC'; % Define Format of Output. Here it is Segmentation CC which corresponds to a full structure obtained by bwconncomp

%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = CC.NumObjects;