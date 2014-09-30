function handles = LoadSegmentedObjects_MPcycle(handles)

% Help for the LOADSEGMENTEDOBJECTS_MPCYCLE module:
% Category: Other
%
% SHORT DESCRIPTION:
% Loads object segmentations from a user defined SEGMENTATION directory. To
% provide path to SEGMENTATION directory use shiftDescriptor.json file (see
% aligncycles.calculateCycleShift.m).
% *************************************************************************
%
% Author:
%    Markus Herrmann <markus.herrmann@imls.uzh.ch>
%


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the objects that you want to measure later?
%defaultVAR01 = Do not use
%infotypeVAR01 = objectgroup indep
ObjectNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Type "Do not use" in unused boxes.
%defaultVAR02 = Do not use
%infotypeVAR02 = objectgroup indep
ObjectNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = 
%defaultVAR03 = Do not use
%infotypeVAR03 = objectgroup indep
ObjectNameList{3} = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 =
%defaultVAR04 = Do not use
%infotypeVAR04 = objectgroup indep
ObjectNameList{4} = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 =
%defaultVAR05 = Do not use
%infotypeVAR05 = objectgroup indep
ObjectNameList{5} = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 =
%defaultVAR06 = Do not use
%infotypeVAR06 = objectgroup indep
ObjectNameList{6} = char(handles.Settings.VariableValues{CurrentModuleNum,6});




%%%%%%%%%%%%%%%%
%%% ANALYSIS %%%
%%%%%%%%%%%%%%%%

% get the filename of the object setgementation as stored by
% SaveSegmentedCells
strOrigImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});

% format the segmentation file name
matDotIndices = strfind(strOrigImageName,'.');
% new CP apparently removes file extensions from image names
if ~isempty(matDotIndices)
    strOrigImageName = strOrigImageName(1,1:matDotIndices(end)-1);
end

% load shift descriptors file only once and store to handles
if ~isfield(handles,'shiftDescriptor')
    % define path to json file
    % If output dir is BATCH directory, assume ALIGNCYCLES directory
    strAlignCyclesDir = handles.Current.DefaultOutputDirectory;
    if strcmp(getlastdir(strAlignCyclesDir),'BATCH')
        strAlignCyclesDir = strrep(strAlignCyclesDir, [filesep,'BATCH'],[filesep,'ALIGNCYCLES']);
    end
    jsonDescriptorFile = fullfile(strAlignCyclesDir,'shiftDescriptor.json');
    % load json file
    handles.shiftDescriptor = loadjson(jsonDescriptorFile);
end

% get filename trunk from handles and modify filename accordingly
strSegmentationFileNameTrunk = handles.shiftDescriptor.SegmentationFileNameTrunk;

% obtain filename and load segmented image
SegmentationImages = cell(1,length(ObjectNameList));
for i = 1:length(ObjectNameList)
    ObjectName = ObjectNameList{i};
    if strcmpi(ObjectName,'Do not use')
        continue
    end
    
    % built absolute SEGMENTATION path from relative path stored in handles
    strSegmentationDir = [strrep(handles.Current.DefaultOutputDirectory, [filesep,'BATCH'], filesep), handles.shiftDescriptor.SegmentationDirectory];
    
    % built full segmentation image filename from filename trunk stored in handles
    strSegmentationFileName = [regexprep(strOrigImageName,'.+(_\w{1}\d{2}_)',sprintf('%s$1',strSegmentationFileNameTrunk)),'_Segmented',ObjectName,'.png'];
    
    % load segmentation image
    strFilePath = fullfile(strSegmentationDir,strSegmentationFileName);
    if fileattrib(strFilePath)
        SegmentationImages{i} = double(imread(fullfile(strSegmentationDir,strSegmentationFileName)));
    else
        error('%s: looking for segmentation file ''%s''. Does not exist!',mfilename,strFilePath)
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    
    for i = 1:length(ObjectNameList)
        ObjectName = ObjectNameList{i};
        if strcmpi(ObjectName,'Do not use')
            continue
        end
        % RGB color
        subplot(2,3,i);
        ColoredLabelMatrixImage = CPlabel2rgb(handles,SegmentationImages{i});
        CPimagesc(ColoredLabelMatrixImage,handles);
        title(sprintf('Loaded %s segmentation , cycle # %d',ObjectName,handles.Current.SetBeingAnalyzed));
    end
end



%%%%%%%%%%%%%%%%%%%%%
%%% STORE RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%


% These fields are stored in identifyprimautomatic, might be that other
% segmentation/object detection modules store more fields... check
% (secondary, and identifyprimlog...)

for i = 1:length(ObjectNameList)
    ObjectName = ObjectNameList{i};
    if strcmpi(ObjectName,'Do not use')
        continue
    end
    
    %%% Saves the segmented image, not edited for objects along the edges or
    %%% for size, to the handles structure.
    fieldname = ['UneditedSegmented',ObjectName];
    handles.Pipeline.(fieldname) = SegmentationImages{i};
    
    %%% Saves the segmented image, only edited for small objects, to the
    %%% handles structure.
    fieldname = ['SmallRemovedSegmented',ObjectName];
    handles.Pipeline.(fieldname) = SegmentationImages{i};
    
    %%% Saves the final segmented label matrix image to the handles structure.
    fieldname = ['Segmented',ObjectName];
    handles.Pipeline.(fieldname) = SegmentationImages{i};
    
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
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(SegmentationImages{i}(:));
    
    %%% Saves the location of each segmented object
    handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
    tmp = regionprops(SegmentationImages{i},'Centroid');
    Centroid = cat(1,tmp.Centroid);
    handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
end

end
