function handles = LoadSegmentedObjects_MPcycle(handles)

% Help for the LOADSEGMENTEDOBJECTS_MPCYCLE module:
% Category: MPcycle
%
% SHORT DESCRIPTION:
% Loads object segmentations from a user defined SEGMENTATION directory. To
% provide path to SEGMENTATION directory use shiftDescriptor.json file (see
% aligncycles.calculateCycleShift.m).
% *************************************************************************
%
% Author:
%    Markus Herrmann <markus.herrmann@imls.uzh.ch>

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

%textVAR07 = Which segmentation do you want to load?
%choiceVAR07 = Original
%choiceVAR07 = Aligned
OriginalOrAligned = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu




%%%%%%%%%%%%%%%%
%%% ANALYSIS %%%
%%%%%%%%%%%%%%%%

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

% built absolute SEGMENTATION path from relative path stored in handles
strSegmentationDir = [strrep(handles.Current.DefaultOutputDirectory, [filesep,'BATCH'], filesep), handles.shiftDescriptor.SegmentationDirectory];
if strcmp(OriginalOrAligned,'Aligned')
    strSegmentationDir = strrep(strSegmentationDir,'SEGMENTATION','SEGMENTATION_ALIGNED');
end
if exist(strSegmentationDir,'dir')==2
    error('%s: segmentation directory ''%s'' does not exist!',mfilename,strSegmentationDir)
end

% obtain filename and load segmented image
SegmentationImages = cell(1,length(ObjectNameList));
for i = 1:length(ObjectNameList)
    ObjectName = ObjectNameList{i};
    if strcmpi(ObjectName,'Do not use')
        continue
    end
    
    % read segmentation folder content
    structSegDir = CPdir(fullfile(strSegmentationDir));
    cellSegFiles = {structSegDir(~[structSegDir.isdir]).name}; 
    % build filenames
    cellFieldnames = flatten(regexp(fieldnames(handles.Pipeline),'Filename.*','match'));
    cellLookup = handles.Pipeline.(cellFieldnames{1});
    strLookup = cellLookup{end};
    [~, strMicroscopeType] = check_image_position(strLookup);
    if strcmpi(strMicroscopeType, 'CV7K') 
        strExpr = char(flatten(regexp(strLookup,'_(\w{1}\d{2}_T\d{4}F\d{3})L\d{2}A\d{2}Z\d{2}C\d{2}','tokens')));
    elseif strcmpi(strMicroscopeType, 'Visi') 
        strExpr = char(flatten(regexp(strLookup,'_(s\d{04}_r\d{02}_c\d{02})_[A-Z]+_C\d{02}','tokens')));
    else
        error('Microscope type could not be determined. Currently implemented are ''Visi'' and ''CV7K''')
    end
    ixImage = cellfun(@(x) ~isempty(x),regexp(cellSegFiles,sprintf('%s.*Segmented%s',strExpr,ObjectName),'once'));
    strSegmentationFileName = char(cellSegFiles(ixImage));
    
    % load segmentation image
    strFilePath = fullfile(strSegmentationDir,strSegmentationFileName);
    if fileattrib(strFilePath)
        SegmentationImages{i} = double(imread(fullfile(strSegmentationDir,strSegmentationFileName)));
        fprintf('%s: loaded segmentation image: %s\n',mfilename,strSegmentationFileName);
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

function C = flatten(A)
% 
% C1 = flatten({{1 {2 3}} {4 5} 6})
% C2 = flatten({{'a' {'b','c'}} {'d' 'e'} 'f'})
% 
% Outputs:
% C1 = 
%     [1]    [2]    [3]    [4]    [5]    [6]
% C2 = 
%     'a'    'b'    'c'    'd'    'e'    'f'
%
% Copyright 2010  The MathWorks, Inc.
C = {};
for i=1:numel(A)  
    if(~iscell(A{i}))
        C = [C,A{i}];
    else
       Ctemp = flatten(A{i});
       C = [C,Ctemp{:}];
       
    end
end
