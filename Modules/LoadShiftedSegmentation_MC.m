function handles = LoadShiftedSegmentation_MC(handles)

% Help for the LoadShiftedSegmentation_MC module:
% Category: Other
%
% SHORT DESCRIPTION:
% Module, which loads Segmentation according to another multiplexing cycle.
% In contrast to the _MP modules, there is no cropping. Regions, which are
% not covered, are substituted with background
%
% Module creates an output, which allows to map back individual object IDs
% to their original object IDs
%
% Note that no shifting of intensity images is required or intended.
%
% Identical labelling of different segmentations (such as Nuclei and Cells)
% can be obtained by specifying both in same module - Multiple instances of
% LoadShiftedSegmentation_MC module act independently of each other
% *************************************************************************
%
% Authors:
%   Thomas Stoeger
%   original template (AlignObjects_MPcycle) by Markus Herrmann
%
% $Revision: 1879 $


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object?
%defaultVAR01 = Cells
%infotypeVAR01 = objectgroup indep
InputObjectName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = For which object? (ignore by using / )
%defaultVAR02 = /
%infotypeVAR02 = objectgroup indep
InputObjectName{2}  = char(handles.Settings.VariableValues{CurrentModuleNum,2});


%%%VariableRevisionNumber = 12



%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZATION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Objects To Process
f = ~cell2mat(cellfun(@(x) strcmp(x, '/'), InputObjectName,'UniformOutput',false));
NamesOfObjectsToShift = InputObjectName(f);
numNamesOfObjectsToShift = sum(f);

if numNamesOfObjectsToShift == 0
   error('No object specified'); 
end

% get the file name of the object segmentation as stored by
% SaveSegmentedCells
strReferenceImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});
currentTIFFfolder = handles.Measurements.Image.PathNames{1}{1,1};
currentPlatePath = getbasedir(currentTIFFfolder);

% Load Shift descriptor
FullShiftDescriptor = mcyc.loadShiftDescriptor(currentPlatePath);

% get shift information for current site
fprintf('Processing %s ', strReferenceImageName);
[SiteSpecificShiftDescriptor] = mcyc.getShiftForSingleSite(FullShiftDescriptor,strReferenceImageName);

xShift = SiteSpecificShiftDescriptor.xShift;
yShift = SiteSpecificShiftDescriptor.yShift;
maxShift = SiteSpecificShiftDescriptor.maxShift;

fprintf('X: %d Y: %d \n', xShift, yShift);

% Check, if any Segmentation should be used
suppressShiftAccordingToSettingFile =  SiteSpecificShiftDescriptor.noShiftIndex;
outsideOfMaximalShift = (abs(xShift) > maxShift) || (abs(yShift) > maxShift);
useNoSegmentation  = suppressShiftAccordingToSettingFile || outsideOfMaximalShift;


%%% Create shifted segmentation image or suppress segmenation (as defined
%%% above)
OrigSegmentationImagePerObject = cell(1,numNamesOfObjectsToShift);
ShiftedSegmentationsPerObject = cell(1,numNamesOfObjectsToShift);
for j=1:numNamesOfObjectsToShift
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% GET SEGMENTATION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    ObjectName = NamesOfObjectsToShift{j};
    currSegmentation = ['_Segmented' ObjectName];
    
    [path_cells_segmentation] = ...
        mcyc.getPathToSegmentationFileFromOtherCycle(SiteSpecificShiftDescriptor,strReferenceImageName,currentPlatePath,currSegmentation);
    
    if ~any(fileattrib(path_cells_segmentation))
        error('Segmentation file does not exist')
    else
        OrigSegmentationImage = double(imread(path_cells_segmentation));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% SHIFT IMAGES     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    if useNoSegmentation == false
        
        x_correction = -xShift; % segmentation has to be shifted opposite to what an image would have to be shifted
        y_correction = -yShift;
        
        ShiftedSegmentationImage = mcyc.shiftImage(OrigSegmentationImage,x_correction,y_correction);
        
    else  % Define case where no segmentation should be used
        ShiftedSegmentationImage = zeros(size(OrigSegmentationImage),'double');
    end
    
    OrigSegmentationImagePerObject{j} = OrigSegmentationImage;
    ShiftedSegmentationsPerObject{j} = ShiftedSegmentationImage;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SYNCHRONIZE SEGMENTATIONS OF MULTIPLE OBJECTS     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Explanation: remove objects, which are not present in at least one other
% segmentation image, this means that upon unique with 'sort' in latter
% part of code, all objects, which are present in same site will have same
% remapping ; Note: In contrast to the original procdure in alignObjects, 
% which was discarding full sites (approx 1/4 of all data), only few single
% cells are removed

% Safety check, if segementations, which should be mapped have same size
if any(var(cell2mat(cellfun(@(x) size(x), ShiftedSegmentationsPerObject,'UniformOutput',false)'),[],1) ~= 0)
    error('Image dimensions differ');
else
    ImageDimensions = size(ShiftedSegmentationsPerObject{1});
end


% set objects, which only occur in one of the reference segmentations
ReferenceIds = cell(1,numNamesOfObjectsToShift);
for j=1:numNamesOfObjectsToShift
    ReferenceIds{j} = unique(ShiftedSegmentationsPerObject{j}(:));
end

for j=1:numNamesOfObjectsToShift
    for jj = 1:numNamesOfObjectsToShift
        if j == jj
            continue;
        else
            isInOtherSegmentation = ismember(ShiftedSegmentationsPerObject{j}(:),ReferenceIds{jj});
            if any(~isInOtherSegmentation)
                ShiftedSegmentationsPerObject{j}(~isInOtherSegmentation) = 0;   % Set to background
            end
        end
    end
end


%%%%%%%%%%%%%%%
%%% RELABEL  %%
%%%%%%%%%%%%%%%

RelabeledSegmentationPerObject = cell(1,numNamesOfObjectsToShift);
OrigLabelsPerObject = cell(1,numNamesOfObjectsToShift);
for j=1:numNamesOfObjectsToShift
    
    %%% Relabel image
    ShiftedSegmentationImage = ShiftedSegmentationsPerObject{j};
    [uniqueOriginalPixelIds, posOfUnique, newLabelForSinglePixels] = unique(ShiftedSegmentationImage(:),'sorted'); % sort ensures same match depending upon original ID!
    
    % Correct new labels for background pixels, if present
    if any(ShiftedSegmentationImage(:) < 0)      % Must be fulfilled for the following logic
        error('Currently only supports segmentations with non-negative values')
    else  % default assumption for background pixels
        hasSegmentationValueOfZero = uniqueOriginalPixelIds(1) == 0;
        if any(hasSegmentationValueOfZero)
            newLabelForSinglePixels = newLabelForSinglePixels-1;
        end
    end
    
    % Create relabelled image
    RelabeledSegmentation = zeros(ImageDimensions);
    RelabeledSegmentation(:) = newLabelForSinglePixels;
    
    % Trace back to original ID of segmentation
    origLabels = ShiftedSegmentationImage(posOfUnique);
    if any(hasSegmentationValueOfZero)
        if length(origLabels) > 1
            origLabels = origLabels(2:end);
        else
            origLabels = NaN;
        end
    end
    
    RelabeledSegmentationPerObject{j} = RelabeledSegmentation;
    OrigLabelsPerObject{j} = origLabels;
end

%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

for j=1:numNamesOfObjectsToShift
    ObjectName = NamesOfObjectsToShift{j};
    origLabels = OrigLabelsPerObject{j};
    RelabeledSegmentation = RelabeledSegmentationPerObject{j};
    
    %%% Saves the segmented image, not edited for objects along the edges or
    %%% for size, to the handles structure.
    fieldname = ['UneditedSegmented',ObjectName];
    handles.Pipeline.(fieldname) = RelabeledSegmentation;
    
    %%% Saves the segmented image, only edited for small objects, to the
    %%% handles structure.
    fieldname = ['SmallRemovedSegmented',ObjectName];
    handles.Pipeline.(fieldname) = RelabeledSegmentation;
    
    %%% Saves the final segmented label matrix image to the handles structure.
    fieldname = ['Segmented',ObjectName];
    handles.Pipeline.(fieldname) = RelabeledSegmentation;
    
    
    %%% Saves the ObjectCount, i.e., the number of segmented objects.
    %%% See comments for the Threshold saving above
    if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
        handles.Measurements.Image.ObjectCountFeatures = {};
        handles.Measurements.Image.ObjectCount = {};
    end
    column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,ObjectName));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    ObjCount = max(RelabeledSegmentation(:));
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;
    
    
    %%% Saves the location of each segmented object
    % Follow CP convention for empty images (e.g.: as in IdentifySecondary
    % module)
    handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
    Centroid = [0 0];		
    if ObjCount ~= 0 
        tmp = regionprops(RelabeledSegmentation,'Centroid');
        Centroid = cat(1,tmp.Centroid);
    end
    handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
    %%% save relation to original object ID to handles
    handles.Measurements.(ObjectName).UnshiftedObjectIdFeatures{handles.Current.SetBeingAnalyzed} = 'UnshiftedObjectId';
    handles.Measurements.(ObjectName).UnshiftedObjectId{handles.Current.SetBeingAnalyzed} = origLabels;
    
end


%%%%%%%%%%%%%%%
%%% DISPLAY %%%
%%%%%%%%%%%%%%%


drawnow

if ~CPisHeadless()
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber)
        
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        
        subplot(2,2,1);
        CPimagesc(OrigSegmentationImagePerObject{1},handles);
        title(sprintf('Original Segmentation on ''%s'' Image', NamesOfObjectsToShift{1}));
        colormap('jet')
        
        subplot(2,2,2);
        CPimagesc(RelabeledSegmentationPerObject{1},handles);
        title(sprintf('Relabelled Segmentation on ''%s'' Image', NamesOfObjectsToShift{1}));
        colormap('jet')
        
        if numNamesOfObjectsToShift > 1 % only display second object, if specified
            
            subplot(2,2,3);
            CPimagesc(OrigSegmentationImagePerObject{2},handles);
            title(sprintf('Original Segmentation on ''%s''', NamesOfObjectsToShift{2}));
            colormap('jet')
            
            subplot(2,2,4);
            CPimagesc(RelabeledSegmentationPerObject{2},handles);
            title(sprintf('Relabelled Segmentation on ''%s''', NamesOfObjectsToShift{2}));
            colormap('jet')
            
        end
        
        drawnow
    end
end


end