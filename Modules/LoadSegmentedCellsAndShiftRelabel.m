function handles = LoadSegmentedCellsAndShiftRelabel(handles)

% Help for the LoadSegmentedCellsAndShiftRelabel module:
% Category: Other
%
% SHORT DESCRIPTION:
% Loads object segmentation from iBRAIN's SEGMENTATION directory and shifts
% it by a global offset and relabels individual objects: eg.: useful if
% channel with spots, which will be linked later, is shifted
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1727 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object?
%defaultVAR01 = Cells
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = x correction in pixels.
%defaultVAR02 = 0
x_correction = str2num(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = y correction in pixels.
%defaultVAR03 = 0
y_correction = str2num(handles.Settings.VariableValues{CurrentModuleNum,3});


%%%VariableRevisionNumber = 2


%%%%%%%%%%%%%%%%
%%% ANALYSIS %%%
%%%%%%%%%%%%%%%%

% get the filename of the object setgementation as stored by
% SaveSegmentedCells
strOrigImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});

% format the segmentation file name
matDotIndices=strfind(strOrigImageName,'.');
% new CP apparently removes file extensions from image names
if ~isempty(matDotIndices)
    strOrigImageName = strOrigImageName(1,1:matDotIndices(end)-1);
end
strSegmentationFileName = [strOrigImageName,'_Segmented',ObjectName,'.png'];

% get the SEGMENTATION directory
% If output dir is BATCH directory, assume SEGMENTATION directory
strSegmentationDir = handles.Current.DefaultOutputDirectory;
if strcmp(getlastdir(strSegmentationDir),'BATCH')
    strSegmentationDir = strrep(strSegmentationDir, [filesep,'BATCH'],[filesep,'SEGMENTATION']);
end

% load segmentation image
strFilePath = fullfile(strSegmentationDir,strSegmentationFileName);
if fileattrib(strFilePath)
    matSegmentationImage = double(imread(fullfile(strSegmentationDir,strSegmentationFileName)));
else
    error('%s: looking for segmentation file ''%s''. Does not exist!',mfilename,strFilePath)
end

%%%%%%%%%%%% Now Shift the Image  %%%%%%%%%%%%%%%%%%%%%%%%%

%%% get the output image
OrigImage = matSegmentationImage;
ShiftedImage = zeros(size(OrigImage));

% define the change
if y_correction > 0
    final_index_one = abs(y_correction)+1:size(OrigImage,1);
    orig_index_one = 1:size(OrigImage,1)-abs(y_correction);
elseif y_correction < 0
    final_index_one = 1:size(OrigImage,1)-abs(y_correction);
    orig_index_one = abs(y_correction-1):size(OrigImage,1);%GG changed +1 to -1
else
    final_index_one = 1:size(OrigImage,1);
    orig_index_one = 1:size(OrigImage,1);
end
    
if x_correction > 0
    final_index_two = abs(x_correction+1):size(OrigImage,2);
    orig_index_two = 1:size(OrigImage,2)-abs(x_correction);
elseif x_correction < 0
    final_index_two = 1:size(OrigImage,2)-abs(x_correction);
    orig_index_two = abs(x_correction-1):size(OrigImage,2);%GG changed +1 to -1
else
    final_index_two = 1:size(OrigImage,2);
    orig_index_two = 1:size(OrigImage,2);
end
    
% Create Shifted Image
ShiftedImage(final_index_one,final_index_two)=OrigImage(orig_index_one,orig_index_two);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
 
RelabeledSeg = zeros(size(ShiftedImage));

fNoZero = ShiftedImage ~= 0;
if any(fNoZero(:))
    pixelsWoZero = ShiftedImage(fNoZero);
    [~, ~, reMapped] = unique(pixelsWoZero);
    RelabeledSeg(fNoZero) = reMapped;
end


%%%%%%%%%%%%% reintegrate into logic of module %%%%%%%%%%%

matSegmentationImage = double(RelabeledSeg);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);

    % RGB color
    ColoredLabelMatrixImage = CPlabel2rgb(handles,matSegmentationImage);
    CPimagesc(ColoredLabelMatrixImage,handles);
    %imagesc(ColoredLabelMatrixImage);
    
    title(sprintf('Loaded %s segmentation , cycle # %d',ObjectName,handles.Current.SetBeingAnalyzed));
end



%%%%%%%%%%%%%%%%%%%%%
%%% STORE RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%


% These fields are stored in identifyprimautomatic, might be that other
% segmentation/object detection modules store more fields... check
% (secondary, and identifyprimlog...)

%%% Saves the segmented image, not edited for objects along the edges or
%%% for size, to the handles structure.
fieldname = ['UneditedSegmented',ObjectName];
handles.Pipeline.(fieldname) = matSegmentationImage;

%%% Saves the segmented image, only edited for small objects, to the
%%% handles structure.
fieldname = ['SmallRemovedSegmented',ObjectName];
handles.Pipeline.(fieldname) = matSegmentationImage;

%%% Saves the final segmented label matrix image to the handles structure.
fieldname = ['Segmented',ObjectName];
handles.Pipeline.(fieldname) = matSegmentationImage;

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
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(matSegmentationImage(:));

%%% Saves the location of each segmented object
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(matSegmentationImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

end
