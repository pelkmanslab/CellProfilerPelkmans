function handles = LoadSegmentedObjectsCP3D(handles)

% Help for the LoadSegmentedCells module:
% Category: Other
%
% SHORT DESCRIPTION:
% Loads object segmentation from iBRAIN's SEGMENTATION directory; this is
% module is based upon LoadSegmentedCells
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object?
%defaultVAR01 = Spots
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});


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
    %%%% inserted to distinguish t
    bnHasCP3DSegmentation = false;
else
    strSegmentationFileName = [strOrigImageName,'_SegmentedCP3D',ObjectName,'.mat'];
    if fileattrib(strFilePath)
        strSegmentationCP3D = load(fullfile(strSegmentationDir,strSegmentationFileName));
        bnHasCP3DSegmentation = true;
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
    switch bnHasCP3DSegmentation
        case false     % 1:1 copy of code of Load SegmentedCells
            % RGB color
            ColoredLabelMatrixImage = CPlabel2rgb(handles,matSegmentationImage);
            CPimagesc(ColoredLabelMatrixImage,handles);
            %imagesc(ColoredLabelMatrixImage);
        case true
            % here the code for displaying stacks of segementation might be inserted
    end
    title(sprintf('Loaded %s segmentation , cycle # %d',ObjectName,handles.Current.SetBeingAnalyzed));
end



%%%%%%%%%%%%%%%%%%%%%
%%% STORE RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%

switch bnHasCP3DSegmentation
    case false      % 1:1 copy of code of Load SegmentedCells
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
        
    case true
        % note that in principle code could have been written that can be
        % applied on CP1 and CP3D, however introduce clear separation to
        % allow potentially independent handling
        fieldname = ['Segmented',ObjectName];
        handles.Pipeline.(fieldname) = strSegmentationCP3D;
end


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

switch bnHasCP3DSegmentation
    case false
        MaxObject =  max(matSegmentationImage(:));
    case true
        if isfield(strSegmentationCP3D,'Format')
            switch strSegmentationCP3D.Format
                case 'SegmentationCC'
                    
                    SegmentationCC = strSegmentationCP3D.Label;
                    MaxObject = SegmentationCC.NumObjects;
                otherwise
                    error(['Image processing was canceled in the ', ModuleName, ' module because the format ', strSegmentationCP3D.Format , 'is not supported.'])
            end
        else
            error(['Image processing was canceled in the ', ModuleName, ' module because no Format was specified within the Segmentation.'])
        end
end

handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = MaxObject;


%%% Saves the location of each segmented object
switch bnHasCP3DSegmentation
    case false
        handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
        tmp = regionprops(matSegmentationImage,'Centroid');
        Centroid = cat(1,tmp.Centroid);
        handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    case true
        
        %%% Saves the location of each segmented object
        switch size(SegmentationCC.ImageSize,2)
            case 2
                handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
            case 3
                handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY','CenterZ'};
            otherwise
                error(['Image processing was canceled in the ', ModuleName, ' module. Currently only centroids of 2D and 3D are supported. '])
        end
        
        if SegmentationCC.NumObjects ~= 0 % determine centroid, if at least one object
            tmp = regionprops(SegmentationCC,'Centroid');
            Centroid = cat(1,tmp.Centroid);
            if isempty(Centroid)   % keep the resettign to 0 0 found in other modules to remain consistent
                Centroid = [0 0];
            end
            handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
        end 
end
