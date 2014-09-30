function handles = BackgroundImageCalculate(handles)

% Help for the BackgroundImageCalculation module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Calculates a background image exploiting the cell segmentations
% *************************************************************************
%
%
% Authors:
%   Nico Battich
%
%
%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the first object in which the mask should be based on?
%choiceVAR01 = Do not use
%infotypeVAR01 = objectgroup
ObjectNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the second object in which the mask should be based on?
%choiceVAR02 = Do not use
%infotypeVAR02 = objectgroup
ObjectNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What did you call the first object which sould be substracted from the mask?
%choiceVAR03 = Do not use
%infotypeVAR03 = objectgroup
ObjectNameList{3} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = What did you call the second object which sould be substracted from the mask?
%choiceVAR04 = Do not use
%infotypeVAR04 = objectgroup
ObjectNameList{4} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = What is te quantile value which should set as background?
%defaultVAR05 = 0.95
numQuanThres = str2num(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Smoothing is done by a Gaussian filter to both illumination functions. Enter the size of the filter.
%defaultVAR06 = 100
SmoothingSize = str2num(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Enter the factor 'x' for the sigma calculation of the Gaussian filter. where sigma = x*size.
%defaultVAR07 = 0.5
SmoothingSigma = str2num(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 =  What did you call the image which mark the segmentation file base name?
%infotypeVAR08 = imagegroup
BaseImageName = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 =What did you call the image you whant to base the background on?
%infotypeVAR09 = imagegroup
OrigImageName = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = What do you want to call the background function?
%defaultVAR10 = BGImageBlue
%infotypeVAR10 = imagegroup indep
BackgroundImage = char(handles.Settings.VariableValues{CurrentModuleNum,10});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% The background calculation is done in the first cycle only.
%BaseImageName = 'OrigBlue';
if handles.Current.SetBeingAnalyzed == 1
    
    %get a list of all images to be used
    OrigImageList = cellfun(@(x) fullfile(handles.Pipeline.(strcat('Pathname',OrigImageName)),x),...
        handles.Pipeline.(strcat('FileList',OrigImageName)), 'uniformoutput',false)';
    
    
    % get the SEGMENTATION directory
    % If output dir is BATCH directory, assume SEGMENTATION directory
    strSegmentationDir = handles.Current.DefaultOutputDirectory;
    if strcmp(getlastdir(strSegmentationDir),'BATCH')
        strSegmentationDir = strrep(strSegmentationDir, [filesep,'BATCH'],[filesep,'SEGMENTATION']);
    end
    
    SegmentationImageList = cell(1,4);
    
    %%% get the segmentation base name
    BaseImageList = cellfun(@(x) fullfile(handles.Pipeline.(strcat('Pathname',BaseImageName)),x),...
        handles.Pipeline.(strcat('FileList',BaseImageName)), 'uniformoutput',false)';
    
    
    
    %get the segmentation file list
    for iObjects = 1:length(ObjectNameList)
        
        if ~strcmpi(ObjectNameList{iObjects},'Do not use')
            
            %segmentation file names
            SegmentationImageList{iObjects} = cellfun(@(x) fullfile(strSegmentationDir,[x(1:end-4),'_Segmented',ObjectNameList{iObjects},'.png']),...
                handles.Pipeline.(strcat('FileList',BaseImageName)), 'uniformoutput',false)';
            
        end
        
    end
    
    % Loop all sites is the image set is less than 2000 else randomly sample 2000 images (bad bad bad bad)
    
    % init image counter
    intImageCounter = 1;
    
    % i think we should randomly sample images...
    RandImageIX = randperm(length(OrigImageList));
    if length(RandImageIX) > 2000
        RandImageIX = RandImageIX(1:2000);
        numIterations = 100;
    else %make sure the image number is divisable by 20
        numIterations = floor(length(RandImageIX)/20);
        RandImageIX = RandImageIX(1:numIterations*20);
    end
    
    
    for iX = 1:numIterations-1
        
        intImageCounter = intImageCounter + 20;
        %intCounter2 = intCounter2 + 1;
        
        %get the original images
        tic
        cellTempOrigImages = arrayfun(@(x) double(imread(OrigImageList{x})),RandImageIX(intImageCounter-20:intImageCounter-1),'uniformoutput',false);
        %cellTempImages = arrayfun(@(x) CPimread(ImageList{x},handles),matRandIndx,'uniformoutput',false);
        %cellTempImages = cellfun(@log10,cellTempImages,'uniformoutput',false);
        %matTempImages = cell2mat(cellfun(@(x) x(:),cellTempImages,'uniformoutput',false)');
        if iX == 1
            %initialize some parameters
            QuantilesImages = nan(size(cellTempOrigImages{1}(:)),numIterations);
        end
        
        
        
        cellTempSegImages = cell(1,4);
        for iObjects = 1:length(ObjectNameList)
            if ~isempty(SegmentationImageList{iObjects})
                %segmentation file names
                cellTempSegImages{iObjects} = arrayfun(@(x) (imread(SegmentationImageList{iObjects}{x})),RandImageIX(intImageCounter-20:intImageCounter-1),'uniformoutput',false);
            end
        end
      
        
        
        % make the base mask
        if isempty(SegmentationImageList{2})
            cellMaskBase = cellTempSegImages{1};
        elseif isempty(SegmentationImageList{1})
            cellMaskBase = cellTempSegImages{2};
        elseif isempty(SegmentationImageList{1}) && isempty(SegmentationImageList{2})
            error('you IDIOT: can you enter a mask object please!!!!!  (someting on the first or second field)')
        else
            cellMaskBase = cellfun(@(a,b) or(a,b), cellTempSegImages{1},cellTempSegImages{2},'uniformoutput',false);
        end
        
        
        % make the mask modifier
        if isempty(SegmentationImageList{4})
            cellMaskModifier = cellTempSegImages{3};
        elseif isempty(SegmentationImageList{3})
            cellMaskModifier = cellTempSegImages{4};
        elseif isempty(SegmentationImageList{1}) && isempty(SegmentationImageList{2})
            cellMaskModifier = [];
        else
            cellMaskModifier = cellfun(@(a,b) or(a,b), cellTempSegImages{3},cellTempSegImages{4},'uniformoutput',false);
        end
        
        %make final mask
        if isempty(cellMaskModifier)
            cellFinaMask = cellMaskBase;
        else
            cellFinaMask = cellfun(@(a,b) ~xor(a,b),cellMaskBase,cellMaskModifier,'uniformoutput',false);
        end
        
        %linearize the matrices
        matTempImages = cell2mat(cellfun(@(x) x(:),cellTempOrigImages,'uniformoutput',false));
        matTempMask = cell2mat(cellfun(@(x) x(:),cellFinaMask,'uniformoutput',false));
            
        
        matMaskedImages = matTempImages;
        matMaskedImages(matTempMask) = nan;
        matMaskedImages = matMaskedImages+1;
        %QuantilesImages(:,iX) = quantile(matMaskedImages,numQuanThres,2);
        matMaskedImages(isnan(matMaskedImages))=0;
        QuantilesImages(:,iX) = getcolumn(sort(matMaskedImages,2),round(size(matMaskedImages,2)*numQuanThres));
        QuantilesImages(QuantilesImages(:,iX)==0,iX) = nan;
        QuantilesImages(:,iX) = QuantilesImages(:,iX)-1;
        toc
    end
    
    
    FinalImage = nanmean(QuantilesImages,2);
    FinalImage = reshape(FinalImage,size(cellTempOrigImages{1},1),size(cellTempOrigImages{1},2));
    
    %smooth the image!
    H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize*SmoothingSigma);
    FinalImageFiltered =  imfilter(FinalImage,H,'replicate');
    
    %save to handle structure
    handles.Pipeline.(BackgroundImage) = FinalImageFiltered/65535;
    handles.Pipeline.([BackgroundImage,'_NotFiltered']) = FinalImage/65535;
    
    drawnow
    
end

%display the results
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
CPfigure(handles,'Image',ThisModuleFigureNumber);
imagesc(handles.Pipeline.(BackgroundImage));
colormap('JET')
colorbar

drawnow
