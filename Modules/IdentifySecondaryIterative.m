function handles = IdentifySecondaryIterative(handles)

% Help for the IdentifySecondaryIterative module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Identifies objects (e.g. cell edges) using "seed" objects identified by
% an Identify Primary module (e.g. nuclei). See the help of Cell Profiler's
% original IdentifySecondary module for more details.
%
% In contrast to the orignal module, sequential rounds of watershedding are
% done. The outcome will be a very precise cell outline segmentation, which
% does not need a lot of human supervision (and thus greatly reduces
% working time). However, since this modulce can take up to 30 min on a
% single image, you might not want to use it, if you do not have access to
% massive parallel computing facilities.
%
% This module identifies secondary objects by sequential watershedding.
% This allows to combine the advantage of a high threshold correction
% factor (correct allocation of pixels to cells within crowded regions)
% with the advantage of a low threshold correction factor (detection of
% cellular periphery in sparse regions).
%
% To prevent false positives, if very few primary objects are present,
% limits for threshold can be used.
%
%
% SUGGESTION TO FIND MINIMAL BACKGROUND - ROBUST CAMERA / IMAGING
% CONDITIONS
% If using sCMOS cameras, the intensity value, which corresponds to
% background signal from the chip is ~100 greyscale values. Distinctive
% cell-specific signal from the cell outline stain can usually be detected
% at ~125 greyscale values. The minimal background value thus corresponds
% to 0.0019 (125/2^16)
%
% SUGGESTION TO FIND MINIMAL BACKGROUND - VIA IDENTIFYSECONDARY
% Get lowest threshold, which does not yet recognize background as
% cells.
% Make a CP pipeline with several IdentifySecondary modules. In each one
% select a different threshold correction value. Use OTSU GLOBAL and
% WATERSHEDDING
% Now start the pipeline. Of all modules, which do not recognize the
% background, use the one with the smallest threshold correction value
% (for us this is frequently around 0.5). Then manually write down the
% exact threshold value of this module. It will be displayed in the
% window opened by this module.
%
% Do not bother, whether the cells are
% correclty segmented. The only important point is that the outline of the
% cell has to be detected fully. Make sure that the test image is
% representative of your assay. Usually spreading cells do require a much
% lower threshold value than cells in crowded environemnts.
%
% SETTING UP THE IDENTIFYSECONDARYITERATIVE MODULE
%
% THRESHOLD CORRECTION FACTORS. IN DESCENDING RANKING.
% Should indicate many
% different thresholds. It starts with the most stringent and starts with
% the lowest. The lowest one should be so low that it would recognize the
% background as an object.
%
% Note that you do not care about:
% x the number of thresholds. The more, the better. Use supercomputing.
% Applying around 20 different ones usually gives very robust results. You
% can not select too many thresholds. You can save days of manual work by
% not selecting (a single) individual threshold(s).
% x the specific value of the lowest threshold: The lowest value has to be
% lower than the lowest threshold correction, which you tested previously.
% It should be a threshold value that recognizes the background. The
% separation from the background will be done by a later option. If the
% last threshold corrction value is too high, the periphery of spreading
% cells might be missed.
%
% An example for the range would be 1.1 1.05 1 0.95 0.9 0.85 0.8 0.75 0.7
% 0.6 0.58 0.55 0.50 0.45 0.4 0.35 0.3 0.25
%
% Usually the best gain in segmentation quality per threshold is achived
% with threshold correction factors close to the threshold correction
% factor, which you would choose in the normal IdentifySecondary module
%
% LOWER AND UPPER BONDS ON THRESHOLD.
% These values correspond to the minimal and maximal values that a
% threshold is allowed to have. They have the format
% SmallestThreshold,HighestThreshold
% For SmallestThreshold you should use the value obtained in c). Leaving
% the maximal value at 1 has worked fine for us all the time. Setting a
% minimal value will prevent recognition of the background as an object
%
%
% *************************************************************************
%
%
% How it works>
% This is a heavily adjusted version of the orignal IdentifySecondary module
% It also has lower memory requiremnts and fixes some bugs of the original
% module (shrinkage of nuclei, expansion of surrounding objects discarded by
% DiscardSinglePixel)
%
% 0) Obtain masks with proper foreground objects.
%
% 1)sequential watershedding
%
% 2)Then one label image is constructed. If a pixel is part of different
% objects at given threshold (which is likely in cell rich regions), it will be
% allocated to the threshold which was defined prior. eg. if thresholds
% specifications were 1 and 0.5 it would be attributed to the object
% identified with a threshold of 1. If threshold specification was 0.5 1
% it would be attributed to the object identified at 0.5
%
% 3) Cleaning up step. It could happen that an object would end up
% separated into multiple fragments (which in most cases would be
% biologically meaningless). Thus all fragments except the one, which
% includes the primary object, are set to background
%
% Authors:
%   Thomas Stoeger
%   Nico Battich
%   Lucas Pelkmans
%
%
% $Revision: 1808 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the primary objects you want to create secondary objects around?
%infotypeVAR01 = objectgroup
PrimaryObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = Cells
%infotypeVAR02 = objectgroup indep
SecondaryObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What did you call the images to be used to find the edges of the secondary objects?
%infotypeVAR03 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Lower and upper bounds on threshold, in the range [0,1]. Pixels with an intensity below the lower bound will never be part of a secondary object.
%defaultVAR04 = 0.0019,1
ThresholdRange = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Threshold correction factors. In descending ranking. eg 0.9 0.8 0.7
%defaultVAR05 = 2 1.5 1.3 0.9 0.7 0.6 0.58 0.55 0.50 0.45 0.4 0.35 0.3 0.25
iThresholdCorrection = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%%%VariableRevisionNumber = 3


import jtlib.SecondarySegmentation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

%%% Retrieves the preliminary label matrix image that contains the primary
%%% segmented objects which have only been edited to discard objects
%%% that are smaller than a certain size.  This image
%%% will be used as markers to segment the secondary objects with this
%%% module.  Checks first to see whether the appropriate image exists.
PrelimPrimaryLabelMatrixImage = CPretrieveimage(handles,['SmallRemovedSegmented', PrimaryObjectName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));

UseAsLabelInCaseNoBackgroundPresent = PrelimPrimaryLabelMatrixImage;
if any(PrelimPrimaryLabelMatrixImage(:) == 0)
    originalSegmentationHasBackground = true;
else
    originalSegmentationHasBackground = false;
end


%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects which will be used to weed out which objects are
%%% real - not on the edges and not below or above the specified size
%%% limits. Checks first to see whether the appropriate image exists.
EditedPrimaryLabelMatrixImage = CPretrieveimage(handles,['Segmented', PrimaryObjectName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));


%%% Checks that the Min and Max threshold bounds have valid values
index = strfind(ThresholdRange,',');
if isempty(index)
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max threshold bounds are invalid.'])
end

MinimumThreshold = str2double(ThresholdRange(1:index-1));
MaximumThreshold = str2double(ThresholdRange(index+1:end));

% Create vector containing the thresholds that should be tested
[isSafe, iThresholdCorrection]= inputVectorsForEvalCP3D(iThresholdCorrection,false);
if isSafe == false
    error(['Image processing was canceled in the ', ModuleName, ' module because input of threshold contained forbidden characters'])
else
    ThresholdCorrection = eval(iThresholdCorrection);
end


%======================== handles free zone - start ===========================

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%

[FinalLabelMatrixImage, EditedPrimaryBinaryImage, ThresholdArray] = SecondarySegmentation(OrigImage, PrelimPrimaryLabelMatrixImage, EditedPrimaryLabelMatrixImage, ...
                                                                                          ThresholdCorrection, MinimumThreshold, MaximumThreshold);
    

%======================== handles free zone - end =============================


if ~isfield(handles.Measurements,SecondaryObjectName)
    handles.Measurements.(SecondaryObjectName) = {};
end

if ~isfield(handles.Measurements,PrimaryObjectName)
    handles.Measurements.(PrimaryObjectName) = {};
end

handles = CPrelateobjects(handles,SecondaryObjectName,PrimaryObjectName,FinalLabelMatrixImage,EditedPrimaryLabelMatrixImage,ModuleName);
%%% [PLab] hack. save memory
clear EditedPrimaryLabelMatrixImage;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the final, segmented label matrix image of secondary objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented',SecondaryObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

%[PLab removed propagation]

%%% Saves the ObjectCount, i.e. the number of segmented objects.
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,SecondaryObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {SecondaryObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(FinalLabelMatrixImage(:));

%%% Saves the location of each segmented object
handles.Measurements.(SecondaryObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(FinalLabelMatrixImage,'Centroid');
%%% [PLab] hack. save memory.
Centroid = cat(1,tmp.Centroid);
if isempty(Centroid)
    Centroid = [0 0];   % follow CP's convention to save 0s if no object
end
handles.Measurements.(SecondaryObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

% [PLab] note that the following CP code would require additional
% calculations, which as default were always done and also used for
% visualization. If it should be included again, the code either has to be
% arranged back or , better, a check included, whether the outline should
% be saved
% %%% Saves images to the handles structure so they can be saved to the hard
% %%% drive, if the user requested.
% try
%     if ~strcmpi(SaveOutlines,'Do not save')
%         handles.Pipeline.(SaveOutlines) = LogicalOutlines;
%     end
% catch dummyError %[PLab] bugfix for error message
%     error(['The object outlines were not calculated by the ', ModuleName, ' module, so these images were not saved to the handles structure. The Save Images module will therefore not function on these images. This is just for your information - image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
% end


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

numThresholdsToTest = length(ThresholdCorrection);

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    
    if CPisHeadless == false
        
        %%%% [PLab] %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Code for visualization  %%%%%%%%%%%
        %%%%%%% Rearranged: Inculde visualization into a conditional statement starting
        %%%%%%% only on local machine, but not CPCluster
        
        
        %%% Calculates the ColoredLabelMatrixImage for displaying in the figure
        %%% window in subplot(2,2,2).
        ColoredLabelMatrixImage = CPlabel2rgb(handles,FinalLabelMatrixImage);
        %%% Calculates OutlinesOnOrigImage for displaying in the figure
        %%% window in subplot(2,2,3).
        %%% Note: these outlines are not perfectly accurate; for some reason it
        %%% produces more objects than in the original image.  But it is OK for
        %%% display purposes.
        %%% Maximum filters the image with a 3x3 neighborhood.
        MaxFilteredImage = ordfilt2(FinalLabelMatrixImage,9,ones(3,3),'symmetric');
        %%% Determines the outlines.
        IntensityOutlines = FinalLabelMatrixImage - MaxFilteredImage;
        %%% [PLab] hack.s ave memory.
        clear MaxFilteredImage;
        %%% Converts to logical.
        warning off MATLAB:conversionToLogical
        LogicalOutlines = logical(IntensityOutlines);
        %%% [PLab] hack.s ave memory.
        clear IntensityOutlines;
        warning on MATLAB:conversionToLogical
        
        % Determines the grayscale intensity to use for the cell outlines.
        %[PLab-HACK] so that images are not so dim!!!!
        ObjectOutlinesOnOrigImage = OrigImage;
        ObjectOutlinesOnOrigImage=ObjectOutlinesOnOrigImage-quantile(OrigImage(:), 0.025);
        ObjectOutlinesOnOrigImage(ObjectOutlinesOnOrigImage<0)=0;
        ObjectOutlinesOnOrigImage(ObjectOutlinesOnOrigImage>quantile(ObjectOutlinesOnOrigImage(:), 0.95))=quantile(ObjectOutlinesOnOrigImage(:), 0.95);
        LineIntensity = quantile(ObjectOutlinesOnOrigImage(:), 0.99);
        
        
        ObjectOutlinesOnOrigImage(LogicalOutlines) = LineIntensity;
        %%% Calculates BothOutlinesOnOrigImage for displaying in the figure
        %%% window in subplot(2,2,4).
        %%% Creates the structuring element that will be used for dilation.
        StructuringElement = strel('square',3);
        %%% Dilates the Primary Binary Image by one pixel (8 neighborhood).
        DilatedPrimaryBinaryImage = imdilate(EditedPrimaryBinaryImage, StructuringElement);
        %%% Subtracts the PrelimPrimaryBinaryImage from the DilatedPrimaryBinaryImage,
        %%% which leaves the PrimaryObjectOutlines.
        PrimaryObjectOutlines = DilatedPrimaryBinaryImage - EditedPrimaryBinaryImage;
        %%% [PLab] hack. save memory.
        clear DilatedPrimaryBinaryImage EditedPrimaryBinaryImage;
        BothOutlinesOnOrigImage = ObjectOutlinesOnOrigImage;
        BothOutlinesOnOrigImage(PrimaryObjectOutlines == 1) = LineIntensity;
        %%% [PLab] hack. save memory.
        clear PrimaryObjectOutlines LineIntensity;
        
        %%%%%%%%%%%%%%%%%%%%%%%% END OF INITIATION OF VISUALIZATION
        %%%%%%%%%%%%%%%%%%%%%%%% (rearrangement)
        
        
        %%% Activates the appropriate figure window.
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
        end
        ObjectCoverage = 100*sum(sum(FinalLabelMatrixImage > 0))/numel(FinalLabelMatrixImage);
        
        %[PLab] display range of thresholds. Which is useful if limits for treshold
        %should be used
        %     uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[0.25 0.01 .6 0.04],...
        %         'BackgroundColor',[.7 .7 .9],'HorizontalAlignment','Left','String',sprintf('Threshold:  %0.3f               %0.1f%% of image consists of objects',Threshold,ObjectCoverage),'FontSize',handles.Preferences.FontSize);
        
        
        ThresholdFirst  = ThresholdArray{1};
        ThresholdLast = ThresholdArray{numThresholdsToTest};
        uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[0.25 0.01 .6 0.04],...
            'BackgroundColor',[.7 .7 .9],'HorizontalAlignment','Left','String',sprintf('Threshold: Start %0.5f End %0.5f                %0.1f%% of image consists of objects',ThresholdFirst,ThresholdLast,ObjectCoverage),'FontSize',handles.Preferences.FontSize);
        
        %%% A subplot of the figure window is set to display the original image.
        subplot(2,2,1);
        CPimagesc(OrigImage,handles);
        title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        %%% A subplot of the figure window is set to display the colored label
        %%% matrix image.
        subplot(2,2,2);
        CPimagesc(ColoredLabelMatrixImage,handles);
        clear ColoredLabelMatrixImage
        title(['Outlined ',SecondaryObjectName]);
        %%% A subplot of the figure window is set to display the original image
        %%% with secondary object outlines drawn on top.
        subplot(2,2,3);
        CPimagesc(ObjectOutlinesOnOrigImage,handles);
        clear ObjectOutlinesOnOrigImage
        title([SecondaryObjectName, ' Outlines on Input Image']);
        %%% A subplot of the figure window is set to display the original
        %%% image with outlines drawn for both the primary and secondary
        %%% objects.
        subplot(2,2,4);
        CPimagesc(BothOutlinesOnOrigImage,handles);
        clear BothOutlinesOnOrigImage;
        title(['Outlines of ', PrimaryObjectName, ' and ', SecondaryObjectName, ' on Input Image']);
    end
end


end
