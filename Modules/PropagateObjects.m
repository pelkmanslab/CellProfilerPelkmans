function handles = PropagateObjects(handles)

% Help for the Propagate Objects module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Propagates identified objects using IdentifySecPropagateSubfunction.
% *************************************************************************
%
% See also IdentifySecondary.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by Pelkmans group.
% Copyright 2014.
%
% Authors:
%   Markus Herrmann
%
% Website: http://www.pelkmanslab.org
%
% $Revision: 1718 $


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the objects that you want to propagate?
%infotypeVAR01 = objectgroup
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the propagated objects?
%defaultVAR02 = propagatedNuclei
%infotypeVAR02 = objectgroup indep
PropagatedObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What did you call the itensity image that you want to use for propagation?
%infotypeVAR03 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Propagation factor, value between 0 and 1 (0=Intensity, 1=Distance)
%defaultVAR04 = 0.03
PropagFactor = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,4}));

%textVAR05 = Threshold correction factor (for thresholding of intensity image)
%defaultVAR05 = 0.8
correctionFactor = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = Select an automatic thresholding method.
%choiceVAR06 = Otsu Global
%choiceVAR06 = Otsu Adaptive
%choiceVAR06 = Otsu PerObject
%choiceVAR06 = MoG Global
%choiceVAR06 = MoG Adaptive
%choiceVAR06 = MoG PerObject
%choiceVAR06 = Background Global
%choiceVAR06 = Background Adaptive
%choiceVAR06 = Background PerObject
%choiceVAR06 = RobustBackground Global
%choiceVAR06 = RobustBackground Adaptive
%choiceVAR06 = RobustBackground PerObject
%choiceVAR06 = RidlerCalvard Global
%choiceVAR06 = RidlerCalvard Adaptive
%choiceVAR06 = RidlerCalvard PerObject
Threshold = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu custom

%textVAR07 = Lower and upper bounds on threshold, in the range [0,1]
%defaultVAR07 = 0,1
ThresholdRange = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = For MoG thresholding, what is the approximate percentage of image covered by objects?
%choiceVAR08 = 10%
%choiceVAR08 = 20%
%choiceVAR08 = 30%
%choiceVAR08 = 40%
%choiceVAR08 = 50%
%choiceVAR08 = 60%
%choiceVAR08 = 70%
%choiceVAR08 = 80%
%choiceVAR08 = 90%
pObject = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%%%VariableRevisionNumber = 2


%%% Checks that the Min and Max threshold bounds have valid values
index = strfind(ThresholdRange,',');
if isempty(index)
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max threshold bounds are invalid.'])
end
MinimumThreshold = ThresholdRange(1:index-1);
MaximumThreshold = ThresholdRange(index+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% Load images from handles
SegmentedImage = double(CPretrieveimage(handles,['Segmented', ObjectName],ModuleName));
SegmentedImageLabel = bwlabel(logical(SegmentedImage));
IntensityImage = double(CPretrieveimage(handles,ImageName,ModuleName,'DontCheckColor','DontCheckScale'));%double(CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale'));


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Create forground mask for propagation step
ThresOtsu = CPthreshold_mask(Threshold,pObject,MinimumThreshold,MaximumThreshold,correctionFactor,IntensityImage,SegmentedImage,ModuleName);
%ThresOtsu = multithresh(IntensityImage.*65536); % reverse CP rescaling !!!
IntensityImageThres = IntensityImage > ThresOtsu;%./65536;

%%% Re-dilate objects up to border of the original object mask while retaining newly gained objects (see IdentifySecondary.m)
PropagatedImage = IdentifySecPropagateSubfunction(SegmentedImageLabel,IntensityImage,logical(IntensityImageThres),PropagFactor);
PropagatedImage = bwlabel(PropagatedImage);


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Calculates the OriginalColoredLabelMatrixImage for displaying in the figure
    %%% window in subplot(2,1,1).
    OriginalColoredLabelMatrixImage = CPlabel2rgb(handles,SegmentedImage);
    %%% Calculates the ShrunkenColoredLabelMatrixImage for displaying in the figure
    %%% window in subplot(2,1,2).
    PropagatedColoredLabelMatrixImage = CPlabel2rgb(handles,PropagatedImage);

    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OriginalColoredLabelMatrixImage,'TwoByOne',ThisModuleFigureNumber)
    end%%% A subplot of the figure window is set to display the original image.
    subplot(2,1,1);
    CPimagesc(OriginalColoredLabelMatrixImage,handles);
    title([ObjectName, ' cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,1,2);
    CPimagesc(PropagatedColoredLabelMatrixImage,handles);
    title([PropagatedObjectName,' cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the final segmented label matrix image to the handles structure.
fieldname = ['Segmented',PropagatedObjectName];
handles.Pipeline.(fieldname) = PropagatedImage;

%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,PropagatedObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {['ObjectCount ' PropagatedObjectName]};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(PropagatedImage(:));

%%% Saves the location of each segmented object
handles.Measurements.(PropagatedObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(PropagatedImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(PropagatedObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};