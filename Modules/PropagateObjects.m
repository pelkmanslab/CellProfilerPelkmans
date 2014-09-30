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

%%%VariableRevisionNumber = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% Load images from handles
SegmentedImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName);
SegmentedImageLabel = bwlabel(logical(SegmentedImage));
IntensityImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%ThresOtsu = multithresh(IntensityImage); %ultimately, one may use CPthreshold !!!
[~,ThresOtsu] = CPthreshold(handles,'Otsu Global','10%','Do not use','Do not use',correctionFactor,IntensityImage,ImageName,ModuleName);

ThresCorrected = ThresOtsu * correctionFactor;
IntensityImageThres = IntensityImage > ThresCorrected;

%%% Re-dilate objects up to border of the original object mask while retaining newly gained objects (see IdentifySecondary.m)
PropagatedImage = IdentifySecPropagateSubfunction(SegmentedImageLabel,IntensityImage,logical(IntensityImageThres),PropagFactor);

%%% Create nicely outlined objects (see MeasureLocalizationOfSpots.m)
if ~max(max(PropagatedImage))==0
    DilatedImage = DilateBackground(PropagatedImage);
else
    DilatedImage = PropagatedImage;
end


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
    PropagatedColoredLabelMatrixImage = CPlabel2rgb(handles,DilatedImage);

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
handles.Pipeline.(fieldname) = DilatedImage;

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
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(DilatedImage(:));

%%% Saves the location of each segmented object
handles.Measurements.(PropagatedObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(DilatedImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(PropagatedObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};