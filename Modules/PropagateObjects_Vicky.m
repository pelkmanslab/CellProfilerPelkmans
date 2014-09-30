function handles = PropagateObjects_Vicky(handles)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% Load images from handles
SegmentedInputImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName);
% SegmentedCellsImage = CPretrieveimage(handles,['Segmented', 'Cells'],ModuleName);

% % Propagate only objects that lie within cells
% SegmentedImage = zeros(size(SegmentedInputImage));
% SegmentedImage((logical(SegmentedInputImage) + logical(SegmentedCellsImage))==2) = 1;

% Label objects that were selected for propagation
SegmentedImageLabel = bwlabel(logical(SegmentedInputImage));
IntensityImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');


%%%%%%%%%%%%%%%%%%%%%
%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Create threshold image for input into subsequent propagation function
% Calculate threshold using Otsu's method 
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

%%% Smooth image
SmoothDisk = getnhood(strel('disk',1,0));%minimum that has to be done to avoid problems with bwtraceboundary

%%% Fill holes
FinalImage = imfill(DilatedImage);

% %%% Label objects
% FinalImage = imdilate(imerode(FinalImage,SmoothDisk),SmoothDisk);

% ObjectAreas = cell2mat(struct2cell(regionprops(FinalImage,'Area')))';
% ValidObjectIndices = find(ObjectAreas>10);
% FinalImage = bwlabel(ismember(FinalImage,ValidObjectIndices));%relabel to get continuous indices


%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

if ~CPisHeadless()
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber)
        %%% Calculates the OriginalColoredLabelMatrixImage for displaying in the figure
        %%% window in subplot(2,1,1).
        OriginalColoredLabelMatrixImage = CPlabel2rgb(handles,SegmentedInputImage);
        %%% Calculates the ShrunkenColoredLabelMatrixImage for displaying in the figure
        %%% window in subplot(2,1,2).
        PropagatedColoredLabelMatrixImage = CPlabel2rgb(handles,FinalImage);

        %%% Activates the appropriate figure window.
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            CPresizefigure(OriginalColoredLabelMatrixImage,'TwoByOne',ThisModuleFigureNumber)
        end%%% A subplot of the figure window is set to display the original image.
        subplot(2,1,1);
        CPimagesc(OriginalColoredLabelMatrixImage,handles);
        title([ObjectName, ' cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        subplot(2,1,2);
        CPimagesc(IntensityImage,handles);
        title([PropagatedObjectName,' cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        B = bwboundaries(FinalImage,'holes');
        hold on
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
        end
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the final segmented label matrix image to the handles structure.
fieldname = ['Segmented',PropagatedObjectName];
handles.Pipeline.(fieldname) = FinalImage;

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
tmp = regionprops(FinalImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(PropagatedObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};