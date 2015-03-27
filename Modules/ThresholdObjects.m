function handles = ThresholdObjects(handles)

% Help for the Threshold Objects module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Thresholds objects using a modified version of CPthreshold.
% *************************************************************************
%
% See also CPthreshold.m subfunction.

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

%textVAR01 = What did you call the objects that you want to threshold?
%infotypeVAR01 = objectgroup
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the intentified objects?
%defaultVAR02 = thresOrganelles
%infotypeVAR02 = objectgroup indep
ThresObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What did you call the itensity image that you want to use for thresholding?
%infotypeVAR03 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Do you want to calculate threshold only within objects?
%defaultVAR04 = Yes
UseMask = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Threshold correction factor (manually fine tune thresholding)
%defaultVAR05 = 1
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

%textVAR09 = Do you want to use a fixed threshold? Enter value you want to use.
%defaultVAR09 = /
ThresholdFixed = char(handles.Settings.VariableValues{CurrentModuleNum,9});

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
IntensityImage = double(CPretrieveimage(handles,ImageName,ModuleName,'DontCheckColor','DontCheckScale'));%double(CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale'));


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

if handles.Current.SetBeingAnalyzed == 1
    
    %%% Calculate threshold only once and then apply this threshold to all
    %%% images to ensure the same threshold is applied to each object
    
    if ~strcmp(ThresholdFixed,'/')
        
        if strcmp(ThresholdFixed,'Pre') 
            %%% Load precalculated threshold values
            tmp = load(fullfile(handles.Current.DefaultOutputDirectory,['Measurements_Image_Threshold_',ImageName,'.mat']));
            tmp = eval(sprintf('tmp.handles.Measurements.Image.Threshold_%s',ImageName));
            tmp = cell2mat(tmp');
            
            %%% Take median of threshold values per image
            ThresholdValue = median(tmp);
            fprintf('%s: median of precalculated thresholds for objects ''%s'' is %d\n',mfilename,ThresObjectName,ThresholdValue);
 
        else
            %%% Use fixed value
            ThresholdValue = str2num(ThresholdFixed);
        end
        
    else
        
        %%% determine whether a mask should be used for threshold calculation
        if strcmp(UseMask,'Yes')
            MaskImage = logical(SegmentedImage);
        else
            MaskImage = [];
        end
        
        %%% Automatically calculate threshold
        ThresholdValue = CPthreshold_mask(Threshold,pObject,MinimumThreshold,MaximumThreshold,correctionFactor,IntensityImage,MaskImage,ModuleName);
        fprintf('%s: automatically calculated threshold for objects ''%s'' is %d\n',mfilename,ThresObjectName,ThresholdValue);
        
    end
    
    %%% store calculated threshold in handles
    handles.Threshold.(ThresObjectName) = ThresholdValue;
    
else    
    
    %%% load previously calculated threshold from handles
    ThresholdValue = handles.Threshold.(ThresObjectName); 
end

%%% Apply threshold
ThresholdImage = IntensityImage > ThresholdValue;
ThresholdImage = bwlabel(ThresholdImage);


%%% Smooth objects a bit to get nicer object outlines
%%% (this may also prevent weird things with area/shape measurements later)
SmoothDisk = getnhood(strel('disk',1,0));
ThresholdImage = bwlabel(imdilate(imerode(logical(ThresholdImage),SmoothDisk),SmoothDisk));


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
    PropagatedColoredLabelMatrixImage = CPlabel2rgb(handles,ThresholdImage);

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
    title([ThresObjectName,' cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the final segmented label matrix image to the handles structure.
fieldname = ['Segmented',ThresObjectName];
handles.Pipeline.(fieldname) = ThresholdImage;

%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ThresObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ThresObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(ThresholdImage(:));

%%% Saves the location of each segmented object
handles.Measurements.(ThresObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(ThresholdImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(ThresObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};