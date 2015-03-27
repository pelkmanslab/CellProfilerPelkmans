function handles = MeasureObjectThreshold(handles)

% Help for the Measure Object Threshold module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Measures several intensity features for identified objects.
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

%textVAR02 = What did you call the itensity image that you want to use for thresholding?
%infotypeVAR02 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = Threshold correction factor (manually fine tune thresholding)
%defaultVAR03 = 1
correctionFactor = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,3}));

%textVAR04 = Select an automatic thresholding method.
%choiceVAR04 = Otsu Global
%choiceVAR04 = Otsu PerObject
%choiceVAR04 = MoG Global
%choiceVAR04 = MoG PerObject
%choiceVAR04 = Background Global
%choiceVAR04 = Background PerObject
%choiceVAR04 = RobustBackground Global
%choiceVAR04 = RobustBackground PerObject
%choiceVAR04 = RidlerCalvard Global
%choiceVAR04 = RidlerCalvard PerObject
Threshold = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu custom

%textVAR05 = Lower and upper bounds on threshold, in the range [0,1]
%defaultVAR05 = 0,1
ThresholdRange = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = For MoG thresholding, what is the approximate percentage of image covered by objects?
%choiceVAR06 = 10%
%choiceVAR06 = 20%
%choiceVAR06 = 30%
%choiceVAR06 = 40%
%choiceVAR06 = 50%
%choiceVAR06 = 60%
%choiceVAR06 = 70%
%choiceVAR06 = 80%
%choiceVAR06 = 90%
pObject = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

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

%%% determine whether a mask should be used for threshold calculation
MaskImage = logical(SegmentedImage);

%%% Automatically calculate threshold
ThresholdValue = CPthreshold_mask(Threshold,pObject,MinimumThreshold,MaximumThreshold,correctionFactor,IntensityImage,MaskImage,ModuleName);
fprintf('%s: automatically calculated threshold for objects ''%s'' in image ''%s'' is %d\n',mfilename,ObjectName,ImageName,ThresholdValue);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

handles.Measurements.Image.(['Threshold_',ImageName,'Features']) = sprintf('%s_%s_min%s_max%s_corr%d',Threshold,pObject,MinimumThreshold,MaximumThreshold,correctionFactor);
handles.Measurements.Image.(['Threshold_',ImageName])(handles.Current.SetBeingAnalyzed) = {ThresholdValue};

