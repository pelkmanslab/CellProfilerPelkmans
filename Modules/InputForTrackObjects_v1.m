function handles = InputForTrackObjects_v1(handles,varargin)

%%%%WARNING:CHANGES MADE HERE BY MAT 110728 to fit with CP v.104553%%%
% Help for the Track Objects module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Allows tracking objects throughout sequential frames of a movie, so that
% each object maintains a unique identity in the output measurements.
% *************************************************************************
% This module must be run after the object to be tracked has been 
% identified using an Identification module (e.g., IdentifyPrimAutomatic).
%
% Settings:
%
% Tracking method:
% Choose between the methods based on which is most consistent from frame
% to frame of your movie. For each, the maximum search distance that a 
% tracked object will looked for is specified with the Neighborhood setting
% below:
%
%   Overlap - Compare the amount of overlaps between identified objects in 
%   the previous frame with those in the current frame. The object with the
%   greatest amount of overlap will be assigned the same label. Recommended
%   for movies with high frame rates as compared to object motion.
%       
%   Distance - Compare the distance between the centroid of each identified
%   object in the previous frame with that of the current frame. The 
%   closest objects to each other will be assigned the same label.
%   Distances are measured from the perimeter of each object. Recommended
%   for movies with lower frame rates as compared to object motion, but
%   the objects are clearly separable.
%
%   Measurement - Compare the specified measurement of each object in the 
%   current frame with that of objects in the previous frame. The object 
%   with the closest measurement will be selected as a match and will be 
%   assigned the same label. This selection requires that you run the 
%   specified Measurement module previous to this module in the pipeline so
%   that the measurement values can be used to track the objects. 
%
% Catagory/Feature Name or Number/Image/Scale:
% Specifies which type of measurement (catagory) and which feature from the
% Measure module will be used for tracking. Select the feature name from 
% the popup box or see each Measure module's help for the numbered list of
% the features measured by that module. Additional details such as the 
% image that the measurements originated from and the scale used as
% specified below if neccesary.
%
% Neighborhood:
% This indicates the region (in pixels) within which objects in the
% next frame are to be compared. To determine pixel distances, you can look
% at the markings on the side of each image (shown in pixel units) or
% using the ShowOrHidePixelData Image tool (under the Image Tools menu of 
% any CellProfiler figure window)
%
% How do you want to display the tracked objects?
% The objects can be displayed as a color image, in which an object with a 
% unique label is assigned a unique color. This same color is maintained
% throughout the object's lifetime. If desired, a number identifiying the 
% object is superimposed on the object.
%
% What number do you want displayed?
% The displayed number is the unique label assigned to the object or the
% progeny identifier.
%
% Do you want to calculate statistics:
% Select whether you want statistics on the tracked objects to be added to
% the measurements for that object. The current statistics are collected:
%
% Features measured:       Feature Number:
% TrajectoryX           |       1
% TrajectoryY           |       2
% DistanceTraveled      |       3
% IntegratedDistance    |       4
% Linearity             |       5
% LostObjectCount       |       6
% NewObjectCount        |       7
%
% In addition to these, the following features are also recorded: Label, 
% Lifetime as a per-object measurement, and the number of unique objects 
% that have appeared and dissappeared in each frame.
%
% Desscription of each feature:
%   Label: Each tracked object is assigned a unique identifier (label). 
%   Results of splits or merges are seen as new objects and assigned a new
%   label.
%
%   Trajectory: The direction of motion (in x and y coordinates) of the 
%   object from the previous frame to the curent frame.
%
%   Distance traveled: The distance traveled by the object from the 
%   previous frame to the curent frame (calculated as the magnititude of 
%   the distance traveled vector).
%
%   Lifetime: The duration (in frames) of the object. The lifetime begins 
%   at the frame when an object appears and is ouput as a measurement when
%   the object disappears. At the final frame of the image set/movie, the 
%   lifetimes of all remaining objects are ouput.
%
%   Integrated distance: The total distance traveled by the object during
%   the lifetime of the object
%   
%   Linearity: A measure of how linear the object trajectity is during the
%   object lifetime. Calculated as (distance from initial to final 
%   location)/(integrated object distance). Value is in range of [0,1].
%
%   LostObjectCount: Number of objects that appear in the previous frame
%   but have no identifiable child in the current frame
%
%   NewObjectCount: Number of objects that appear in the current frame but
%   have no identifiable parent in the previous frame 
%
% What do you want to call the image with the tracked objects?
% Specify a name to give the image showing the tracked objects. This image
% can be saved with a SaveImages module placed after this module.
%
% Additional notes:
%
% In the figure window, a popupmenu allows you to display the objects as a
% solid color or as an outline with the current objects in color and the
% previous objects in white.
%
% Since the movie is processed sequentially by frame, it cannot be broken
% up into batches for execution on a distributed cluster.
%
% If running on a cluster and saving the colored image with text labels,
% the labels will not show up in the final result. This is a limitation of
% using MATLAB's hardcopy command.
%
% See also: Any of the Measure* modules, IdentifyPrimAutomatic

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003--2008.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 8988 $

% MBray 2009_03_25: Comments on variables for pyCP upgrade
%
% Recommended variable order (setting, followed by current variable in MATLAB CP)
% (1) What tracking method would you like to use? (TrackingMethod)
% (2) What did you call the objects you want to track? (ObjectName)
% (3a) What category of measurement do you want to use (MeasurementCategory)
% (3b) What feature do you want to use? (MeasurementFeature)
% (3c) (If the answer to (3b) involves a scale) What scale was used to 
%      calculate the feature? (SizeScale) 
%      (If the answer to (3b) involves an image) What image was used to
%      calculate the feature? (ImageName)
% (4) Within what pixel distance will objects be considered to find a
%   potential match? (PixelRadius)
% (5) How do you want to display the tracked objects? (DisplayType)
% (6) What do you want to call the resulting image with tracked,
%   color-coded objects? Type "Do not use" to ignore. (DataImage)
%
% (i) Option (3) should appear only if the user selects "Measuremement" in (2)
% (ii) The Measurement category/feature/image/scale settings in (3a,b,c) should only be shown if
% the measurement hierarchy requires it.
% (iii) The setting to collect statistics (CollectStatistics) should be
% removed as it should be collected by default.

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Choose a tracking method
%choiceVAR01 = Overlap
%choiceVAR01 = Distance
%choiceVAR01 = Measurements
%inputtypeVAR01 = popupmenu
TrackingMethod = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What did you call the objects you want to track?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = (Tracking by Measurements only) What category of measurement you want to use?
%inputtypeVAR03 = popupmenu category
MeasurementCategory = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = (Tracking by Measurements only) What feature you want to use? (Enter the feature number or name - see help for details)
%defaultVAR04 = Do not use
%inputtypeVAR04 = popupmenu measurement
MeasurementFeature = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = (Tracking by Measurements only) For INTENSITY, AREAOCCUPIED or TEXTURE features, which image's measurements do you want to use?
%infotypeVAR05 = imagegroup
%inputtypeVAR05 = popupmenu
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = (Tracking by Measurements only) For TEXTURE, RADIAL DISTRIBUTION, OR NEIGHBORS features, what previously measured size scale (TEXTURE OR NEIGHBORS) or previously used number of bins (RADIALDISTRIBUTION) do you want to use?
%defaultVAR06 = 1
%inputtypeVAR06 = popupmenu scale
SizeScale = str2double(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Choose the neighborhood (in pixels) within which objects will be evaluated to find a potential match.
%defaultVAR07 = 50
PixelRadius = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = How do you want to display the tracked objects?
%choiceVAR08 = Color and Number
%choiceVAR08 = Color Only
%inputtypeVAR08 = popupmenu
DisplayType = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = If you chose an option with "Number" above, select the number you want displayed (ProgenyID is not yet working)
%choiceVAR09 = Object ID
%choiceVAR09 = Progeny ID
%inputtypeVAR09 = popupmenu
LabelMethod = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%textVAR10 = Do you want to calculate statistics?
%choiceVAR10 = No
%choiceVAR10 = Yes
%inputtypeVAR10 = popupmenu
CollectStatistics = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = What do you want to call the resulting image with tracked, color-coded objects? Type "Do not use" to ignore.
%defaultVAR11 = Do not use
%infotypeVAR11 = imagegroup indep
DataImage = char(handles.Settings.VariableValues{CurrentModuleNum,11});

%textVAR12 = How many sites are composing your movie? .
%defaultVAR12 = 25
NumOfSites = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,12}));

%textVAR13 = How many timepoints are composing your movie? .
%defaultVAR13 = 200
NumOftimepoints = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,13}));

