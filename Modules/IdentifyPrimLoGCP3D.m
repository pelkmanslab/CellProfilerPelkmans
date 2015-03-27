function handles = IdentifyPrimLoGCP3D(handles)
% Help for the IdentifyPrimLGCP3D module:
% Category: Object Processing
%

% SHORT DESCRIPTION:
% Will Determine Spots in 3D Image stacks after Laplacian Of Gaussian (LoG)
% enhancing of spots. Many of the input arguments are optional. Note that
% while an external script has to be run in order to choose robust values,
% manual selection of the parameters can often yield good estimates, if
% the signal is clear enough.
%
% WHAT DID YOU CALL THE IMAGES YOU WANT TO PROCESS?
% Object detection should be done on this image.
%
% HOW DO YOU WANT TO CALL THE OBJECTS IDENTIFIED BY THIS MODULE?
% This is the name of the the spots identified after thresholding the LoG
% image.
%
% WHICH SPOT ENHANCEMENT DO YOU WANT TO USE?
% You can either enhance spots only within their plane or alternatively
% with the including information from adjacent planes, as described by Raj
% et al. 2009.
%
% OBJECTSIZE
% This value corresponds to the approximate size of you spots in the 2D plane. It should
% be their diameter in pixels. The LoG will use a mask of this size to
% enhance radial signal of that size. Note that in practice the specific value
% does not affect the number of spots, if spots are bright (eg. pixel size 5
% or 6).
%
% INTENSITY QUANTA PER IMAGE
% Prior to spot detection the images are rescaled according to their
% intensity. Since the specific value of minimal and maximal intensities
% are frequently not robust across multiple images, intensity quantile are
% used instead. [0 1] would correspond to using the single dimmest pixel
% for minimal intensity and the single brightest pixel for maximal
% intensity. [0.01 0.90] would mean that the minimum intensity is derived
% from the pixel, which is the 1% brightest pixel of all and that the
% maximum intensity is derived from the pixel, which is the 90% brightest
% pixel .
%
% INTENSITY BORERS FOR INTENSITY RESCALING OF IMAGES
% Most extreme values that the image intensity minimum and image intensity
% maximum (as defined by the quanta) are allowed to have
% [LowestPossibleGreyscaleValueForImageMinimum
% HighestPossibleGreyscaleValueForImageMinimum
% LowestPossibleGreyscaleValueForImageMaximum
% HighestPossibleGreyscaleValueForImageMaximum]
% To ignore individual values, place a NaN.
% Note that these parameters very strongly depend upon the variability of
% your illumination source. When using a robust confocal microscope you can
% set the lowest and highest possible values to values,  which are very
% close (or even identical). If your light source is variable during the
% acquisition (which can be the case with Halogen lamps) you might choose 
% less strict borders to detect spots of varying intensites.
%
% THRESHOLD OF SPOT DETECTION
% This is the threshold value for spot detection. The higher it is the more
% stringent your spot detection is. Use external script to determine a
% threshold where the spot number is robust against small variations in the
% threshold.
%
% WHAT IS THE MINIMAL INTENSITY OF A VOXEL WITHIN A SPOT?
% Minimal greyscale value of a voxel, which a voxel has to have in order to
% be recognized to be within a spot. Opitonal argument to make spot
% detection even more robust against very dim spots. In practice, we have
% never observed that this parameter would have any influence on the spot
% detection. However, you might include it as an additional safety measure.
%
%
% [TS]
% *************************************************************************

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images you want to process?
%infotypeVAR01 = imagegroup
StackName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = Spots
%infotypeVAR02 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which Spot enhancement do you want to use?
%choiceVAR03 = 3D LoG, Raj
%choiceVAR03 = 2D LoG
iFilterType = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = ObjectSize
%defaultVAR04 = 6
iHsize = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Intensity Quanta Per Image
%defaultVAR05 = [0.01 0.99]
iImgLimes = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Intensity borders for intensity rescaling of images
%[MinOfMinintens MaxOfMinintens MinOfMaxintens MaxOfMaxintens]
%defaultVAR06 = [NaN 120 500 NaN]
iRescaleThr = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Threshold of Spot Detection
%defaultVAR07 = 0.01
iDetectionThr = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = What is the minimal intensity of a voxel within a spot?
%defaultVAR08 = /
iObjIntensityThr = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CHECK INPUT   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter Size
try
    iHsize = str2double(iHsize);
catch errFilterSize
    error(['Image processing was canceled in the ', ModuleName, ' module because the object size could not be converted to a number.'])
end

if iHsize<=2
    error(['Image processing was canceled in the ', ModuleName, ' module because the object size was too small. Has to be at least 3'])
end

% Intensity Quanta Of Image
[isSafe iImgLimes]= inputVectorsForEvalCP3D(iImgLimes,true);
if isSafe ==false
    error(['Image processing was canceled in the ', ModuleName, ' module because Intensity Quanta per Image contain forbidden characters.'])
end

% Rescale Thresholds
[isSafe iRescaleThr]= inputVectorsForEvalCP3D(iRescaleThr,true);
if isSafe ==false
    error(['Image processing was canceled in the ', ModuleName, ' module because Rescaling Boundaries contain forbidden characters.'])
end

% Detection Threshold
try
    iDetectionThr = str2double(iDetectionThr);
catch errDetectionThr
    error(['Image processing was canceled in the ', ModuleName, ' module because the Detection Threshold could not be converted to a number.'])
end

% Object intensity Threshold
if iObjIntensityThr == '/'
    iObjIntensityThr = [];
else
    try
        iObjIntensityThr = str2double(iObjIntensityThr);
    catch errObjIntensityThr
        error(['Image processing was canceled in the ', ModuleName, ' module because the Stepsize for deblending could not be converted to a number.'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%

Image = handles.Pipeline.(StackName);


op = fspecialCP3D(iFilterType,iHsize);
% Detect objects, note that input vectors are eval'ed
[ObjCount SegmentationCC] = ObjByFilter(Image,op,iDetectionThr,eval(iImgLimes),eval(iRescaleThr),iObjIntensityThr,false,[],[]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Save Segmentation to Pipeline and define format
fieldname = ['Segmented', ObjectName];
handles.Pipeline.(fieldname).Label = SegmentationCC;
handles.Pipeline.(fieldname).Format = 'SegmentationCC';

%%% Saves the ObjectCount, i.e. the number of segmented objects:
% This is saved in the .Image measurement so that different objects are
% stored together independend of whether they were initially derived
% from a 2D image or a 3D stack
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;


%%% Saves the location of each segmented object
switch size(SegmentationCC.ImageSize,2)
    case 2
        handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
    case 3
        handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY','CenterZ'};
    otherwise
        error(['Image processing was canceled in the ', ModuleName, ' module. Currently only centroids of 2D and 3D are supported. '])
        % note that it would be very easy to add more dimensions by removing this induced error.
        % The only reason for inducing this crash is to prevent
        % amiguity resulting from not having the dimensions not named
        % clearly. Thus if someone wants to add more, one should
        % conciously make a simple adaptation 
end

if SegmentationCC.NumObjects ~= 0 % determine centroid, if at least one object
    tmp = regionprops(SegmentationCC,'Centroid');
    Centroid = cat(1,tmp.Centroid);
    if isempty(Centroid)   % keep the resettign to 0 0 found in other modules to remain consistent
        Centroid = [0 0];
    end
    handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
end


%%%%%%%%%%%%%%%%%%%
%%% DISPLAY %%%%%%%
%%%%%%%%%%%%%%%%%%%


% Create Occupancy image, which shows how manz Z planes have an object
occImg = createOccupancyImage(SegmentationCC,'XY');
    

drawnow
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    
    % Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
                CPresizefigure(handles.Pipeline.(StackName),'TwoByTwo',ThisModuleFigureNumber);
    end
    
    % Make heatmap showing occupancy of individual pixels with objects
    imagesc(occImg);
    colormap('JET')
    colorbar
    
    %CPimagesc(occImg,handles);
    title(sprintf('Z planes with object. Total Amount of Objects is %d', ObjCount))
end



end
