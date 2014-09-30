function handles = IdentifyPrimAutomatic_Region(handles)

% Help for the Identify Primary Automatic Region module (modified Identify
% Primary Automatic module):
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Identifies objects given only an image as input.
% *************************************************************************
%
% This module identifies primary objects (e.g. nuclei) in grayscale images
% that show bright objects on a dark background. The module has many
% options which vary in terms of speed and sophistication.
%
% Requirements for the images to be fed into this module:
% * If the objects are dark on a light background, they must first be
% inverted using the Invert Intensity module.
% * If you are working with color images, they must first be converted to
% grayscale using the Color To Gray module.
%
% Overview of the strategy ('Settings' below has more details):
%   Properly identifying primary objects (nuclei) that are well-dispersed,
% non-confluent, and bright relative to the background is straightforward
% by applying a simple threshold to the image. This is fast but usually
% fails when nuclei are touching. In CellProfiler, several automatic
% thresholding methods are available, including global and adaptive, using
% Otsu's (Otsu, 1979) and our own version of a Mixture of Gaussians
% algorithm (O. Friman, unpublished). 
%  Some identified nuclei are discarded or merged together if
% the user chooses. Incomplete nuclei touching the border of the image can
% be discarded. Objects smaller than a user-specified size range, which are
% likely to be fragments of real nuclei, can be discarded. Alternately, any
% of these small objects that touch a valid nucleus can be merged together
% based on a set of heuristic rules; for example similarity in intensity
% and statistics of the two objects.
%
% For more details, see the Settings section below and also the notation
% within the code itself (Developer's version).
%
% Malpica, N., de Solorzano, C. O., Vaquero, J. J., Santos, A., Vallcorba,
% I., Garcia-Sagredo, J. M., and del Pozo, F. (1997). Applying watershed
% algorithms to the segmentation of clustered nuclei. Cytometry 28,
% 289-297.
% Meyer, F., and Beucher, S. (1990). Morphological segmentation. J Visual
% Communication and Image Representation 1, 21-46.
% Ortiz de Solorzano, C., Rodriguez, E. G., Jones, A., Pinkel, D., Gray, J.
% W., Sudar, D., and Lockett, S. J. (1999). Segmentation of confocal
% microscope images of cell nuclei in thick tissue sections. Journal of
% Microscopy-Oxford 193, 212-226.
% Wahlby, C. (2003) Algorithms for applied digital image cytometry, Ph.D.,
% Uppsala University, Uppsala.
% Wahlby, C., Sintorn, I. M., Erlandsson, F., Borgefors, G., and Bengtsson,
% E. (2004). Combining intensity, edge and shape information for 2D and 3D
% segmentation of cell nuclei in tissue sections. J Microsc 215, 67-76.
%
% Settings:
%
% Typical diameter of objects, in pixel units (Min,Max):
% This is a very important parameter which tells the module what you are
% looking for. Most options within this module use this estimate of the
% size range of the objects in order to distinguish them from noise in the
% image. For example, for some of the identification methods, the smoothing
% applied to the image is based on the minimum size of the objects. A comma
% should be placed between the minimum and the maximum diameters. The units
% here are pixels so that it is easy to zoom in on objects and determine
% typical diameters. To measure distances easily, use the CellProfiler
% Image Tool, 'ShowOrHidePixelData', in any open window. Once this tool is
% activated, you can draw a line across objects in your image and the
% length of the line will be shown in pixel units. Note that for non-round
% objects, the diameter here is actually the 'equivalent diameter', meaning
% the diameter of a circle with the same area as the object.
%
% Discard objects outside the diameter range:
% You can choose to discard objects outside the specified range of
% diameters. This allows you to exclude small objects (e.g. dust, noise,
% and debris) or large objects (e.g. clumps) if desired. See also the
% FilterByObjectMeasurement module to further discard objects based on some
% other measurement. During processing, the window for this module will
% show that objects outlined in green were acceptable, objects outlined in
% red were discarded based on their size, and objects outlined in yellow
% were discarded because they touch the border.
%
% Try to merge 'too small' objects with nearby larger objects:
% Use caution when choosing 'Yes' for this option! This is an experimental
% functionality that takes objects that were discarded because they were
% smaller than the specified Minimum diameter and tries to merge them with
% other surrounding objects. This is helpful in cases when an object was
% incorrectly split into two objects, one of which is actually just a tiny
% piece of the larger object. However, this could be dangerous if you have
% selected poor settings which produce many tiny objects - the module
% will take a very long time and you will not realize that it is because
% the tiny objects are being merged. It is therefore a good idea to run the
% module first without merging objects to make sure the settings are
% reasonably effective.
%
% Select automatic thresholding method:
%    The threshold affects the stringency of the lines between the objects
% and the background. You can have the threshold automatically calculated
% using several methods, or you can enter an absolute number between 0 and
% 1 for the threshold (to see the pixel intensities for your images in the
% appropriate range of 0 to 1, use the CellProfiler Image Tool,
% 'ShowOrHidePixelData', in a window showing your image). There are
% advantages either way. An absolute number treats every image identically,
% but is not robust to slight changes in lighting/staining conditions
% between images. An automatically calculated threshold adapts to changes
% in lighting/staining conditions between images and is usually more
% robust/accurate, but it can occasionally produce a poor threshold for
% unusual/artifactual images. It also takes a small amount of time to
% calculate.
%    The threshold which is used for each image is recorded as a
% measurement in the output file, so if you find unusual measurements from
% one of your images, you might check whether the automatically calculated
% threshold was unusually high or low compared to the other images.
%    There are five methods for finding thresholds automatically, Otsu's
% method, the Mixture of Gaussian (MoG) method, the Background method, the
% Robust Background method and the Ridler-Calvard method. The Otsu method
% uses our version of the Matlab function graythresh (the code is in the
% CellProfiler subfunction CPthreshold). Our modifications include taking
% into account the max and min values in the image and log-transforming the
% image prior to calculating the threshold. Otsu's method is probably best
% if you don't know anything about the image, or if the percent of the
% image covered by objects varies substantially from image to image. If you
% know the object coverage percentage and it does not vary much from image
% to image, the MoG can be better, especially if the coverage percentage is
% not near 50%. Note, however, that the MoG function is experimental and
% has not been thoroughly validated. The Background method is very simple
% and is appropriate for images in which most of the image is background.
% It finds the mode of the histogram of the image, which is assumed to be
% the background of the image, and chooses a threshold at twice that value
% (which you can adjust with a Threshold Correction Factor, see below).
% This can be very helpful, for example, if your images vary in overall
% brightness but the objects of interest are always twice (or actually, any
% constant) as bright as the background of the image. The Robust background
% method trims the brightest and dimmest 5% of pixel intensities off first
% in the hopes that the remaining pixels represent a gaussian of intensity
% values that are mostly background pixels. It then calculates the mean and
% standard deviation of the remaining pixels and calculates the threshold
% as the mean + 2 times the standard deviation. The Ridler-Calvard method
% is simple and its results are often very similar to Otsu's. It chooses an
% initial threshold, and then iteratively calculates the next one by taking
% the mean of the average intensities of the background and foreground
% pixels determined by the first threshold, repeating this until the
% threshold converges.
%    You can also choose between Global, Adaptive, and Per object
% thresholding:
% Global: one threshold is used for the entire image (fast).
% Adaptive: the threshold varies across the image - a bit slower but
% provides more accurate edge determination which may help to separate
% clumps, especially if you are not using a clump-separation method (see
% below).
% Per object: if you are using this module to find child objects located
% *within* parent objects, the per object method will calculate a distinct
% threshold for each parent object. This is especially helpful, for
% example, when the background brightness varies substantially among the
% parent objects. Important: the per object method requires that you run an
% IdentifyPrim module to identify the parent objects upstream in the
% pipeline. After the parent objects are identified in the pipeline, you
% must then also run a Crop module as follows: the shape in which to crop
% is the name of the parent objects, and the image to be cropped is the one
% that you will want to use within this module to identify the children
% objects (e.g., ChildrenStainedImage). Then, set this
% IdentifyPrimAutomatic module to identify objects within the
% CroppedChildrenStainedImage.
%
% Threshold correction factor:
% When the threshold is calculated automatically, it may consistently be
% too stringent or too lenient. You may need to enter an adjustment factor
% which you empirically determine is suitable for your images. The number 1
% means no adjustment, 0 to 1 makes the threshold more lenient and greater
% than 1 (e.g. 1.3) makes the threshold more stringent. For example, the
% Otsu automatic thresholding inherently assumes that 50% of the image is
% covered by objects. If a larger percentage of the image is covered, the
% Otsu method will give a slightly biased threshold that may have to be
% corrected using a threshold correction factor.
%
% Lower and upper bounds on threshold:
% Can be used as a safety precaution when the threshold is calculated
% automatically. For example, if there are no objects in the field of view,
% the automatic threshold will be unreasonably low. In such cases, the
% lower bound you enter here will override the automatic threshold.
%
% Approximate percentage of image covered by objects:
% An estimate of how much of the image is covered with objects. This
% information is currently only used in the MoG (Mixture of Gaussian)
% thresholding but may be used for other thresholding methods in the future
% (see below).
%
% Special note on saving images: Using the settings in this module, object
% outlines can be passed along to the module OverlayOutlines and then saved
% with the SaveImages module. Objects themselves can be passed along to the
% object processing module ConvertToImage and then saved with the
% SaveImages module. This module produces several additional types of
% objects with names that are automatically passed along with the following
% naming structure: (1) The unedited segmented image, which includes
% objects that are outside the size range, can be saved using the name:
% UneditedSegmented + whatever you called the objects (e.g.
% UneditedSegmentedNuclei). (2) The segmented image which excludes objects
% smaller than your selected size range can be saved using the name:
% SmallRemovedSegmented + whatever you called the objects (e.g.
% SmallRemovedSegmented Nuclei).
%
% See also IdentifyPrimManual, IdentifySecondary.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Authors:
%   Anne E. Carpenter
%   Thouis Ray Jones
%   In Han Kang
%   Ola Friman
%   Steve Lowe
%   Joo Han Chang
%   Colin Clarke
%   Mike Lamprecht
%   Peter Swire
%   Rodrigo Ipince
%   Vicky Lay
%   Jun Liu
%   Chris Gang
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1879 $
%
%
% Modified: Markus Herrmann

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%%% NOTE: We cannot indent the variables or they will not be read
%%% properly.

%textVAR01 = What did you call the images you want to process?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = RegionObjects
%infotypeVAR02 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Typical diameter of objects, in pixel units (Min,Max):
%defaultVAR03 = 40,1000000
SizeRange = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Discard objects outside the diameter range?
%choiceVAR04 = Yes
%choiceVAR04 = No
ExcludeSize = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = Try to merge too small objects with nearby larger objects?
%choiceVAR05 = No
%choiceVAR05 = Yes
MergeChoice = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu

%textVAR06 = Select an automatic thresholding method or enter an absolute threshold in the range [0,1].  To choose a binary image, select "Other" and type its name.  Choosing 'All' will use the Otsu Global method to calculate a single threshold for the entire image group. The other methods calculate a threshold for each image individually. "Set interactively" will allow you to manually adjust the threshold during the first cycle to determine what will work well.
%choiceVAR07 = Otsu Global
%choiceVAR07 = Otsu Adaptive
%choiceVAR07 = Otsu PerObject
%choiceVAR07 = MoG Global
%choiceVAR07 = MoG Adaptive
%choiceVAR07 = MoG PerObject
%choiceVAR07 = Background Global
%choiceVAR07 = Background Adaptive
%choiceVAR07 = Background PerObject
%choiceVAR07 = RobustBackground Global
%choiceVAR07 = RobustBackground Adaptive
%choiceVAR07 = RobustBackground PerObject
%choiceVAR07 = RidlerCalvard Global
%choiceVAR07 = RidlerCalvard Adaptive
%choiceVAR07 = RidlerCalvard PerObject
%choiceVAR07 = All
%choiceVAR07 = Set interactively
Threshold = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu custom

%textVAR07 = Threshold correction factor
%defaultVAR07 = 1
ThresholdCorrection = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,7})); %#ok Ignore MLint

%textVAR08 = Lower and upper bounds on threshold, in the range [0,1]
%defaultVAR08 = 0,1
ThresholdRange = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = For MoG thresholding, what is the approximate percentage of image covered by objects?
%choiceVAR09 = 10%
%choiceVAR09 = 20%
%choiceVAR09 = 30%
%choiceVAR09 = 40%
%choiceVAR09 = 50%
%choiceVAR09 = 60%
%choiceVAR09 = 70%
%choiceVAR09 = 80%
%choiceVAR09 = 90%
pObject = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = Size of smoothing filter, in pixel units (if you are distinguishing between clumped objects). Enter 0 for low resolution images with small objects (~< 5 pixel diameter) to prevent any image smoothing.
%defaultVAR10 = Automatic
SizeOfSmoothingFilter = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = What do you want to call the outlines of the identified objects (optional)?
%defaultVAR11 = RegionOutlines
%infotypeVAR11 = outlinegroup indep
SaveOutlines = char(handles.Settings.VariableValues{CurrentModuleNum,11});

%textVAR12 = Do you want to fill holes in identified objects?
%choiceVAR12 = Yes
%choiceVAR12 = No
FillHolesOption = char(handles.Settings.VariableValues{CurrentModuleNum,12});
%inputtypeVAR12 = popupmenu

%textVAR13 = What do you want to call the region mask? 
%defaultVAR13 = RegionMask
%infotypeVAR13 = imagegroup indep
RegionMaskImageName = char(handles.Settings.VariableValues{CurrentModuleNum,13});

%%%VariableRevisionNumber = 12


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a variable,
%%% "OrigImage".
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

if islogical(OrigImage)
    error(['Image processing was canceled in the ', ModuleName, ' module because the input image is binary (black/white). The input image must be grayscale.']);
end

%%% Chooses the first word of the method name (removing 'Global' or 'Adaptive').
ThresholdMethod = strtok(Threshold);
%%% Checks if a custom entry was selected for Threshold, which means we are using an incoming binary image rather than calculating a threshold.
if isempty(strmatch(ThresholdMethod,{'Otsu','MoG','Background','RobustBackground','RidlerCalvard','All','Set'},'exact'))
    %if ~(strncmp(Threshold,'Otsu',4) || strncmp(Threshold,'MoG',3) || strfind(Threshold,'Background') ||strncmp(Threshold,'RidlerCalvard',13) || strcmp(Threshold,'All') || strcmp(Threshold,'Set interactively'))
    if isnan(str2double(Threshold))
        GetThreshold = 0;
        BinaryInputImage = CPretrieveimage(handles,Threshold,ModuleName,'MustBeGray','CheckScale');
    else
        GetThreshold = 1;
    end
else
    GetThreshold = 1;
end
%%% I don't think this should be here, Anne 12-4-06
% GetThreshold = 1;

%%% Checks that the Min and Max diameter parameters have valid values
index = strfind(SizeRange,',');
if isempty(index),
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max size entry is invalid.'])
end
MinDiameter = SizeRange(1:index-1);
MaxDiameter = SizeRange(index+1:end);

MinDiameter = str2double(MinDiameter);
if isnan(MinDiameter) | MinDiameter < 0 %#ok Ignore MLint
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min diameter entry is invalid.'])
end
if strcmpi(MaxDiameter,'Inf')
    MaxDiameter = Inf;
else
    MaxDiameter = str2double(MaxDiameter);
    if isnan(MaxDiameter) | MaxDiameter < 0 %#ok Ignore MLint
        error(['Image processing was canceled in the ', ModuleName, ' module because the Max diameter entry is invalid.'])
    end
end
if MinDiameter > MaxDiameter
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min Diameter is larger than the Max Diameter.'])
end

%%% Checks that the Min and Max threshold bounds have valid values
index = strfind(ThresholdRange,',');
if isempty(index)
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max threshold bounds are invalid.'])
end
MinimumThreshold = ThresholdRange(1:index-1);
MaximumThreshold = ThresholdRange(index+1:end);

%%% Check the smoothing filter size parameter
if ~strcmpi(SizeOfSmoothingFilter,'Automatic')
    SizeOfSmoothingFilter = str2double(SizeOfSmoothingFilter);
    if isnan(SizeOfSmoothingFilter) | isempty(SizeOfSmoothingFilter) | SizeOfSmoothingFilter < 0 %| SizeOfSmoothingFilter > min(size(OrigImage)) %#ok Ignore MLint
        %%% I commented out the part where we check that the size of
        %%% smoothing filter is greater than the image, because I think it
        %%% does not yield errors when that is the case, and in some
        %%% dynamic pipelines, the size of the image may be very small in
        %%% some cases and large in others.
        error(['Image processing was canceled in the ', ModuleName, ' module because the specified size of the smoothing filter is not valid or unreasonable.'])
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

if GetThreshold
    [handles,OrigThreshold] = CPthreshold(handles,Threshold,pObject,MinimumThreshold,MaximumThreshold,ThresholdCorrection,OrigImage,ImageName,ModuleName,ObjectName);
else
    OrigThreshold = 0;
end
Threshold = OrigThreshold;


    
%%% Apply a slight smoothing before thresholding to remove
%%% 1-pixel objects and to smooth the edges of the objects.
%%% Note that this smoothing is hard-coded, and not controlled
%%% by the user, but it is omitted if the user selected 0 for
%%% the size of the smoothing filter.
if SizeOfSmoothingFilter == 0
    %%% No blurring is done.
    BlurredImage = OrigImage;
else
    sigma = 1;
    FiltLength = 2*sigma;                                              % Determine filter size, min 3 pixels, max 61
    [x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);   % Filter kernel grid
    f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                 % Gaussian filter kernel
    %                BlurredImage = conv2(OrigImage,f,'same');                             % Blur original image
    %%% This adjustment prevents the outer borders of the image from being
    %%% darker (due to padding with zeros), which causes some objects on the
    %%% edge of the image to not  be identified all the way to the edge of the
    %%% image and therefore not be thrown out properly.
    BlurredImage = conv2(OrigImage,f,'same') ./ conv2(ones(size(OrigImage)),f,'same');
end
if GetThreshold
    Objects = BlurredImage > Threshold;                                   % Threshold image
else
    Objects = BinaryInputImage;
end
fieldname = ['CropMask', ImageName];
if isfield(handles.Pipeline,fieldname)
    %%% Retrieves previously selected cropping mask from handles
    %%% structure.
    BinaryCropImage = handles.Pipeline.(fieldname);
    try Objects = Objects & BinaryCropImage;
    catch
        error('The image in which you want to identify objects has been cropped, but there was a problem recognizing the cropping pattern.');
    end
end
Threshold = mean(Threshold(:));                                       % Use average threshold downstreams
if strcmp(FillHolesOption,'Yes')
    Objects = imfill(double(Objects),'holes');                            % Fill holes
end
drawnow
    
%%% Label the objects
Objects = bwlabel(Objects);

%%% Label the objects
Objects = bwlabel(Objects);

%%% Merge small objects
if strcmp(MergeChoice,'Yes')
    NumberOfObjectsBeforeMerge = max(Objects(:));
    Objects = MergeObjects(Objects,OrigImage,[MinDiameter MaxDiameter]);
    NumberOfObjectsAfterMerge = max(Objects(:));
    NumberOfMergedObjects = NumberOfObjectsBeforeMerge-NumberOfObjectsAfterMerge;
end

%%% Will be stored to the handles structure
UneditedLabelMatrixImage = Objects;

%%% Get diameters of objects and calculate the interval
%%% that contains 90% of the objects
tmp = regionprops(Objects,'EquivDiameter');
Diameters = [0;cat(1,tmp.EquivDiameter)];
SortedDiameters = sort(Diameters);
NbrInTails = max(round(0.05*length(Diameters)),1);
Lower90Limit = SortedDiameters(NbrInTails);
Upper90Limit = SortedDiameters(end-NbrInTails+1);

%%% Locate objects with diameter outside the specified range
tmp = Objects;
if strcmp(ExcludeSize,'Yes')
    %%% Create image with object intensity equal to the diameter
    DiameterMap = Diameters(Objects+1);
    %%% Remove objects that are too small
    Objects(DiameterMap < MinDiameter) = 0;
    %%% Will be stored to the handles structure
    SmallRemovedLabelMatrixImage = Objects;
    %%% Remove objects that are too big
    Objects(DiameterMap > MaxDiameter) = 0;
else
    %%% Will be stored to the handles structure even if it's unedited.
    SmallRemovedLabelMatrixImage = Objects;
end

%%% Store objects that fall outside diameter range for display
DiameterExcludedObjects = tmp - Objects;

%%% Relabel the objects
[Objects,NumOfObjects] = bwlabel(Objects > 0);
FinalLabelMatrixImage = Objects;

     
%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Indicate objects in original image and color excluded objects in red
tmp = OrigImage/max(OrigImage(:));
OutlinedObjectsR = tmp;
OutlinedObjectsG = tmp;
OutlinedObjectsB = tmp;
PerimObjects = bwperim(Objects > 0);
PerimDiameter = bwperim(DiameterExcludedObjects > 0);
OutlinedObjectsR(PerimObjects) = 0; OutlinedObjectsG(PerimObjects) = 1; OutlinedObjectsB(PerimObjects) = 0;
OutlinedObjectsR(PerimDiameter) = 1; OutlinedObjectsG(PerimDiameter)   = 0; OutlinedObjectsB(PerimDiameter)   = 0;

FinalOutline = false(size(OrigImage,1),size(OrigImage,2));
FinalOutline(PerimObjects) = 1;
FinalOutline(PerimDiameter) = 0;

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber)
    end
    subplot(2,2,1)
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    hx = subplot(2,2,2);
    im = CPlabel2rgb(handles,Objects);
    CPimagesc(im,handles);
    title(['Identified ',ObjectName]);
    hy = subplot(2,2,3);
    OutlinedObjects = cat(3,OutlinedObjectsR,OutlinedObjectsG,OutlinedObjectsB);
    CPimagesc(OutlinedObjects,handles);
    title(['Outlined ', ObjectName]);
    
    %%% Report numbers
    posx = get(hx,'Position');
    posy = get(hy,'Position');
    bgcolor = get(ThisModuleFigureNumber,'Color');
    uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[posx(1)-0.05 posy(2)+posy(4)-0.04 posx(3)+0.1 0.04],...
        'BackgroundColor',bgcolor,'HorizontalAlignment','Left','String',sprintf('Threshold:  %0.3f',Threshold),'FontSize',handles.Preferences.FontSize);
    uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[posx(1)-0.05 posy(2)+posy(4)-0.08 posx(3)+0.1 0.04],...
        'BackgroundColor',bgcolor,'HorizontalAlignment','Left','String',sprintf('Number of identified objects: %d',NumOfObjects),'FontSize',handles.Preferences.FontSize);
    uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[posx(1)-0.05 posy(2)+posy(4)-0.16 posx(3)+0.1 0.08],...
        'BackgroundColor',bgcolor,'HorizontalAlignment','Left','String',sprintf('90%% of objects within diameter range [%0.1f, %0.1f] pixels',...
        Lower90Limit,Upper90Limit),'FontSize',handles.Preferences.FontSize);
    ObjectCoverage = 100*sum(sum(Objects > 0))/numel(Objects);
    uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[posx(1)-0.05 posy(2)+posy(4)-0.20 posx(3)+0.1 0.04],...
        'BackgroundColor',bgcolor,'HorizontalAlignment','Left','String',sprintf('%0.1f%% of image consists of objects',ObjectCoverage),'FontSize',handles.Preferences.FontSize);
    if strcmp(MergeChoice,'Yes')
        uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[posx(1)-0.05 posy(2)+posy(4)-0.32 posx(3)+0.1 0.04],...
            'BackgroundColor',bgcolor,'HorizontalAlignment','Left','String',sprintf('Number of Merged Objects:  %d',NumberOfMergedObjects),'FontSize',handles.Preferences.FontSize);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the segmented image, not edited for objects along the edges or
%%% for size, to the handles structure.
fieldname = ['UneditedSegmented',ObjectName];
handles.Pipeline.(fieldname) = UneditedLabelMatrixImage;

%%% Saves the segmented image, only edited for small objects, to the
%%% handles structure.
fieldname = ['SmallRemovedSegmented',ObjectName];
handles.Pipeline.(fieldname) = SmallRemovedLabelMatrixImage;

%%% Saves the final segmented label matrix image to the handles structure.
fieldname = ['Segmented',ObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

%%% Saves images to the handles structure so they can be saved to the hard
%%% drive, if the user requested.
if ~strcmpi(SaveOutlines,'Do not save')
    try    handles.Pipeline.(SaveOutlines) = FinalOutline;
    catch error(['The object outlines were not calculated by the ', ModuleName, ' module, so these images were not saved to the handles structure. The Save Images module will therefore not function on these images. This is just for your information - image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
    end
end

if strcmp(MergeChoice,'Yes')
    %%% Saves the NumberOfMergedObjects to the handles structure.
    %%% See comments for the Threshold saving below
    if ~isfield(handles.Measurements.Image,'NumberOfMergedObjectsFeatures')
        handles.Measurements.Image.NumberOfMergedObjectsFeatures = {};
        handles.Measurements.Image.NumberOfMergedObjects = {};
    end
    column = find(~cellfun('isempty',strfind(handles.Measurements.Image.NumberOfMergedObjectsFeatures,ObjectName)));
    if isempty(column)
        handles.Measurements.Image.NumberOfMergedObjectsFeatures(end+1) = {ObjectName};
        column = length(handles.Measurements.Image.NumberOfMergedObjectsFeatures);
    end
    handles.Measurements.Image.NumberOfMergedObjects{handles.Current.SetBeingAnalyzed}(1,column) = NumberOfMergedObjects;
end

%%% Saves the Threshold value to the handles structure. Storing
%%% the threshold is a little more complicated than storing
%%% other measurements because several different modules will
%%% write to the handles.Measurements.Image.Threshold
%%% structure, and we should therefore probably append the
%%% current threshold to an existing structure.
% First, if the Threshold fields don't exist, initialize them
if ~isfield(handles.Measurements.Image,'ThresholdFeatures')
    handles.Measurements.Image.ThresholdFeatures = {};
    handles.Measurements.Image.Threshold = {};
end
% Search the ThresholdFeatures to find the column for this object type
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ThresholdFeatures,ObjectName)));
% If column is empty it means that this particular object has
% not been segmented before. This will typically happen for the
% first cycle. Append the feature name in the
% handles.Measurements.Image.ThresholdFeatures matrix
if isempty(column)
    handles.Measurements.Image.ThresholdFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ThresholdFeatures);
end
handles.Measurements.Image.Threshold{handles.Current.SetBeingAnalyzed}(1,column) = Threshold;

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
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(FinalLabelMatrixImage(:));


%%% Saves mask of identified regions
RegionMaskImage = im2bw(FinalLabelMatrixImage,0);

fieldname = RegionMaskImageName;
handles.Pipeline.(fieldname) = RegionMaskImage;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%[BS, 2006-04-19] Add Upper90DiameterLimit to stored output,
%%% might be  an indication for out-of-focus images.
%%% Saves some more Object Analysis stuff.
%%% See comments for the Threshold saving above

%%% BS HACK TO GET NUMBERS OF EXCLUDED OBJECTS
[tmpXXX,tmpExcludedDiameterObjects] = bwlabel(DiameterExcludedObjects > 0);
clear tmpXXX;

ObjectCoverage = 100*sum(sum(Objects > 0))/numel(Objects);

if ~isfield(handles.Measurements.Image,'ObjectDiameterFeatures')
    handles.Measurements.Image.ObjectDiameterFeatures = {};
    handles.Measurements.Image.ObjectDiameter = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectDiameterFeatures,ObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectDiameterFeatures(end+1) = {['Upper90Limit ' ObjectName]};
    column = length(handles.Measurements.Image.ObjectDiameterFeatures);
    handles.Measurements.Image.ObjectDiameterFeatures(end+1) = {['MaxDiameter ' ObjectName]};
    handles.Measurements.Image.ObjectDiameterFeatures(end+1) = {['ObjectCoverage ' ObjectName]};
    handles.Measurements.Image.ObjectDiameterFeatures(end+1) = {['DiameterExcludedObjects ' ObjectName]};
end
handles.Measurements.Image.ObjectDiameter{handles.Current.SetBeingAnalyzed}(1,column) = Upper90Limit;
handles.Measurements.Image.ObjectDiameter{handles.Current.SetBeingAnalyzed}(1,column+1) = MaxDiameter;
handles.Measurements.Image.ObjectDiameter{handles.Current.SetBeingAnalyzed}(1,column+2) = ObjectCoverage;
handles.Measurements.Image.ObjectDiameter{handles.Current.SetBeingAnalyzed}(1,column+3) = tmpExcludedDiameterObjects;

%%%[BS, END OF HACK]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Saves the location of each segmented object
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(FinalLabelMatrixImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};



%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

function Objects = MergeObjects(Objects,OrigImage,Diameters)

%%% Find the object that we should try to merge with other objects. The object
%%% numbers of these objects are stored in the variable 'MergeIndex'. The objects
%%% that we will try to merge are either the ones that fall below the specified
%%% MinDiameter threshold, or relatively small objects that are above the MaxEccentricity
%%% threshold. These latter objects are likely to be cells where two maxima have been
%%% found and the watershed transform has divided cells into two parts.
MinDiameter = Diameters(1);
MaxDiameter = Diameters(2);
MaxEccentricity = 0.75;      % Empirically determined value
props = regionprops(Objects,'EquivDiameter','PixelIdxList','Eccentricity');   % Get diameters of the objects
EquivDiameters = cat(1,props.EquivDiameter);
Eccentricities = cat(1,props.Eccentricity);
IndexEccentricity = intersect(find(Eccentricities > MaxEccentricity),find(EquivDiameters < (MinDiameter + (MaxDiameter - MinDiameter)/4)));
IndexDiameter = find(EquivDiameters < MinDiameter);
MergeIndex = unique([IndexDiameter;IndexEccentricity]);

% Try to merge until there are no objects left in the 'MergeIndex' list.
[sr,sc] = size(OrigImage);
while ~isempty(MergeIndex)

    % Get next object to merge
    CurrentObjectNbr = MergeIndex(1);

    %%% Identify neighbors and put their label numbers in a list 'NeighborsNbr'
    %%% Cut a patch so we don't have to work with the entire image
    [r,c] = ind2sub([sr sc],props(CurrentObjectNbr).PixelIdxList);
    rmax = min(sr,max(r) + 3);
    rmin = max(1,min(r) - 3);
    cmax = min(sc,max(c) + 3);
    cmin = max(1,min(c) - 3);
    ObjectsPatch = Objects(rmin:rmax,cmin:cmax);
    BinaryPatch = double(ObjectsPatch == CurrentObjectNbr);
    GrownBinaryPatch = conv2(BinaryPatch,double(getnhood(strel('disk',2))),'same') > 0;
    Neighbors = ObjectsPatch .*GrownBinaryPatch;
    NeighborsNbr = setdiff(unique(Neighbors(:)),[0 CurrentObjectNbr]);


    %%% For each neighbor, calculate a set of criteria based on which we decide if to merge.
    %%% Currently, two criteria are used. The first is a Likelihood ratio that indicates whether
    %%% the interface pixels between the object to merge and its neighbor belong to a background
    %%% class or to an object class. The background class and object class are modeled as Gaussian
    %%% distributions with mean and variance estimated from the image. The Likelihood ratio determines
    %%% to which of the distributions the interface voxels most likely belong to. The second criterion
    %%% is the eccentrity of the object resulting from a merge. The more circular, i.e., the lower the
    %%% eccentricity, the better.
    LikelihoodRatio    = zeros(length(NeighborsNbr),1);
    MergedEccentricity = zeros(length(NeighborsNbr),1);
    for j = 1:length(NeighborsNbr)

        %%% Get Neigbor number
        CurrentNeighborNbr = NeighborsNbr(j);

        %%% Cut patch which contains both original object and the current neighbor
        [r,c] = ind2sub([sr sc],[props(CurrentObjectNbr).PixelIdxList;props(CurrentNeighborNbr).PixelIdxList]);
        rmax = min(sr,max(r) + 3);
        rmin = max(1,min(r) - 3);
        cmax = min(sc,max(c) + 3);
        cmin = max(1,min(c) - 3);
        ObjectsPatch = Objects(rmin:rmax,cmin:cmax);
        OrigImagePatch = OrigImage(rmin:rmax,cmin:cmax);

        %%% Identify object interiors, background and interface voxels
        BinaryNeighborPatch      = double(ObjectsPatch == CurrentNeighborNbr);
        BinaryObjectPatch        = double(ObjectsPatch == CurrentObjectNbr);
        GrownBinaryNeighborPatch = conv2(BinaryNeighborPatch,ones(3),'same') > 0;
        GrownBinaryObjectPatch   = conv2(BinaryObjectPatch,ones(3),'same') > 0;
        Interface                = GrownBinaryNeighborPatch.*GrownBinaryObjectPatch;
        Background               = ((GrownBinaryNeighborPatch + GrownBinaryObjectPatch) > 0) - BinaryNeighborPatch - BinaryObjectPatch - Interface;
        WithinObjectIndex        = find(BinaryNeighborPatch + BinaryObjectPatch);
        InterfaceIndex           = find(Interface);
        BackgroundIndex          = find(Background);

        %%% Calculate likelihood of the interface belonging to the background or to an object.
        WithinObjectClassMean   = mean(OrigImagePatch(WithinObjectIndex));
        WithinObjectClassStd    = std(OrigImagePatch(WithinObjectIndex)) + sqrt(eps);
        BackgroundClassMean     = mean(OrigImagePatch(BackgroundIndex));
        BackgroundClassStd      = std(OrigImagePatch(BackgroundIndex)) + sqrt(eps);
        InterfaceMean           = mean(OrigImagePatch(InterfaceIndex)); %#ok Ignore MLint
        LogLikelihoodObject     = -log(WithinObjectClassStd^2) - (InterfaceMean - WithinObjectClassMean)^2/(2*WithinObjectClassStd^2);
        LogLikelihoodBackground = -log(BackgroundClassStd^2) - (InterfaceMean - BackgroundClassMean)^2/(2*BackgroundClassStd^2);
        LikelihoodRatio(j)      =  LogLikelihoodObject - LogLikelihoodBackground;

        %%% Calculate the eccentrity of the object obtained if we merge the current object
        %%% with the current neighbor.
        MergedObject =  double((BinaryNeighborPatch + BinaryObjectPatch + Interface) > 0);
        tmp = regionprops(MergedObject,'Eccentricity');
        MergedEccentricity(j) = tmp(1).Eccentricity;

        %%% Get indexes for the interface pixels in original image.
        %%% These indexes are required if we need to merge the object with
        %%% the current neighbor.
        tmp = zeros(size(OrigImage));
        tmp(rmin:rmax,cmin:cmax) = Interface;
        tmp = regionprops(double(tmp),'PixelIdxList');
        OrigInterfaceIndex{j} = cat(1,tmp.PixelIdxList); %#ok Ignore MLint
    end

    %%% Let each feature rank which neighbor to merge with. Then calculate
    %%% a score for each neighbor. If the neighbors is ranked 1st, it will get
    %%% 1 point; 2nd, it will get 2 points; and so on. The lower score the better.
    [ignore,LikelihoodRank]   = sort(LikelihoodRatio,'descend'); %#ok Ignore MLint % The higher the LikelihoodRatio the better
    [ignore,EccentricityRank] = sort(MergedEccentricity,'ascend'); %#ok Ignore MLint % The lower the eccentricity the better
    NeighborScore = zeros(length(NeighborsNbr),1);
    for j = 1:length(NeighborsNbr)
        NeighborScore(j) = find(LikelihoodRank == j) +  find(EccentricityRank == j);
    end

    %%% Go through the neighbors, starting with the highest ranked, and merge
    %%% with the first neighbor for which certain basic criteria are fulfilled.
    %%% If no neighbor fulfil the basic criteria, there will be no merge.
    [ignore,TotalRank] = sort(NeighborScore); %#ok Ignore MLint
    for j = 1:length(NeighborsNbr)
        CurrentNeighborNbr = NeighborsNbr(TotalRank(j));

        %%% To merge, the interface between objects must be more likely to belong to the object class
        %%% than the background class. The eccentricity of the merged object must also be lower than
        %%% for the original object.
        if LikelihoodRatio(TotalRank(j)) > 0 && MergedEccentricity(TotalRank(j)) < Eccentricities(CurrentObjectNbr)

            %%% OK, let's merge!
            %%% Assign the neighbor number to the current object
            Objects(props(CurrentObjectNbr).PixelIdxList) = CurrentNeighborNbr;

            %%% Assign the neighbor number to the interface pixels between the current object and the neigbor
            Objects(OrigInterfaceIndex{TotalRank(j)}) = CurrentNeighborNbr;

            %%% Add the pixel indexes to the neigbor index list
            props(CurrentNeighborNbr).PixelIdxList = cat(1,...
                props(CurrentNeighborNbr).PixelIdxList,...
                props(CurrentObjectNbr).PixelIdxList,...
                OrigInterfaceIndex{TotalRank(j)});

            %%% Remove the neighbor from the list of objects to be merged (if it's there).
            MergeIndex = setdiff(MergeIndex,CurrentNeighborNbr);
        end
    end

    %%% OK, we are done with the current object, let's go to the next
    MergeIndex = MergeIndex(2:end-1);
end

%%% Finally, relabel the objects
Objects = bwlabel(Objects > 0);