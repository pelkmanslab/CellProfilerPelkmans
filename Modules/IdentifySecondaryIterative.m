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
% cellular periphery in sparese regions).
%
% To prevent false positives, if very few primary objects are present,
% limits for threshold can be used. Reference threshold boundaries are
% displayed when using the module locally (note that more digits are shown
% than in standard module to allow finer grained settings for cv7k images
% (where dim pixels might already be separateble from background)
%
% In addition to standard IdentifySeondary module , this module is not
% affected by the presence of small false positive nuclei, that were
% filtered out by the DiscardSinglePixel... module (which is frequent in
% pipelines with confocal images)
%
% HOW TO USE IT
% a) change IdentifySecondary module 
% go to the code of the intial IdentifySecondary module and look for
% 'Threshold:  %0.3f               %0.1f%% of image consists of objects'
% change it to
% 'Threshold:  %0.5f               %0.1f%% of image consists of objects'
% This way the graphical user interface will display a very precise value
% of the chosen threshold.
%
% b) Get lowest threshold, which does not yet recognize background as
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
% c) Now make a new CP pipeline with the IdentifySecondaryIterative module.
% THRESHOLD CORRECTION FACTORS. IN DESCENDING RANKING. Should indicate many
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
% d) Select LOWER AND UPPER BONDS ON THRESHOLD. 
% These values correspond to the minimal and maximal values that a
% threshold is allowed to have. They have the format
% SmallestThreshold,HighestThreshold 
% For SmallestThreshold you should use the value obtained in c). Leaving
% the maximal value at 1 has worked fine for us all the time. Setting a
% minimal value will prevent recognition of the background as an object
%
% e) Go to a vending machine, get some Swiss chocolate and have a 
% small break. Later check, if segmentation has worked.
%
% f) If everything is fine, run on parallel computing.
% Consider to plan around 1min per tested threshold. Since this module will
% be slow it usually determines the amount of image cycles, which can be
% processed in one batch. Batch sizes should be small enough to allow jobs
% to finish within the time that your parallel computing facility grants to 
% individual jobs. On the other hand, you can seriously slow down
% file systems (such as lustre) if the amount of output files 
% exceeds 100.000 files and you would have parallel jobs requesting the 
% listing of this directory.
%
% THRESHOLD SELECTION METHOD
% NOTE THAT this module was only tested with the Otsu Global method. Other
% options have been disabled, but might be reactivated easily by adjusting
% the code.
%
% *************************************************************************
%
%
% How it works>
% This is a heavily adjusted version of the orignal IdentifySecondary module
% It also has lower memory requiremnts.
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
% [TS]
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

%textVAR04 = Select an automatic thresholding method or enter an absolute threshold in the range [0,1]. To choose a binary image, select "Other" and type its name.  Choosing 'All' will use the Otsu Global method to calculate a single threshold for the entire image group. The other methods calculate a threshold for each image individually. Set interactively will allow you to manually adjust the threshold during the first cycle to determine what will work well.
%choiceVAR04 = Otsu Global
iThreshold = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu custom

%textVAR05 = Threshold correction factors. In descending ranking. eg 0.9 0.8 0.7
%defaultVAR05 = 0.9 0.7 0.6 0.58 0.55 0.50 0.45 0.4 0.35 0.3 0.25
iThresholdCorrection = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Lower and upper bounds on threshold, in the range [0,1].
%defaultVAR06 = 0.00300,1
ThresholdRange = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = For MoG thresholding, what is the approximate percentage of image covered by objects?
%choiceVAR07 = 10%
%choiceVAR07 = 20%
%choiceVAR07 = 30%
%choiceVAR07 = 40%
%choiceVAR07 = 50%
%choiceVAR07 = 60%
%choiceVAR07 = 70%
%choiceVAR07 = 80%
%choiceVAR07 = 90%
pObject = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu


%%%VariableRevisionNumber = 3

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

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects which will be used to weed out which objects are
%%% real - not on the edges and not below or above the specified size
%%% limits. Checks first to see whether the appropriate image exists.
EditedPrimaryLabelMatrixImage = CPretrieveimage(handles,['Segmented', PrimaryObjectName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));

%%% Chooses the first word of the method name (removing 'Global' or 'Adaptive').
ThresholdMethod = strtok(iThreshold);
%%% Checks if a custom entry was selected for Threshold, which means we are using an incoming binary image rather than calculating a threshold.
if isempty(strmatch(ThresholdMethod,{'Otsu','MoG','Background','RobustBackground','RidlerCalvard','All','Set'},'exact'))
    %if ~(strncmp(Threshold,'Otsu',4) || strncmp(Threshold,'MoG',3) || strfind(Threshold,'Background') ||strncmp(Threshold,'RidlerCalvard',13) || strcmp(Threshold,'All') || strcmp(Threshold,'Set interactively'))
    if isnan(str2double(iThreshold))
        GetThreshold = 0;
        BinaryInputImage = CPretrieveimage(handles,iThreshold,ModuleName,'MustBeGray','CheckScale');
    else
        GetThreshold = 1;
    end
else
    GetThreshold = 1;
end

%%% Checks that the Min and Max threshold bounds have valid values
index = strfind(ThresholdRange,',');
if isempty(index)
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max threshold bounds are invalid.'])
end

MinimumThreshold = str2double(ThresholdRange(1:index-1));
MaximumThreshold = str2double(ThresholdRange(index+1:end));

% Create vector containing the thresholds that should be tested
[isSafe iThresholdCorrection]= inputVectorsForEvalCP3D(iThresholdCorrection,false);
if isSafe == false
    error(['Image processing was canceled in the ', ModuleName, ' module because input of threshold contained forbidden characters'])
else
    ThresholdCorrection = eval(iThresholdCorrection);
end



%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow
numThresholdsToTest = length(ThresholdCorrection);
ThresholdArray = cell(numThresholdsToTest,1); % [modified by TS to include multiple thrsholds]

%obtain first threshold via CPthreshold. note that this will also create a
%threshold measuremnt, this measuremnt will not be produced for sequential
%thresholds to prevent over/writing and conflicts that arise if
%numThreshold>=3
if GetThreshold
    % [TS] force to use minimal treshold value of 0 and maximum of 1, to ensure
    % equal thresholding for all tested tresholds
    [handles,ThresholdArray{1}] = CPthreshold(handles,iThreshold,pObject,'0','1',ThresholdCorrection(1),OrigImage,ImageName,ModuleName,SecondaryObjectName);
else
    ThresholdArray{1} = 0; % should never be used
end

%%%% [TS] start modification for obtaining multiple thresholds  %%%%%%%%%%%%%%%%%
if numThresholdsToTest>1
    for k=2:numThresholdsToTest
        %%% STEP 1a: Marks at least some of the background
        if GetThreshold
            refThreshold = ThresholdArray{1};
            ThresholdArray{k} = refThreshold .* ThresholdCorrection(k) ./ThresholdCorrection(1);
        else
            ThresholdArray{k} = 0; % should never be used
        end
    end
end

% now fix thresholds outside of range. Could be made nicer by direcly
% calling a function for fixing thresholds for CP standard case (k=1) and
% [TS] modification for k>=2
for k=1:numThresholdsToTest
    % note that CP adresses the threshold in such a way that it could be
    % either a number or a matrix.-> the internally generated threshold
    % might be either of it. The following lines should support both.
    reconstituteThresholdImage = ThresholdArray{k};
    bnSomethingOutsidRange = false;
    
    f = reconstituteThresholdImage(:) < MinimumThreshold;
    if any(f)
        reconstituteThresholdImage(f) = MinimumThreshold;
        bnSomethingOutsidRange = true;
    end
    
    f = reconstituteThresholdImage(:) > MaximumThreshold;
    if any(f)
        reconstituteThresholdImage(f) = MaximumThreshold;
        bnSomethingOutsidRange = true;
    end
    
    if bnSomethingOutsidRange == true
        ThresholdArray{k} = reconstituteThresholdImage;
    end
end

%%%% [TS] end modification for obtaining multiple thresholds  %%%%%%%%%%%%%%%%%


%%%% [TS] Start modification> DISMISS only border  %%%%%%%%%%%%%%%%%%%%%%
%%% Preliminary objects, which were not identified as object proper, still
%%% serve as seeds for allocating pixels to secondary object. While this
%%% makes sense for nuclei, which were discared in the primary module due to
%%% their location at the image border (and have a surrounding cytoplasm),
%%% it can lead to wrong segmenations, if a false positive nucleus, that was
%%% filtered away , eg. by the DiscardSinglePixel... module , was present

%%% corrsponds to one line from STEP 10, moved up. Allows proper
%%% initialzing for reconstitution
%%% Converts the EditedPrimaryBinaryImage to binary.
EditedPrimaryBinaryImage = im2bw(EditedPrimaryLabelMatrixImage,.5);

% Replace the way the mask PrelimPrimaryBinaryImage is generated
%%% Use a shared line from STEP 0. This will allow proper initializing for reconstitution.
%%% Converts the PrelimPrimaryLabelMatrixImage to binary.
%%% OLD> PrelimPrimaryBinaryImage = im2bw(PrelimPrimaryLabelMatrixImage,.5);

%%% Get IDs of objects at image border
R= PrelimPrimaryLabelMatrixImage([1 end],:);
C= PrelimPrimaryLabelMatrixImage(:,[1 end]);
BoderObjIDs = unique([R C']);
while any(BoderObjIDs==0)
    BoderObjIDs = BoderObjIDs(2:end);
end
clear R; clear C;

PrelimPrimaryBinaryImage = false(size(EditedPrimaryBinaryImage));

f =     ismember(PrelimPrimaryLabelMatrixImage,BoderObjIDs) | ... % objects at border
    EditedPrimaryBinaryImage;            % proper objects

PrelimPrimaryBinaryImage(f) = true;


%%%% [TS] End modification> DISMISS only border  %%%%%%%%%%%%%%%%%%%%%




%%%% [TS] %%%%%%%%%%%%% Start of SHARED code for precalculations %%%%%%%%%%%
% note that fragments of original function were replaced by TS to prevent
% redundant calculations

drawnow

%%% Creates the structuring element that will be used for dilation.
StructuringElement = strel('square',3);
%%% Dilates the Primary Binary Image by one pixel (8 neighborhood).
DilatedPrimaryBinaryImage = imdilate(PrelimPrimaryBinaryImage, StructuringElement);
%%% Subtracts the PrelimPrimaryBinaryImage from the DilatedPrimaryBinaryImage,
%%% which leaves the PrimaryObjectOutlines.
PrimaryObjectOutlines = DilatedPrimaryBinaryImage - PrelimPrimaryBinaryImage;


%%% STEP 4: Calculate the Sobel image, which reflects gradients, which will
%%% be used for the watershedding function.
drawnow
%%% Calculates the 2 sobel filters.  The sobel filter is directional, so it
%%% is used in both the horizontal & vertical directions and then the
%%% results are combined.
filter1 = fspecial('sobel');
filter2 = filter1';
%%% Applies each of the sobel filters to the original image.
I1 = imfilter(OrigImage, filter1);
I2 = imfilter(OrigImage, filter2);
%%% Adds the two images.
%%% The Sobel operator results in negative values, so the absolute values
%%% are calculated to prevent errors in future steps.
AbsSobeledImage = abs(I1) + abs(I2);
clear I1; clear I2;                  %%% [NB] hack. save memory

%%%% [TS] %%%%%%%%%%%%% End of SHARED code for precalculations %%%%%%%%%%%


%%%%%% [TS] %%%%%%%%%%%%%%%  ITERATION CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intialize output

cellFinalLabelMatrixImage = cell(numThresholdsToTest,1);


for k=1:numThresholdsToTest
    
    % STEP 0
    %%% Thresholds the original image.
    if GetThreshold
        ThresholdedOrigImage = OrigImage > ThresholdArray{k};
    else
        ThresholdedOrigImage = logical(BinaryInputImage);
    end
    
    
    %%% STEP 1b: Marks at least some of the background
    
    %%% Inverts the image.
    InvertedThresholdedOrigImage = imcomplement(ThresholdedOrigImage);
    clear ThresholdedOrigImage;             %%% [NB] hack. save memory.
    
    %%% STEP 3: Produce the marker image which will be used for the first
    %%% watershed.
    drawnow
    %%% Combines the foreground markers and the background markers.
    BinaryMarkerImagePre = PrelimPrimaryBinaryImage | InvertedThresholdedOrigImage;
    %%% Overlays the PrimaryObjectOutlines to maintain distinctions between each
    %%% primary object and the background.
    BinaryMarkerImage = BinaryMarkerImagePre;
    clear BinaryMarkerImagePre;             %%% [NB] hack. save memory.
    BinaryMarkerImage(PrimaryObjectOutlines == 1) = 0;
    
    
    %%% STEP 5: Perform the first watershed.
    drawnow
    
    %%% Overlays the foreground and background markers
    Overlaid = imimposemin(AbsSobeledImage, BinaryMarkerImage);
    clear BinaryMarkerImage;  %%% [NB] hack. save memory.
    
    %%% Perform the watershed on the marked absolute-value Sobel Image.
    BlackWatershedLinesPre = watershed(Overlaid);
    clear Overlaid;                 %%% [NB] hack. save memory.
    
    %%% Bug workaround (see step 9).
    %%% [NB, WATERSHED BUG IN VERSION 2011A OR HIGHER HAS BEEN FIXED. SO CHECK VERSION FIRST]
    if verLessThan('matlab', '7.12.0') && ~isunix
        BlackWatershedLinesPre2 = im2bw(BlackWatershedLinesPre,.5);
        BlackWatershedLines = bwlabel(BlackWatershedLinesPre2);
        %%% [NB] hack. save memory.
        clear BlackWatershedLinesPre2 BlackWatershedLinesPre;
    else
        %%% [BS, QUICK AND DIRTY HACK FROM PEKLMANS]
        BlackWatershedLines = double(BlackWatershedLinesPre);
        %%% [NB] hack. save memory.
        clear BlackWatershedLinesPre;
        %%% END OF BS-HACK BUGFIX FOR VERSION 2011 AND LATER?
    end
    
    %%% STEP 6: Identify and extract the secondary objects, using the watershed
    %%% lines.
    drawnow
    %%% The BlackWatershedLines image is a label matrix where the watershed
    %%% lines = 0 and each distinct object is assigned a number starting at 1.
    %%% This image is converted to a binary image where all the objects = 1.
    SecondaryObjects1 = im2bw(BlackWatershedLines,.5);
    %%% [NB] hack. save memory.
    clear BlackWatershedLines;
    %%% Identifies objects in the binary image using bwlabel.
    %%% Note: Matlab suggests that in some circumstances bwlabeln is faster
    %%% than bwlabel, even for 2D images.  I found that in this case it is
    %%% about 10 times slower.
    LabelMatrixImage1 = bwlabel(SecondaryObjects1,4);
    %%% [NB] hack. save memory.
    clear SecondaryObjects1;
    drawnow
    
    %%% STEP 7: Discarding background "objects".  The first watershed function
    %%% simply divides up the image into regions.  Most of these regions
    %%% correspond to actual objects, but there are big blocks of background
    %%% that are recognized as objects. These can be distinguished from actual
    %%% objects because they do not overlap a primary object.
    
    %%% The following changes all the labels in LabelMatrixImage1 to match the
    %%% centers they enclose (from PrelimPrimaryBinaryImage), and marks as background
    %%% any labeled regions that don't overlap a center. This function assumes
    %%% that every center is entirely contained in one labeled area.  The
    %%% results if otherwise may not be well-defined. The non-background labels
    %%% will be renumbered according to the center they enclose.
    
    %%% Finds the locations and labels for different regions.
    area_locations = find(LabelMatrixImage1);
    area_labels = LabelMatrixImage1(area_locations);
    %%% Creates a sparse matrix with column as label and row as location,
    %%% with the value of the center at (I,J) if location I has label J.
    %%% Taking the maximum of this matrix gives the largest valued center
    %%% overlapping a particular label.  Tacking on a zero and pushing
    %%% labels through the resulting map removes any background regions.
    map = [0 full(max(sparse(area_locations, area_labels, PrelimPrimaryBinaryImage(area_locations))))];
    
    ActualObjectsBinaryImage = map(LabelMatrixImage1 + 1);
    clear area_labels area_locations map;              %%% [NB] hack. save memory.
    
    
    %%% STEP 8: Produce the marker image which will be used for the second
    %%% watershed.
    drawnow
    %%% The module has now produced a binary image of actual secondary
    %%% objects.  The gradient (Sobel) image was used for watershedding, which
    %%% produces very nice divisions between objects that are clumped, but it
    %%% is too stringent at the edges of objects that are isolated, and at the
    %%% edges of clumps of objects. Therefore, the stringently identified
    %%% secondary objects are used as markers for a second round of
    %%% watershedding, this time based on the original (intensity) image rather
    %%% than the gradient image.
    
    %%% Creates the structuring element that will be used for dilation.
    StructuringElement = strel('square',3);
    %%% Dilates the Primary Binary Image by one pixel (8 neighborhood).
    DilatedActualObjectsBinaryImage = imdilate(ActualObjectsBinaryImage, StructuringElement);
    %%% Subtracts the PrelimPrimaryBinaryImage from the DilatedPrimaryBinaryImage,
    %%% which leaves the PrimaryObjectOutlines.
    ActualObjectOutlines = DilatedActualObjectsBinaryImage - ActualObjectsBinaryImage;
    %%% [NB] hack. save memory.
    clear DilatedActualObjectsBinaryImage;
    %%% Produces the marker image which will be used for the watershed. The
    %%% foreground markers are taken from the ActualObjectsBinaryImage; the
    %%% background markers are taken from the same image as used in the first
    %%% round of watershedding: InvertedThresholdedOrigImage.
    BinaryMarkerImagePre2 = ActualObjectsBinaryImage | InvertedThresholdedOrigImage;
    %%% [NB] hack. save memory.
    clear InvertedThresholdedOrigImage ActualObjectsBinaryImage;
    %%% Overlays the ActualObjectOutlines to maintain distinctions between each
    %%% secondary object and the background.
    BinaryMarkerImage2 = BinaryMarkerImagePre2;
    %%% [NB] hack. save memory.
    clear BinaryMarkerImagePre2;
    
    BinaryMarkerImage2(ActualObjectOutlines == 1) = 0;
    
    %%% STEP 9: Perform the second watershed.
    %%% As described above, the second watershed is performed on the original
    %%% intensity image rather than on a gradient (Sobel) image.
    drawnow
    %%% Inverts the original image.
    InvertedOrigImage = imcomplement(OrigImage);
    %%% Overlays the foreground and background markers onto the
    %%% InvertedOrigImage, so there are black secondary object markers on top
    %%% of each dark secondary object, with black background.
    MarkedInvertedOrigImage = imimposemin(InvertedOrigImage, BinaryMarkerImage2);
    %%% [NB] hack. save memory.
    clear BinaryMarkerImage2 BinaryMarkerImage2;
    
    %%% Performs the watershed on the MarkedInvertedOrigImage.
    SecondWatershedPre = watershed(MarkedInvertedOrigImage);
    %%% [NB] hack.save memory
    clear MarkedInvertedOrigImage;
    %%% BUG WORKAROUND:
    %%% There is a bug in the watershed function of Matlab that often results in
    %%% the label matrix result having two objects labeled with the same label.
    %%% I am not sure whether it is a bug in how the watershed image is
    %%% produced (it seems so: the resulting objects often are nowhere near the
    %%% regional minima) or whether it is simply a problem in the final label
    %%% matrix calculation. Matlab has been informed of this issue and has
    %%% confirmed that it is a bug (February 2004). I think that it is a
    %%% reasonable fix to convert the result of the watershed to binary and
    %%% remake the label matrix so that each label is used only once. In later
    %%% steps, inappropriate regions are weeded out anyway.
    
    %%% [NB, WATERSHED BUG IN VERSION 2011A OR HIGHER HAS BEEN FIXED. SO CHECK VERSION FIRST]
    if verLessThan('matlab', '7.12.0') && ~isunix
        SecondWatershedPre2 = im2bw(SecondWatershedPre,.5);
        SecondWatershed = bwlabel(SecondWatershedPre2);
        %%% [NB] hack.save memory
        clear SecondWatershedPre2;
    else
        %%% [BS, QUICK AND DIRTY HACK FROM PEKLMANS]
        SecondWatershed = double(SecondWatershedPre);
        %%% END OF BS-HACK BUGFIX FOR VERSION 2011 AND LATER?
    end
    %%% [NB] hack.save memory
    clear SecondWatershedPre;
    drawnow
    
    %%% STEP 10: As in step 7, remove objects that are actually background
    %%% objects.  See step 7 for description. This time, the edited primary object image is
    %%% used rather than the preliminary one, so that objects whose nuclei are
    %%% on the edge of the image and who are larger or smaller than the
    %%% specified size are discarded.
    
    %%% Finds the locations and labels for different regions.
    area_locations2 = find(SecondWatershed);
    area_labels2 = SecondWatershed(area_locations2);
    %%% Creates a sparse matrix with column as label and row as location,
    %%% with the value of the center at (I,J) if location I has label J.
    %%% Taking the maximum of this matrix gives the largest valued center
    %%% overlapping a particular label.  Tacking on a zero and pushing
    %%% labels through the resulting map removes any background regions.
    map2 = [0 full(max(sparse(area_locations2, area_labels2, EditedPrimaryBinaryImage(area_locations2))))];
    FinalBinaryImagePre = map2(SecondWatershed + 1);
    %%% [NB] hack. save memory
    clear SecondWatershed area_labels2 map2;
    
    %%% Fills holes in the FinalBinaryPre image.
    FinalBinaryImage = imfill(FinalBinaryImagePre, 'holes');
    %%% [NB] hack. save memory
    clear FinalBinaryImagePre;
    %%% Converts the image to label matrix format. Even if the above step
    %%% is excluded (filling holes), it is still necessary to do this in order
    %%% to "compact" the label matrix: this way, each number corresponds to an
    %%% object, with no numbers skipped.
    ActualObjectsLabelMatrixImage3 = bwlabel(FinalBinaryImage);
    %%% [NB] hack. save memory
    clear FinalBinaryImage;
    %%% The final objects are relabeled so that their numbers
    %%% correspond to the numbers used for nuclei.
    %%% For each object, one label and one label location is acquired and
    %%% stored.
    [LabelsUsed,LabelLocations] = unique(EditedPrimaryLabelMatrixImage);
    %%% The +1 increment accounts for the fact that there are zeros in the
    %%% image, while the LabelsUsed starts at 1.
    LabelsUsed(ActualObjectsLabelMatrixImage3(LabelLocations(2:end))+1) = EditedPrimaryLabelMatrixImage(LabelLocations(2:end));
    FinalLabelMatrixImagePre = LabelsUsed(ActualObjectsLabelMatrixImage3+1);
    %%% [NB] hack. save memory
    clear FinalBinaryImage LabelsUsed LabelLocations;
    %%% The following is a workaround for what seems to be a bug in the
    %%% watershed function: very very rarely two nuclei end up sharing one
    %%% "cell" object, so that one of the nuclei ends up without a
    %%% corresponding cell.  I am trying to determine why this happens exactly.
    %%% When the cell is measured, the area (and other
    %%% measurements) are recorded as [], which causes problems when dependent
    %%% measurements (e.g. perimeter/area) are attempted.  It results in divide
    %%% by zero errors and the mean area = NaN and so on.  So, the Primary
    %%% label matrix image (where it is nonzero) is written onto the Final cell
    %%% label matrix image pre so that every primary object has at least some
    %%% pixels of secondary object.
    FinalLabelMatrixImage = FinalLabelMatrixImagePre;
    %%% [NB] hack. save memory
    clear FinalLabelMatrixImagePre;
    FinalLabelMatrixImage(EditedPrimaryLabelMatrixImage ~= 0) = EditedPrimaryLabelMatrixImage(EditedPrimaryLabelMatrixImage ~= 0);
    
    %[TS] insert to allow easy collecition of segmentations at all
    %different thresholds
    if max(FinalLabelMatrixImage(:))<intmax('uint16')
        cellFinalLabelMatrixImage{k} = uint16(FinalLabelMatrixImage); % if used for cells, few objects, reduce memory load
    else
        cellFinalLabelMatrixImage{k} = FinalLabelMatrixImage;
    end
    
    clear FinalLabelMatrixImage; % memory==low
    
end
%%%% [TS] %%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of iteration  %%%%%%%%%%%

clear AbsSobeledImage;
clear PrelimPrimaryBinaryImage;



%%%% [TS] %%%%%%%%%% ABSOLUTE SEGEMENTATION  Start  %%%%%%%%%%%

% this code combines knowledge of about the segementation at individual
% thresholds to one common segmentation, which will be superior and
% combines the advantage of high threshold (less/no false allocation to
% wrong cell) with the advantage of low thresholds (inclusion of cell
% boundaries)


% A) Reverse projection
FinalLabelMatrixImage  = zeros(size(cellFinalLabelMatrixImage{1}),'double');
for k=numThresholdsToTest:-1:1
    f = cellFinalLabelMatrixImage{k} ~=0;
    FinalLabelMatrixImage(f) = cellFinalLabelMatrixImage{k}(f);
end

% make a second loop, which creates second best object ID
% FinalLabelMatrixImageSurrogate  = zeros(size(FinalLabelMatrixImage),'double');
% if numThresholdsToTest > 1
%     for k=numThresholdsToTest:-1:1
%         f = cellFinalLabelMatrixImage{k} ~=0 && cellFinalLabelMatrixImage{k} ~= FinalLabelMatrixImage;
%         FinalLabelMatrixImageSurrogate(f) = cellFinalLabelMatrixImage{k}(f);
%     end
% end


% B) Make sure objects are separted

% Dilate segmentation by one pixel and reassign IDs. This is necessary
% because edge detection is done in next step to create 0 intensity pixels
% between IDa-IDb. However, without dilation to background, background-IDa
% boundaries would become extended in next step

% use code from spot qualtiy control showSpotsInControl.m
DistanceToDilate = 1;
%%% Creates the structuring element using the user-specified size.
StructuringElementMini = strel('disk', DistanceToDilate);
%%% Dilates the preliminary label matrix image (edited for small only).
DilatedPrelimSecObjectLabelMatrixImageMini = imdilate(FinalLabelMatrixImage, StructuringElementMini);
%%% Converts to binary.
DilatedPrelimSecObjectBinaryImageMini = im2bw(DilatedPrelimSecObjectLabelMatrixImageMini,.5);
%%% Computes nearest neighbor image of nuclei centers so that the dividing
%%% line between secondary objects is halfway between them rather than
%%% favoring the primary object with the greater label number.
[~, Labels] = bwdist(full(FinalLabelMatrixImage>0)); % We want to ignore MLint error checking for this line.
%%% Remaps labels in Labels to labels in FinalLabelMatrixImage.
if max(Labels(:)) == 0,
    Labels = ones(size(Labels));
end
ExpandedRelabeledDilatedPrelimSecObjectImageMini = FinalLabelMatrixImage(Labels);
RelabeledDilatedPrelimSecObjectImageMini = zeros(size(ExpandedRelabeledDilatedPrelimSecObjectImageMini));
RelabeledDilatedPrelimSecObjectImageMini(DilatedPrelimSecObjectBinaryImageMini) = ExpandedRelabeledDilatedPrelimSecObjectImageMini(DilatedPrelimSecObjectBinaryImageMini);
% Stop using code from showSpotsInControl.m
clear ExpandedRelabeledDilatedPrelimSecObjectImageMini;
% Create Boundaries

I1 = imfilter(RelabeledDilatedPrelimSecObjectImageMini, filter1);   % [TS] reuse sobel filters from above
I2 = imfilter(RelabeledDilatedPrelimSecObjectImageMini, filter2);
AbsSobeledImage = abs(I1) + abs(I2);
clear I1; clear I2;                  %%% [NB] hack. save memory
edgeImage = AbsSobeledImage>0;    % detect edges
FinalLabelMatrixImage = RelabeledDilatedPrelimSecObjectImageMini .* ~edgeImage;   % set edges in Labelmatrix to zero
clear Labels; clear ExpandedRelabeledDilatedPrelimSecObjectImageMini;
clear edgeImage;

if max(FinalLabelMatrixImage(:)) ~= 0       % check if an object is present Empty Image Handling
    
    % C) Remove regions no longer connected to the primary object
    % Take code from Neighbour module
    distanceToObjectMax = 3;
    loadedImage = FinalLabelMatrixImage;
    props = regionprops(loadedImage,'BoundingBox');
    BoxPerObj = cat(1,props.BoundingBox);
    
    N = floor(BoxPerObj(:,2)-distanceToObjectMax-1);                    f = N < 1;                      N(f) = 1;
    S = ceil(BoxPerObj(:,2)+BoxPerObj(:,4)+distanceToObjectMax+1);   	f = S > size(loadedImage,1);    S(f) = size(loadedImage,1);
    W = floor(BoxPerObj(:,1)-distanceToObjectMax-1);                    f = W < 1;                      W(f) = 1;
    E = ceil(BoxPerObj(:,1)+BoxPerObj(:,3)+distanceToObjectMax+1);      f = E > size(loadedImage,2);    E(f) = size(loadedImage,2);
    
    % create empty output
    FinalLabelMatrixImage2  = zeros(size(FinalLabelMatrixImage));
    numObjects =size(BoxPerObj,1);
    if numObjects>=1  % if objects present
        for k=1: numObjects  % loop through individual objects to safe computation
            miniImage = FinalLabelMatrixImage(N(k):S(k),W(k):E(k));
            bwminiImage = miniImage>0;
            labelmini = bwlabel(bwminiImage);
            
            miniImageNuclei = EditedPrimaryLabelMatrixImage(N(k):S(k),W(k):E(k));
            bwParentOfInterest = miniImageNuclei == k;
            
            % now find the most frequent value. note that preobject will not be
            % completely within child at border of image
            
            NewChildID = labelmini(bwParentOfInterest);
            NewChildID = NewChildID(NewChildID>0);
            WithParentIX = mode(NewChildID);
            bwOutCellBody = labelmini == WithParentIX;
            
            
            % now map back the linear indices
            [r c] = find(bwOutCellBody);
            
            % get indices for final image (note that mini image might have
            % permitted regions of other cells).
            r = r-1+N(k);
            c = c-1+W(k);
            w = sub2ind(size(FinalLabelMatrixImage2),r,c);
            
            % Update Working copy of Final Segmentation image based on linear indices.
            FinalLabelMatrixImage2(w) = k;
            
        end
        
    end
    % Now mimik standard outupt of calculations of standard module
    FinalLabelMatrixImage = FinalLabelMatrixImage2;
end

% duplicate penultimate row and column. Thus pixels at border will carry 
% an object ID (and are detected by iBrain function to discard border cells); 
FinalLabelMatrixImage(:,1)= FinalLabelMatrixImage(:,2);
FinalLabelMatrixImage(:,end)= FinalLabelMatrixImage(:,(end-1));
FinalLabelMatrixImage(1,:)= FinalLabelMatrixImage(2,:);
FinalLabelMatrixImage(end,:)= FinalLabelMatrixImage((end-1),:);


%%%% [TS] %%%%%%%%%% ABSOLUTE SEGEMENTATION  End  %%%%%%%%%%%

if ~isfield(handles.Measurements,SecondaryObjectName)
    handles.Measurements.(SecondaryObjectName) = {};
end

if ~isfield(handles.Measurements,PrimaryObjectName)
    handles.Measurements.(PrimaryObjectName) = {};
end

handles = CPrelateobjects(handles,SecondaryObjectName,PrimaryObjectName,FinalLabelMatrixImage,EditedPrimaryLabelMatrixImage,ModuleName);
%%% [NB] hack. save memory
clear EditedPrimaryLabelMatrixImage;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the final, segmented label matrix image of secondary objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented',SecondaryObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

%[TS removed propagation]

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
%%% [NB] hack. save memory.
Centroid = cat(1,tmp.Centroid);
if isempty(Centroid)
    Centroid = [0 0];
end
handles.Measurements.(SecondaryObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

% [TS] note that the following CP code would require additional
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
% catch dummyError %[TS] bugfix for error message
%     error(['The object outlines were not calculated by the ', ModuleName, ' module, so these images were not saved to the handles structure. The Save Images module will therefore not function on these images. This is just for your information - image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
% end


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    
    %%%% [TS] %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Code for visualization  %%%%%%%%%%%
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
    %%% [NB] hack.s ave memory.
    clear MaxFilteredImage;
    %%% Converts to logical.
    warning off MATLAB:conversionToLogical
    LogicalOutlines = logical(IntensityOutlines);
    %%% [NB] hack.s ave memory.
    clear IntensityOutlines;
    warning on MATLAB:conversionToLogical
    %%% Determines the grayscale intensity to use for the cell outlines.
    %[NB-HACK] so that images are not so dim!!!!
    LineIntensity = quantile(OrigImage(:), 0.99);
    %%% Overlays the outlines on the original image.
    ObjectOutlinesOnOrigImage = OrigImage;
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
    %%% [NB] hack. save memory.
    clear DilatedPrimaryBinaryImage EditedPrimaryBinaryImage;
    BothOutlinesOnOrigImage = ObjectOutlinesOnOrigImage;
    BothOutlinesOnOrigImage(PrimaryObjectOutlines == 1) = LineIntensity;
    %%% [NB] hack. save memory.
    clear PrimaryObjectOutlines LineIntensity;
    
    %%%%%%%%%%%%%%%%%%%%%%%% END OF INITIATION OF VISUALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%% (rearrangement)
    
    
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
    end
    ObjectCoverage = 100*sum(sum(FinalLabelMatrixImage > 0))/numel(FinalLabelMatrixImage);
    
    %[TS] display range of thresholds. Which is useful if limits for treshold
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