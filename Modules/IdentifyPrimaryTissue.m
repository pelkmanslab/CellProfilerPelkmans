function handles = IdentifyPrimaryTissue(handles)

% Help for IdentifyPrimaryTissue
% Category: Object Processing
%
%
% DESCRIPTION:
% Multistage algorithm for automated segmentation of tissue-section images.
% Performs segmentation of primary objects (nuclei) in a multistep process:
% Step 1: 
% Integration of nuclear and membranous staining images 
% Step 2: 
% Separation of foreground objects from background using LoG filter.
% Step 3:
% Propagation of objects up to outlines in original nuclear staining
% intensity image.

%
% PARAMETERS:
% Laplacion of Gaussian/shape based segmentation for low intensity images.
% For detailed information on the LoG see LowIntensityLoG.m.
%
% LoG Sigma: as low as possible without subsegmenting/eroding objects or 
% missing low intensity objects. Cells should apper as uniform objects
% (no holes) in the LoG output image. Low=more detail,rather oversegmentation,
% worse detection for low intensity. High=less detail, undersegmentation,
% better detection of low intensity objects.
%
% Propagation Factor: Finetunes IdentifySecPropagateSubfunction function.
% For help with subfunction see IdentifySecondary.m.
% Value between 0 and 1. 
% 0 = Intensity only; 1 = Distance only
%
% Threshold Correction Factor: Finetunes threshold for creation of
% foreground mask (required for propagate subfunction). Prevent object
% erosion; thus select threshold rather less sensitive. Values between 0
% and 1. Smaller = less sensitive.
%
% LoG Cutoff: As high as possible to avoid excluding excessive object parts,
% but still low enough separate clumped objects. Pick lowest values from
% LoG image regions connecting clumped objects. For subsequent shape based
% segmentation Log Cutoff can be set zero to completely disable it or to a 
% rather high value in order not to exclude low intensity objects or erode
% objects.
%
% LoG Marker Cutoff: To prevent object erosion but still allow background
% separation, object markers (eroded objects) are created via a
% cuttoff (may be much more stringent than LoG Cutoff). Only objects
% containing a minimum fraction of marker overlay are later considered as
% true objects. In general as low as possible without losing objects. Value
% should be lower than the lowest positive background value, but higher
% than the lowest value of the object with the highest values. Smaller than
% LoG Cutoff (if not disabled).
%
% LoG Minimum Marker to Object Fraction: determines the minimum fraction of
% marker area/object area to still detect an object. This allows to
% distinguish background (wich may also contain very few marker pixels)
% from true objects (containing a larger fraction of marker pixels). Low
% Marker Cutoff requires low marker fraction, high Marker Cutoff allows
% high marker fraction.
%
% LoG Minimum Object Area: exclued very small objects eg dirt.
%
% Erosion/Dilation Smoothing: size of smoothing filter.
%
% Stack Threshold: percentage of stacks that must include the object
% 
% Maximal Background Value: if the 99th quantile intensity of an image is
% below this value, the image is regarded empty and is not processed by the
% module.
%
%
% DEPENDENCIES: 
% IdentifySecPropagateSubfunction.mexmaci64
% DilateBackground.m
% imnorm.m
% LowIntensityLoG.m
%
%
% [Markus Herrmann]


%Todo: 
%-maybe save perimeter freatures to handles?
%-thresholding correction factor as variable, create test mode to check thresholing
%-consider as a final step to merge objects if they together create a
% "nice" nucleus (e.g. optimizing for solidity)


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

drawnow
[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);
%textVAR01 = What did you call the NUCLEI staining images you want to process?
%infotypeVAR01 = imagegroup
NucleiStackName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the MEMBRANE staining images you want to process?
%infotypeVAR02 = imagegroup
MembraneStackName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the objects identified by this module?
%defaultVAR03 = Nuclei
%infotypeVAR03 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = LoG Sigma
%defaultVAR04 = 7
LoGSigma = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,4}));

%textVAR05 = Propagation factor, value between 0 and 1 (0=Intensity, 1=Distance)
%defaultVAR05 = 0.03
PropagFactor = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = Threshold correction factor 
%defaultVAR06 = 0.8
correctionFactor = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,6}));

%textVAR07 = LoG Marker Cutoff
%defaultVAR07 = 9
LoGMarkerCutoff = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = LoG minimum marker to object fraction
%defaultVAR08 = 0.6
LoGObjMarkerRatio = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,8}));

%textVAR09 = LoG minimum object area
%defaultVAR09 = 80
LoGObjSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,9}));

%textVAR10 = Erosion/dilation smoothing (remove sharp edges), disk radius (0=off)
%defaultVAR10 = 4
SmoothingSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,10}));

%textVAR11 = Membrane segmentation: percentage of stacks that must include the object
%defaultVAR11 = 0.9
stackSegmentationThreshold = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,11}));

%textVAR12 = Discard border objects
%choiceVAR12 = No
%choiceVAR12 = Yes
DiscardBorder = char(handles.Settings.VariableValues{CurrentModuleNum,12});
%inputtypeVAR12 = popupmenu

%textVAR13 = Maximal background value
%defaultVAR13 = 200
MaxBackgroundValue = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,13}));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stackDapi = handles.Pipeline.(NucleiStackName);
stackDapi = double(stackDapi);
stackMembrane = handles.Pipeline.(MembraneStackName);
stackMembrane = double(stackMembrane);

% imDapi = CPretrieveimage(handles,NucleiStackName,ModuleName,'MustBeGray','CheckScale');
%imDapi = imDapi .* 65535;%reverse rescaling done by CP!

LoGSize = 8 * LoGSigma;

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%

%% Load images %%

%%% Load DAPI stacks
i = size(stackDapi,3);
% Add MIP to stack
stackDapi(:,:,i+1) = max(stackDapi,[],3);

%%% Load Membrane stacks
j = size(stackMembrane,3);
% Add MIP to stack
stackMembrane(:,:,j+1) = max(stackMembrane,[],3);

%%% Preallocate matrices
SegmentationPerStack(:,:,size(stackDapi,3)) = zeros(size(stackDapi(:,:,1)));
ThresholdMaskPerStack(:,:,size(stackDapi,3)) = zeros(size(stackDapi(:,:,1)));

% Only process image, if image is not empty (i.e. if maximum pixel intensity lies above background)
imMIPDapi = stackDapi(:,:,end);
imMIPMembrane = stackMembrane(:,:,end);
if quantile(imMIPDapi(:),0.99) < MaxBackgroundValue
    
    imThresholdStackOverlay = zeros(size(imMIPDapi));
    imSegmentationStackOverlay = zeros(size(imMIPDapi));
    imProjSegmentation = zeros(size(imMIPDapi));
    imFinalObjectsIntensityBorder = zeros(size(imMIPDapi));
    ImageProcessed = false;
        
    %%% Normalize images
    imNormDapi = imnorm(imMIPDapi,0.01,0.99);
    imNormMembrane = imnorm(imMIPMembrane,0.01,0.99);
    
    %%% Subtract membrane image from DAPI image
    imSub = immultiply(imNormDapi,imcomplement(imNormMembrane));
    
else
    ImageProcessed = true;
    
    for i = 1:size(stackDapi,3)
        
        
        %% Integrate nuclear and nuclear membrane images into one image %%
        
        %%% Load images
        imDapi = stackDapi(:,:,i);
        imMembrane = stackMembrane(:,:,i);
        
        %%% Normalize images
        imNormDapi = imnorm(imDapi,0.01,0.99);
        imNormMembrane = imnorm(imMembrane,0.01,0.99);
        
        %%% 'Subtract' membrane image from DAPI image (weight)
        imSub = immultiply(imNormDapi,imcomplement(imNormMembrane));
        
        
        %% Seperate foreground (objects) from background using LoG filter %%
        
        % Parts of the code are derived from IdentifyPrimLoGShape.m
        
        %%% Apply LoG (Laplacion of Gaussian) filter
        LoGCutoff = 0;
        [imLoGObjectMask,~,~] = LowIntensityLoG(imSub,LoGSigma,[LoGSize,LoGSize],LoGCutoff,LoGMarkerCutoff,LoGObjMarkerRatio,LoGObjSize);
        imLoGObjectsLabel = bwlabel(imLoGObjectMask);
        
        %%% Apply a slight smoothing before thresholding (copied from IdentifyPrimAutomatic.m)
        sigma = 1;
        FiltLength = 2*sigma;                                              % Determine filter size, min 3 pixels, max 61
        [x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);   % Filter kernel grid
        f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                 % Gaussian filter kernel
        imDapiBlurred = conv2(imDapi,f,'same') ./ conv2(ones(size(imDapi)),f,'same');
        
        %% Perform thresholding
        ThresOtsu = multithresh(imDapiBlurred); %ultimately, one may use CPthreshold !!!
        ThresCorrected = ThresOtsu * correctionFactor;
        imDapiThres = imDapi > ThresCorrected;
        
        %%% Re-dilate objects up to border of the original object mask while retaining newly gained objects (see IdentifySecondary.m)
        imPropagatedMembraneMask = IdentifySecPropagateSubfunction(imLoGObjectsLabel,imDapi,logical(imDapiThres),PropagFactor);
        
        %%% Create nicely outlined objects (see MeasureLocalizationOfSpots.m)
        if ~max(max(imPropagatedMembraneMask))==0
            imDilatedMembraneMask = DilateBackground(imPropagatedMembraneMask);
        else
            imDilatedMembraneMask = imPropagatedMembraneMask;
        end
        
        %%% Smooth image
        if SmoothingSize == 0
            SmoothDisk = getnhood(strel('disk',1,0));%minimum that has to be done to avoid problems with bwtraceboundary
        else
            SmoothDisk = getnhood(strel('disk',SmoothingSize,0));%helps to avoid unnatural sharp edges/high curvatures for further cutting
        end
        
        %%% Fill holes
        imFinalMembraneMask = imfill(imDilatedMembraneMask);
        
        %%% Label objects
        imFinalMembraneMask = imdilate(imerode(imFinalMembraneMask,SmoothDisk),SmoothDisk);
        
        %%% Store segmentation for each stack
        %LoGMask(:,:,i) = logical(imLoGObjectMask);
        ThresholdMaskPerStack(:,:,i) = logical(imDapiThres);
        SegmentationPerStack(:,:,i) = logical(imFinalMembraneMask);
        
    end
    
    
    %% Combine segmentations of individual z-stacks into one final segmenation image %%
    
    %%% Combine segmentations
    imThresholdStackOverlay = sum(ThresholdMaskPerStack,3);
    imSegmentationStackOverlay = sum(SegmentationPerStack,3);
    % Use 3D information of z-stacks to create gaps between clumped nuclei
    imProjThresholdMask = imThresholdStackOverlay > floor(size(stackDapi,3)*0.6);
    imProjThresholdMask = imfill(imProjThresholdMask,'holes');
    imProjSegmentation = imSegmentationStackOverlay > floor(size(stackDapi,3)*stackSegmentationThreshold);
    imOverlay = imProjSegmentation + imProjThresholdMask;
    imProjSegmentation = imOverlay==2;

    % imLower = imSegmentationStackOverlay < floor(size(stackDapi,3)*stackSegmentationThreshold) & imSegmentationStackOverlay > 0;%floor(size(stackDapi,3)*0.35);
    % imLowerSolidity = rplabel(imLower,[],'Solidity');
    % imOff = imLowerSolidity < 0.8 & imLower > 0;
    % imOffDil = imdilate(imOff,strel('disk',2));
    % imProjectedSegmentation = logical(imSegmentationStackOverlay) .* ~imOffDil;% > floor(size(stackDapi,3)*stackSegmentationThreshold);
    
    % %%% Recover objects that fall below the 'stackSegmentationThreshold', but have typical nuclear features (specifically, high solidity values)
    % imLowerEroded = imerode(imLower,strel('disk',2));
    % imLowerErodedSolidity = rplabel(imLowerEroded,[],'Solidity');
    % imMissed = imLowerErodedSolidity > 0.8;% & imLeftArea > 100;
    % imRecoveredSegmentation = logical(imProjectedSegmentation + imMissed);
    
    %%% Recover objects that are missed by the LoG filter (as a consequence of nuclear SUN2 staining)
    imThresRecover = imOverlay==1;
    % Only recover objects within a certain range of area and solidity
    ThresAreas = cell2mat(struct2cell(regionprops(imThresRecover,'Area')))';
    ValidThresAreas = find(ThresAreas<2000 & ThresAreas>100);
    imThresArea = ismember(bwlabel(imThresRecover),ValidThresAreas);%figure,imagesc(rplabel(imThresRecover,[],'Area'))
    ThresSolidity = cell2mat(struct2cell(regionprops(imThresRecover,'Solidity')))';
    ValidThresSolidity = find(ThresSolidity>0.8);
    imThresSolidity = ismember(bwlabel(imThresRecover),ValidThresSolidity);%figure,imagesc(rplabel(imThresRecover,[],'Solidity'))
    imMissed = imfill((imThresArea + imThresSolidity)==2,'holes');
    imRecovSegmentation = logical(imProjSegmentation + imMissed);
    
    %%% Discard small objects
    propsObj = regionprops(imRecovSegmentation,'Area','PixelIdxList');
    ProjectedSegmentationArea = cat(1,propsObj.Area);
    Objects = 1:length(ProjectedSegmentationArea);
    ObjectsArea = Objects(ProjectedSegmentationArea>100);
    PixelIndex = cat(1,propsObj(ObjectsArea).PixelIdxList);
    imSegmentedObjects = zeros(size(imRecovSegmentation));
    imSegmentedObjects(PixelIndex) = 1;
    imSegmentedObjects = imfill(imSegmentedObjects);
    
    %%% Propagate object outlines
    imSegmentedObjectsLabel = bwlabel(imSegmentedObjects);
    imSegmentedObjects = IdentifySecPropagateSubfunction(imSegmentedObjectsLabel,imDapi,logical(imProjThresholdMask),0.03);%0
    if ~max(imSegmentedObjects(:))==0
        imSegmentedObjects = DilateBackground(imSegmentedObjects);
    end
    
    %%% Smooth image
    imFinalObjects = imdilate(imerode(imSegmentedObjects,SmoothDisk),SmoothDisk);
    imFinalObjects = imfill(imFinalObjects);
    
    
    %% Final clean-up %%
    
    %%% Size filter: discard small objects
    ObjectAreas = cell2mat(struct2cell(regionprops(imFinalObjects,'Area')))';
    ValidObjectIndices = find(ObjectAreas>100);
    imFinalObjectsArea = bwlabel(ismember(imFinalObjects,ValidObjectIndices));%relabel to get continuous indices
    %%% Border filter: disard border objects
    if strcmp(DiscardBorder,'Yes')
        BorderIndices = setdiff(unique([imFinalObjectsArea(1,1:end),imFinalObjectsArea(end,1:end),imFinalObjectsArea(1:end,1)',imFinalObjectsArea(1:end,end)']),0);
        BorderObjectmask = ismember(imFinalObjectsArea,BorderIndices);
        imFinalObjectsIntensityBorder = bwlabel(imFinalObjectsArea.*~BorderObjectmask);%relabel to get continuous indices
    else
        imFinalObjectsIntensityBorder = imFinalObjectsArea;
    end
    
end
    




%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%

drawnow
% %Save colormap
% ColormapBackup = handles.Preferences.IntensityColorMap;
% %Set custom colormap
% handles.Preferences.IntensityColorMap = jet;

if ImageProcessed
    
    % Create overlay image
    cellB = bwboundaries(imFinalObjectsIntensityBorder,'holes');
    matB = cell2mat(cellB);
    linInd = sub2ind(size(imNormDapi),matB(:,1),matB(:,2));
    imSegOver = zeros(size(imNormDapi));
    imSegOver(linInd) = 1;
    imSegOver = imdilate(logical(imSegOver),strel('disk',1));
    imOutlineShapeSeparatedOverlay = imoverlay(imNormDapi,imSegOver,[1,1,0]);
    % Create object label image
    imObjectsLabel = label2rgb(bwlabel(imFinalObjectsIntensityBorder),'jet',[1 1 1],'shuffle');
    
else
    
    imSegmentationStackOverlay = zeros(size(imMIPDapi));
    imOutlineShapeSeparatedOverlay = zeros(size(imMIPDapi));
    imObjectsLabel = zeros(size(imMIPDapi));
    imFinalObjects = zeros(size(imMIPDapi));
    
end

    
%GUI output
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    CPfigure(handles,'PrimObj: Membrane segmentation',ThisModuleFigureNumber);
    subplot(2,3,1), CPimagesc(imSub,handles), colormap('default'),
    title(['Subtracted image (MIP), cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,3,2), CPimagesc(imSegmentationStackOverlay,handles), colormap('default'),
    title(['Stack segmentation overlay, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,3,3), CPimagesc(imThresholdStackOverlay,handles), colormap('default'),
    title(['Stack threshold overlay, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,3,4), CPimagesc(imProjSegmentation,handles), colormap('default'),
    title(['Stack projected segmentation overlay, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,3,5), CPimagesc(imOutlineShapeSeparatedOverlay,handles),
    title(['Final objects outline, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,3,6), CPimagesc(imObjectsLabel,handles),
    title(['Final objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
end
       

% %Restore old colormap
% handles.Preferences.IntensityColorMap = ColormapBackup;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%todo:
%-consider saving more data to handle structure as done in IdentifyPrimAutormatic.m

fieldname = ['UneditedSegmented',ObjectName];%not edited for size or edge
handles.Pipeline.(fieldname) = imFinalObjects;

fieldname = ['SmallRemovedSegmented',ObjectName];%for IdentifySecondary.m
handles.Pipeline.(fieldname) = imFinalObjects;

fieldname = ['Segmented',ObjectName];%final label image
handles.Pipeline.(fieldname) = imFinalObjectsIntensityBorder;

%%% Saves location of each segmented object
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(imFinalObjectsIntensityBorder,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves ObjectCount, i.e. number of segmented objects.
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(imFinalObjectsIntensityBorder(:));
end