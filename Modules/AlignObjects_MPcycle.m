function handles = AlignObjects_MPcycle(handles)

% Help for the ALIGNOBJECTS_MPCYCLE module:
% Category: Other
%
% SHORT DESCRIPTION: 
% Loads shift descriptors (structure stored as .json file) from iBRAIN's
% ALIGNCYCLES directory and obtains intensity as well as segmentation image
% from the handles. Then it aligns the intensity and the segmentation
% images. To this end, both images are cropped according to the loaded
% shift descriptors.
%
% To segment different object groups, which can have different amount of
% objects (like cells and single vesicles), use separate instances of the
% AlignObjects_MPcycle module. Individual objects within an object group
% require a 1:1 matching
% *************************************************************************
%
% Author:
%    Markus Herrmann <markus.herrmann@imls.uzh.ch>
%
% $Revision: 1879 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the intensity images that you want to shift?
%infotypeVAR01 = imagegroup
IntImNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = 
%choiceVAR02 = Do not use
%infotypeVAR02 = imagegroup
IntImNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = 
%choiceVAR03 = Do not use
%infotypeVAR03 = imagegroup
IntImNameList{3} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = 
%choiceVAR04 = Do not use
%infotypeVAR04 = imagegroup
IntImNameList{4} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = What did you call the objects that you want to measure later?
%choiceVAR05 = Do not use
%infotypeVAR05 = objectgroup
ObjectNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu

%textVAR06 = 
%choiceVAR06 = Do not use
%infotypeVAR06 = objectgroup
ObjectNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%textVAR07 =
%choiceVAR07 = Do not use
%infotypeVAR07 = objectgroup
ObjectNameList{3} = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu

%textVAR08 =
%choiceVAR08 = Do not use
%infotypeVAR08 = objectgroup
ObjectNameList{4} = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 =
%choiceVAR09 = Do not use
%infotypeVAR09 = objectgroup
ObjectNameList{5} = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 =
%choiceVAR10 = Do not use
%infotypeVAR10 = objectgroup
ObjectNameList{6} = char(handles.Settings.VariableValues{CurrentModuleNum,10});
%inputtypeVAR10 = popupmenu

%textVAR11 = How do you want to call the shifted intensity images?
%defaultVAR11 = AlignedGreen
%infotypeVAR11 = imagegroup indep
IntImOutputNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,11});

%textVAR12 = 
%defaultVAR12 = /
%infotypeVAR12 = imagegroup indep
IntImOutputNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,12});

%textVAR13 = 
%defaultVAR13 = /
%infotypeVAR13 = imagegroup indep
IntImOutputNameList{3} = char(handles.Settings.VariableValues{CurrentModuleNum,13});

%textVAR14 = 
%defaultVAR14 = /
%infotypeVAR14 = imagegroup indep
IntImOutputNameList{4} = char(handles.Settings.VariableValues{CurrentModuleNum,14});

%%%VariableRevisionNumber = 12



%%%%%%%%%%%%%%%%%%
%%% PROCESSING %%%
%%%%%%%%%%%%%%%%%%

%%% retrieve shift descriptors from handles
shift = handles.shiftDescriptor;

f = isnan(shift.xShift);   % [TS] : assume NaN (no shift calculated) to be no shift
shift.xShift(f) = 0;

f = isnan(shift.yShift);
shift.yShift(f) = 0;

%%% get index of current image
strOrigImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});
strLookup = regexprep(strOrigImageName,'A\d{2}Z\d{2}C\d{2}','A\\d{2}Z\\d{2}C\\d{2}');
index = find(cell2mat(regexp(cellstr(shift.fileName),strLookup)));

IntensityImages = cell(1,length(IntImNameList));
IntensityOutputImages = cell(1,length(IntImNameList));
for j = 1:length(IntImNameList)
    ImageName = IntImNameList{j};
    if strcmpi(ImageName,'Do not use')
        continue
    end
    
    %%% retrieve intensity images
    IntensityImages{j} = CPretrieveimage(handles,IntImNameList{j},ModuleName,'MustBeGray','CheckScale');
    
    %%% shift/crop intensity images
    if abs(shift.yShift(index))>shift.maxShift || abs(shift.xShift(index))>shift.maxShift % don't shift images if shift values are very high (reflects empty images)
        IntensityOutputImages{j} = zeros(size(IntensityImages{j}(1+shift.lowerOverlap : end-shift.upperOverlap, 1+shift.rightOverlap : end-shift.leftOverlap)));
    else
        IntensityOutputImages{j} = IntensityImages{j}(1+shift.lowerOverlap-shift.yShift(index) : end-(shift.upperOverlap+shift.yShift(index)), 1+shift.rightOverlap-shift.xShift(index) : end-(shift.leftOverlap+shift.xShift(index)));
    end
end

% do the same for segmentation images
SegmentationImages = cell(1,sum(~strcmp(ObjectNameList,'Do not use')));
SegmentationOutputImages = cell(1,sum(~strcmp(ObjectNameList,'Do not use')));
for i = 1:sum(~strcmp(ObjectNameList,'Do not use'))
    ObjectName = ObjectNameList{i};
    if strcmpi(ObjectName,'Do not use')
        continue
    end
    
    %%% load segmentation images
    SegmentationImages{i} = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName,'MustBeGray','DontCheckScale');
    
    %%% crop segmenation images
    if abs(shift.yShift(index))>shift.maxShift || abs(shift.xShift(index))>shift.maxShift
        SegmentationOutputImages{i} = zeros(size(SegmentationImages{i}(1+shift.lowerOverlap : end-shift.upperOverlap, 1+shift.rightOverlap : end-shift.leftOverlap)));
    else
        SegmentationOutputImages{i} = SegmentationImages{i}(1+shift.lowerOverlap : end-shift.upperOverlap, 1+shift.rightOverlap : end-shift.leftOverlap);
    end
end

% track, wheher all segementations should be cleared (e.g.: if object count 
% differs or the shift is above the maximally allowed tolerance for shifts);
doEraseSegementation = false; 


%%% if more than one object, make sure that all objects have identical object counts
if length(SegmentationOutputImages)>1
    LabeledSegmentationOutputImages = SegmentationOutputImages;
    % get object count
    objectNum = cell2mat(cellfun(@(x) length(unique(x(:))),SegmentationOutputImages,'Uniformoutput',false));  
    [~,IX] = sort(objectNum,'ascend');
    
    for k = IX(2:end)
        aIx = unique(SegmentationOutputImages{k}(:));
        aIx(aIx==0)=[];
        bIx = unique(SegmentationOutputImages{IX(1)}(:));
        bIx(bIx==0)=[];
        % get common objects
        [b,Locb] = ismember(aIx,bIx);
        [~,Loca] = ismember(bIx,aIx(b));

        % give objects new labels (discard objects that are not part of both object classes)
        for m = 1:length(Locb)
            LabeledSegmentationOutputImages{k}(SegmentationOutputImages{k}==aIx(m)) = Locb(m);
        end
        for n = 1:length(Loca)
            LabeledSegmentationOutputImages{IX(1)}(SegmentationOutputImages{IX(1)}==bIx(n)) = Loca(n);
        end
    end

    % check that object count is finally really identical
    checkNum = cell2mat(cellfun(@(x) max(unique(x(:))),LabeledSegmentationOutputImages,'Uniformoutput',false));
    if length(unique(checkNum))>1
        % [TS 14-05-04] : changed default behaviour: instead of error
        % (which will break analysis of plate-wide projects), there will be no
        % segmentation used
        fprintf('object count must not be different between segmented objects: Assume absent segmentation \n');
        doEraseSegementation = true;
    end
    
else
    %%% only one object, so business as usual
    LabeledSegmentationOutputImages{1} = bwlabel(SegmentationOutputImages{1});  
end

%%% if shift exeeds maximally tolerated valued, exclude site from analysis
if logical(shift.noShiftIndex(index))
    doEraseSegementation = true;
end

% Replace segmentation by empty image, in case there was a problem
if doEraseSegementation == true;
    for k = 1:length(ObjectNameList)
        LabeledSegmentationOutputImages{k} = bwlabel(zeros(size(IntensityOutputImages{1})));
    end
end


%%%%%%%%%%%%%%%
%%% DISPLAY %%%
%%%%%%%%%%%%%%%


drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
       
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    
    %%% Calculates the object outlines, which are overlaid on the intensity
    %%% image.
    %%% Creates the structuring element that will be used for dilation.
    StructuringElement = strel('square',3);
    %%% Converts the FinalLabelMatrixImage to binary.
    FinalBinaryImage = im2bw(SegmentationImages{1},.5);
    %%% Dilates the FinalBinaryImage by one pixel (8 neighborhood).
    DilatedBinaryImage = imdilate(FinalBinaryImage, StructuringElement);
    %%% Subtracts the FinalBinaryImage from the DilatedBinaryImage,
    %%% which leaves the PrimaryObjectOutlines.
    PrimaryObjectOutlines = DilatedBinaryImage - FinalBinaryImage;
    %%% Overlays the object outlines on the original image.
    ObjectOutlinesOnIntImage = IntensityImages{1};
    %%% Determines the grayscale intensity to use for the cell outlines.
    LineIntensity = max(IntensityOutputImages{1}(:));
    ObjectOutlinesOnIntImage(PrimaryObjectOutlines == 1) = LineIntensity;
    %%% display original images
    subplot(2,2,1); 
    CPimagesc(ObjectOutlinesOnIntImage,handles);
    title([ObjectNameList{1}, sprintf(' Outlines on ''%s'' Image', IntImNameList{1})]);
    
    %%% Calculates the object outlines, which are overlaid on the intensity
    %%% image.
    %%% Creates the structuring element that will be used for dilation.
    StructuringElement2 = strel('square',3);
    %%% Converts the FinalLabelMatrixImage to binary.
    FinalBinaryImage2 = im2bw(LabeledSegmentationOutputImages{1},.5);
    %%% Dilates the FinalBinaryImage by one pixel (8 neighborhood).
    DilatedBinaryImage2 = imdilate(FinalBinaryImage2, StructuringElement2);
    %%% Subtracts the FinalBinaryImage from the DilatedBinaryImage,
    %%% which leaves the PrimaryObjectOutlines.
    PrimaryObjectOutlines2 = DilatedBinaryImage2 - FinalBinaryImage2;
    %%% Overlays the object outlines on the original image.
    ObjectOutlinesOnIntImage2 = IntensityOutputImages{1};
    %%% Determines the grayscale intensity to use for the cell outlines.
    LineIntensity2 = max(IntensityOutputImages{1}(:));
    ObjectOutlinesOnIntImage2(PrimaryObjectOutlines2 == 1) = LineIntensity2;
    %%% display aligned images
    subplot(2,2,2); 
    CPimagesc(ObjectOutlinesOnIntImage2,handles);
    title([ObjectNameList{1}, sprintf(' Outlines on ''%s'' Image', IntImOutputNameList{1})]);
    
    if sum(~cellfun(@isempty,IntensityOutputImages))>1
        
        %%% Calculates the object outlines, which are overlaid on the intensity
        %%% image.
        %%% Creates the structuring element that will be used for dilation.
        StructuringElement = strel('square',3);
        %%% Converts the FinalLabelMatrixImage to binary.
        FinalBinaryImage = im2bw(SegmentationImages{1},.5);
        %%% Dilates the FinalBinaryImage by one pixel (8 neighborhood).
        DilatedBinaryImage = imdilate(FinalBinaryImage, StructuringElement);
        %%% Subtracts the FinalBinaryImage from the DilatedBinaryImage,
        %%% which leaves the PrimaryObjectOutlines.
        PrimaryObjectOutlines = DilatedBinaryImage - FinalBinaryImage;
        %%% Overlays the object outlines on the original image.
        ObjectOutlinesOnIntImage = IntensityImages{2};
        %%% Determines the grayscale intensity to use for the cell outlines.
        LineIntensity = max(IntensityOutputImages{2}(:));
        ObjectOutlinesOnIntImage(PrimaryObjectOutlines == 1) = LineIntensity;
        %%% display original images
        subplot(2,2,3);
        CPimagesc(ObjectOutlinesOnIntImage,handles);
        title([ObjectNameList{1}, sprintf(' Outlines on ''%s'' Image', IntImNameList{2})]);
        
        %%% Calculates the object outlines, which are overlaid on the intensity
        %%% image.
        %%% Creates the structuring element that will be used for dilation.
        StructuringElement2 = strel('square',3);
        %%% Converts the FinalLabelMatrixImage to binary.
        FinalBinaryImage2 = im2bw(LabeledSegmentationOutputImages{1},.5);
        %%% Dilates the FinalBinaryImage by one pixel (8 neighborhood).
        DilatedBinaryImage2 = imdilate(FinalBinaryImage2, StructuringElement2);
        %%% Subtracts the FinalBinaryImage from the DilatedBinaryImage,
        %%% which leaves the PrimaryObjectOutlines.
        PrimaryObjectOutlines2 = DilatedBinaryImage2 - FinalBinaryImage2;
        %%% Overlays the object outlines on the original image.
        ObjectOutlinesOnIntImage2 = IntensityOutputImages{2};
        %%% Determines the grayscale intensity to use for the cell outlines.
        LineIntensity2 = max(IntensityOutputImages{2}(:));
        ObjectOutlinesOnIntImage2(PrimaryObjectOutlines2 == 1) = LineIntensity2;
        %%% display aligned images
        subplot(2,2,4);
        CPimagesc(ObjectOutlinesOnIntImage2,handles);
        title([ObjectNameList{1}, sprintf(' Outlines on ''%s'' Image', IntImOutputNameList{2})]);
        
    end
    
    drawnow
end



%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%


%%% save shifted images to handles
for j = 1:length(IntImNameList)
    ImageName = IntImNameList{j};
    if strcmpi(ImageName,'Do not use')
        continue
    end
    if strcmpi(IntImOutputNameList{j},'/')
        error('name for output intensity image not specified')
    else
        handles.Pipeline.(IntImOutputNameList{j}) = IntensityOutputImages{j};
    end
end

for i = 1:length(ObjectNameList)
    ObjectName = ObjectNameList{i};
    if strcmpi(ObjectName,'Do not use')
        continue
    end
    
    %%% Saves the segmented image, not edited for objects along the edges or
    %%% for size, to the handles structure.
    fieldname = ['UneditedSegmented',ObjectName];
    handles.Pipeline.(fieldname) = LabeledSegmentationOutputImages{i};
    
    %%% Saves the segmented image, only edited for small objects, to the
    %%% handles structure.
    fieldname = ['SmallRemovedSegmented',ObjectName];
    handles.Pipeline.(fieldname) = LabeledSegmentationOutputImages{i};
    
    %%% Saves the final segmented label matrix image to the handles structure.
    fieldname = ['Segmented',ObjectName];
    handles.Pipeline.(fieldname) = LabeledSegmentationOutputImages{i};
    
    %%% Saves the ObjectCount, i.e., the number of segmented objects.
    %%% See comments for the Threshold saving above
    if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
        handles.Measurements.Image.ObjectCountFeatures = {};
        handles.Measurements.Image.ObjectCount = {};
    end
    column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,ObjectName));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(LabeledSegmentationOutputImages{i}(:));
    
    %%% Saves the location of each segmented object
    handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
    tmp = regionprops(LabeledSegmentationOutputImages{i},'Centroid');
    Centroid = cat(1,tmp.Centroid);
    handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
    %%% Saves the parent-children relations of each segmented object
    % for parent objects
    if isfield(handles.Measurements.(ObjectName),'ChildrenFeatures')
        childrenCount = [];%not optimal, but saver
        for obj = 1:length(ObjectNameList)
            if strcmpi(ObjectNameList{obj},'Do not use')
                continue
            end
            % identify parent objects
            ixParent = find(strcmp(handles.Measurements.(ObjectName).ChildrenFeatures,ObjectNameList{obj}));
            if ~isempty(ixParent)
                % calculate new object counts for the children objects
                parentObj = max(unique(LabeledSegmentationOutputImages{ixParent}));
                for subobj = 1:parentObj
                    childObj = unique(LabeledSegmentationOutputImages{obj}(LabeledSegmentationOutputImages{ixParent}==subobj));
                    childObj(childObj==0) = [];
                    childrenCount(subobj,ixParent) = length(childObj);
                end
                % save new object counts to handles
                handles.Measurements.(ObjectName).Children{handles.Current.SetBeingAnalyzed} = childrenCount;
            end
        end
        
    end
    % for children objects
    if isfield(handles.Measurements.(ObjectName),'ParentFeatures')
        parentCount = [];
        for obj = 1:length(ObjectNameList)
            if strcmpi(ObjectNameList{obj},'Do not use')
                continue
            end
            % identify children object
            ixChild = find(strcmp(handles.Measurements.(ObjectName).ParentFeatures,ObjectNameList{obj}));
            if ~isempty(ixChild)
                % calculate new parent ids
                childObj = max(unique(LabeledSegmentationOutputImages{ixChild}));%minus 1 to get rid of background 0
                for subobj = 1:childObj
                    parentObj = unique(LabeledSegmentationOutputImages{obj}(LabeledSegmentationOutputImages{ixChild}==subobj));
                    parentObj(parentObj==0) = [];
                    parentCount(subobj,ixChild) = parentObj;
                end
                % save new object counts to handles
                handles.Measurements.(ObjectName).Parent{handles.Current.SetBeingAnalyzed} = parentCount;
            end
        end
    end
    
    
end

end