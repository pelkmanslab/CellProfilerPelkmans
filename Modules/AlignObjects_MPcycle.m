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

%textVAR15 = Type "Pre" (Pre) to load shift descriptor file from previous multiplexing cycle. Type period (.) for using the one stored in the handles.
%defaultVAR15 = .
AlternativeShiftDescriptor = char(handles.Settings.VariableValues{CurrentModuleNum,15});

%%%VariableRevisionNumber = 12



%%%%%%%%%%%%%%%%%%
%%% PROCESSING %%%
%%%%%%%%%%%%%%%%%%

%%% retrieve shift descriptors from handles
if strncmp(AlternativeShiftDescriptor,'.',1)
    if length(AlternativeShiftDescriptor) == 1
        shift = handles.shiftDescriptor;
    else
        if strncmp(AlternativeShiftDescriptor,'./../',5)
            strAlignCyclesDir = fullfile(handles.Current.DefaultOutputDirectory,AlternativeShiftDescriptor);
            jsonDescriptorFile = fullfile(strAlignCyclesDir,'shiftDescriptor.json');
            fprintf(['=================================================' ...
                '=================================================' ...
                '=================================================' ...
                '\n\n%s: You loaded an alternative shift descriptor file:\n%s\n\n' ...
                '=================================================' ...
                '=================================================' ...
                '=================================================' ...
                '\n'],mfilename,jsonDescriptorFile);
            % load json file
            shift = loadjson(jsonDescriptorFile);
        else
            error('Please specify relative path (relative from default output directory) to ALIGNCYCLES folder that contains the shift descriptor file')
        end
    end
elseif strcmp(AlternativeShiftDescriptor,'Pre')
    % define path to json file
    % If output dir is BATCH directory, assume ALIGNCYCLES directory
    strAlignCyclesDir = handles.Current.DefaultOutputDirectory;
    if strcmp(getlastdir(strAlignCyclesDir),'BATCH')
        strAlignCyclesDir = strrep(strAlignCyclesDir, [filesep,'BATCH'],[filesep,'ALIGNCYCLES']);
    end
    % use shift descriptor of previous cycle
    cycleNum = str2double(cell2mat(flatten(regexp(strAlignCyclesDir,'.*?([0-9]+).ALIGNCYCLES','tokens'))));
    strAlignCyclesDir = strrep(strAlignCyclesDir,[num2str(cycleNum),filesep,'ALIGNCYCLES'],[num2str(cycleNum-1),filesep,'ALIGNCYCLES']);
    jsonDescriptorFile = fullfile(strAlignCyclesDir,'shiftDescriptor.json');
    fprintf(['=================================================' ...
                    '=================================================' ...
                    '=================================================' ...
                    '\n\n%s: You loaded an alternative shift descriptor file:\n%s\n\n' ...
                    '=================================================' ...
                    '=================================================' ...
                    '=================================================' ...
                    '\n'],mfilename,jsonDescriptorFile);
    % load json file
    shift = loadjson(jsonDescriptorFile);
elseif strcmp(AlternativeShiftDescriptor,'.')
    shift = handles.shiftDescriptor;
else
    error('Please specify which shift descriptor file should be used!')
end

f = isnan(shift.xShift);   % [TS] : assume NaN (no shift calculated) to be no shift
shift.xShift(f) = 0;

f = isnan(shift.yShift);
shift.yShift(f) = 0;

%%% get index of current image
strOrigImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});
% do this with a subfunction that allows using strings from different microscopes
% VisiScope: \d+_r\d{02}_c\d{02}_[A-Z]_C\d{02}
% Yokogawa: T\d{04}F\d{03}L\d{02}A\d{2}Z\d{2}C\d{2}
[~, strMicroscopeType] = check_image_position(strOrigImageName);
if strcmpi(strMicroscopeType, 'CV7K')
    strLookup = regexprep(strOrigImageName,'A\d{2}Z\d{2}C\d{2}','A\\d{2}Z\\d{2}C\\d{2}');
    strLookupShort = regexprep(strLookup,'.*(_[A-Z]\d{2}_)','$1','once');
elseif strcmpi(strMicroscopeType, 'Visi') %'_s\d{04}_r\d{02}_c\d{02}_[A-Z]+_C\d{02}'
    strLookup = regexprep(strOrigImageName,'_[A-Z]+_C\d{02}','_[A-Z]+_C\\d{02}');
    strLookupShort = regexprep(strLookup,'.*(_[A-Z]+_C\d{02})','$1','once');
else
    error('Microscope type could not be determined. Currently implemented are ''Visi'' and ''CV7K''')
end

index = find(cellfun(@(x) ~isempty(x),regexp(cellstr(shift.fileName),strLookupShort)));

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
    if shift.noShiftIndex(index)
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
    
    %%% crop segmentation images
    if shift.noShiftIndex(index)
        SegmentationOutputImages{i} = bwlabel(zeros(size(SegmentationImages{i}(1+shift.lowerOverlap : end-shift.upperOverlap, 1+shift.rightOverlap : end-shift.leftOverlap))));
    else
        SegmentationOutputImages{i} = SegmentationImages{i}(1+shift.lowerOverlap : end-shift.upperOverlap, 1+shift.rightOverlap : end-shift.leftOverlap);
    end
end

if sum(~strcmp(ObjectNameList,'Do not use'))>0
    
    %%% if more than one object, make sure that all objects have identical object counts
    if length(SegmentationOutputImages)>1
        LabeledSegmentationOutputImages = cellfun(@(x) zeros(size(x)),SegmentationOutputImages,'Uniformoutput',false);
        % get object counts
        woZeros = cellfun(@(x) x(x~=0),SegmentationOutputImages,'Uniformoutput',false);
        objectNum = cell2mat(cellfun(@(x) length(unique(x(:)) ~= 0),woZeros,'Uniformoutput',false));
        [~,IX] = sort(objectNum,'ascend');
        
        % take objects with lowest and second lowest counts
        aIx = unique(SegmentationOutputImages{IX(2)}(:));
        aIx(aIx==0)=[];
        bIx = unique(SegmentationOutputImages{IX(1)}(:));
        bIx(bIx==0)=[];
        % get common objects
        b = ismember(aIx,bIx);
        a = ismember(bIx,aIx(b));
        if ~(isequal(aIx(b),bIx(a)))
            error('object count must not be different between segmented objects')
        end
        
        % assign objects common, continuous labels
        aIx = aIx(b);
        bIx = bIx(a);
        for m = 1:length(aIx)
            for k = IX(2:end)
                LabeledSegmentationOutputImages{k}(SegmentationOutputImages{k}==aIx(m)) = m;
            end
            LabeledSegmentationOutputImages{IX(1)}(SegmentationOutputImages{IX(1)}==bIx(m)) = m;
        end
        
        %%% Note: New objects labels have to be assigned to make labels
        %%% continuous starting from 1 (necessary for downstream CP
        %%% modules). This is cycle specific, however, and could in
        %%% principle lead to different object Ids between cycles. To check
        %%% this, we relate the new object ids to the original ones, so
        %%% that we can later map the measurements back to the correct
        %%% objects.
        
        % keep track of original object labels (of the uncropped images)
        comIx = aIx;
                    
        % check that object count is finally really identical
        checkNum = cell2mat(cellfun(@(x) length(unique(x(:))),LabeledSegmentationOutputImages,'Uniformoutput',false));
        if length(unique(checkNum))>1
            error('object count must not be different between segmented objects');
        end
        
    elseif length(SegmentationOutputImages)==1
        comIx = unique(SegmentationOutputImages{1});
        comIx(comIx==0) = []; % remove background 
        LabeledSegmentationOutputImages{1} = bwlabel(SegmentationOutputImages{1});
    end

end


%%%%%%%%%%%%%%%
%%% DISPLAY %%%
%%%%%%%%%%%%%%%


drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    
    IntensityImages4Display = IntensityOutputImages(cellfun(@(x) ~isempty(x),IntensityOutputImages));
       
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    
    %%% Display aligned image
    subplot(1,2,1); 
    CPimagesc(IntensityImages4Display{1},handles);
    title([ObjectNameList{1}, sprintf(' Outlines on ''%s'' Image', IntImOutputNameList{1})]);
    if sum(~strcmp(ObjectNameList,'Do not use'))>0
        SegmentationImage4Display = LabeledSegmentationOutputImages{1};
    elseif isfield(handles.Pipeline,'SegmentedCells')
        SegmentationImage4Display = handles.Pipeline.('SegmentedCells');
    else
        SegmentationImage4Display = [];
    end
    %%% Calculates object outlines, which are overlaid on the intensity image.
    B = bwboundaries(SegmentationImage4Display,'holes');
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
    end
    hold off
    
    if sum(~cellfun(@isempty,IntensityImages4Display))>1
        
        %%% Display aligned image
        subplot(1,2,2);
        CPimagesc(IntensityImages4Display{2},handles);
        title([ObjectNameList{2}, sprintf(' Outlines on ''%s'' Image', IntImOutputNameList{2})]);
        %%% Calculates object outlines, which are overlaid on the intensity image.
        B = bwboundaries(SegmentationImage4Display,'holes');
        hold on
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
        end
        hold off
        
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
                handles.Measurements.(ObjectName).Children(handles.Current.SetBeingAnalyzed) = {childrenCount};
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
                handles.Measurements.(ObjectName).Parent(handles.Current.SetBeingAnalyzed) = {parentCount};
            end
        end
    end
    
    % save relation to original object ID to handles
    handles.Measurements.(ObjectName).OrigObjectIDFeatures{handles.Current.SetBeingAnalyzed} = 'OriginalObjectId';
    handles.Measurements.(ObjectName).OrigObjectID(handles.Current.SetBeingAnalyzed) = {comIx};
    
    
end

end
