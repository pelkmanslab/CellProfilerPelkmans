function handles = RelateAndJoinSegmentation(handles)

% Help for the Relate module:
% Category: Object Processing
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What objects are the children objects (subobjects)?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
SubObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What are the parent objects?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
ParentName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What do you want to call the objects identified by this module?
%defaultVAR03 = Organelle
%infotypeVAR03 = objectgroup indep
SegmentedObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What do you want to call the subobjects identified by this module?
%defaultVAR04 = SubOrganelle
%infotypeVAR04 = objectgroup indep
SegmentedSubObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = How did you call the corresponding intensity image?
%infotypeVAR05 = imagegroup
%inputtypeVAR05 = popupmenu
IntImageName = char(handles.Settings.VariableValues{CurrentModuleNum,5});


%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
SubObjectLabelMatrix = CPretrieveimage(handles,['Segmented', SubObjectName],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
ParentObjectLabelMatrix = CPretrieveimage(handles,['Segmented', ParentName],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieve intensity image
IntensityImage = CPretrieveimage(handles,IntImageName,ModuleName,'MustBeGray','DontCheckScale');


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

[handles,ChildList,FinalParentList] = CPrelateobjects(handles,SubObjectName,ParentName,SubObjectLabelMatrix,ParentObjectLabelMatrix,ModuleName);

%%% Since the label matrix starts at zero, we must include this value in
%%% the list to produce a label matrix image with children re-labeled to
%%% their parents values. This does not get saved and is only for display.
if ~isempty(FinalParentList)
    FinalParentListLM = [0;FinalParentList];
    NewObjectParentLabelMatrix = FinalParentListLM(SubObjectLabelMatrix+1);
    CurrentObjNhood = bwmorph(NewObjectParentLabelMatrix,'dilate',1);
    CurrentObjNhood = DilateBackground(bwlabel(CurrentObjNhood));%bwmorph(CurrentObjNhood,'erode',1);
    %CurrentObjNhood(edge(CurrentObjNhood))=0;
    CurrentObjNhood = ParentObjectLabelMatrix&CurrentObjNhood;
    NewObjectParentLabelMatrix = zeros(size(ParentObjectLabelMatrix));
    NewObjectParentLabelMatrix(CurrentObjNhood) = ParentObjectLabelMatrix(CurrentObjNhood);
    
else
    NewObjectParentLabelMatrix = SubObjectLabelMatrix;
end

NewSubobjectLabelMatrix = bwlabel(NewObjectParentLabelMatrix);

%merge the segments of a same cell

%[CurrentObjNhood,CurrentObjLabels] = bwdist(CurrentObjNhood);
%CurrentObjNhood = (CurrentObjNhood < 2).*NewObjectParentLabelMatrix(CurrentObjLabels);
% CurrentObjNhood = bwmorph(NewObjectParentLabelMatrix,'dilate',1);
% CurrentObjNhood(edge(CurrentObjNhood))=0;          
% CurrentObjNhood = ParentObjectLabelMatrix&CurrentObjNhood;
% NewObjectParentLabelMatrix = zeros(size(ParentObjectLabelMatrix));
% NewObjectParentLabelMatrix(CurrentObjNhood) = ParentObjectLabelMatrix(CurrentObjNhood);
%  
% CurrentObjNhood(CurrentObjNhood) = ParentObjectLabelMatrix(CurrentObjNhood);
% figure;imagesc(mattest)


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%

drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(ParentObjectLabelMatrix,'TwoByTwo',ThisModuleFigureNumber);
    end
    subplot(2,2,1);
    ColoredParentLabelMatrixImage = CPlabel2rgb(handles,ParentObjectLabelMatrix);
    CPimagesc(ColoredParentLabelMatrixImage,handles);
    title(['Parent Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,2,2);
    ColoredSubObjectLabelMatrixImage = CPlabel2rgb(handles,SubObjectLabelMatrix);
    CPimagesc(ColoredSubObjectLabelMatrixImage,handles);
    title('Original Sub Objects');
    subplot(2,2,3);
    ColoredNewObjectParentLabelMatrix = CPlabel2rgb(handles,NewObjectParentLabelMatrix);
    CPimagesc(ColoredNewObjectParentLabelMatrix,handles);
    title('New Sub Objects');
    OutlineOverlay = IntensityImage;
    B = bwboundaries(NewObjectParentLabelMatrix,'holes');
    subplot(2,2,4), CPimagesc(OutlineOverlay,handles),
    title(['Outlines of New Sub Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
    end
    hold off
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE TO HANDLES %%%
%%%%%%%%%%%%%%%%%%%%%%%

%save the object segmentation
fieldname = ['Segmented',SegmentedObjectName];
handles.Pipeline.(fieldname) = NewObjectParentLabelMatrix;%FinalLabelMatrixImage;

%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end

% ======= hack 2014/02/17 [MH] ==============================================================================================
% Use of 'strfind' is very dangerous, since it is not specific and finds any occurance (also substring)! I thus
% replaced it by 'strcmp'.

% save input objects to handles
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,SubObjectName));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {SubObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
elseif ~strcmp(handles.Measurements.Image.ObjectCountFeatures(column),SubObjectName)
    handles.Measurements.Image.ObjectCountFeatures(column) = {SubObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(SubObjectLabelMatrix(:));
% save output subobjects to handles
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,SegmentedSubObjectName));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {SegmentedSubObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
elseif ~strcmp(handles.Measurements.Image.ObjectCountFeatures(column),SegmentedSubObjectName)
    handles.Measurements.Image.ObjectCountFeatures(column) = {SegmentedSubObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(NewSubobjectLabelMatrix(:));
% ======= hack end ==========================================================================================================

% save output object to handles
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,SegmentedObjectName));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {SegmentedObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
elseif ~strcmp(handles.Measurements.Image.ObjectCountFeatures(column),SegmentedObjectName)
    handles.Measurements.Image.ObjectCountFeatures(column) = {SegmentedObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(NewObjectParentLabelMatrix(:));

%%% Saves the location of each segmented object
handles.Measurements.(SegmentedObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(NewObjectParentLabelMatrix,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(SegmentedObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves the location of each segmented subobject
handles.Measurements.(SegmentedSubObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(NewSubobjectLabelMatrix,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(SegmentedSubObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves parent-child relationship for each segmented object
[handles,ChildList_obj,FinalParentList_obj] = CPrelateobjects(handles,SegmentedObjectName,ParentName,NewObjectParentLabelMatrix,ParentObjectLabelMatrix,ModuleName);
[handles,ChildList_subobj,FinalParentList_subobj] = CPrelateobjects(handles,SegmentedSubObjectName,ParentName,NewSubobjectLabelMatrix,ParentObjectLabelMatrix,ModuleName);

%%% Saves the final, segmented label matrix images of objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented', SegmentedObjectName];
handles.Pipeline.(fieldname) = NewObjectParentLabelMatrix;
fieldname = ['Segmented', SegmentedSubObjectName];
handles.Pipeline.(fieldname) = NewSubobjectLabelMatrix;

end