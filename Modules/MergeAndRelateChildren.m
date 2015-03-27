function handles = MergeAndRelateChildren(handles)

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
%defaultVAR03 = SubOrganelle
%infotypeVAR03 = objectgroup indep
SegmentedObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = How did you call the corresponding intensity image?
%infotypeVAR04 = imagegroup
%inputtypeVAR04 = popupmenu
IntImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});


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

[~,~,ParentList] = CPrelateobjects(handles,SubObjectName,ParentName,SubObjectLabelMatrix,ParentObjectLabelMatrix,ModuleName);

%%% Since the label matrix starts at zero, we must include this value in
%%% the list to produce a label matrix image with children re-labeled to
%%% their parents values. This does not get saved and is only for display.
if ~isempty(ParentList)
    FinalParentListLM = [0;ParentList];
    NewObjectParentLabelMatrix = FinalParentListLM(SubObjectLabelMatrix+1);
    CurrentObjNhood = bwlabel(NewObjectParentLabelMatrix); 
    CurrentObjNhood = ParentObjectLabelMatrix&CurrentObjNhood;
    NewObjectParentLabelMatrix = zeros(size(ParentObjectLabelMatrix));
    NewObjectParentLabelMatrix(CurrentObjNhood) = ParentObjectLabelMatrix(CurrentObjNhood);
    
else
    NewObjectParentLabelMatrix = SubObjectLabelMatrix;
end

SegmentedObjectLabelMatrix = bwlabel(NewObjectParentLabelMatrix);



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
    subplot(2,2,4);
    CPimagesc(OutlineOverlay,handles);
    title(['Outlines of New Sub Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
    end
    hold off
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE TO HANDLES %%%
%%%%%%%%%%%%%%%%%%%%%%%

%save the object segmentation
fieldname = ['Segmented',SegmentedObjectName];
handles.Pipeline.(fieldname) = SegmentedObjectLabelMatrix;%FinalLabelMatrixImage;

%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end

%%% Saves output subobjects to handles
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,SegmentedObjectName));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {SegmentedObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
elseif ~strcmp(handles.Measurements.Image.ObjectCountFeatures(column),SegmentedObjectName)
    handles.Measurements.Image.ObjectCountFeatures(column) = {SegmentedObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(SegmentedObjectLabelMatrix(:));

%%% Saves the location of each segmented subobject
handles.Measurements.(SegmentedObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(SegmentedObjectLabelMatrix,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(SegmentedObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves parent-child relationship for each segmented object
[handles,~,~] = CPrelateobjects(handles,SegmentedObjectName,ParentName,SegmentedObjectLabelMatrix,ParentObjectLabelMatrix,ModuleName);

% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % this would be an alternative way of doing it
%
% [numChildren, ParentID]= relateobjectsCP3D(NewSubobjectLabelMatrix,ParentObjectLabelMatrix,'Centroid');
% 
% % Save Results
% handles = CPaddmeasurements(handles,ParentName,'Children',SegmentedObjectName,numChildren);
% handles = CPaddmeasurements(handles,SegmentedObjectName,'Parent',ParentName,ParentID);
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%%% Saves the final, segmented label matrix images of objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented', SegmentedObjectName];
handles.Pipeline.(fieldname) = SegmentedObjectLabelMatrix;

end