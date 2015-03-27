function handles = JoinObjectSegmentations(handles)

% Help for the Relate module:
% Category: Object Processing
%
% This module join two object segmentations,to create new objects. 
%
% NB
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the first objects to join?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
SubObjectName1 = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What did you call the second object to join?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
SubObjectName2 = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What do you want to call the objects identified by this module?
%defaultVAR03 = Organelle
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
SubObjectLabelMatrix1 = CPretrieveimage(handles,['Segmented', SubObjectName1],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
SubObjectLabelMatrix2 = CPretrieveimage(handles,['Segmented', SubObjectName2],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieve intensity image
IntensityImage = CPretrieveimage(handles,IntImageName,ModuleName,'MustBeGray','DontCheckScale');


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%join segmentation
NewObjectMatrix = SubObjectLabelMatrix1>0 | SubObjectLabelMatrix2>0;
NewObjectMatrix = imfill(double(NewObjectMatrix),'holes');
NewSubobjectLabelMatrix = bwlabel(NewObjectMatrix);

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%

drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(NewSubobjectLabelMatrix,'TwoByTwo',ThisModuleFigureNumber);
    end
    subplot(2,2,1);
    ColoredLabelMatrixImage1 = CPlabel2rgb(handles,SubObjectLabelMatrix1);
    CPimagesc(ColoredLabelMatrixImage1,handles);
    title(['Original Sub Objects 1, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,2,2);
    ColoredLabelMatrixImage2 = CPlabel2rgb(handles,SubObjectLabelMatrix2);
    CPimagesc(ColoredLabelMatrixImage2,handles);
    title('Original Sub Objects 2');
    subplot(2,2,3);
    ColoredNewSubobjectLabel = CPlabel2rgb(handles,NewSubobjectLabelMatrix);
    CPimagesc(ColoredNewSubobjectLabel,handles);
    title('New Objects');
    OutlineOverlay = IntensityImage;
    B = bwboundaries(NewSubobjectLabelMatrix,'noholes');
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
handles.Pipeline.(fieldname) = NewSubobjectLabelMatrix;%FinalLabelMatrixImage;

%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end

% ======= hack 2014/02/17 [MH] ==============================================================================================
% Use of 'strfind' is very dangerous, since it is not specific and finds any occurance (also substring)! I thus
% replaced it by 'strcmp'.

% save output subobjects to handles
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,SegmentedObjectName));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {SegmentedObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
elseif ~strcmp(handles.Measurements.Image.ObjectCountFeatures(column),SegmentedObjectName)
    handles.Measurements.Image.ObjectCountFeatures(column) = {SegmentedObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(NewSubobjectLabelMatrix(:));
% ======= hack end ==========================================================================================================


%%% Saves the location of each segmented object
handles.Measurements.(SegmentedObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(NewSubobjectLabelMatrix,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(SegmentedObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves the final, segmented label matrix images of objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented', SegmentedObjectName];
handles.Pipeline.(fieldname) = NewSubobjectLabelMatrix;

%%% Save the SmallRemovedSegmented image in pipeline to allow usage of
%%% identify secondary
fieldname = ['SmallRemovedSegmented', SegmentedObjectName];
handles.Pipeline.(fieldname) = NewSubobjectLabelMatrix;







end