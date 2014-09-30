function handles = DynamicObjectFilter(handles)
% Help for DynamicObjectFilter
% Category: Object Processing
%
% Module to filter objects using interpolated time dependent values. 
%
% The current timepoint is extracted from the filename. Object properties
% may be loaded from the handles structure via the eval field (relevant 
% variable names are given in the labels). Properties could originate from
% Shape, Area or Intensity measurements. This requires knowledge how the
% data is stored in the handles structure. The easiest way to find out is
% to look at the code of the module creating the measurement data.
% Values have to be loaded such to result in a vector containing one value
% for each object.
%
% Example for retrieving the mean object intensity created by
% MeasureObjectIntensity.m:
% (['Intensity_',MeasureImageName]){handles.Current.SetBeingAnalyzed}(:,2)
% Example for retrieving the object area created by MeasureObjectAreaShape.m
% AreaShape{handles.Current.SetBeingAnalyzed}(:,1)
%
% Objects are subsequntly filtered by comparing the property to an
% interpolated value. The value is interpolated between two timepoints, no
% exapolation is performed but the first/last value assumed. Objects are
% relabeled and saved back to the handles structure
% ATTENTION: Measurements are now invalid and have to be reperformed!
% (Object labels have changed and do not match the measurement anymore)
%
% Recomended pipeline useage:
%   -IdentifyPrimary...
%   -MeasureObject...
%   -DynamicObjectFilter
%   -MeasureObject...
%
% [Anatol Schwab 28.11.12]
%
% $Revision: 1809 $


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow
[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images containing the timepoint in the name? ('ImageName')
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the images measurements were created from (if applicable)?('MeasureImageName')
%infotypeVAR02 = imagegroup
MeasureImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What did you call the objects you want to filter? ('ObjectName')
%choiceVAR03 = Do not use
%infotypeVAR03 = objectgroup
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Timepoint identifier in filename
%defaultVAR04 = T
TIdentifier = handles.Settings.VariableValues{CurrentModuleNum,4};

%textVAR05 = First timepoint
%defaultVAR05 = 1
TFirst = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = Last timepoint
%defaultVAR06 = 1
TLast = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,6}));

%textVAR07 = retrieve variable (eval!): handles.Measurements.(ObjectName).*
%defaultVAR07 = (['Intensity_',MeasureImageName]){handles.Current.SetBeingAnalyzed}(:,2)
VarIdentifier=handles.Settings.VariableValues{CurrentModuleNum,7};

%textVAR08 = Compare method
%choiceVAR08 = >
%choiceVAR08 = >=
%choiceVAR08 = ==
%choiceVAR08 = <=
%choiceVAR08 = <
CompareMethod = str2func(handles.Settings.VariableValues{CurrentModuleNum,8});%fancy trick!
%inputtypeVAR08 = popupmenu

%textVAR09 = First value (at first timepoint)
%defaultVAR09 = 0
ValFirst = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,9}));

%textVAR10 = Last value (at last timepoint)
%defaultVAR10 = 1
ValLast = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,10}));

%textVAR11 = Interpolation of values
%choiceVAR11 = linear
InterpolateMethod = handles.Settings.VariableValues{CurrentModuleNum,11};
%inputtypeVAR11 = popupmenu

%%%VariableRevisionNumber = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LabelMatrixImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName,'MustBeGray','DontCheckScale');
ObjectLabels=setdiff(unique(LabelMatrixImage),0);
%retrieve data to compare from handles
%the result should be an array with properties for each object!
CompareValue=[];
EvalStr=strcat('CompareValue=','handles.Measurements.(ObjectName).',VarIdentifier,';');
try
    eval(EvalStr);
catch caughtError
    disp(['can not execute ',EvalStr]);
end
if ~isempty(CompareValue)
    %get timepoint
    ImageFileName=handles.Pipeline.(['Filename' ImageName]){end};
    TimePointMatch=regexp(ImageFileName,[TIdentifier,'\d{1,5}[A-Za-z.]'],'match');
    if length(TimePointMatch)==1&&~isempty(TimePointMatch{1})
        TimePoint=str2num(TimePointMatch{1}(2:end-1));
        TimeFraction=(TimePoint-TFirst)/(TLast-TFirst);%relative position 0 to 1 in interpolation range
    end
    %interpolate
    switch(InterpolateMethod)
        case{'linear'}
            if TimePoint<=TFirst
                InterpolatedValue=ValFirst;
            elseif TimePoint>=TLast
                InterpolatedValue=ValLast;
            else
                InterpolatedValue=ValFirst+(ValLast-ValFirst)*TimeFraction;
            end
    end

    %perform comparisson
    ValidLabels=[];
    for i=1:size(ObjectLabels,1)
        CurrentLabel=ObjectLabels(i);
        if CompareMethod(CompareValue(i),InterpolatedValue)
            ValidLabels=[ValidLabels,CurrentLabel];
        end
    end
    %print results
    drawnow
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber);
        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            CPresizefigure('','NarrowText',ThisModuleFigureNumber)
        end
        currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);
        TextString=['Timepoint: ', num2str(TimePoint),' Interpolated: ',num2str(InterpolatedValue)];
        uicontrol(currentfig,'style','text','units','normalized','fontsize',handles.Preferences.FontSize,'HorizontalAlignment','left','string',TextString,'position',[.05 .85 .95 .1],'BackgroundColor',[.7 .7 .9])
        TextString=['Loaded objects: ',num2str(size(ObjectLabels,1)),' Saved objects: ',num2str(size(ValidLabels,2))];
        uicontrol(currentfig,'style','text','units','normalized','fontsize',handles.Preferences.FontSize,'HorizontalAlignment','left','string',TextString,'position',[.05 .7 .95 .1],'BackgroundColor',[.7 .7 .9])
    end

    NewLabelImage=bwlabel(ismember(LabelMatrixImage,ValidLabels));
    %filter objects, save modified results
    fieldname = ['Segmented',ObjectName];%final label image
    handles.Pipeline.(fieldname) =NewLabelImage;
    
    %%% Saves the location of each segmented object
    handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
    tmp = regionprops(NewLabelImage,'Centroid');
    Centroid = cat(1,tmp.Centroid);
    handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

    %%% Saves the ObjectCount, i.e., the number of segmented objects.
    if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
        handles.Measurements.Image.ObjectCountFeatures = {};
        handles.Measurements.Image.ObjectCount = {};
    end
    column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ObjectName)));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(NewLabelImage(:));

end
end

