function handles = MeasureObjectIntensityPercentiles(handles)

% Help for the Measure ObjectIntensity Percentiles module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Measures several intensity features for identified objects.
% *************************************************************************
%
% Given an image with objects identified (e.g. nuclei or cells), this
% module extracts intensity features for each object based on a
% corresponding grayscale image. Measurements are recorded for each object.
%
% Measurement:             Feature Number:
% MedianIntensity         |       1
% 25PercentileIntensity   |       2
% 75PercentileIntensity   |       3
% 5PercentileIntensity    |       4
% 95PercentileIntensity   |       5
%
% How it works:
% Is derived from original CP module Measure Object Intensity. Except of
% using mean, min and max, robust intensity measurements (Median,
% percentiles) are computed
%
% Website: http://www.cellprofiler.org
%
% $Revision: 4526 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the greyscale images you want to measure?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the objects that you want to measure?
%choiceVAR02 = Do not use
%infotypeVAR02 = objectgroup
ObjectNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = Type "Do not use" in unused boxes.
%choiceVAR03 = Do not use
%infotypeVAR03 = objectgroup
ObjectNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 =
%choiceVAR04 = Do not use
%infotypeVAR04 = objectgroup
ObjectNameList{3} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 =
%choiceVAR05 = Do not use
%infotypeVAR05 = objectgroup
ObjectNameList{4} = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu

%textVAR06 =
%choiceVAR06 = Do not use
%infotypeVAR06 = objectgroup
ObjectNameList{5} = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%textVAR07 =
%choiceVAR07 = Do not use
%infotypeVAR07 = objectgroup
ObjectNameList{6} = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu

%%%VariableRevisionNumber = 2

%%% Set up the window for displaying the results
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber);
    CPfigure(handles,'Text',ThisModuleFigureNumber);
    columns = 1;
end

%%% START LOOP THROUGH ALL THE OBJECTS
for i = 1:length(ObjectNameList)
    ObjectName = ObjectNameList{i};
    if strcmpi(ObjectName,'Do not use')
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    drawnow
    
    %%% Reads (opens) the image you want to analyze and assigns it to a variable,
    %%% "OrigImage".
    OrigImage = CPretrieveimage(handles,ImageName,ModuleName);
    
    %%% If the image is three dimensional (i.e. color), the three channels
    %%% are added together in order to measure object intensity.
    if ~ismatrix(OrigImage)        
        error(['Image processing was canceled in the ', ModuleName, ' module. There was a problem with the dimensions. The image must be grayscale.'])
    end
    
    %%% Retrieves the label matrix image that contains the segmented objects which
    %%% will be measured with this module.
    LabelMatrixImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName,'MustBeGray','DontCheckScale');
    
    %%% For the cases where the label matrix was produced from a cropped
    %%% image, the sizes of the images will not be equal. So, we crop the
    %%% LabelMatrix and try again to see if the matrices are then the
    %%% proper size. Removes Rows and Columns that are completely blank.
    if any(size(OrigImage) < size(LabelMatrixImage))
        ColumnTotals = sum(LabelMatrixImage,1);
        RowTotals = sum(LabelMatrixImage,2)';
        warning off all
        ColumnsToDelete = ~logical(ColumnTotals);
        RowsToDelete = ~logical(RowTotals);
        warning on all
        drawnow
        CroppedLabelMatrix = LabelMatrixImage;
        CroppedLabelMatrix(:,ColumnsToDelete,:) = [];
        CroppedLabelMatrix(RowsToDelete,:,:) = [];
        clear LabelMatrixImage
        LabelMatrixImage = CroppedLabelMatrix;
        %%% In case the entire image has been cropped away, we store a single
        %%% zero pixel for the variable.
        if isempty(LabelMatrixImage)
            LabelMatrixImage = 0;
        end
    end
    
    if any(size(OrigImage) ~= size(LabelMatrixImage))
        error(['Image processing was canceled in the ', ModuleName, ' module. The size of the image you want to measure is not the same as the size of the image from which the ',ObjectName,' objects were identified.'])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MAKE MEASUREMENTS & SAVE TO HANDLES STRUCTURE %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    drawnow
    
    numMeasurements = 5;
    
    %%% Initialize measurement structure
    BasicFeatures    = {'MedianIntensity',...
                        '25pIntensity',...
                        '75pIntensity',...
                        '05pIntensity',...
                        '95pIntensity'};
    
    %%% Get pixel indexes (fastest way), and count objects
    props = regionprops(LabelMatrixImage,'PixelIdxList');
    ObjectCount = length(props);
    if ObjectCount > 0        
        Basic = zeros(ObjectCount,numMeasurements);       % note that original CP module initializes insufficiently
        for Object = 1:ObjectCount            
            %%% Measure basic set of Intensity features
            Basic(Object,1) = median(OrigImage(props(Object).PixelIdxList));
            Basic(Object,2) = quantile(OrigImage(props(Object).PixelIdxList),0.25);
            Basic(Object,3) = quantile(OrigImage(props(Object).PixelIdxList),0.75);
            Basic(Object,4) = quantile(OrigImage(props(Object).PixelIdxList),0.05);
            Basic(Object,5) = quantile(OrigImage(props(Object).PixelIdxList),0.95);
        end
    else
        Basic(1,1:numMeasurements) = 0;
    end
    %%% Save measurements
    handles.Measurements.(ObjectName).(['RobustIntensity_',ImageName,'Features']) = BasicFeatures;
    handles.Measurements.(ObjectName).(['RobustIntensity_',ImageName])(handles.Current.SetBeingAnalyzed) = {Basic};
    
    %%% Report measurements
    if any(findobj == ThisModuleFigureNumber);
        FontSize = handles.Preferences.FontSize;
        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            delete(findobj('parent',ThisModuleFigureNumber,'string','R'));
            delete(findobj('parent',ThisModuleFigureNumber,'string','G'));
            delete(findobj('parent',ThisModuleFigureNumber,'string','B'));
        end
        %%%% This first block writes the same text several times
        %%% Header
        
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0 0.95 1 0.04],...
            'HorizontalAlignment','center','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'fontweight','bold','string',sprintf(['Average intensity features for ', ImageName,', cycle #%d'],handles.Current.SetBeingAnalyzed));
        
        %%% Number of objects
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.05 0.85 0.3 0.03],...
            'HorizontalAlignment','left','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'fontweight','bold','string','Number of objects:');
        
        %%% Text for Basic features
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.05 0.8 0.3 0.03],...
            'HorizontalAlignment','left','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'fontweight','bold','string','Intensity feature:');
        for k = 1:numMeasurements
            uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.05 0.8-0.04*k 0.3 0.03],...
                'HorizontalAlignment','left','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
                'fontsize',FontSize,'string',BasicFeatures{k});
        end
        
        %%% The name of the object image
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.35+0.1*(columns-1) 0.9 0.1 0.03],...
            'HorizontalAlignment','center','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'fontweight','bold','string',ObjectName);
        
        %%% Number of objects
        uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.35+0.1*(columns-1) 0.85 0.1 0.03],...
            'HorizontalAlignment','center','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
            'fontsize',FontSize,'string',num2str(ObjectCount));
        
        if ObjectCount > 0
            %%% Basic features
            for k = 1:numMeasurements
                uicontrol(ThisModuleFigureNumber,'style','text','units','normalized', 'position', [0.35+0.1*(columns-1) 0.8-0.04*k 0.1 0.03],...
                    'HorizontalAlignment','center','BackgroundColor',[.7 .7 .9],'fontname','Helvetica',...
                    'fontsize',FontSize,'string',sprintf('%0.5f',mean(Basic(:,k))));
            end
        end
        %%% This variable is used to write results in the correct column
        %%% and to determine the correct window size
        columns = columns + 1;
    end
end