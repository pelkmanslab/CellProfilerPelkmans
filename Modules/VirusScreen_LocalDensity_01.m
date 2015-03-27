function handles = VirusScreen_LocalDensity_01(handles)
warning off all

% Help for the VirusScreen_Cluster_01 module:
% Category: Other
%
% SHORT DESCRIPTION:
% *************************************************************************
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images you want to process?
%infotypeVAR01 = imagegroup
strImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = Which object would you like to use for the grid density measurement?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Horizontal grid size (48 for non binned, 24 for binned cellWoRx images)
%defaultVAR03 = 24
intXGridSize = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,3})); %#ok Ignore MLint

%textVAR04 = Vertical grid size (52 for non binned, 26 for binned cellWoRx images)
%defaultVAR04 = 26
intYGridSize = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,4})); %#ok Ignore MLint

%%%VariableRevisionNumber = 1



%%% Reads (opens) the image you want to analyze
OrigImage = CPretrieveimage(handles,strImageName,ModuleName,'MustBeGray','CheckScale');
intImageWidth = size(OrigImage,2);
intImageHeight = size(OrigImage,1);

clear OrigImage;

xStep = intXGridSize;
yStep = intYGridSize;

if ~isfield(handles.Measurements.(strObjectName),'GridNucleiCount')
    handles.Measurements.(strObjectName).GridNucleiCount = {};
end
if ~isfield(handles.Measurements.(strObjectName),'GridNucleiEdges')
    handles.Measurements.(strObjectName).GridNucleiEdges = {};
end
if ~isfield(handles.Measurements.(strObjectName),'GridNucleiIntegratedDensity')
    handles.Measurements.(strObjectName).GridNucleiIntegratedDensity = {};
end


% disp(sprintf('after\tLoading\t%.2f',toc))
%%%%%%%%%%%%%%%%%%%%%
%%% GRID FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%

intSetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
NucleiLocations = handles.Measurements.(strObjectName).Location{intSetBeingAnalyzed};

GridNucleiIntegratedDensity = zeros(size(NucleiLocations,1),1);
GridNucleiCount = [];
GridNucleiEdges = [];
GridNucleiNonEdges = [];

GridIndexes = {};
GridLengths = [];
GridBW = [];

%disp(sprintf('ImageWidth %g ImageHeight %g xStep %g yStep %g\n', intImageWidth, intImageHeight, xStep, yStep));

if size(NucleiLocations,1) > 1

    GridNucleiCount = zeros(length(NucleiLocations(:,1)),1);
    GridNucleiEdges = zeros(length(NucleiLocations(:,1)),1);
    GridNucleiNonEdges = ones(length(NucleiLocations(:,1)),1);

    GridX = 0;
    for xPos = 0:xStep:intImageWidth-1 %-1 to skip last class
        GridX = GridX + 1;

        GridY = 0;
        for yPos = 0:yStep:intImageHeight-1 %-1 to skip last class
            GridY = GridY + 1;

            %disp(sprintf('xPos %g yPos %g GridX %g GridY %g\n', xPos, yPos, GridX, GridY));
            
            GridIndexes{GridY,GridX} = find((NucleiLocations(:,1) >= xPos) & (NucleiLocations(:,1) < xPos+xStep) & (NucleiLocations(:,2) >= yPos) & (NucleiLocations(:,2) < yPos+yStep));
            GridLengths(GridY,GridX) = length(GridIndexes{GridY,GridX});

            if GridLengths(GridY,GridX) > 0
                GridBW(GridY,GridX) = 1;
            else
                GridBW(GridY,GridX) = 0;                
            end
            
            GridNucleiCount(GridIndexes{GridY,GridX},1) = GridLengths(GridY,GridX);
        end
    end

% disp(sprintf('after\tGrids \t%.2f',toc))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINDING COLONY EDGES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PaddedGridBW = padarray(GridBW,[1 1],'replicate');
    %find edges, where all 8 surrounding squares are checked for being empty
    PaddedGridEdges = bwperim(PaddedGridBW,8);
    GridEdges = PaddedGridEdges(2:end-1,2:end-1);

    GridNucleiEdges(cell2mat(GridIndexes(GridEdges)),1) = 1;
    GridNucleiNonEdges(cell2mat(GridIndexes(GridEdges)),1) = 0;

% disp(sprintf('after\tEdges. \t%.2f',toc))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SMOOTHING LOCAL DENSITY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%     PaddedGridLengths = padarray(GridLengths,[2 2],'replicate');
%     PaddedGridLengthsResized = imresize(PaddedGridLengths,[(intImageHeight + 4*yStep) (intImageWidth + 4*xStep)],'bicubic');
%     PSF = fspecial('gaussian',[(2*yStep),(2*xStep)],(xStep+yStep)/2);
%     PaddedGridLengthsResizedBlurred = imfilter(PaddedGridLengthsResized,PSF,'conv');   
%     GridLengthsResizedBlurred = PaddedGridLengthsResizedBlurred((2*yStep)+1:end-(2*yStep),(2*xStep)+1:end-(2*xStep));

    smoothingsize = round( sqrt( size(GridLengths,1)^2 + size(GridLengths,2)^2 ) / 10) + 1;

    H = fspecial('disk',smoothingsize); % is 5 optimal?!
    GridLengthsBlurred = imfilter(GridLengths,H,'replicate');    
    GridLengthsBlurredResized = imresize(GridLengthsBlurred,[intImageHeight intImageWidth],'bicubic');

%     [min(GridLengths(:)) max(GridLengths(:))]    

    for i = 1:size(NucleiLocations,1)
        GridNucleiIntegratedDensity(i,1) = GridLengthsBlurredResized(round(NucleiLocations(i,2)),round(NucleiLocations(i,1)));
    end

%     [min(GridNucleiIntegratedDensity(:)) max(GridNucleiIntegratedDensity(:))]
    
    
end %if there are nuclei

% disp(sprintf('after\tSmooth\t%.2f',toc))
%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber) && size(NucleiLocations,1) > 1
    drawnow

    % load image again?! :(
    OrigImage = CPretrieveimage(handles,strImageName,ModuleName,'MustBeGray','CheckScale');
    LabelMatrixImage = CPretrieveimage(handles,['Segmented' strObjectName],ModuleName,'MustBeGray','DontCheckScale');

    Filter = find(GridNucleiNonEdges);
    FinalLabelMatrixImage = LabelMatrixImage;
    for i=1:numel(Filter)
        FinalLabelMatrixImage(FinalLabelMatrixImage == Filter(i)) = 0;
    end

    x = sortrows(unique([LabelMatrixImage(:) FinalLabelMatrixImage(:)],'rows'),1);
    x(x(:,2)>0,2)=1:sum(x(:,2)>0);
    LookUpColumn = x(:,2);

    FinalLabelMatrixImage = LookUpColumn(FinalLabelMatrixImage+1);

    %%% Note: these outlines are not perfectly accurate; for some reason it
    %%% produces more objects than in the original image.  But it is OK for
    %%% display purposes.
    %%% Maximum filters the image with a 3x3 neighborhood.
    MaxFilteredImage = ordfilt2(FinalLabelMatrixImage,9,ones(3,3),'symmetric');
    %%% Determines the outlines.
    IntensityOutlines = FinalLabelMatrixImage - MaxFilteredImage;
    %%% Converts to logical.
    warning off MATLAB:conversionToLogical
    LogicalOutlines = logical(IntensityOutlines);
    warning on MATLAB:conversionToLogical
    %%% Determines the grayscale intensity to use for the cell outlines.
    LineIntensity = max(OrigImage(:));
    %%% Overlays the outlines on the original image.
    ObjectOutlinesOnOrigImage = OrigImage;
    ObjectOutlinesOnOrigImage(LogicalOutlines) = LineIntensity;

    drawnow   
    
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
    end
    %%% A subplot of the figure window is set to display the original
    %%% image.
    subplot(2,3,1);
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the label
    %%% matrix image.
    subplot(2,3,2);
    
    im = CPlabel2rgb(handles,LabelMatrixImage);
    CPimagesc(im,handles);
    title(['Original ',strObjectName]);    

    %%% A subplot of the figure window is set to display the Overlaid image,
    %%% where the maxima are imposed on the inverted original image
    ColoredLabelMatrixImage = CPlabel2rgb(handles,FinalLabelMatrixImage);

    subplot(2,3,5);
    CPimagesc(ColoredLabelMatrixImage,handles);
    title(['Edge ' strObjectName]);
    
    subplot(2,3,4);
    CPimagesc(ObjectOutlinesOnOrigImage,handles);
    title(['Outlines on Input Image']);

    subplot(2,3,3);
    
% BSHACK IN MY OWN MODULE: SMOOTH GRIDLENGTHS

    CPimagesc(GridLengthsBlurredResized,handles);
    title(['Grid Local Density - Blurred & Resized']);

%    CPimagesc(GridLengths,handles);
%    title(['Grid Local Density']);
% *************************************************        
    
    subplot(2,3,6);
    CPimagesc(GridEdges,handles);
    title(['Grid Edges']);

end

% disp(sprintf('after\tDisplay\t%.2f',toc))
%%%%%%%%%%%%%%%%%%%%
%%% STORING DATA %%%
%%%%%%%%%%%%%%%%%%%%

handles.Measurements.(strObjectName).GridNucleiCount{intSetBeingAnalyzed} = GridNucleiCount;
handles.Measurements.(strObjectName).GridNucleiEdges{intSetBeingAnalyzed} = GridNucleiEdges;
handles.Measurements.(strObjectName).GridNucleiIntegratedDensity{intSetBeingAnalyzed} = GridNucleiIntegratedDensity;
