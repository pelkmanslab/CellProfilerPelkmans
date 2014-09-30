function handles = SeparateRegions(handles)

% Help for the Separate Regions module: Category: Object Processing
%
% SHORT DESCRIPTION: 
% Applies the region mask created with the Identify Region module to create
% separate images of presumed tumour and stroma regions.
% *************************************************************************
%
% The module loads an already created region mask and the corresponding
% image of the nuclei channel. It then applies the mask to the image and
% creates separate images for the different regions.

% CellProfiler is distributed under the GNU General Public License. See the
% accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research. Copyright
% 2003,2004,2005.
%
% Author:
%   Markus Herrmann
%
%


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images in which you want to separate
%regions? 
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the images of the region mask? 
%infotypeVAR02 = imagegroup
RegionMaskName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the tumour region image? 
%defaultVAR03 = TumourRegions
%infotypeVAR03 = imagegroup indep
TumourRegionImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What do you want to call the stroma region image? 
%defaultVAR04 = StromaRegions
%infotypeVAR04 = imagegroup indep
StromaRegionImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});


%%%VariableRevisionNumber = 1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the images you want to analyze and assigns them to
%%% variables.
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');
RegionMask = CPretrieveimage(handles,RegionMaskName,ModuleName,'MustBeGray','CheckScale');


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow


%%% creates region images

TumourRegionMask = RegionMask;
StromaRegionMask = ~RegionMask;

% extracts regions form subtracted Dapi image for subsequent segmentation
TumourRegionImage = OrigImage;
TumourRegionImage(StromaRegionMask) = min(TumourRegionImage(:));
StromaRegionImage = OrigImage;
StromaRegionImage(TumourRegionMask) = min(StromaRegionImage(:));


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% detect edges for outline
ImTumourEdge = edge(RegionMask, 'sobel');
ImTumourEdge = bwmorph(ImTumourEdge, 'thicken',10);

%%% creates image that displays region outline in Dapi image
Rd = OrigImage;
Rd(ImTumourEdge)=max(Rd(:)); 
Gd = OrigImage;
Bd = OrigImage;
ImRegionOutline(:,:,1) = Rd;
ImRegionOutline(:,:,2) = Gd;
ImRegionOutline(:,:,3) = Bd;


ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber)
    end
    %%% A subplot of the figure window is set to display the original
    %%% image.
    subplot(2,2,1); 
    CPimagesc(OrigImage,handles); 
    title(['Original Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the outlined
    %%% regions.
    subplot(2,2,2); 
    CPimagesc(ImRegionOutline,handles); 
    title(['Regions Outlined in Original Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the tumour region. 
    subplot(2,2,3); 
    CPimagesc(TumourRegionImage,handles); 
    title('Tumour Regions');
    %%% A subplot of the figure window is set to display the stroma region.
    subplot(2,2,4); 
    CPimagesc(StromaRegionImage,handles); 
    title('Stroma Regions');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the Inverted image to the handles structure so it can be used by
%%% subsequent modules.
fieldname = TumourRegionImageName;
handles.Pipeline.(fieldname) = TumourRegionImage;
fieldname = StromaRegionImageName;
handles.Pipeline.(fieldname) = StromaRegionImage;