function handles = IdentifyRegions(handles)

% Help for the Identify Regions module: Category: Object Processing
%
% SHORT DESCRIPTION: Creates a region mask that separates presumed tumour
% from and stroma regions based on combined information from 2 channels
% (DAPI and cell outline marker).
% *************************************************************************
%
% The module loads images of two corresponding channels, normalizes
% both images, and creates an image overlay. On the resulting overlay image
% it then performs edge detection, closes the image to create connected
% objects, fills holes of connected objects, and opens the image to remove
% small objects.

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

%textVAR01 = What did you call the images of the first channel in which you want to identify regions?
%infotypeVAR01 = imagegroup
ImageNameChannel1 = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the images of the second channel in which you want to identify regions?
%infotypeVAR02 = imagegroup
ImageNameChannel2 = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the region mask? 
%defaultVAR03 = RegionMask 
%infotypeVAR03 = imagegroup indep
RegionMaskImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%%%VariableRevisionNumber = 1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the images you want to analyze and assigns them to
%%% variables.
OrigImageChannel1 = CPretrieveimage(handles,ImageNameChannel1,ModuleName,'MustBeGray','CheckScale');
OrigImageChannel2 = CPretrieveimage(handles,ImageNameChannel2,ModuleName,'MustBeGray','CheckScale');



%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% normalizes images

% normalizes channel 1
lowerQuantileChannel1 = quantile(OrigImageChannel1(:),0.01);
upperQuantileChannel1 = quantile(OrigImageChannel1(:),0.99);
if lowerQuantileChannel1 ~= upperQuantileChannel1
    tempim = OrigImageChannel1;
    tempim(tempim<lowerQuantileChannel1) = lowerQuantileChannel1;
    tempim(tempim>upperQuantileChannel1) = upperQuantileChannel1;
    tempim = tempim - lowerQuantileChannel1;
    tempim = double(tempim);
    ImNormChannel1 = tempim ./ max(tempim(:));
else
    ImNormChannel1 = zeros(size(OrigImageChannel1));
end

% normalizes channel 2
lowerQuantileChannel2 = quantile(OrigImageChannel2(:),0.01);
upperQuantileChannel2 = quantile(OrigImageChannel2(:),0.99);
if lowerQuantileChannel2 ~= upperQuantileChannel2
    tempim = OrigImageChannel2;
    tempim(tempim<lowerQuantileChannel2) = lowerQuantileChannel2;
    tempim(tempim>upperQuantileChannel2) = upperQuantileChannel2;
    tempim = tempim - lowerQuantileChannel2;
    tempim = double(tempim);
    ImNormChannel2 = tempim ./ max(tempim(:));
else
    ImNormChannel2 = zeros(size(OrigImageChannel2));
end
clear tempim

%%% detecs regions and creates region mask

% creates overlay of Dapi and Ecad image
ImOver = ImNormChannel1 + ImNormChannel2;

% creates black-and-white image using threshold (Otsu's)
level = graythresh(ImOver);
ImBwOverlay = im2bw(ImOver, level);

% merges signals (morphological closing)
ImBwClose = bwmorph(ImBwOverlay,'close');

% mirrors borders of the image
padsize = 100;
ImBwClose = padarray(ImBwClose, [padsize,padsize], 'symmetric');

% fills holes in regions
ImBwFill = imfill(ImBwClose, 'holes');
%%%%ImBwCloseFill=bwfillholes(ImBwClose);

% resizes image to original size
ImBwFill = ImBwFill((padsize+1):(end-padsize),(padsize+1):(end-padsize));

% defines tumour region by removing non-merged signals (morphological
% opening)
se2 = strel('disk',20);
RegionMaskImage = imopen(ImBwFill,se2);
RegionMaskImage = bwmorph(RegionMaskImage, 'thicken', 5);



%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% detect edges for outline
ImTumourEdge = edge(RegionMaskImage, 'sobel');
ImTumourEdge = bwmorph(ImTumourEdge, 'thicken',10);

%%% creates image that displays region outline in Dapi image
R1 = ImNormChannel1;
R1(ImTumourEdge)=max(R1(:)); 
G1 = ImNormChannel1;
B1 = ImNormChannel1;
ImRegionOutlineChannel1(:,:,1) = R1;
ImRegionOutlineChannel1(:,:,2) = G1;
ImRegionOutlineChannel1(:,:,3) = B1;

%%% creates image that displays region outline in Ecad image
R2 = ImNormChannel2;
R2(ImTumourEdge)=max(R2(:)); 
G2 = ImNormChannel2;
B2 = ImNormChannel2;
ImRegionOutlineChannel2(:,:,1) = R2;
ImRegionOutlineChannel2(:,:,2) = G2;
ImRegionOutlineChannel2(:,:,3) = B2;

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImageChannel1,'TwoByTwo',ThisModuleFigureNumber)
    end
    %%% A subplot of the figure window is set to display the original image
    %%% of channel 1.
    subplot(2,2,1); 
    CPimagesc(OrigImageChannel1,handles); 
    title(['Channel 1, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the original image
    %%% of channel 2.
    subplot(2,2,2); 
    CPimagesc(OrigImageChannel2,handles); 
    title(['Channel 2, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the region mask
    %%% outlined on the original image of channel 1.
    subplot(2,2,3); 
    CPimagesc(ImRegionOutlineChannel1,handles); 
    title('Outlined Regions Channel 1');
    %%% A subplot of the figure window is set to display the region mask
    %%% outlined on the original image of channel 2.
    subplot(2,2,4); 
    CPimagesc(ImRegionOutlineChannel2,handles); 
    title('Outlined Regions Channel 2');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the Inverted image to the handles structure so it can be used by subsequent modules.
fieldname = RegionMaskImageName;
handles.Pipeline.(fieldname) = RegionMaskImage;