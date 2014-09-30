function handles = PresegmentNuclei(handles)

% Help for the Presegmentation Nuclei module: Category: Object Processing
%
% SHORT DESCRIPTION: Makes use of information derived from a second channel
% (cell outline marker) in the first channel (nuclear marker) to facilitate
% subesequent segmentation of the first channel.
% *************************************************************************
%
% The module loads images of two channels (nuclear marker) and a
% corresponding second channel (cell outline marker), normalizes both
% images, and multiplies the image of channel one with the inversed image
% of channel two.

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

%textVAR01 = What did you call the images of channel 1? 
%infotypeVAR01 = imagegroup
ImageNameChannel1 = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the images of channel 2? 
%infotypeVAR02 = imagegroup
ImageNameChannel2 = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the presegmented image? 
%defaultVAR03 = PresegBlue
%infotypeVAR03 = imagegroup indep
PresegmentedImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

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
    tempim=OrigImageChannel2;
    tempim(tempim<lowerQuantileChannel2) = lowerQuantileChannel2;
    tempim(tempim>upperQuantileChannel2) = upperQuantileChannel2;
    tempim = tempim - lowerQuantileChannel2;
    tempim = double(tempim);
    ImNormChannel2 = tempim ./ max(tempim(:));
else
    ImNormChannel2 = zeros(size(OrigImageChannel2));
end
clear tempim

%%% multiplies the normalized image of channel 1 with the inversed,
%%% normalized image of channel 2)
PresegmentedImage = ImNormChannel1.*(1-ImNormChannel2);



%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow


ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImageChannel1,'TwoByTwo',ThisModuleFigureNumber)
    end
    %%% A subplot of the figure window is set to display the original
    %%% image.
    subplot(2,2,1); 
    CPimagesc(OrigImageChannel1,handles); 
    title(['Channel 1, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the Inverted
    %%% Image.
    subplot(2,2,2); 
    CPimagesc(OrigImageChannel2,handles); 
    title(['Channel 2, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the original
    %%% image.
    subplot(2,2,3); 
    CPimagesc(PresegmentedImage,handles); 
    title('Presegmented Channel 1');
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the Inverted image to the handles structure so it can be used by subsequent modules.
fieldname = PresegmentedImageName;
handles.Pipeline.(fieldname) = PresegmentedImage;