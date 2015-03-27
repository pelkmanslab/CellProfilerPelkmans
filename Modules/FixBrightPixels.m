function handles = FixBrightPixels(handles)

% Help for the FixBrightPixels module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Removes extremly bright pixels by substituting them with a local median.
% This prevents bright pixels from yokogawa to suppress spot detection in
% adjacent regions
% *************************************************************************
%
% Threshold intensity (note that this value will be divided by 65535 as in
% ordinary CP modules)
%
% This module will use local median to substitute individual pixels
%


% Website: http://www.cellprofiler.org
%
% $Revision: 4433 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the image where bright pixels should be removed?
%infotypeVAR01 = imagegroup
OrigImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the image without these bright pixels?
%defaultVAR02 = FixedGreen
%infotypeVAR02 = imagegroup indep
FixedImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What is the highest allowed intensity (will be divided by 65535)
%defaultVAR03 = 9000
PixelIntensity = str2double(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Pixels within this distance will be used for fixing.
%defaultVAR04 = 1
PixelDistance = str2double(handles.Settings.VariableValues{CurrentModuleNum,4});

%%%VariableRevisionNumber = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow


%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
try
    OrigImage = CPretrieveimage(handles,OrigImageName,ModuleName,'MustBeGray','CheckScale');
catch
    ErrorMessage = lasterr;
    error(['Image processing was canceled in the ' ModuleName ' module because: ' ErrorMessage(33:end)]);
end

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Creates Threshold for intensity 
PixelIntensity = PixelIntensity ./65535; % note that 65535 is hardcoded at multiple positions in CP and raw image intensities are rescaled by this value

% get Pixels above Threhold
f = OrigImage(:)>PixelIntensity;
[PixelRow PixelColumn] =ind2sub(size(OrigImage), find(f));

% create Output
OutputImage = OrigImage;

for k=1:size(PixelColumn,1)
    N = max([(PixelRow(k)-PixelDistance) 1]);
    S = min([(PixelRow(k)+PixelDistance) size(OrigImage,1)]);
    W = max([(PixelColumn(k)-PixelDistance) 1]);
    E = min([(PixelColumn(k)+PixelDistance) size(OrigImage,2)]);

    Pix = OrigImage(N:S,W:E);
    Pix = sort(Pix);
    Pix = Pix(1:(end-1)); % remove highest pixel 
    MedianPix = nanmedian(Pix);  % get median value of neighbours
    
    OutputImage(PixelRow(k),PixelColumn(k)) = MedianPix;
    
end



%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByOne',ThisModuleFigureNumber)
    end
    %%% A subplot of the figure window is set to display the original
    %%% image and the smoothed image.
    subplot(2,1,1);
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,1,2);
    CPimagesc(OutputImage,handles);
    title(['Fixed Image. Substituted ' num2str(sum(f)) ' pixels']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the processed image to the handles structure so it can be used by
%%% subsequent modules.
handles.Pipeline.(FixedImageName) = OutputImage;
