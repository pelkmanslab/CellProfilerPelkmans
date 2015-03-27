function handles = SubtractBackground2(handles)

% Help for the Subtract Background module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Calculates the minimum pixel intensity value for the entire set of images
% and subtracts this value from every pixel in every image.
% *************************************************************************
%
% Note that this is not an illumination correction module. It subtracts a
% single value from every pixel across the image.
%
% The intensity due to camera or illumination or antibody background
% (intensity where no cells are sitting) can in good conscience be
% subtracted from the images, but it must be subtracted from every pixel,
% not just the pixels where cells actually are sitting.  This is because we
% assume that this staining is additive with real staining. This module
% calculates the lowest possible pixel intensity across the entire image
% set and subtracts this background value from every pixel in every image.
% This module is identical to the Apply Threshold module (in shift mode),
% except in the SubtractBackground module, the threshold is automatically
% calculated as the 10th lowest pixel value. This will not push any values
% below zero (therefore, we aren't losing any information). It moves the
% baseline up and looks prettier (improves signal to noise) without any
% 'ethical' concerns.
%
% If images have already been quantified and you want to apply the concept
% of this module without reprocessing your images, then multiply the
% background threshold calculated by this module during the first image
% cycle by the number of pixels in the image to get the number that should
% be subtracted from the intensity measurements.
%
% If you want to run this module only to calculate the proper threshold to
% use, simply run the module as usual and use the button on the Status
% window to stop processing after the first image cycle.
%
% How it works:
% Sort each image's pixel values and pick the 10th lowest pixel value as
% the minimum. Typical images have a million pixels. The lowest pixel value
% is chosen because it might be zero if it is a stuck pixel. It is quite
% certain that there will not be 10 stuck pixels so this should be safe.
% Then, take the minimum of these values from all the images. This scalar
% value should be subtracted from every pixel in the image. CellProfiler is
% not calculating a different value for each pixel position in the image
% because in a small image set, that position may always be occupied by
% real staining.
%
% See also ApplyThreshold.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Authors:
%   Anne E. Carpenter
%   Thouis Ray Jones
%   In Han Kang
%   Ola Friman
%   Steve Lowe
%   Joo Han Chang
%   Colin Clarke
%   Mike Lamprecht
%   Peter Swire
%   Rodrigo Ipince
%   Vicky Lay
%   Jun Liu
%   Chris Gang
%
% Website: http://www.cellprofiler.org
%
% $Revision: 4129 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the image to be corrected?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = Type the text that one type of image has in common (for TEXT options), or their position in each group (for ORDER option):
%defaultVAR02 = DAPI
TextToFind = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What do you want to call the background subtracted image?
%defaultVAR03 = CorrBlue
%infotypeVAR03 = imagegroup indep
BGSImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%%%VariableRevisionNumber = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable and loads the BackgroundImage
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

strBatchDir = handles.Current.DefaultOutputDirectory;
TmpBackground = load(fullfile(strBatchDir,sprintf('BackgroundEstimation_%02s.mat',TextToFind)));
TmpBackground = TmpBackground.Background;
TmpBackground = TmpBackground/65535;
%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% The first time the module is run, the threshold shifting value must be
%%% calculated.
if handles.Current.SetBeingAnalyzed == 1
    try
        drawnow
        %%% Retrieves the path where the images are stored from the handles
        %%% structure.
        fieldname = ['Pathname', ImageName];
        try Pathname = handles.Pipeline.(fieldname);
        catch error(['Image processing was canceled in the ', ModuleName, ' module because it must be run using images straight from a load images module (i.e. the images cannot have been altered by other image processing modules). This is because the Subtract Background module calculates an illumination correction image based on all of the images before correcting each individual image as CellProfiler cycles through them. One solution is to process the entire batch of images using the image analysis modules preceding this module and save the resulting images to the hard drive, then start a new stage of processing from the ', ModuleName,' module onward.'])
        end
        %%% Retrieves the list of filenames where the images are stored from the
        %%% handles structure.
        fieldname = ['FileList', ImageName];
        FileList = handles.Pipeline.(fieldname);
        if size(FileList,1) == 2
            error(['Image processing was canceled in the ', ModuleName, ' module because it cannot function on movies.']);
        end
        
        %%Subtract the Background from the actual image
        BGSImage = OrigImage-TmpBackground;
        BGSImage(BGSImage < 0) = 0;
        
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
        %%% image, some intermediate images, and the final corrected image.
        subplot(2,1,1);
        CPimagesc(OrigImage,handles);
        title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        %%% The mean image does not absolutely have to be present in order to
        %%% carry out the calculations if the illumination image is provided,
        %%% so the following subplot is only shown if MeanImage exists in the
        %%% workspace.
        subplot(2,1,2);
        CPimagesc(BGSImage,handles);
        title('Corrected Image');
        %%% Displays the text.
    end
else BGSImage = OrigImage;
end % This end goes with the if MinimumTenthMinimumPixelValue ~= 0 line above.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the corrected image to the handles structure so it can be used by
%%% subsequent modules.
handles.Pipeline.(BGSImageName) = BGSImage;
