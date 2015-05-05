function handles = IlluminationCorrectionZScoring(handles)
% Help for the "ILLUMINATIONCORRECTIONZSCORING" module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Uses std and mean illumination function calculated by iBrain to correct
% uneven illumination/lighting/shading via Z-scoring.
% *************************************************************************
%
% Illumination correction of raw images is essential for subsequent steps
% in the image analysis pipeline. It ensures correct object detection and
% accurate measurements of intensity features, reducing biases due to
% uneven illumination of the sample as well as positional differences in
% the signal gain resulting from the detection system. This illumination
% correction algorithm exploit the statistical power of the large number of
% images acquired per channel to learn pixel-wise illumination and signal
% gain biases. Briefly, This module uses the standard deviation and mean
% intensity values per pixel in a channel calculated by iBrain to correct
% the illumination bias using per-pixel z-scoring. For more information see
% Battich et al. 2013 (Supplementary Info.), and Stoeger et al. 2015.
%
% Authors:
%   Nico Battich
%   Berend Snijder
%   Yauhen Yakimovich
%   Thomas Stoeger
%
% $Revision: 1808 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images to be used to calculate the illumination functions (input)?
%infotypeVAR01 = imagegroup
InputName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = How do you want to call the corrected images (output)?
%defaultVAR02 = CorrBlue
%infotypeVAR02 = imagegroup indep
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,2});


%textVAR03 = Do you also want to smooth mean/std statistics values before using them for correction?
%choiceVAR03 = No
%choiceVAR03 = Yes
DoMeanStdSmoothing = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Smoothing is done by a Gaussian filter to both illumination functions. Enter the size of the filter.
%defaultVAR04 = 50
SmoothingSize = str2double(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Enter the factor 'x' for the sigma calculation of the Gaussian filter. where sigma = x*size.
%defaultVAR05 = 0.5
SmoothingSigma = str2double(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Enter the relative path name to the folder where the illumination correction files are located (starting with "./../"). Type "Pre" (Pre) to load files from previous multiplexing cycle. Type period (.) for default directory.
%defaultVAR06 = .
AlternativeIllCorrFolder = char(handles.Settings.VariableValues{CurrentModuleNum,6});


%%%VariableRevisionNumber = 4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% The illumination correction image was calculated using all the incoming
% images by iBRAIN. Load the mean and std statistics assuming iBrain file
% structure.


% Get the channel of the image.
[intChannelNumber] = check_image_channel(handles.Pipeline.(strcat('FileList',InputName)){handles.Current.SetBeingAnalyzed});
intZstackNumber = 0;
fprintf('Applying illumination correction on channel number: %d\n', intChannelNumber)

% Store the stat in the Measurement field, with channel and szstack
% specific fieldnames.
if strcmp(AlternativeIllCorrFolder,'Pre')
    strStatFieldName = sprintf('illcor_ch%03dz%03d_Pre',intChannelNumber,intZstackNumber);
else
    strStatFieldName = sprintf('illcor_ch%03dz%03d',intChannelNumber,intZstackNumber);
end

if handles.Current.SetBeingAnalyzed == 1
    % Get the BATCH directory and load illcor statistics.
    
    % Get path to folder where files are located
    if strncmp(AlternativeIllCorrFolder,'.',1)
        if length(AlternativeIllCorrFolder) == 1
            VarPathname = handles.Current.DefaultOutputDirectory;
        else
            VarPathname = fullfile(handles.Current.DefaultOutputDirectory,AlternativeIllCorrFolder(2:end));
            fprintf(['=================================================' ...
                '=================================================' ...
                '=================================================' ...
                '\n\n%s: You specified an alternative path:\n%s\n\n' ...
                '=================================================' ...
                '=================================================' ...
                '=================================================' ...
                '\n'],mfilename,VarPathname);
        end
    elseif strcmp(AlternativeIllCorrFolder,'Pre')
        if isfield(handles.Pipeline,['Pathname',InputName])
            VarPathname = fullfile(strrep(handles.Pipeline.(['Pathname',InputName]),'TIFF','BATCH'));
            fprintf(['=================================================' ...
                '=================================================' ...
                '=================================================' ...
                '\n\n%s: You specified an alternative path in the LoadImages module:\n%s\n\n' ...
                '=================================================' ...
                '=================================================' ...
                '=================================================' ...
                '\n'],mfilename,VarPathname);
        else
            error(['Image processing was canceled in the ', ModuleName, ' module because the fieldname "',['Pathname',InputName],'" does not exist within handles.Pipeline. Be sure that you have saved it correctly using the LoadImages module'])
        end
    end
    strBatchDir = VarPathname;
    if ~exist(strBatchDir,'dir')
        error(['Image processing was canceled in the ', ModuleName, ' module because the directory "',strBatchDir,'" does not exist. Be sure that no spaces or unusual characters exist in your typed entry and that the pathname of the directory begins with /.'])
    end
    TempStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',intChannelNumber,intZstackNumber)));
    
    TempStats.stat_values.mean = double(TempStats.stat_values.mean);
    TempStats.stat_values.std = double(TempStats.stat_values.std);
    
    % Optional smoothing of illumination statistics computed by iBRAIN.
    switch DoMeanStdSmoothing
        case 'Yes'
            % Create gaussian filter handle for for smoothing
            H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize*SmoothingSigma);
            handles.Measurements.Image.([strStatFieldName,'_mean']) =  imfilter(TempStats.stat_values.mean,H,'symmetric');
            handles.Measurements.Image.([strStatFieldName,'_std']) = imfilter(TempStats.stat_values.std,H,'symmetric');
        case 'No'
            handles.Measurements.Image.([strStatFieldName,'_mean']) = TempStats.stat_values.mean;
            handles.Measurements.Image.([strStatFieldName,'_std']) = TempStats.stat_values.std;
    end
    clear TempStats;
end


%%%%%%%%%%%%%%%%%%
%%% LOAD IMAGE %%%
%%%%%%%%%%%%%%%%%%
strImageToImport = fullfile( ...
    handles.Pipeline.(strcat('Pathname',InputName)), ...
    handles.Pipeline.(strcat('FileList',InputName)){handles.Current.SetBeingAnalyzed});
OrigImage = double(imread(strImageToImport));


%%%%%%%%%%%%%%%%%%
%%% CORRECTION %%%
%%%%%%%%%%%%%%%%%%
IllumFilt_Mean = handles.Measurements.Image.([strStatFieldName,'_mean']);
IllumFilt_STD = handles.Measurements.Image.([strStatFieldName,'_std']);

% do illumination correction
ImageOutput = IllumCorrect(OrigImage,IllumFilt_Mean,IllumFilt_STD,1);

% store non-scaled for visualization
ImageOutputPlot = ImageOutput;
% Rescale from 0  to 1
ImageOutput = ImageOutput/65535;

% Check if there are pixels with a value higher than one, fix them and give
% warning.
f = ImageOutput(:)>1;
HighPixels = sum(f);
if HighPixels > 0
    fprintf('%s: %.3d pixels with values higher than 1 were found.\nAll these will be set to 1.\n', mfilename, HighPixels)
    ImageOutput(f) = 1;
end

%%%%%%%%%%%%%%%%%%%
%%% SAVE OUTPUT %%%
%%%%%%%%%%%%%%%%%%%

%save to handle structure
handles.Pipeline.(OutputName) = ImageOutput;
fieldname = ['Filename', OutputName];
handles.Pipeline.(fieldname) = handles.Pipeline.(strcat('FileList',InputName))(handles.Current.SetBeingAnalyzed);
fieldname = ['Pathname', OutputName];
handles.Pipeline.(fieldname) = handles.Pipeline.(strcat('Pathname',InputName));


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);

if any(findobj == ThisModuleFigureNumber)
    if CPisHeadless ==  false
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
        subplot(2,2,1);
        
        CPimagesc(IllumFilt_Mean,handles);
        colormap('JET')
        colorbar
        title('Mean Intensity Filter [Log10(intensity)]')
        
        subplot(2,2,3);
        CPimagesc(IllumFilt_STD,handles);
        colormap('JET')
        colorbar
        title('STD Intensity Filter [Log10(intensity)]')
        
        subplot(2,4,3);
        CPimagesc(OrigImage,handles);
        colormap('JET')
        title('Original Image')
        
        subplot(2,4,4);
        CPimagesc(ImageOutputPlot,handles);
        colormap('JET')
        title('Corrected Image')
        
        subplot(2,4,7);
        hold on
        hist(OrigImage(:),round(length(OrigImage(:))/20))
        hold off
        set(gca,'xlim',[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.95)])
        ylabel('Pixel Count')
        xlabel('Intensity')
        title('Original Image Histogram')
        
        subplot(2,4,8);
        hist(ImageOutputPlot(:),round(length(ImageOutputPlot(:))/20))
        set(gca,'xlim',[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.95)])
        title('Corrected Image Histogram')
        ylabel('Pixel Count')
        xlabel('Intensity')
        
        drawnow
    end
end
