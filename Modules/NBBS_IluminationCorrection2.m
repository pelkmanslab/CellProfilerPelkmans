function handles = NBBS_IluminationCorrection2(handles)

% Help for the NBBS_IluminationCorrection_Calc module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Calculates a mean illumination function and the std of the illumination
% mean, used to correct uneven illumination/lighting/shading via Z-scoring.
% *************************************************************************
%
% This module calculates the mean and std values for the intensity of each
% pixel position in a set of images from the same wave length. Such
% fuctions can later be used to to correct the intensity values of each
% pixel by substracting the mean function and normalizing with the std
% function. Such proceedure sould not only normalize for the background
% intensities but also for amplitude differences whic result from optical
% aberrances or differences in illumination.
% Note: the images for this module should be directed loaded by the
% LoadImages module. Alternatively a precalculated illumination correction
% files (from iBrain) can be placed in the default output folder.
%
% Authors:
%   Nico Battich
%   Berend Snijder
%
% $Revision: 1808 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Have the Illumination functions been calculated by iBrain and are they in the BATCH directory?
%choiceVAR01 = Yes
%choiceVAR01 = No
iBrainTakeOver = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu


%textVAR02 = What did you call the images to be used to calculate the illumination functions?
%infotypeVAR02 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the mean illumination function?
%defaultVAR03 = MeanIllumBlue
%infotypeVAR03 = imagegroup indep
Mean_IlluminationImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What do you want to call the std illumination function?
%defaultVAR04 = StdIllumBlue
%infotypeVAR04 = imagegroup indep
Std_IlluminationImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Smoothing is done by a Gaussian filter to both illumination functions. Enter the size of the filter.
%defaultVAR05 = 100
SmoothingSize = str2num(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Enter the factor 'x' for the sigma calculation of the Gaussian filter. where sigma = x*size.
%defaultVAR06 = 0.5
SmoothingSigma = str2num(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Enter the number of iterations runned to approximate the std value per pixel.
%defaultVAR07 = 20
IterationNum = str2num(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Do you want to apply the illumination correction?
%choiceVAR08 = Yes
%choiceVAR08 = No
AplyIllumCorr = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 = Do you also want to apply background subtraction? (iBRAIN only)
%choiceVAR09 = Yes
%choiceVAR09 = No
AplyBackgroundSubtraction = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = Optionally enter a background-subtraction correction factor (<1 = less subtraction, >1 = more subtraction
%defaultVAR10 = 1
BackgroundCorrectionCorrectionFactor = str2num(handles.Settings.VariableValues{CurrentModuleNum,10}); %#ok<*ST2NM>

%textVAR11 = How do you want to call the corrected images?
%defaultVAR11 = CorrBlue
%infotypeVAR11 = imagegroup indep
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,11});

%textVAR12 = Do you also want to smoothing the iBrain derived correction function?
%choiceVAR12 = No
%choiceVAR12 = Yes
iDoBrainSmoothing = char(handles.Settings.VariableValues{CurrentModuleNum,12});
%inputtypeVAR12 = popupmenu

%%%VariableRevisionNumber = 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow


% [TS] Initialize variables, which are independently iBrainTakeOver
if handles.Current.SetBeingAnalyzed == 1
    H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize*SmoothingSigma); % Filter for smoothing
end

%%% The illumination correction image is calculated using all the incoming
%%% images. This module is therefore only runned in the first cycle. This
%%% also allows for batch processing

if  ~strncmpi(iBrainTakeOver,'y',1)
    if handles.Current.SetBeingAnalyzed == 1
        

        %get a list of all images to be used
        % [TS1, 2012-07-07: reorganized code for initializing images so that only
        % essential information is created. Will bypass slowest step of module and reduce time per call of module by 1-2 min]
        ImageList = cellfun(@(x) fullfile(handles.Pipeline.(strcat('Pathname',ImageName)),x),...
        handles.Pipeline.(strcat('FileList',ImageName)), 'uniformoutput',false)';
        
        % Initialize the Mean Image
        LogExampleImage = double(imread(ImageList{1}));

        % Initialize the STD Image
        LogStdPixel = nan(length(LogExampleImage(:)),IterationNum);
        LogMeanPixel = nan(length(LogExampleImage(:)),IterationNum);
        NumImagesToLoad = 25;
        if length(ImageList)<25*2
            NumImagesToLoad = floor(length(ImageList)/2);
        end

        for i = 1:IterationNum
            fprintf('%d out of %d... \n',i,IterationNum)
            tic
            matRandIndx = randperm(length(ImageList));
            matRandIndx = matRandIndx(1:NumImagesToLoad)';
            cellTempImages = arrayfun(@(x) double(imread(ImageList{x})),matRandIndx,'uniformoutput',false);
            %cellTempImages = arrayfun(@(x) CPimread(ImageList{x},handles),matRandIndx,'uniformoutput',false);
            %cellTempImages = arrayfun(@(x) x(x == 0) = 1;,matRandIndx,'uniformoutput',false);
            cellTempImages = cellfun(@log10,cellTempImages,'uniformoutput',false);
            matTempImages = cell2mat(cellfun(@(x) x(:),cellTempImages,'uniformoutput',false)');
            matTempImages(matTempImages == -inf) = 0;
            LogStdPixel(:,i) = nanstd(matTempImages')';
            LogMeanPixel(:,i) = nanmean(matTempImages')';
            toc
        end


        LogStdPixel_Mean = nanmean(LogStdPixel')';
        LogStdImage = reshape(LogStdPixel_Mean,size(LogExampleImage,1),size(LogExampleImage,2));

        LogMeanPixel_Mean = nanmean(LogMeanPixel')';
        LogMeanImage = reshape(LogMeanPixel_Mean,size(LogExampleImage,1),size(LogExampleImage,2));


        %[TS] deactivated -> calculate prior to make accesible to iBrain
        %loaded H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize*SmoothingSigma);

        IllumFilt_STD = imfilter(LogStdImage,H,'replicate');
        IllumFilt_Mean = imfilter(LogMeanImage,H,'replicate');
        BackgroundSubtractionVal = (10^mean(IllumFilt_Mean(:))) * BackgroundCorrectionCorrectionFactor;
        
        %save to handle structure
        handles.Pipeline.(Mean_IlluminationImageName) = (10.^IllumFilt_Mean)/65535;
        handles.Pipeline.(Std_IlluminationImageName) = (10.^IllumFilt_STD)/65535;
        handles.Pipeline.(strcat(Mean_IlluminationImageName,'_init')) = IllumFilt_Mean;
        handles.Pipeline.(strcat(Std_IlluminationImageName,'_init')) = IllumFilt_STD;
        handles.Pipeline.(strcat(Mean_IlluminationImageName,'_background')) = BackgroundSubtractionVal;

    else
        IllumFilt_STD = handles.Pipeline.(strcat(Std_IlluminationImageName,'_init'));
        IllumFilt_Mean = handles.Pipeline.(strcat(Mean_IlluminationImageName,'_init'));
        BackgroundSubtractionVal = handles.Pipeline.(strcat(Mean_IlluminationImageName,'_background'));
    end
    
else
    
    % load the illum functions assuming iBrain file structure % note: load
    % the functions every cycle because of z-stacks which might make the batch-data file too big.
    
    %get the channel of the image
    [intChannelNumber] = check_image_channel(handles.Pipeline.(strcat('FileList',ImageName)){handles.Current.SetBeingAnalyzed})
    intZstackNumber = 0;

    % [BS] store the stat in the Measurement field, with channel and
    % szstack specific fieldnames 
    strStatFieldName = sprintf('illcor_ch%03dz%03d',intChannelNumber,intZstackNumber);
    
    if handles.Current.SetBeingAnalyzed == 1
        %get the BATCH directory and load illcor data
        strBatchDir = handles.Current.DefaultOutputDirectory;
        TempStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',intChannelNumber,intZstackNumber)));
        
        %[TS] adjusted to include option for smoothing iBrain prepared correction function
        switch iDoBrainSmoothing
            case 'Yes'
                ValueMean = imfilter(double(TempStats.stat_values.mean),H,'replicate');
                ValueStd = imfilter(double(TempStats.stat_values.std),H,'replicate');
            case 'No'
                ValueMean = double(TempStats.stat_values.mean);
                ValueStd = double(TempStats.stat_values.std);
        end
        
        handles.Measurements.Image.([strStatFieldName,'_mean']) = ValueMean;
        handles.Measurements.Image.([strStatFieldName,'_std']) = ValueStd;
        clear TempStats
    end
    
    IllumFilt_Mean = handles.Measurements.Image.([strStatFieldName,'_mean']);
    IllumFilt_STD = handles.Measurements.Image.([strStatFieldName,'_std']);
    BackgroundSubtractionVal = (10^mean(IllumFilt_Mean(:))) * BackgroundCorrectionCorrectionFactor;

end


%do correction
if strncmpi(AplyIllumCorr,'y',1)
    %    [handles, ImageOutputPlot, OrigImage] = DoCorrectionAndSave(handles,ImageList,IllumFilt_Mean,IllumFilt_STD,OutputName);
    % OrigImage = double(imread(ImageList{handles.Current.SetBeingAnalyzed})); %Deactivated by TS
  
    % Replacement by  [TS], see comment [TS1], allows to speed up module
    strImageToImport = fullfile(...
        handles.Pipeline.(strcat('Pathname',ImageName)),...
        handles.Pipeline.(strcat('FileList',ImageName)){handles.Current.SetBeingAnalyzed});
    OrigImage  = double(imread(strImageToImport));
    % End of replacement
    
    OrigImage(OrigImage == 0) = 1;
    ImageOutput = (log10(OrigImage)-IllumFilt_Mean)./IllumFilt_STD;
    ImageOutput = (ImageOutput.*mean(IllumFilt_STD(:)))+mean(IllumFilt_Mean(:));
    ImageOutput = 10.^ImageOutput;
    
    if strncmpi(AplyBackgroundSubtraction,'y',1)
        % [BS, 2012-06-25] Added optional background subtraction
        ImageOutput = ImageOutput-BackgroundSubtractionVal;
        ImageOutput(ImageOutput<0) = 0;
    end    
    
    % store non-scaled for visualization
    ImageOutputPlot = ImageOutput;    
    % rescale 0 / 1
    ImageOutput = ImageOutput/65535;
    
    % fix potentially bad pixels (eg. dead pixels: std of 0)
    ImageOutput = fixNonNumericalValueInImage(ImageOutput);
    
    %save to handle structure
    handles.Pipeline.(OutputName) = ImageOutput;
    
    
end

%display the results
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
CPfigure(handles,'Image',ThisModuleFigureNumber);
CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
subplot(2,2,1);
try
    CPimagesc(IllumFilt_Mean,handles);
    colormap('JET')
    colorbar
    title('Mean Intensity Filter [Log10(intensity)]')
end
subplot(2,2,3);
try
    CPimagesc(IllumFilt_STD,handles);
    colormap('JET')
    colorbar
    title('STD Intensity Filter [Log10(intensity)]')
end

subplot(2,4,3);
try
    CPimagesc(OrigImage,handles);
    colormap('JET')
    title('Original Image')
end
subplot(2,4,4);
try
    CPimagesc(ImageOutputPlot,handles);
    colormap('JET')
    if strncmpi(AplyBackgroundSubtraction,'y',1)
        title(sprintf('Corrected Image\n+Background Subtraction'),'FontSize',7)
    else
        title('Corrected Image')
    end
end

subplot(2,4,7);
try
    hold on
    hist(OrigImage(:),round(length(OrigImage(:))/20))
    vline(BackgroundSubtractionVal,'--r')
    hold off
    if ~strncmpi(AplyBackgroundSubtraction,'y',1)
        set(gca,'xlim',[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.95)])
    end
    ylabel('Pixel Count')
    xlabel('Intensity')
    title('Original Image Histogram')
end
subplot(2,4,8);
try
    if ~strncmpi(AplyBackgroundSubtraction,'y',1)
        hist(ImageOutputPlot(:),round(length(ImageOutputPlot(:))/20))
        set(gca,'xlim',[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.95)])
        title('Corrected Image Histogram')
    else
        hist(ImageOutputPlot(ImageOutputPlot>0),round(length(ImageOutputPlot(ImageOutputPlot>0))/20))
        title('Corrected Image Histogram (values > 0)')
    end
    ylabel('Pixel Count')
    xlabel('Intensity')
    
end
drawnow


%
% function [handles, ImageOutputPlot, OrigImage] = DoCorrectionAndSave(handles,ImageList,IllumFilt_Mean,IllumFilt_STD,OutputName)
%
% % do the image correction
%
% OrigImage = double(imread(ImageList{handles.Current.SetBeingAnalyzed}));
% OrigImage(OrigImage == 0) = 1;
% ImageOutput = (log10(OrigImage)-IllumFilt_Mean)./IllumFilt_STD;
% ImageOutput = (ImageOutput.*mean(IllumFilt_STD(:)))+mean(IllumFilt_Mean(:));
% ImageOutput = 10.^ImageOutput;
% ImageOutputPlot = ImageOutput;
% ImageOutput = ImageOutput/65535;
%
% %save to handle structure
% handles.Pipeline.(OutputName) = ImageOutput;
