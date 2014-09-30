function handles = NBBS_IluminationCorrection(handles)

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
% Note: the images for tis module sould be directed loaded by the
% LoadImages module
%
% Authors:
%   Nico Battich
%   Berend Snijder
%
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

%textVAR08 = Do you wnat to apply the illumination correction?
%choiceVAR08 = Yes
%choiceVAR08 = No
AplyIllumCorr = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 = If Yes, how do you want to call the corrected images?
%defaultVAR09 = CorrBlue
%infotypeVAR09 = imagegroup indep
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,9});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%get a list of all images to be used
ImageList = cellfun(@(x) fullfile(handles.Pipeline.(strcat('Pathname',ImageName)),x),...
    handles.Pipeline.(strcat('FileList',ImageName)), 'uniformoutput',false)';



%%% The illumination correction image is calculated using all the incoming
%%% images. This module is therefore only runned in the first cycle. This
%%% also allows for batch processing

if  ~strncmpi(iBrainTakeOver,'y',1)
    if handles.Current.SetBeingAnalyzed == 1
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


        H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize*SmoothingSigma);

        IllumFilt_STD = imfilter(LogStdImage,H,'replicate');
        IllumFilt_Mean = imfilter(LogMeanImage,H,'replicate');

        %save to handle structure
        handles.Pipeline.(Mean_IlluminationImageName) = (10.^IllumFilt_Mean)/65535;
        handles.Pipeline.(Std_IlluminationImageName) = (10.^IllumFilt_STD)/65535;
        handles.Pipeline.(strcat(Mean_IlluminationImageName,'_init')) = IllumFilt_Mean;
        handles.Pipeline.(strcat(Std_IlluminationImageName,'_init')) = IllumFilt_STD;

    else
        
        IllumFilt_STD = handles.Pipeline.(strcat(Std_IlluminationImageName,'_init'));
        IllumFilt_Mean = handles.Pipeline.(strcat(Mean_IlluminationImageName,'_init'));
    end
    
else
    
    % load the illum functions assuming iBrain file structure % note: load
    % the functions every cycle because of z-stacks which might make the batch-data file too big.
    
    %get the channel of the image
    [intChannelNumber] = check_image_channel(handles.Pipeline.(strcat('FileList',ImageName)){handles.Current.SetBeingAnalyzed});
    intZstackNumber = 0;
    
    %get the BATCH directory
    strBatchDir = handles.Current.DefaultOutputDirectory;
    
    %get the images that match the Z-stack and the channel
    TempStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',intChannelNumber,intZstackNumber)));
    
    IllumFilt_Mean = TempStats.stat_values.mean;
    IllumFilt_STD = TempStats.stat_values.std;

end

%do correction
if strncmpi(AplyIllumCorr,'y',1)
    %    [handles, ImageOutputPlot, OrigImage] = DoCorrectionAndSave(handles,ImageList,IllumFilt_Mean,IllumFilt_STD,OutputName);
    OrigImage = double(imread(ImageList{handles.Current.SetBeingAnalyzed}));
    OrigImage(OrigImage == 0) = 1;
    ImageOutput = (log10(OrigImage)-IllumFilt_Mean)./IllumFilt_STD;
    ImageOutput = (ImageOutput.*mean(IllumFilt_STD(:)))+mean(IllumFilt_Mean(:));
    ImageOutput = 10.^ImageOutput;
    ImageOutputPlot = ImageOutput;
    ImageOutput = ImageOutput/65535;
    
    %save to handle structure
    handles.Pipeline.(OutputName) = ImageOutput;
    
    
end

%display the results
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
CPfigure(handles,'Image',ThisModuleFigureNumber);
CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
subplot(2,2,1);
imagesc(IllumFilt_Mean);
colormap('JET')
colorbar
title('Mean Intensity Filter [Log10(intensity)]')
subplot(2,2,3);
imagesc(IllumFilt_STD);
colormap('JET')
colorbar
title('STD Intensity Filter [Log10(intensity)]')


subplot(2,4,3);
imagesc(OrigImage,[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.999)]);
colormap('JET')
title('Original Image')
subplot(2,4,4);
imagesc(ImageOutputPlot,[quantile(ImageOutputPlot(:), 0.001) quantile(ImageOutputPlot(:), 0.999)]);
colormap('JET')
title('Corrected Image')

subplot(2,4,7);
hist(OrigImage(:),round(length(OrigImage(:))/20))
set(gca,'xlim',[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.95)])
ylabel('Pixel Count')
xlabel('Intensity')
title('Original Image Histogram')
subplot(2,4,8);
hist(ImageOutputPlot(:),round(length(ImageOutputPlot(:))/20))
set(gca,'xlim',[quantile(OrigImage(:), 0.001) quantile(OrigImage(:), 0.95)])
ylabel('Pixel Count')
xlabel('Intensity')
title('Corrected Image Histogram')
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
