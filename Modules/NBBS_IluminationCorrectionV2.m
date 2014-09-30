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

%textVAR02 = If yes to previous option do you want to runn the illumination correction per site or over all sites 
%choiceVAR02 = PerSite
%choiceVAR02 = OverAllSites
PerSiteCorrection = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What did you call the images to be used to calculate the illumination functions?
%infotypeVAR03 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = What do you want to call the mean illumination function?
%defaultVAR04 = MeanIllumBlue
%infotypeVAR04 = imagegroup indep
Mean_IlluminationImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = What do you want to call the std illumination function?
%defaultVAR05 = StdIllumBlue
%infotypeVAR05 = imagegroup indep
Std_IlluminationImageName = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Smoothing is done by a Gaussian filter to both illumination functions. Enter the size of the filter.
%defaultVAR06 = 100
SmoothingSize = str2num(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Enter the factor 'x' for the sigma calculation of the Gaussian filter. where sigma = x*size.
%defaultVAR07 = 0.5
SmoothingSigma = str2num(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Enter the number of iterations runned to approximate the std value per pixel.
%defaultVAR08 = 20
IterationNum = str2num(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = Do you wnat to apply the illumination correction?
%choiceVAR09 = Yes
%choiceVAR09 = No
AplyIllumCorr = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = If Yes, how do you want to call the corrected images?
%defaultVAR10 = CorrBlue
%infotypeVAR10 = imagegroup indep
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = Do you wnat to apply Background substraction after Illumination correction?
%choiceVAR11 = Yes
%choiceVAR11 = No
ApplyBackSubs = char(handles.Settings.VariableValues{CurrentModuleNum,11});
%inputtypeVAR11 = popupmenu

%textVAR12 = If Yes, What is the value to substract?
%defaultVAR12 = 0.0014
SubstractParam = str2num(handles.Settings.VariableValues{CurrentModuleNum,12});

%textVAR13 = If Yes, How do you want to call the substracted images?
%defaultVAR13 = SubsCorrBlue
%infotypeVAR13 = imagegroup indep
OutputBackName = char(handles.Settings.VariableValues{CurrentModuleNum,13});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%get a list of all images to be used
ImageList = cellfun(@(x) fullfile(handles.Pipeline.(strcat('Pathname',ImageName)),x),...
    handles.Pipeline.(strcat('FileList',ImageName)), 'uniformoutput',false)';
ImageNameList = handles.Pipeline.(strcat('FileList',ImageName));


%%% The illumination correction image is calculated using all the incoming
%%% images. This module is therefore only runned in the first cycle. This
%%% also allows for batch processing

if handles.Current.SetBeingAnalyzed == 1 && ~strncmpi(iBrainTakeOver,'y',1)
    
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
    
    
elseif strncmpi(iBrainTakeOver,'y',1) %&& handles.Current.SetBeingAnalyzed == 1
    
    %if 1st cycle get the maximum image position 
    if handles.Current.SetBeingAnalyzed == 1
        
        handles.Pipeline.intMaxImagePosition = max(cellfun(@check_image_position,ImageNameList));        
        handles.Pipeline.intMaxImagePosition = 9;
    end
    
    % load the illum functions assuming iBrain file structure % note: load
    % the functions every cycle because of z-stacks which might make the batch-data file too big.
    strImageName = handles.Pipeline.(strcat('FileList',ImageName)){handles.Current.SetBeingAnalyzed};
    
    %get the channel of the image
    intChannelNumber = check_image_channel(strImageName);
    %intZstackNumber = 0;
    
    %get the site information 
    %[intRow, intColumn, strWellName, intTimepoint] = filterimagenamedata(strImageName);
    
    %get microscope type and image position
    [intImagePosition,strMicroscopeType] = check_image_position(strImageName);
    
    %get position of site in the well
    matImageSnake = get_image_snake(handles.Pipeline.intMaxImagePosition, strMicroscopeType);
    RowIdx = matImageSnake(2,intImagePosition)+1;
    ColumnIdx = matImageSnake(1,intImagePosition)+1;
    %matImageSnake = matImageSnake'+1; matImageSnake(intImagePosition,1) = column index
    
    %get the BATCH directory
    strBatchDir = handles.Current.DefaultOutputDirectory;
    
    %Note that the loading of the correction functions is the most time
    %consuming step so try to cache it 
    
    %First look for the file in the local cache
    if isfield(handles.Pipeline, 'PathIllumCached')
        %[NB] this will create a problem in the file is deleted for
        %specially in the cluster
        if exist(handles.Pipeline.PathIllumCached)
            load(handles.Pipeline.PathIllumCached);
        else
            %if the file does not exist loaded again
            warning('%s: Caching failure. File not found. Loading file from BATCH directory.',mfilename)
            TempStats = load(fullfile(strBatchDir,sprintf('illcor2.mat')));
        end
    else
        TempStats = load(fullfile(strBatchDir,sprintf('illcor2.mat')));
        %get the caching directory
        if isunix && strncmpi(getenv('HOSTNAME'),'brutus',6)
            strCachePath = fullfile('/cluster/work/biol/bsnijder/','IllumCorrCaching');
        else
            strCachePath = fullfile(tempdir,'IllumCorrCaching');
        end
        %if the caching directory does not exist then make it
        if ~isdir(strCachePath)
            mkdir(strCachePath);
        end
        %generate temp name
        handles.Pipeline.PathIllumCached = fullfile(strCachePath,strcat(getlastdir(tempname),'.mat'));
        % save the file
        save(handles.Pipeline.PathIllumCached,'TempStats');
    end
    
    IllumFilt_Mean = TempStats.cellBeforeCorrectionStats{RowIdx,ColumnIdx,intChannelNumber}.mean;
    IllumFilt_STD = TempStats.cellBeforeCorrectionStats{RowIdx,ColumnIdx,intChannelNumber}.std;
    
else
    return
end

%do correction
if strncmpi(AplyIllumCorr,'y',1)
    
    OrigImage = double(imread(ImageList{handles.Current.SetBeingAnalyzed}));
    CorrectedImage = IllumCorrect(OrigImage,IllumFilt_Mean,IllumFilt_STD,1);
    ImageOutputPlot = CorrectedImage;
    
    % division required for compatibility with cellprofiler
    ImageOutput = CorrectedImage./65535;%intmax('uint16');
    ImageOutput(ImageOutput>1) = 1; 
    %save to handle structure
    handles.Pipeline.(OutputName) = ImageOutput;
    
    
end

if strncmpi(ApplyBackSubs,'y',1)

    BackSubsImage = ImageOutput-SubstractParam;
    BackSubsImage(BackSubsImage<0) = 0;
        
    %save to handle structure
    handles.Pipeline.(OutputBackName) = BackSubsImage;
    
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


