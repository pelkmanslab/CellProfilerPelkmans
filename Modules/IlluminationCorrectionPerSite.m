function handles = IlluminationCorrectionPerSite(handles)

% Help for the IluminationCorrectionPerSite module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Calculates a mean illumination function and the std of the illumination
% mean, on both a per site and plate-wise basis, which are then used to correct uneven illumination/lighting/shading via Z-scoring.
% *************************************************************************
%
% This module is based on NBBS_IlluminationCorrectonPerSite2. It calculates the mean and std values for the intensity of each
% pixel position in a set of images from the same wave length and same site. Such
% functions can later be used to to correct the intensity values of each
% pixel by substracting the mean function and normalizing with the std
% function. Such procedure should not only normalize for the background
% intensities but also for amplitude differences which result from optical
% aberrances or differences in illumination. The module also calculates the
% mean and std values for the intensity of each pixel position in the whole
% image set, irrespective of site, which are used to transform the images
% that were normalised per site. This avoids differences in the corrected
% images on a per site basis, which arises if the amount of empty space differs
% per site. 

% Masking can be performed such that only defined pixels are subject to
% illumination correction. This mask may be predefined and stored in the
% BATCH folder as a logical matrix the same size as the images, with 1
% indicating a masked pixel and 0 a non-masked pixel. (These files must be
% stored in the form 'illcorr_mask_site%d.mat'). Alternatively, the mask may be calculated with pixels above a certain threshold, defined from the
% std intensity per site, corrected, with the others masked from
% correction. Masking is useful to avoid over correction of empty
% space. The default of zero results in no masking. You can also
% choose to expand the mask.
% Note: the images for this module should be directly loaded by the
% LoadImages module. In addition, at present precalculated illumination correction
% files (from iBrain) are required, both plate-wise and per site
% (see prepare_batch_measure_illcor_stats.m).
%
% Authors:
%   Victoria Green
%   Nico Battich
%   Berend Snijder
%   
%
% $Revision: 1809 $

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

%textVAR08 = Do you want to apply a mask?
%choiceVAR08 = Yes
%choiceVAR08 = No
ApplyMaskExp = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 = Are the mask(s) pre-defined and located in the BATCH directory?
%choiceVAR09 = Yes
%choiceVAR09 = No
MaskTakeOver = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = Otionally enter the std intensity threshold to define mask.
%defaultVAR10 = 0
intThresh = str2num(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = Optionally enter the size of the mask expansion.
%defaultVAR11 = 100
intExpMask = str2num(handles.Settings.VariableValues{CurrentModuleNum,11});

%textVAR12 = Do you want to apply the illumination correction?
%choiceVAR12 = Yes
%choiceVAR12 = No
AplyIllumCorr = char(handles.Settings.VariableValues{CurrentModuleNum,12});
%inputtypeVAR12 = popupmenu

%textVAR13 = Do you also want to apply background subtraction? (iBRAIN only)
%choiceVAR13 = Yes
%choiceVAR13 = No
AplyBackgroundSubtraction = char(handles.Settings.VariableValues{CurrentModuleNum,13});
%inputtypeVAR13 = popupmenu

%textVAR14 = Optionally enter a background-subtraction correction factor (<1 = less subtraction, >1 = more subtraction
%defaultVAR14 = 1
BackgroundCorrectionCorrectionFactor = str2num(handles.Settings.VariableValues{CurrentModuleNum,14}); %#ok<*ST2NM>

%textVAR15 = How do you want to call the corrected images?
%defaultVAR15 = CorrBlue
%infotypeVAR15 = imagegroup indep
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,15});

%textVAR16 = Do you also want to smoothing the iBrain derived correction function?
%choiceVAR16 = No
%choiceVAR16 = Yes
iDoBrainSmoothing = char(handles.Settings.VariableValues{CurrentModuleNum,16});
%inputtypeVAR16 = popupmenu

%%%VariableRevisionNumber = 4


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
%%% also allows for batch processing. At present, only calculated
%%% plate-wise per channel! Therefore require iBRIAN calculated files for
%%% per site.

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
    [intChannelNumber] = check_image_channel(handles.Pipeline.(strcat('FileList',ImageName)){handles.Current.SetBeingAnalyzed});
    intZstackNumber = 0;
        
    %get microscope type and image position
    [intImagePosition,strMicroscopeType] = check_image_position(handles.Pipeline.(strcat('FileList',ImageName)){handles.Current.SetBeingAnalyzed});
    
    % [BS] store the stat in the Measurement field, with site, channel and
    % szstack specific fieldnames 
    strStatFieldName = sprintf('illcor_s%03dch%03dz%03d',intImagePosition,intChannelNumber,intZstackNumber);
    
    %if handles.Current.SetBeingAnalyzed == 1
    if ~isfield(handles.Measurements.Image, [strStatFieldName '_mean'])
        % Unload previous illcor measurements if found.
        handles.Measurements.Image = removeNames(handles.Measurements.Image, [strStatFieldName '.*']);
        
        % Load illcor measurements for the next site
        % get the BATCH directory and load illcor data
        strBatchDir = handles.Current.DefaultOutputDirectory;
        TempStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',intChannelNumber,intZstackNumber)));
        TempSiteStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_site%d_channel%03d_zstack%03d.mat',intImagePosition,intChannelNumber,intZstackNumber)));
        
        % Avoid potential floating type precision mixing.
        TempStats.stat_values.mean = double(TempStats.stat_values.mean);
        TempStats.stat_values.std = double(TempStats.stat_values.std);
        TempSiteStats.stat_values.mean = double(TempSiteStats.stat_values.mean);
        TempSiteStats.stat_values.std = double(TempSiteStats.stat_values.std);
        
        
        %[TS] adjusted to include option for smoothing iBrain prepared correction function
        switch iDoBrainSmoothing
            case 'Yes'
                ValueMean = imfilter(TempStats.stat_values.mean,H,'replicate');
                ValueStd = imfilter(TempStats.stat_values.std,H,'replicate');
                ValueSiteMean = imfilter(TempSiteStats.stat_values.mean,H,'replicate');
                ValueSiteStd = imfilter(TempSiteStats.stat_values.std,H,'replicate');
            case 'No'
                ValueMean = TempStats.stat_values.mean;
                ValueStd = TempStats.stat_values.std;
                ValueSiteMean = TempSiteStats.stat_values.mean;
                ValueSiteStd = TempSiteStats.stat_values.std;
        end
        
        handles.Measurements.Image.([strStatFieldName,'_mean']) = ValueMean;
        handles.Measurements.Image.([strStatFieldName,'_std']) = ValueStd;
        handles.Measurements.Image.([strStatFieldName,'_sitemean']) = ValueSiteMean;
        handles.Measurements.Image.([strStatFieldName,'_sitestd']) = ValueSiteStd;
        clear TempStats
        clear TempSiteStats
    end
    
    IllumFilt_Mean = handles.Measurements.Image.([strStatFieldName,'_mean']);
    IllumFilt_STD = handles.Measurements.Image.([strStatFieldName,'_std']);
    IllumFilt_SiteMean = handles.Measurements.Image.([strStatFieldName,'_sitemean']);
    IllumFilt_SiteSTD = handles.Measurements.Image.([strStatFieldName,'_sitestd']);
    BackgroundSubtractionVal = (10^mean(IllumFilt_Mean(:))) * BackgroundCorrectionCorrectionFactor;

    % define mask
    
    if  strncmpi(ApplyMaskExp,'y',1)
        if strncmpi(MaskTakeOver,'y',1)
            strBatchDir = handles.Current.DefaultOutputDirectory;
            matMaskFilename = fullfile(strBatchDir,sprintf('illcormask_site%d.mat',intImagePosition));
            if ~os.path.exists(matMaskFilename)
                % check in the plate folder
                d = @os.path.dirname;
                matMaskFilename = fullfile(d(d(matMaskFilename)), os.path.basename(matMaskFilename));
                if ~os.path.exists(matMaskFilename)
                    % check in the  folder                    
                    matMaskFilename = fullfile(d(d(matMaskFilename)), os.path.basename(matMaskFilename));
                end
            end
            if ~os.path.exists(matMaskFilename)
                error(['Filed to find mask file: ' matMaskFilename])
            end
            matMask = load(fullfile(matMaskFilename));
            matMask = matMask.matMask;
        else
            matMask = IllumFilt_SiteSTD < intThresh; % at present mask threshold empirically determined and input by user, but could be calculated from IllumFilt_SiteSTD.
            SegMask = bwlabel(matMask);
            M = regionprops(SegMask,'Area','Centroid','PixelIdxList');
            [matArea,Index] = sort([M.Area],'descend');
            matMask = zeros(size(IllumFilt_SiteSTD)); % reinitiate mask
            matMask(M(Index(1),1).PixelIdxList) = 1; % take biggest area only for mask

            % do mask expansion
            matMask = bwmorph(matMask,'thicken',intExpMask);
            SE = strel('square',10); % dilate expanded mask to fill holes
            matMask = imdilate(matMask,SE);
        end
    else
        matMask = zeros(size(IllumFilt_SiteSTD));
    end
    
    % find indices of mask for correction
    MaskIdx = find(matMask==1);
   
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
    ImageOutput = (log10(OrigImage)-IllumFilt_SiteMean)./IllumFilt_SiteSTD;
    ImageOutput = (ImageOutput.*mean(IllumFilt_STD(:)))+mean(IllumFilt_Mean(:));
    ImageOutput = 10.^ImageOutput;
    ImageOutput(MaskIdx) = OrigImage(MaskIdx); % replace corrected values outside mask with original values.
    
    if strncmpi(AplyBackgroundSubtraction,'y',1)
        % [BS, 2012-06-25] Added optional background subtraction
        ImageOutput = ImageOutput-BackgroundSubtractionVal;
        ImageOutput(ImageOutput<0) = 0;
    end    
    
    % store non-scaled for visualization
    ImageOutputPlot = ImageOutput;    
    % rescale 0 / 1
    ImageOutput = ImageOutput/65535;
    
    %save to handle structure
    handles.Pipeline.(OutputName) = ImageOutput;
    
    
end

%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

if ~CPisHeadless()

    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
    subplot(2,5,1);
    try
        CPimagesc(IllumFilt_Mean,handles);
        colormap('JET')
        colorbar
        title('Mean Intensity Filter [Log10(intensity)]')
    end
    subplot(2,5,2);
    try
        CPimagesc(IllumFilt_STD,handles);
        colormap('JET')
        colorbar
        title('STD Intensity Filter [Log10(intensity)]')
    end
    subplot(2,5,6);
    try
        CPimagesc(IllumFilt_SiteMean,handles);
        colormap('JET')
        colorbar
        title('Site Mean Intensity Filter [Log10(intensity)]')
    end
    subplot(2,5,7);
    try
        CPimagesc(IllumFilt_SiteSTD,handles);
        colormap('JET')
        colorbar
        title('Site STD Intensity Filter [Log10(intensity)]')
    end

    subplot(2,5,3);
    try
        CPimagesc(OrigImage,handles);
        colormap('JET')
        title('Original Image')
    end
    subplot(2,5,8);
    try
        CPimagesc(ImageOutputPlot,handles);
        colormap('JET')
        if strncmpi(AplyBackgroundSubtraction,'y',1)
            title(sprintf('Corrected Image\n+Background Subtraction'),'FontSize',7)
        else
            title('Corrected Image')
        end
    end

    subplot(2,5,4);
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
    subplot(2,5,9);
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
    subplot(2,5,5);
    try
        hold on
        hist(IllumFilt_SiteSTD(:),round(length(IllumFilt_SiteSTD(:))/20))
        vline(intThresh,'--r')
        hold off
        ylabel('Pixel Count')
        xlabel('Site STD Intensity Filter [Log10(intensity)]')
        title('Site STD Intensity Filter Histogram')
    end
    subplot(2,5,10);
    try
        CPimagesc(matMask,handles);
        colormap('JET')
        title('Correction Mask')    
    end
    
end 


drawnow

end

function s = removeNames(s, matchStr)

names = fieldnames(s);
for iName = 1:numel(names)
    name = names{iName};
    if regexpi(name, matchStr)
        s = rmfield(s, name);
    end
end

end


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
