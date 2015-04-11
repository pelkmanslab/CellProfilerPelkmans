function handles = IlluminationCorrectionPerSite(handles)

% Help for the IluminationCorrectionPerSite module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Uses std and mean illumination function calculated by iBrain, on both a
% per site and plate-wise basis, which are then used to correct uneven
% illumination/lighting/shading via Z-scoring, with optional masking of
% pixels from the correction. Useful when imaging entire wells such that
% each site comprises differing regions of empty space.
% *************************************************************************
%
% This module is based on IlluminationCorrectonZScoring.m. Briedly, it calculates the mean and std values for the intensity of each
% pixel position in a set of images from the same wave length and same site. Such
% functions can later be used to to correct the intensity values of each
% pixel by substracting the mean function and normalizing with the std
% function. Such a procedure should not only normalize for the background
% intensities, but also for amplitude differences which result from optical
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
% stored in the form 'illcormask_site%d.mat'). Alternatively, the mask may be calculated with pixels above a certain threshold, defined from the
% std intensity per site, corrected, with the others masked from
% correction. Masking is useful to avoid over correction of empty
% space. The default of zero results in no masking. You can also
% choose to expand the mask.

% Note: the images for this module should be directly loaded by the
% LoadImages module. In addition, precalculated illumination correction
% files (from iBrain) are required, both plate-wise and per site
% (see prepare_batch_measure_illcor_stats.m).
%
% Authors:
%   Victoria Green
%   Nico Battich
%   Berend Snijder
%   Yauhen Yakimovich
%   Thomas Stoeger
%   
%
% $Revision: 1809 $


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

%textVAR06 = Enter the relative path name to the folder where the illumination correction files are located (starting with "./../"). Type period (.) for default directory.
%defaultVAR06 = .
AlternativeIllCorrFolder = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Do you want to apply a mask?
%choiceVAR07 = Yes
%choiceVAR07 = No
ApplyMaskExp = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu

%textVAR08 = Are the mask(s) pre-defined and located in the BATCH directory?
%choiceVAR08 = Yes
%choiceVAR08 = No
MaskTakeOver = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 = Otionally enter the std intensity threshold to define mask.
%defaultVAR09 = 0
intThresh = str2double(handles.Settings.VariableValues{CurrentModuleNum,09});

%textVAR10 = Optionally enter the size of the mask expansion.
%defaultVAR10 = 100
intExpMask = str2double(handles.Settings.VariableValues{CurrentModuleNum,10});


%%%VariableRevisionNumber = 5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drawnow

% The illumination correction image was calculated using all the incoming
% images by iBRAIN. Load the mean and std statistics assuming iBrain file 
% structure. 
    
% get the channel of the image
[intChannelNumber] = check_image_channel(handles.Pipeline.(strcat('FileList',InputName)){handles.Current.SetBeingAnalyzed});
intZstackNumber = 0;
        
% get microscope type and image position
[intImagePosition,strMicroscopeType] = check_image_position(handles.Pipeline.(strcat('FileList',InputName)){handles.Current.SetBeingAnalyzed});

% store the stat in the Measurement field, with site, channel and
% szstack specific fieldnames
strStatFieldName = sprintf('illcor_s%03dch%03dz%03d',intImagePosition,intChannelNumber,intZstackNumber);
    
%if handles.Current.SetBeingAnalyzed == 1
if ~isfield(handles.Measurements.Image, [strStatFieldName '_mean'])
    % Unload previous illcor measurements if found.
    handles.Measurements.Image = removeNames(handles.Measurements.Image, [strStatFieldName '.*']);
    
    % Load illcor measurements for the next site
    % get the BATCH directory and load illcor data
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
    end
    strBatchDir = VarPathname;
    
    TempStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',intChannelNumber,intZstackNumber)));
    TempSiteStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_site%d_channel%03d_zstack%03d.mat',intImagePosition,intChannelNumber,intZstackNumber)));
    
    % Avoid potential floating type precision mixing.
    TempStats.stat_values.mean = double(TempStats.stat_values.mean);
    TempStats.stat_values.std = double(TempStats.stat_values.std);
    TempSiteStats.stat_values.mean = double(TempSiteStats.stat_values.mean);
    TempSiteStats.stat_values.std = double(TempSiteStats.stat_values.std);
    
    
    % Optional smoothing of illumination statistics computed by iBRAIN.
    switch DoMeanStdSmoothing
        case 'Yes'
            % Create gaussian filter handle for for smoothing
            H = fspecial('gaussian',[SmoothingSize SmoothingSize],SmoothingSize*SmoothingSigma);
            handles.Measurements.Image.([strStatFieldName,'_mean']) = imfilter(TempStats.stat_values.mean,H,'symmetric');
            handles.Measurements.Image.([strStatFieldName,'_std']) = imfilter(TempStats.stat_values.std,H,'symmetric');
            handles.Measurements.Image.([strStatFieldName,'_sitemean']) = imfilter(TempSiteStats.stat_values.mean,H,'symmetric');
            handles.Measurements.Image.([strStatFieldName,'_sitestd']) = imfilter(TempSiteStats.stat_values.std,H,'symmetric');
        case 'No'
            handles.Measurements.Image.([strStatFieldName,'_mean']) = TempStats.stat_values.mean;
            handles.Measurements.Image.([strStatFieldName,'_std']) = TempStats.stat_values.std;
            handles.Measurements.Image.([strStatFieldName,'_sitemean']) = TempSiteStats.stat_values.mean;
            handles.Measurements.Image.([strStatFieldName,'_sitestd']) = TempSiteStats.stat_values.std;
    end
    clear TempStats
    clear TempSiteStats
    
end

% define correction matrices
IllumFilt_Mean = handles.Measurements.Image.([strStatFieldName,'_mean']);
IllumFilt_STD = handles.Measurements.Image.([strStatFieldName,'_std']);
IllumFilt_SiteMean = handles.Measurements.Image.([strStatFieldName,'_sitemean']);
IllumFilt_SiteSTD = handles.Measurements.Image.([strStatFieldName,'_sitestd']);


%%%%%%%%%%%%%%%%%%%
%%% DEFINE MASK %%%
%%%%%%%%%%%%%%%%%%%

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
            error(['Failed to find mask file: ' matMaskFilename])
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

OrigImage(OrigImage == 0) = 1;
ImageOutput = (log10(OrigImage)-IllumFilt_SiteMean)./IllumFilt_SiteSTD;
ImageOutput = (ImageOutput.*mean(IllumFilt_STD(:)))+mean(IllumFilt_Mean(:));
ImageOutput = 10.^ImageOutput;
ImageOutput(MaskIdx) = OrigImage(MaskIdx); % replace corrected values outside mask with original values.

% fix potentially broken pixels (which are not variable) but only attempt
% to do so if there is no mask (which implies empty space and therefore
% invariable pixels).
if  strncmpi(ApplyMaskExp,'y',1)
    ImageOutput = ImageOutput;
else
    ImageOutput = fixNonNumericalValueInImage(ImageOutput);
end

% store non-scaled for visualization
ImageOutputPlot = ImageOutput;
% Rescale from 0  to 1
ImageOutput = ImageOutput/65535;

% Check if there are pixels with a value higher than one, fix them and give
% warning.
f = ImageOutput(:)>1;
HighPixels = sum(f);
if HighPixels > 0
    warning('%.3d pixels with values higher than 1 were found.\nAll these will be set to 1.\n',HighPixels)
    ImageOutput(f) = 1;
end

%save to handle structure
handles.Pipeline.(OutputName) = ImageOutput;
fieldname = ['Filename', OutputName];
handles.Pipeline.(fieldname) = handles.Pipeline.(strcat('FileList',InputName))(handles.Current.SetBeingAnalyzed);
fieldname = ['Pathname', OutputName];
handles.Pipeline.(fieldname) = handles.Pipeline.(strcat('Pathname',InputName));


%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

if ~CPisHeadless()

    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
    subplot(2,4,1);
    try
        CPimagesc(IllumFilt_Mean,handles);
        colormap('JET')
        colorbar
        title('Mean Intensity Filter [Log10(intensity)]')
    end
    subplot(2,4,2);
    try
        CPimagesc(IllumFilt_STD,handles);
        colormap('JET')
        colorbar
        title('STD Intensity Filter [Log10(intensity)]')
    end
    subplot(2,4,5);
    try
        CPimagesc(IllumFilt_SiteMean,handles);
        colormap('JET')
        colorbar
        title('Site Mean Intensity Filter [Log10(intensity)]')
    end
    subplot(2,4,6);
    try
        CPimagesc(IllumFilt_SiteSTD,handles);
        colormap('JET')
        colorbar
        title('Site STD Intensity Filter [Log10(intensity)]')
    end

    subplot(2,4,3);
    try
        CPimagesc(OrigImage,handles);
        colormap('JET')
        title('Original Image')
    end
    subplot(2,4,7);
    try
        CPimagesc(ImageOutputPlot,handles);
        colormap('JET')
        if strncmpi(AplyBackgroundSubtraction,'y',1)
            title(sprintf('Corrected Image\n+Background Subtraction'),'FontSize',7)
        else
            title('Corrected Image')
        end
    end
    subplot(2,4,4);
    try
        hold on
        hist(IllumFilt_SiteSTD(:),round(length(IllumFilt_SiteSTD(:))/20))
        vline(intThresh,'--r')
        hold off
        ylabel('Pixel Count')
        xlabel('Site STD Intensity Filter [Log10(intensity)]')
        title('Site STD Intensity Filter Histogram')
    end
    subplot(2,4,8);
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

