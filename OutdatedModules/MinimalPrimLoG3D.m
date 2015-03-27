function handles = MinimalPrimLoG3D(handles)

% Help for the Identify Primary LoG module:
% Category: Object Processing

%note> final cropping of padded image may temporally requires dublication
%of memory: rather substract and then remove edges

% some of the variables have ancient names indicating their origin in the
% code for finding the treshold



drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images you want to process?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = Spots
%infotypeVAR02 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});


%textVAR03 = MinQuant
%defaultVAR03 = 0.01
sMinQuant = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = UpQuant
%defaultVAR04 = 0.99
sMaxQuant = char(handles.Settings.VariableValues{CurrentModuleNum,4});


%textVAR05 = Minimum of UpThreshold.
%defaultVAR05 = 400
sImgThresholdUpMin = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = PixelSizeSpot
%defaultVAR06 = 3
sRadius = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = ThresholdSpots
%defaultVAR07 = 0.15
sThresholdForSpots = char(handles.Settings.VariableValues{CurrentModuleNum,7});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  INITIALIZE SETTINGS  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% usually these settings are from the parameter finding function 
settings.General.Radius =str2double(sRadius);
settings.General.ImgThresholdUpMin = str2double(sImgThresholdUpMin);
settings.General.QuantMinIntensityPerImage=str2double(sMinQuant);
settings.General.QuantMaxIntensityPerImage=str2double(sMaxQuant);
settings.General.DownsamplingForIntensityQuanta=10;

% Select List of Files
strFilenames = handles.Pipeline.(ImageName).FileNames;



% numThresholdsToTest = length(settings.AutoSetup.DetermineSpotThreshold.ThresholdsToTest);
numRadius = settings.General.Radius;
numAllowedImgThresholdUpMin= single(settings.General.ImgThresholdUpMin);
strGeneralSettings=fieldnames(settings.General);

%%%%%%%% DIMENSIONS OF IMAGE SET %%%%%%%%%%%%%%
% 
% if size(strFilenames,2)>1
    numImageDimension = 3;
% else
%     numImageDimension = 2;
% end

%%%%%%% DIMENSIONS OF FILTERING %%%%%%%%%%%%%%%%%

hereIsFilterDimension=NaN;
for k=1:size(strGeneralSettings,1)
    if strcmp(strGeneralSettings{k,1},'FilterDimension')
        hereIsFilterDimension =k;
    end
end
if hereIsFilterDimension > 0
    if max(size(settings.General.FilterDimension)) >1
        numFilterDimension = settings.General.FilterDimension(1,1);
        fprintf('Filterdimensions is larger than a 1x1 matrix. Will try to first value %d.', numFilterDimension);
    else
        numFilterDimension=settings.General.FilterDimension;
    end
    
else
    
    numFilterDimension = numImageDimension;
    
end


%%%%%%% SIZE OF MASK FOR CONNECTING PIXELS TO OBJECTS  %%%%%%%

% 
% hereIsConnectivitySize=NaN;
% for k=1:size(strGeneralSettings,1)
%     if strcmp(strGeneralSettings{k,1},'ConnectivitySize')
%         hereIsConnectivitySize =k;
%     end
% end
% 
% if hereIsConnectivitySize>0
%     if max(size(settings.General.ConnectivitySize)) >1
%         numConnectivitySize = settings.General.ConnectivitySize(1,1);
%         fprintf('ConnectivitySize is larger than a 1x1 matrix. Will use first value %d.', numConnectivitySize);
%     else
%         numConnectivitySize=settings.General.ConnectivitySize;
%     end
% else
%     switch numImageDimension   %OTHERWISE BWCONNCOMP STANDARD VALUES OF 8 FOR 2D AND 26 FOR 3D
%         case 2
%             numConnectivitySize = 8;
%         case 3
            numConnectivitySize = 26;
%     end
%     
% end

fprintf('Will apply LoG filtering along %d dimensions. \n Objects are connected along %d dimensions using a mask of %d pixels.\n', numFilterDimension, numImageDimension,numConnectivitySize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  SPOT COUNTING  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
numDataSets=size(strFilenames,1);
% SpotsCount=NaN(numDataSets,numThresholdsToTest);
% PixelIdxList = cell(numDataSets,1);

for j=1:numDataSets
    % obtain quantile minimal and maximal intensity 
    [ImgThresholdMin, ImgThresholdMax] = getQuantiles(strFilenames(j,:),settings);
    ImgThresholdMin = single(ImgThresholdMin);
    ImgThresholdMax = single(ImgThresholdMax);
    
    % load image set in origial resolution
    % (note that also previously loaded temporally within getQuantile,
    % but try to reduce memory during quantile, which is slow in 3D, if not
    % downscaled)
    Images = cell(1,size(strFilenames,2));
    for k = 1:size(strFilenames,2)
        Images{k} = single(imread(strFilenames{j,k}));
    end
    Images=cat(3,Images{:});
    
    %Rescale Image Intensity 
    if ImgThresholdMax < numAllowedImgThresholdUpMin
        ImgThresholdMax = numAllowedImgThresholdUpMin;
    end
    
    
    Images = Images./ImgThresholdMax;
    ImgThresholdMin = ImgThresholdMin/ImgThresholdMax;
    Images = Images - ImgThresholdMin;
    Images = 1-Images;
    
    %create LoG filter for expected spot size
    sigma = (numRadius-1)/3;
    op = fspecial('log',numRadius,sigma);
    op = op - sum(op(:))/numel(op); % make the op to sum to zero
    
    
% Pad image to fix border artifact 
% has been part of IdentifyPrimLoG2, but not of Raj et al. LoG
padsize = ceil(numRadius./2);
Images = padarray(Images,[padsize padsize],'replicate');

    
    %Perform Filtering
    switch numImageDimension
        case 2
            Images = imfilter(single(Images),op,'replicate');
        case 3
            switch numFilterDimension
                case 2
                    for k=1: size(Images,3)
                        Images(:,:,k) =  imfilter(single(Images(:,:,k)),op,'replicate');
                    end
                case 3
                    %If 3 dimensions, follow code of Raj et al., used for
                    % singlemoleculefish: (see his comments & code commented below)
                    % Here, we amplify the signal by making the filter "3-D"
                    % H = 1/3*cat(3,H,H,H);
                    op = 1/3*cat(3,op,op,op);
                    Images = imfilter(single(Images),op,'replicate');
                otherwise
                    fprintf('No valid Dimension for filtering defined (%d). Possible: 2 or 3.',numFilterDimension);
            end
        otherwise
            fprintf('No valid Dimension for Image found (%d). Possible: 2 or 3.',numImageDimension);
    end

    Images = Images(padsize+1:end-padsize,padsize+1:end-padsize,:);

       
    bw=Images>str2double(sThresholdForSpots);
    vislabel = bwconncomp(bw,numConnectivitySize);          % compared to Raj and previous IdentifyPrimLog2 use bwconncomp since fewer memory required

fieldname = ['Segmented', ObjectName];
handles.Pipeline.(fieldname).Label = vislabel;
handles.Pipeline.(fieldname).Format = 'SegmentationCC'; % Define Format of Output. Here it is Segmentation CC which corresponds to a full structure obtained by bwconncomp




end