function handles = SimpleMulti_LoadSegmentedObjects(handles)

% Help for the SmartAlignment_LoadSegmentedObject module:
% Category: Other
%
% SHORT DESCRIPTION:
% Module, which loads the segmentation of an arbitrary different reference
% acquisition (second user input). 
%
% Since multiple acquisition usually have a slight shift, a image-channel 
% (e.g.: DAPI) has to be specified for the current acquistion and the
% reference acquistion. Note that the alignment is very robust against
% illumination artifacts and thus can be done with images that are not
% illumination corrected.
%
% DIFFERENCES TO MPCYCLE MODULES:
%
% a) LABELLING OF OBJECTS CHANGES BETWEEN ACQUISITIONS. This is because 
%    SimpleMulti_LoadSegmentedObjects uses continous labels for objects.
%    This prevents several mistakes in data output, which would be introduced
%    when this CP-default assumption of several modules is broken. This
%    wrong data would be difficult to spot and only in rare occassions can
%    results in a wrong number of objects (which will also cause
%    computational error besides being indicative that measurments are
%    allocated to wrong objects!). Data can still be combined via the 
%    UnshiftedObjectId measuremnent created by this module. See below for
%    detailed example.
% b) Only a single CellProfiler module is required.
% c) Nothing has to be precomputed (e.g. no shiftdescriptor). 
%    Simply adding this module suffices.
% d) no later parts of iBrain has to be changed since every
%    acquisition will be a normal ordinary acquistion that is not different
%    from any ordinary acquistion. This also means that every function used
%    for in-detail data exploration, after iBrain, will work as on any
%    other ordinary experiment.
% e) Acquistions do not have to be named a certain way. (allowing ad-hoc
%    rescanning/combination of experiments)
% f) Images are not cropped. This however does not affect the data output!
% g) Acquistions do not have to be organised in hierarchical manner (e.g.:
%    segmentation can be from first or last acquistion and there is no
%    necessity to wait for other acquistions).
%
%
% If data from multiple different acquistions should be combined, please
% use loading function (e.g.: getRawProbModelData2) separately. Fusion is
% possible by using the Measurements_(whateverObject)_UnshiftedObjectId.mat
% measurement created by this module. For each single object this contains
% the object ID from the corresponding reference acquisition. e.g. along
% the following pseudo-code, where Rows Columns ObjectID and UnshiftedIDs
% refer to the numerical index of the columns that contain the
% corresponding information
% 
% ID_within_reference_acquistion  = matMeta_reference(:,[Rows Columns ObjectID]);
% ID_within_current_acquistion = [matMeta_current(:,[Rows Columns]) matData_current(:, UnshiftedIDs) ]     ;
%
% [is_also_in_reference_acquistion, pos_within_reference_acquistion] = ismember(ID_within_current_acquistion, ID_within_reference_acquistion, 'rows')
% 
% if any(~is_also_in_reference_acquistion)
%   error('some data could not be matched. Please check code, possibly something is wrong in CP module')
% else
%   joinedData = [matData_current matData_reference(ID_within_reference_acquistion,:)];
% end
% *************************************************************************
%
% Authors:
%   Thomas Stoeger
%
% $Revision: 1879 $


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = Which objects do you want to import? (e.g.: Nuclei)
%defaultVAR01 = Nuclei
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%pathnametextVAR02 = From which reference acquistion should objects be imported?
%defaultVAR02 = /BIOL/sonas/biol_uzh_pelkmans_s5/Data/Users/Thomas/151020-AribitraryPlate/
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which images of the current pipeline should be used for aligning the acquisitions?
%infotypeVAR03 = imagegroup
OrigImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Which filter describes the images of the reference acquisition that should be used for aligning the acquisitions?
%defaultVAR04 = C01.
filterForImagesOfTrans = char(handles.Settings.VariableValues{CurrentModuleNum,4});


%%%VariableRevisionNumber = 12



%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZATION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%

strTransPlate = npc(Pathname);
if ~any(fileattrib(strTransPlate))
    error('Could not find reference plate');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Computation  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%% FOR FIRST CYCLE %%%%%%%%%
if handles.Current.SetBeingAnalyzed == 1
      makeExternalBuffer(handles, strTransPlate, filterForImagesOfTrans, ObjectName);
end

%%%%%% FOR ALL CYCLES %%%%%%%%%%
createShiftedSegmentation = true; % depending on results, this switch will possbily be turned (TS151021: creating a segmentation without objects).

[TiffFolder_trans, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate);
[strCorrespondingImage_trans, couldFindSameSite_image, strCorrespondingSegmentation_trans, couldFindSameSite_segmentation] = getFileNameViaReferenceFile(handles, ObjectName, OrigImageName);

% Obtain coordinates of overalap of cis image (this acquisition) and trans image (reference acquistion)
if couldFindSameSite_image == true;
    Image_Cis = CPretrieveimage(handles,OrigImageName,ModuleName);
    readFun = @(x) double(imread(x)) ./ (2^16-1);
    Image_Trans = readFun(fullfile(TiffFolder_trans, strCorrespondingImage_trans));
    [NSWE_Cis, NSWE_Trans]  = getNSWEofOverlappingImageparts(Image_Cis, Image_Trans);
else
    createShiftedSegmentation = false;
end

% Obtain Segmentation from reference acquisition
if couldFindSameSite_segmentation == true;
    OrigSegmentation = double(imread(fullfile(SegmentationFolder_trans, strCorrespondingSegmentation_trans)));
else
    createShiftedSegmentation = false;
    OrigSegmentation = zeros(size(Image_Trans));
end


if createShiftedSegmentation == true  
    ShiftedSegmentationImage = zeros(size(Image_Cis),'double');

    % SHIFT ORIGINAL SEGMENTATION %
    ShiftedSegmentationImage(NSWE_Cis(1):NSWE_Cis(2), NSWE_Cis(3):NSWE_Cis(4)) = ...
        OrigSegmentation(NSWE_Trans(1):NSWE_Trans(2), NSWE_Trans(3):NSWE_Trans(4));
        
    % RELABEL  %
    
    % Since several CellProfiler modules will introduce error, if labelling
    % is not continous, make continuous labels;
    
    uOrigLabel  = unique(OrigSegmentation);
    uShiftedLabel = unique(ShiftedSegmentationImage);
    
    if ~(any(uOrigLabel>0) && any(uShiftedLabel>0))  % in case that at least one segmentation is empty
        createShiftedSegmentation = false;
    else
        uOrigLabel = uOrigLabel(uOrigLabel~=0);   % ignore background
        
        lFun = @(x) x(:);
        if ~isequal(lFun(uOrigLabel), lFun(unique(1:length(uOrigLabel)))) % check if labelling has been continous in original segmentation
            error('Labelling within original segmentation is not continuous. This will cause wrong data (but no compuational errors) in several standard modules of CellProfiler, and is therefore strongly advised not to do! Please reconsider your full analysis, instead of making your pipeline of ignore this discontinuous labelling!');
        else
            uShiftedLabel = uShiftedLabel(uShiftedLabel~=0); % ignore background
            elementsInShiftedLabel = length(uShiftedLabel);
            
            uShiftedLabel = sort(uShiftedLabel,'ascend'); % do relabelling
            relabelledSegmentation = zeros(size(ShiftedSegmentationImage));
            for j=1:elementsInShiftedLabel
                c = uShiftedLabel(j);
                f = ShiftedSegmentationImage == c;
                relabelledSegmentation(f) = j;
            end
            correspondingOrigLabelForShifted = uShiftedLabel;
        end
    end
    
end

if createShiftedSegmentation == false  % create blank image (e.g.: if no object has been present, or some reference could not be found)
    relabelledSegmentation = zeros(size(Image_Cis),'double'); % default values
    correspondingOrigLabelForShifted = 0;
    Centroid = [0 0];
else % create further measurements, if valid objects are present
    tmp = regionprops(relabelledSegmentation,'Centroid');
    Centroid = cat(1,tmp.Centroid);
end


% %%%%%%%%%%%%%%
% %%% OUTPUT %%%
% %%%%%%%%%%%%%%

%% Saves the segmented image, not edited for objects along the edges or
%%% for size, to the handles structure.
fieldname = ['UneditedSegmented',ObjectName];
handles.Pipeline.(fieldname) = relabelledSegmentation;

%%% Saves the segmented image, only edited for small objects, to the
%%% handles structure.
fieldname = ['SmallRemovedSegmented',ObjectName];
handles.Pipeline.(fieldname) = relabelledSegmentation;

%%% Saves the final segmented label matrix image to the handles structure.
fieldname = ['Segmented',ObjectName];
handles.Pipeline.(fieldname) = relabelledSegmentation;


%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,ObjectName));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
ObjCount = max(relabelledSegmentation(:));
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;


%%% Saves the location of each segmented object
% Follow CP convention for empty images (e.g.: as in IdentifySecondary
% module)
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% save relation to original object ID to handles
handles.Measurements.(ObjectName).UnshiftedObjectIdFeatures{handles.Current.SetBeingAnalyzed} = 'UnshiftedObjectId';
handles.Measurements.(ObjectName).UnshiftedObjectId{handles.Current.SetBeingAnalyzed} = correspondingOrigLabelForShifted;



%%%%%%%%%%%%%%%
%%% DISPLAY %%%
%%%%%%%%%%%%%%%


drawnow

if ~CPisHeadless()
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber)
        
        CPfigure(handles,'Image',ThisModuleFigureNumber);        
        
        %%% A subplot of the figure window is set to display the original image.
        subplot(2,2,1);
        CPimagesc(Image_Cis,handles);
        title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        %%% A subplot of the figure window is set to display the colored label
        %%% matrix image.
        subplot(2,2,2);
        if couldFindSameSite_image == true
            im = Image_Trans;
        else
            im = zeros(size(Image_Cis));
        end
        CPimagesc(im,handles);
        title('Image from reference plate');    
        
        subplot(2,2,3);
        ObjectOutlinesOnOrigImage = combineImageAndSegmentation(Image_Cis, relabelledSegmentation);
        CPimagesc(ObjectOutlinesOnOrigImage,handles);
        title([ObjectName, ' Outlines on Input Image']);
        
        
        subplot(2,2,4);
        CPimagesc(relabelledSegmentation,handles);
        colormap('jet');
        clear ColoredLabelMatrixImage
        title(['Segmentation of ',ObjectName]);
        
        drawnow
    end
end


end


function makeExternalBuffer(handles, strTransPlate, filterForImagesOfTrans, ObjectName)
% note that preferentially all data would be stored in pipeline. However
% storing custom fields in handles.pipelines can lead to various errors.
% While storing data in handles.measurements is technically possible there
% would be various mistakes occuring if order of sites becomes rearranged
% in batch files. Since there will always only be one segmentation of a
% given name, there can be a buffer, that links to that name, and which is
% overwritten, if a new pipeline is made


[TiffFolder_trans, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate);



% Images
eB.strImages_trans = getFilesAndDirectories(TiffFolder_trans, filterForImagesOfTrans);
[eB.Row_image_trans, eB.intColumn_image_trans, eB.intImagePosition_image_trans] = cellfun(@(x) MetaFromImageName(x), eB.strImages_trans, 'UniformOutput', true);

% Segmentations 
filterForSegmentationsOfTrans = ['Segmented' ObjectName '\.'];
eB.strSegmentations_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans);
[eB.Row_segmentations_trans, eB.intColumn_segmentations_trans, eB.intImagePosition_segmentations_trans] = cellfun(@(x) MetaFromImageName(x), eB.strSegmentations_trans, 'UniformOutput', true);


strBufferFile = getFileNameOfBuffer(handles, ObjectName);

if any(fileattrib(strBufferFile))
    [~, ex] = fileparts(strBufferFile);
    fprintf([ex ' already exists. It will be overwritten with current one.\n']);
end

save(strBufferFile, 'eB');
end

function strBufferFile = getFileNameOfBuffer(handles, ObjectName)

outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' ObjectName '.mat'];
strBufferFile = fullfile(outDir, ex );

end




function [strCorrespondingImage_trans, couldFindSameSite_image, strCorrespondingSegmentation_trans, couldFindSameSite_segmentation] = getFileNameViaReferenceFile(handles, ObjectName, OrigImageName)
% Avoid repeated, independent access to buffer file by reading segmentation
% and images after same loading event. Note that the code thus becomes more
% ugly and less readable.


SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){SetBeingAnalyzed};

strBufferFile = getFileNameOfBuffer(handles, ObjectName);

if ~any(fileattrib(strBufferFile))
   error('Could not find buffer file for object, that should be created during first cycle');    
else
   load(strBufferFile); 
end

[Row_cis, intColumn_cis, intImagePosition_cis] = MetaFromImageName(ImageName_cis);

% get corresponding image from trans
correspondsToSameSite_image = eB.Row_image_trans == Row_cis & eB.intColumn_image_trans == intColumn_cis & eB.intImagePosition_image_trans == intImagePosition_cis;
if sum(correspondsToSameSite_image) == 0
    couldFindSameSite_image = false;
    strCorrespondingImage_trans = '';
elseif sum(correspondsToSameSite_image) == 1;
    couldFindSameSite_image = true;
    strCorrespondingImage_trans = eB.strImages_trans{correspondsToSameSite_image};
else
    error('Could not unambiguously find corresponding image of other dataset. Please set more stringent filters.');
end

% get corresponding segmentation from trans
correspondsToSameSite_segmentation = eB.Row_segmentations_trans == Row_cis & eB.intColumn_segmentations_trans == intColumn_cis & eB.intImagePosition_segmentations_trans == intImagePosition_cis;
if sum(correspondsToSameSite_segmentation) == 0
    couldFindSameSite_segmentation = false;
    strCorrespondingSegmentation_trans = '';
elseif sum(correspondsToSameSite_segmentation) == 1;
    couldFindSameSite_segmentation = true;
    strCorrespondingSegmentation_trans = eB.strSegmentations_trans{correspondsToSameSite_segmentation};
else
    error('Could not unambiguously find corresponding image of other dataset. Please set more stringent filters.');
end

end



function [TiffFolder_trans, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate)


if any(strfind(strTransPlate, [filesep 'TIFF']))
    error('Reference directory must refer to a plate folder, not the TIFF folder');
end

if any(strfind(strTransPlate, [filesep 'SEGMENTATION']))
    error('Reference directory must refer to a plate folder, not the SEGMENTATION folder');
end

TiffFolder_trans = fullfile(strTransPlate, 'TIFF');
if ~any(fileattrib(TiffFolder_trans))
    error('Could not find TIFF folder of other plate');
end

SegmentationFolder_trans = fullfile(strTransPlate, 'SEGMENTATION');
if ~any(fileattrib(SegmentationFolder_trans))
    error('Could not find SEGMENTATION folder of other plate');
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% COORDINATES OF OVERLAPPING IMAGE   %%%%%%%%%%%%%%%%%%%

function [NSWE_Cis, NSWE_Trans] = getNSWEofOverlappingImageparts(Image_Cis, Image_Trans)
% upon providing two images (Image_Cis, Image_Trans), The top row (N -orth),
% lowest row (S -outh),  left column (W -est) and right column (E -east) of
% the overlapping region are returned;

[imageHeight, imageWidth] = size(Image_Cis);

fourierTrans = fft2(Image_Trans);
fourierCis = fft2(Image_Cis);
[registrationOutput] = dftregistration(fourierTrans,fourierCis,1);

SurplusRows = registrationOutput(3);
SurplusColumns = registrationOutput(4);

% Coordinates of for Rows

if SurplusRows > 0
    
    N_cis = 1;
    S_cis = imageHeight - SurplusRows;
    
    N_trans = 1 + SurplusRows;
    S_trans = imageHeight;
    
elseif SurplusRows == 0
    
    N_cis = 1;
    S_cis = imageHeight;
    N_trans = 1;
    S_trans = imageHeight;
    
elseif SurplusRows < 0
    
    N_cis = -SurplusRows + 1;
    S_cis = imageHeight;
    
    N_trans = 1;
    S_trans = imageHeight + SurplusRows;
    
else
    error('something terribly wrong');
end

% Coordinates of for Columns

if SurplusColumns > 0
    
    W_cis = 1;
    E_cis = imageWidth - SurplusColumns;
    
    W_trans = 1 + SurplusColumns;
    E_trans = imageWidth;
    
elseif SurplusRows == 0
    
    W_cis = 1;
    E_cis = imageWidth;
    
    W_trans = 1;
    E_trans = imageWidth;
    
elseif SurplusRows < 0
    
    W_cis = -SurplusColumns + 1;
    E_cis = imageWidth;
    
    W_trans = 1;
    E_trans = imageWidth + SurplusColumns;
    
else
    error('something terribly wrong');
end

NSWE_Cis = [N_cis, S_cis, W_cis, E_cis];
NSWE_Trans = [N_trans, S_trans, W_trans, E_trans];


end



function ObjectOutlinesOnOrigImage = combineImageAndSegmentation(OrigImage, FinalLabelMatrixImage)
% This code is copy/pasted from original IdentifySecondary


MaxFilteredImage = ordfilt2(FinalLabelMatrixImage,9,ones(3,3),'symmetric');
%%% Determines the outlines.
IntensityOutlines = FinalLabelMatrixImage - MaxFilteredImage;
%%% [PLab] hack.s ave memory.
clear MaxFilteredImage;
%%% Converts to logical.
warning off MATLAB:conversionToLogical
LogicalOutlines = logical(IntensityOutlines);
%%% [PLab] hack.s ave memory.
clear IntensityOutlines;
warning on MATLAB:conversionToLogical

ObjectOutlinesOnOrigImage = OrigImage;

qmin = quantile(ObjectOutlinesOnOrigImage(:),0.01); % [TS] Keep intention of original code by rescaling without outliers
qmax = quantile(ObjectOutlinesOnOrigImage(:),0.99);
ObjectOutlinesOnOrigImage(ObjectOutlinesOnOrigImage < qmin) = qmin;
ObjectOutlinesOnOrigImage(ObjectOutlinesOnOrigImage > qmax) = qmax;
ObjectOutlinesOnOrigImage = (ObjectOutlinesOnOrigImage-qmin) ./ (qmax-qmin);

ObjectOutlinesOnOrigImage = insertRedLine(ObjectOutlinesOnOrigImage, LogicalOutlines);

end


function RGBimage = insertRedLine(image, bwLineImage)

if size(image,3) == 1
    image = repmat(image,[1 1 3]);
end

r = image(:,:,1);
g = image(:,:,2);
b = image(:,:,3);

r(bwLineImage) = 1;
g(bwLineImage) = 0;
b(bwLineImage) = 0;

RGBimage = cat(3, r, g, b);
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% CODE FROM MATLABCENTRAL   %%%%%%%%%%%%%%%%%%%

function [output, Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk
% and James R. Fienup.
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
% "Efficient subpixel image registration algorithms," Opt. Lett. 33,
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image,
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register,
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0,
    CCmax = sum(sum(buf1ft.*conj(buf2ft)));
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    output=[error,diffphase];
    
    % Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
    % peak
elseif usfac == 1,
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc);
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    md2 = fix(m/2);
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
    % Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
    
    % Compute crosscorrelation and locate the peak
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;
    
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac;
        col_shift = round(col_shift*usfac)/usfac;
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid
        [max1,loc1] = max(CC);
        [max2,loc2] = max(max1);
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;
        
        % If upsampling = 2, no additional pixel shift refinement
    else
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end
return
end


function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1)
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end %#ok<*EXIST>
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return

end
