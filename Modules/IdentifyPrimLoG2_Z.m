function handles = IdentifyPrimLoG(handles)

% Help for the Identify Primary LoG module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
%
% Changed by Berend Snijder to try and make iBRAIN / old CellProfiler
% compatible!!!
%
% Identifies the centers of blob-like primary objects.  The result
% consists of only a single pixel per object, located near the center
% of the object.
%
% *************************************************************************
%
% This module identifies the centers of blob-like primary objects
% (e.g. nuclei) in grayscale images that show bright objects on a dark
% background.  When the objects of interest are fairly round and of
% even size, this module may be more sensitive than the methods in
% IdentifyPrimAutomatic and therefore detect objects that would
% otherwise be lost.
% 
% The result consists of only a single pixel per object, located near
% the center of the object; the IdentifySecondary module can be used
% to fill out the object based on this center point.
%
% SETTINGS:
%
% The radius parameter should be set to the approximate radius of the
% objects of interest.  The algorithm is not very sensitive to this
% parameter.
%
% The threshold parameter informs the algorithm how inclusive to be when
% looking for objects.  Internally, each potential object is assigned
% a score that depends on both how bright the object is and how
% blob-like its shape is.  Only objects that score above the threshold
% are returned.  If the thresold is too high, objects will be lost; 
% if it is too low, spurious objects will be found. The threshold 
% can be determined experimentally, but the 'Automatic' setting 
% will make a guess using RobustBackground Global's thresholding 
% method on the transformed image.  RobustBackground is useful because it
% makes little assumption of the intensity histogram, and thus 
% can be protective against out-of-focus or empty images.  If you want the 
% threshold to be consistent across images, then you can use the threshold found by 
% the 'Automatic' setting as a starting point for manual threshold input adjustment.
% Also, if the threshold is consistently high or low, then you can adjust 
% by a multiplicative correction factor by inserting it after a comma, e.g.
% "Automatic,1.5". 
%
% ALGORITHM DETAILS:
%
% The module works by convolving the image with the Laplacian of
% Gaussian (LoG) kernel.  This is equivalent to convolving with the
% Gaussian kernel and then with the Laplace operator.  The regional
% maxima in the filter response that exceed the specificed threshold
% are identified as objects.  The radius parameter specifies the width
% of the kernel.
%
% Ultimately, this module will become an option in
% IdentifyPrimAutomatic, so that its options for maxima suppression
% and finding edges between clumps can be used.
%
% $Revision: 7941 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

% Notes for PyCP
%
% FROM CP ToDo:
% Anne 2008_01_30: Merge IdentifyPrimaryLoG into the regular IdentifyPrimAutomatic module. 
% In particular, be sure that the help describes under what conditions the different options 
% are useful. Think carefully about how to add the variable that is LoG-specific to the module.
%
% Anne 2008_05_12: Describe what the LoG is doing - we think that we are looking for 
% minima (or maxima) of the LoG which makes it a maxima-minima finder of the original image, 
% whereas many people look for zero crossings of the LoG which would be using it as an edge 
% detector. Is that right?
%
% David 2009_04_17:  As the above older comments say, this module should be integrated into
% IDPrimAuto.  However, this ID module is different than other primary segmentation modules
% in that it only finds the center pixel of objects and must utilize a subsequent IDSecondary 
% after to grow objects.  So we either:
% (1) Treat LoG as a special thresholding method that only ouputs single pixel objects
% or
% (2) Treat LoG as a declumping method.  In this case, another thresholding method would define 
% foreground/background, and LoG would find single pixels within the foreground.
% 
% In either case, we would need to decide whether we automatically apply a watershed/propagation
% after the initial single-pixel finding method.  I prefer (2) above and would opt for
% automatically propagating the single-pixels within the foreground objects, since that is almost 
% always done anyway, and would save the step of adding an IDSecondary.
%
% Settings:
% The first two settings map obviously.
% The diameter parameter is only single, so that if LoG is chosen, one of the diameter boxes 
% should gray-out.
% The threshold parameter is very sensitive, and the user was blind to what this should be at first,
% so the 'Automatic' threshold was added recently to use Otsu to guess.  This functionality would be
% equivalent to the 'Automatic' setting in the existing declumping settings.

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images you want to process?
%choiceVAR01 = /
%infotypeVAR01 = imagegroup
ImageName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the images you want to process?
%choiceVAR02 = /
%infotypeVAR02 = imagegroup
ImageName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What did you call the images you want to process?
%choiceVAR03 = /
%infotypeVAR03 = imagegroup
ImageName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = What did you call the images you want to process?
%choiceVAR04 = /
%infotypeVAR04 = imagegroup
ImageName{4} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = What did you call the images you want to process?
%choiceVAR05 = /
%infotypeVAR05 = imagegroup
ImageName{5} = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu

%textVAR06 = What did you call the images you want to process?
%choiceVAR06 = /
%infotypeVAR06 = imagegroup
ImageName{6} = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu


%textVAR07 = What do you want to call the objects identified by this module?
%defaultVAR07 = Nuclei
%infotypeVAR07 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Typical diameter of objects, in pixel units:
%defaultVAR08 = 10
Radius = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = Score threshold for match.  Enter a number, leave as 'Automatic', or adjust the Automatic threshold with a multiplicative correction factor (separating comma is necessary), e.g. 'Automatic,1.2'.
%defaultVAR09 = Automatic
ThresholdStr = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ImageName(1 < cellfun(@length ,ImageName))

OrigImage = cellfun(@(x) CPretrieveimage(handles,x,ModuleName,'MustBeGray','CheckScale'),ImageName(1 < cellfun(@length ,ImageName)), 'UniformOutput',false);
OrigImage = cellfun(@double,OrigImage,'UniformOutput',false);

Radius = str2double(Radius);

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow
%%=========================================================================

[ExpFinalImage, Threshold] = cellfun(@(x) ReturnSpotImage(x,Radius,ThresholdStr,handles,ModuleName),OrigImage,'UniformOutput',false);
ExpFinalImage = cellfun(@double,ExpFinalImage,'UniformOutput',false);
Threshold = Threshold{1};
matImtemp = nan(size(ExpFinalImage{1},1),size(ExpFinalImage{1},2),length(ExpFinalImage));

for i = 1:length(ExpFinalImage)
    matImtemp(:,:,i) = ExpFinalImage{i};    
end

%%%% till here is ok
%%size( ExpFinalImage{2})

matImtemp = sum(matImtemp,3);
matImtemp2 = matImtemp > 1;

% % cellCeantroids = (struct2cell(regionprops(matImtemp2, 'centroid')));
% % matCeantroids = nan(length(cellCeantroids),2);
% % bw = false(size(matImtemp2));
% % for i = 1:length(cellCeantroids)
% % matCeantroids(i,:) = round(cellCeantroids{i});
% % bw(matCeantroids(i,1),matCeantroids(i,2)) = true;
% % end
% % 
% % sum(matImtemp(:)==nan)
% % figure()
% % imagesc(matImtemp2)
% % colorbar
% % 
% % 
% % figure()
% % imagesc(matImtemp(:,:,1))
% % colorbar
% % 
% % 
% % figure()
% % imagesc(bw)
% % colorbar
% % 
% % 
% % ExpFinalImage{i}
% % 
% % 
% % max(sum(matImtemp,3))
% % 



bw = matImtemp2;
FinalLabelMatrixImage = bwlabel(bw);

OrigImage = OrigImage{1};

%%========================================================================= 
% The dilated mask is used only for visualization.
dilated = imdilate(bw, strel('disk', 2));
vislabel = bwlabel(dilated);
% r = OrigImage;
% g = OrigImage;
% b = OrigImage;
r = (OrigImage - min(OrigImage(:))) / max(OrigImage(:));
g = (OrigImage - min(OrigImage(:))) / max(OrigImage(:));
b = (OrigImage - min(OrigImage(:))) / max(OrigImage(:));

r(dilated) = max(r(:));
g(dilated) = 0;
b(dilated) = 0;
visRGB = cat(3, r, g, b);

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
  h_fig = CPfigure(handles,'Image',ThisModuleFigureNumber);
  subplot(2,2,1)
  hImage = CPimagesc(visRGB, handles);
  hAx = gca;
  title(hAx,[ObjectName, ' cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);

  subplot(2,2,2)
  hImage2 = CPimagesc(OrigImage, handles);
  
  subplot(2,2,4)
  hImage2 = CPimagesc(dilated, handles);

  
%   ud(1).img = visRGB;
%   ud(2).img = ac_scaled;
%   ud(3).img = OrigImage;
%   ud(1).title = [ObjectName ' , cycle # ',num2str(handles.Current.SetBeingAnalyzed)];
%   ud(2).title = ['Laplacian of Gaussian transformed ' ObjectName ', cycle # ',num2str(handles.Current.SetBeingAnalyzed)];
%   ud(3).title = ['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)];
%   uicontrol(h_fig, 'Style', 'popup',...
%                     'String', 'Objects on Original Image|Laplacian of Gaussian Transformed|Input Image',...
%                     'UserData',ud,...
%                     'units','normalized',...
%                     'position',[.01 .95 .25 .04],...
%                     'backgroundcolor',[.7 .7 .9],...
%                     'tag','PopupImage',...
%                     'Callback', @CP_ImagePopupmenu_Callback);
  
  text(0.1,-0.08,...
      ['Threshold: ' num2str(mean(Threshold(:))) ', Number of objects: ' num2str(sum(bw(:)))],...
      'Color','black',...
      'fontsize',handles.Preferences.FontSize,...
      'Units','Normalized',...
      'Parent',hAx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

prefixes = {'Segmented', 'SmallRemovedSegmented'};
for i=1:length(prefixes)
  prefix = prefixes{i};
  fieldname = [prefix, ObjectName];
  handles = CPaddimages(handles,fieldname,FinalLabelMatrixImage);
end

handles = CPsaveObjectCount(handles, ObjectName, FinalLabelMatrixImage);
handles = CPsaveObjectLocations(handles, ObjectName, FinalLabelMatrixImage);


function  [ExpFinalImage, Threshold] = ReturnSpotImage(OrigImage,Radius,ThresholdStr,handles,ModuleName)

%% here calculate the LoG of each image and detect the spots

im = double(OrigImage) - double(min(OrigImage(:)));
if any(im(:))
    im = im / max(im(:));
end

% % Set regions outside of CropMasks equal to 0
% MaskFieldname = ['CropMask', ImageName];
% if CPisimageinpipeline(handles, MaskFieldname)
%     %%% Retrieves previously selected cropping mask from handles
%     %%% structure.
%     PreviousCropMask = CPretrieveimage(handles,MaskFieldname,ModuleName);
%     try 
%         im(~PreviousCropMask) = 0;
%     catch
%         error('The image in which you want to identify objects has been cropped, but there was a problem recognizing the cropping pattern.');
%     end
% end

ac = lapofgau(1 - im, 5);
%%figure()
%%imagesc(ac)

% IntermediateImage = ac - min(ac(:));
% %%% The maximum of the image is brought to 1.
% ac_scaled = IntermediateImage ./ max(IntermediateImage(:));

ac_min = min(ac(:));
ac_range = max(ac(:)) - ac_min;
ac_scaled = (ac - ac_min) ./ ac_range;

if ~isnan(str2double(ThresholdStr))
    Threshold = str2double(ThresholdStr);
elseif strcmpi('Automatic',ThresholdStr(1:9))
    if numel(ThresholdStr) > length('Automatic')
        assert(strcmp(ThresholdStr(10),','),'A comma must follow ''Automatic'' if a threshold correction factor is used.')
        Threshold_correction = str2double(ThresholdStr(length('Automatic')+1:end));
    else
        Threshold_correction = 1;
    end
%     [handles,Threshold_scaled] = CPthreshold(handles,'RobustBackground Global',0,'0','1',Threshold_correction,ac_scaled,'LoG',ModuleName);
    disp('BEREND: RidlerCalvard Adaptive')
%     'Otsu','MoG','Background','RobustBackground','RidlerCalvard','Kapur'}
    [handles,Threshold_scaled] = CPthreshold(handles,'RobustBackground Global',0,'0','1',Threshold_correction,ac_scaled,'LoG',ModuleName);

    %% Un-scale threshold
    Threshold = (Threshold_scaled .* ac_range) + ac_min;
end

%%size(ac)
%%size(Threshold)

if numel(Threshold)==1
    ac(ac < Threshold) = Threshold;
else
    ac(ac < Threshold) = Threshold(ac < Threshold);
end
ac = ac - Threshold;

%%size(ac)
%%size(Threshold)


bw = false(size(im));
if any(ac(:))
    indices = find(imregionalmax(ac));
    maxima = sortrows([indices ac(indices)], -2);
    bw(maxima(:,1)) = true;
end



%%figure()
%%imagesc(bw)

%% Mask final outcome
if exist('PreviousCropMask','var')
    bw = bw & PreviousCropMask;
end

ExpFinalImage = bw;
bwf = nan(size(bw,1),size(bw,2),3);
bwf(:,:,1) = double(bw);
bwf(2:end,1:end,2) = bw(1:end-1,1:end);
bwf(1:end,2:end,3) = bw(1:end,1:end-1);
% bwf(1:end-1,1:end,4) = bw(2:end,1:end);
% bwf(1:end,1:end-1,5) = bw(1:end,2:end);
bwf = sum(bwf,3);
%ExpFinalImage = bwmorph(bw, 'thicken', 1);

ExpFinalImage = bwf;


% % FinalLabelMatrixImage = bwlabel(bw);
% % 
% %  
% % ExpFinalLabelMatrixImage = bwlabel(ExpFinalLabelMatrixImage);
% % 
% % figure()
% % imagesc(ExpFinalLabelMatrixImage == ExpFinalLabelMatrixImage)
% % numTempIndex = ExpFinalLabelMatrixImage(:) == 0;
% % ExpFinalLabelMatrixImage(numTempIndex) = nan;
% % matTemIm = ExpFinalLabelMatrixImage == ExpFinalLabelMatrixImage



function f = lapofgau(im, s)
% im: image matrix (2 dimensional)
% s: filter width
% f: filter output.
% Author: Baris Sumengen - sumengen@ece.ucsb.edu

sigma = (s-1)/3;
op = fspecial('log',s,sigma); 
op = op - sum(op(:))/numel(op); % make the op to sum to zero

%% Pad image to fix border artifact
padsize = ceil(s./2);
im = padarray(im,[padsize padsize],'replicate');

%% Use 'single' here, since imfilter is optimized on Intel chips for
%% certain array classes (see Matlab documentation)
f = double(imfilter(single(im),op));

%% Crop pad
f = f(padsize+1:end-padsize,padsize+1:end-padsize);

function handles = CPaddimages(handles, varargin)
% Add images to the handles.Pipeline structure.
% Location will be "handles.Pipeline.ImageName".

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 7247 $

% Parse out varargin. The added data can be numeric, logical or a structure
% (e.g, movie)
if mod(length(varargin),2) ~= 0 || ...
   ~all(cellfun(@ischar,varargin(1:2:end)) & ...
   (cellfun(@isnumeric,varargin(2:2:end)) | cellfun(@islogical,varargin(2:2:end)) | cellfun(@isstruct,varargin(2:2:end))))
    error('The argument list must be of the form: ''ImageName1'', ImageData1, etc');
else
    ImageName = varargin(1:2:end);
    ImageData = varargin(2:2:end);
end

CPvalidfieldname(ImageName);

% Checks have passed, add the data
if ~isfield(handles.Pipeline,'ImageGroupFields')
    % If no image groups, add to the handles.Pipeline structure
    for i = 1:length(ImageName)
        handles.Pipeline.(ImageName{i}) = ImageData{i};
    end
else
    % If no image groups, add to the appropriate
    % handles.Pipeline.GroupFileList structure
    for i = 1:length(ImageName)
        handles.Pipeline.GroupFileList{handles.Pipeline.CurrentImageGroupID}.(ImageName{i}) = ImageData{i};
    end
end


function CPvalidfieldname(fieldname)

% Throw an error if the field name does not start with
% an alphabetic character, if it contains characters other
% than alphanumerics and underbar or if it is more than 63
% characters long.

% $Revision: 5791 $

if length(fieldname) > namelengthmax
    error(['The field name, "',fieldname,'", is more than ',num2str(namelengthmax),' characters long.']);
end
if isempty(regexp(fieldname,'^[A-Za-z]', 'once' ))
    error(['The field name, "',fieldname,'", does not start with an alphabetic character.']);
end
if isempty(regexp(fieldname,'^[A-Za-z][A-Za-z0-9_]{0,62}$', 'once' ))
    error(['The field name, "',fieldname,'", contains characters other than alphanumerics and "_"']);
end

% CPSAVEOBJECTCOUNT Save the count of segmented objects.
%   The function returns a new version of the handles structure, in which
%   the number of segmented objects has been saved.
%
%   Example:
%      handles = CPsaveObjectCount(handles, 'Cells', labelMatrix)
%      creates handles.Measurements.Cells{i}.Count_Cells.
function handles = CPsaveObjectCount(handles, objectName, labels)
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2008.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5777 $
% handles = CPaddmeasurements(handles, 'Image', CPjoinstrings('Count', objectName), CPjoinstrings('Count', objectName), max(labels(:)));
% which('CPaddmeasurements')
handles = CPaddmeasurements(handles, 'Image', CPjoinstrings('Count', objectName), max(labels(:)));

            
            
function handles = CPsaveObjectLocations(handles, objectName, labels)
% CPSAVEOBJECTLOCATIONS Save the location of each segmented object.
%   The function returns a new version of the handles structure, in
%   which the location of each segmented object has been saved.
%
%   Example:
%      handles = CPsaveObjectLocations(handles, 'Cells', cellLabelMatrix)
%      creates handles.Measurements.Cells{1}.Location_Center_X and
%      handles.Measurements.Cells{1}.Location_Center_Y.
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2008.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5777 $
tmp = regionprops(labels, 'Centroid');
centroids = cat(1,tmp.Centroid);
if isempty(centroids)
  centroids = zeros(0,2);
end
% % handles = CPaddmeasurements(handles,Object,   Measure,                Feature,            Data)
% handles = CPaddmeasurements(handles, objectName, 'Location_Center_X', 'Location_Center_X', centroids(:,1));
% handles = CPaddmeasurements(handles, objectName, 'Location_Center_Y', 'Location_Center_Y', centroids(:,2));
% handles = CPaddmeasurements(handles,Object,   Measure,                Feature,            Data)
% which('CPaddmeasurements')
handles = CPaddmeasurements(handles, objectName, 'Location_Center_X', centroids(:,1));
handles = CPaddmeasurements(handles, objectName, 'Location_Center_Y', centroids(:,2));
                        
function string = CPjoinstrings(varargin)
%CPjoinstrings Build underscore-separated string from parts.
%
%   CPjoinstrings(D1,D2, ... ) builds a string from 
%   D1,D2, etc specified.  This is conceptually equivalent to
%
%      F = [D1 '_' D2 '_' ... '_' DN] 
%
%   Care is taken to handle the cases where the directory
%   parts D1, D2, etc. may contain an empty string, in which case there are
%   not two consecutive underscores output.
%
%   Examples
%   See also FILESEP, PATHSEP, FILEPARTS.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2008.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5777 $

error(nargchk(2, Inf, nargin, 'struct'));

sepchar = '_';
string = varargin{1};

for i=2:nargin,
   part = varargin{i};
   if isempty(string) || isempty(part)
      string = [string part];
   else
      % Handle the three possible cases
      if (string(end)==sepchar) && (part(1)==sepchar),
         string = [string part(2:end)];
      elseif (string(end)==sepchar) || (part(1)==sepchar )
         string = [string part];
      else
         string = [string sepchar part];
      end
   end
end



function handles = CPaddmeasurements(handles, ObjectName, FeatureName, Data, ImageSetNumber)



% Add measurements of a feature to the handles.Measurements structure.
% Location will be "handles.Measurements.ObjectName.FeatureName".
% ObjectName can be "Image".  
%
% Data can be multiple doubles, or a single string (only if ObjectName is "Image").

%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 8037 $
% fprintf('Berend was here in %s\n',which('CPaddmeasurements'))

if nargin < 5,
    ImageSetNumber = handles.Current.SetBeingAnalyzed;
end

% Check that either this is a new measurement being added in the first
% set, or an old measurement being appended to in a later set.
if isscalar(ImageSetNumber)
	FirstSet = ImageSetNumber == 1;
elseif isvector(ImageSetNumber)
	FirstSet = ImageSetNumber(1) == 1;
end
OldMeasurement = ...
    isfield(handles.Measurements, ObjectName) && ...
    isfield(handles.Measurements.(ObjectName), FeatureName);
BatchProcessing = isfield(handles.Current, 'BatchInfo');

CPvalidfieldname(FeatureName)

%%% Don't allow overwriting of measurements, *except* in batch processing.
if (FirstSet && OldMeasurement && ~BatchProcessing),
    error(['Image processing was canceled because you are attempting to recreate the same measurements, please remove redundant module (#', handles.Current.CurrentModuleNumber, ').']);
end

if (~FirstSet) && (~OldMeasurement) && (~ strcmp(ObjectName, 'Experiment')),
    error(['This should not happen.  CellProfiler Coding Error.  Attempting to add new measurement ', ObjectName, '.',  FeatureName, ' in set ', int2str(ImageSetNumber) ' that was not added in first set.']);
end

%%% Verify we can add this type of Measurement to this type of object
if ischar(Data) && (~ strcmp(ObjectName, 'Image')),
    error(['This should not happen.  CellProfiler Coding Error.  Attempting to add string measurement to non-image ', ObjectName, '.', FeatureName]);
elseif ~strcmp(ObjectName, 'Image') && ~isvector(Data) && ~isempty(Data)
    error(['This should not happen.  CellProfiler Coding Error.  Attempting to add multidimensional (', int2str(size(Data)), ') measurement ', ObjectName, '.', FeatureName]);
elseif strcmp(ObjectName, 'Image') && isnumeric(Data) && ~isscalar(Data),
    error(['This should not happen.  CellProfiler Coding Error.  Attempting to add non-scalar (', int2str(size(Data)), ') measurement to ', ObjectName, '.', FeatureName]);
end


%%% Checks have passed, add the data.
if strcmp(ObjectName, 'Experiment'),
    handles.Measurements.(ObjectName).(FeatureName) = Data;
else
	if isscalar(ImageSetNumber)
		handles.Measurements.(ObjectName).(FeatureName){ImageSetNumber} = Data;
	elseif isvector(ImageSetNumber)
		ImageSetNumber = ImageSetNumber(:)';
		handles.Measurements.(ObjectName).(FeatureName)(ImageSetNumber) = reshape(Data,size(ImageSetNumber));
	end
end




%%%%%%%%%%%%%%%%%%%%%%
function [handles,Threshold,varargout] = CPthreshold(handles,Threshold,pObject,MinimumThreshold,MaximumThreshold,ThresholdCorrection,OrigImage,ImageName,ModuleName,ObjectVar)
%
% Returns an automatically computed threshold, and if requested in
% varargout, the Otsu and Kapur measures of thresholding quality
% (weighted variance and sum of entropies, resp.).
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 8068 $

if nargin == 9
    ObjectVar = [];
end

isImageGroups = isfield(handles.Pipeline,'ImageGroupFields');
if ~isImageGroups
    SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
else
    SetBeingAnalyzed = handles.Pipeline.GroupFileList{handles.Pipeline.CurrentImageGroupID}.SetBeingAnalyzed;
end

% Make sure the input image is double precision
OrigImage = double(OrigImage);

%%% If we are running the Histogram data tool we do not want to limit the
%%% threshold with a maximum of 1 or minimum of 0; otherwise we check for
%%% values outside the range here.
if ~strcmpi('Histogram Data tool',ModuleName)
    %%% Check the MinimumThreshold entry. If no minimum threshold has been
    %%% set, set it to zero. Otherwise make sure that the user gave a valid
    %%% input.
    if strcmp(MinimumThreshold,'Do not use')
        MinimumThreshold = 0;
    else
        MinimumThreshold = str2double(MinimumThreshold);
        if isnan(MinimumThreshold) |  MinimumThreshold < 0 | MinimumThreshold > 1 %#ok Ignore MLint
            error(['The Minimum threshold entry in the ', ModuleName, ' module is invalid.'])
        end
    end

    if strcmp(MaximumThreshold,'Do not use')
        MaximumThreshold = 1;
    else
        MaximumThreshold = str2double(MaximumThreshold);
        if isnan(MaximumThreshold) | MaximumThreshold < 0 | MaximumThreshold > 1 %#ok Ignore MLint
            error(['The Maximum bound on the threshold in the ', ModuleName, ' module is invalid.'])
        end

        if MinimumThreshold > MaximumThreshold,
            error(['Min bound on the threshold larger than the Max bound on the threshold in the ', ModuleName, ' module.'])
        end
    end
end

%%% If the image was produced using a cropping mask, we do not want to
%%% include the Masked part in the calculation of the proper threshold,
%%% because there will be many zeros in the image.  So, we check to see
%%% whether there is a crop mask image in the handles structure that goes
%%% along with the image of interest.

%%% However, we do not need to retrieve the binary crop mask in a few
%%% cases: if the user is going to set the threshold interactively, if they
%%% are using All images together to calculate the threshold, or if they
%%% have manually entered a numerical value for the threshold.
if strcmp(Threshold,'Set interactively') || strcmp(Threshold,'All') || ~isempty(str2num(Threshold)) 
    %%% In these cases, don't do anything.
else
    fieldname = ['CropMask', ImageName];
    if CPisimageinpipeline(handles, fieldname)
        %%% Retrieves crop mask from handles structure. In some cases, it
        %%% might be a label matrix, if we are cropping based on objects
        %%% being present, so we make sure the resulting image is of the
        %%% logical class. This yields a warning at times, because the
        %%% image has values above one, so we temporarily turn the warning
        %%% off.
        OrigWarnState = warning('off','MATLAB:conversionToLogical');
        RetrievedCropMask = CPretrieveimage(handles,fieldname,ModuleName);
        RetrievedBinaryCropMask = logical(RetrievedCropMask);
        warning(OrigWarnState);
        %%% Handle the case where there are no pixels on in the mask, in
        %%% which case the threshold should be set to a numnber higher than
        %%% any pixels in the original image. In this case, the automatic
        %%% calculations are aborted below, because now Threshold = a
        %%% numerical value rather than a method name like 'Otsu', etc. So
        %%% from this point forward, in this case, the threshold is as if
        %%% entered manually.
        if (~any(RetrievedBinaryCropMask)),
            Threshold = 1;
        end
        %%% Checks whether the size of the RetrievedBinaryCropMask matches
        %%% the size of the OrigImage.
        if numel(OrigImage) == numel(RetrievedBinaryCropMask)
            BinaryCropMask = RetrievedBinaryCropMask;
            %%% Masks the image based on its BinaryCropMask and
            %%% simultaneously makes it a linear set of numbers.
            LinearMaskedImage = OrigImage(BinaryCropMask~=0);
        else
            Warning = CPwarndlg(['In CPthreshold, within the ',ModuleName,' module, the retrieved binary crop mask image (handles.Pipeline.',fieldname,') is not being used because it does not match the size of the original image(',ImageName,').']);
            %%% If the sizes do not match, then it is as if the crop mask
            %%% does not exist. I don't think this should ever actually
            %%% happen, but it might be needed for debugging.
        end
    end
    %%% If we have not masked the image for some reason, we need to create
    %%% the LinearMaskedImage variable, and simultaneously make it a linear
    %%% set of numbers.
    if ~exist('LinearMaskedImage','var')
        LinearMaskedImage = OrigImage(:);
    end
end

%%% STEP 1. Find threshold and apply to image
if ~isempty(strfind(Threshold,'Global')) || ~isempty(strfind(Threshold,'Adaptive')) || ~isempty(strfind(Threshold,'PerObject'))
    if ~isempty(strfind(Threshold,'Global'))
        MethodFlag = 0;
    elseif ~isempty(strfind(Threshold,'Adaptive'))
        MethodFlag = 1;
    elseif ~isempty(strfind(Threshold,'PerObject'))
        MethodFlag = 2;
    end
    %%% Chooses the first word of the method name (removing 'Global' or 'Adaptive' or 'PerObject').
    ThresholdMethod = strtok(Threshold);
    %%% Makes sure we are using an existing thresholding method.
    if isempty(strmatch(ThresholdMethod,{'Otsu','MoG','Background','RobustBackground','RidlerCalvard','Kapur'},'exact'))
        error(['The method chosen for thresholding, ',Threshold,', in the ',ModuleName,' module was not recognized by the CPthreshold subfunction. Adjustment to the code of CellProfiler is needed; sorry for the inconvenience.'])
    end
    
    %%% For all methods, Global or Adaptive or PerObject, we want to
    %%% calculate the global threshold. Sends the linear masked image to
    %%% the appropriate thresholding subfunction.
    eval(['Threshold = ',ThresholdMethod,'(LinearMaskedImage,handles,ImageName,pObject);']);

    %%% This evaluates to something like: Threshold =
    %%% Otsu(LinearMaskedImage,handles,ImageName,pObject);

    %%% The global threshold is used to constrain the Adaptive or PerObject
    %%% thresholds.
    GlobalThreshold = Threshold;

    %%% For Global, we are done. There are more steps involved for Adaptive
    %%% and PerObject methods.

    if MethodFlag == 1 %%% The Adaptive method.
        %%% Choose the block size that best covers the original image in
        %%% the sense that the number of extra rows and columns is minimal.
        %%% Get size of image
        [m,n] = size(OrigImage);
        %%% Deduce a suitable block size based on the image size and the
        %%% percentage of image covered by objects. We want blocks to be
        %%% big enough to contain both background and objects. The more
        %%% uneven the ratio between background pixels and object pixels
        %%% the larger the block size need to be. The minimum block size is
        %%% about 50x50 pixels. The line below divides the image in 10x10
        %%% blocks, and makes sure that the block size is at least 50x50
        %%% pixels.
        BlockSize = max(50,min(round(m/10),round(n/10)));
        %%% Calculates a range of acceptable block sizes as plus-minus 10%
        %%% of the suggested block size.
        BlockSizeRange = floor(1.1*BlockSize):-1:ceil(0.9*BlockSize);
        [ignore,index] = min(ceil(m./BlockSizeRange).*BlockSizeRange-m + ceil(n./BlockSizeRange).*BlockSizeRange-n); %#ok Ignore MLint
        BestBlockSize = BlockSizeRange(index);
        %%% Pads the image so that the blocks fit properly.
        RowsToAdd = BestBlockSize*ceil(m/BestBlockSize) - m;
        ColumnsToAdd = BestBlockSize*ceil(n/BestBlockSize) - n;
        RowsToAddPre = round(RowsToAdd/2);
        RowsToAddPost = RowsToAdd - RowsToAddPre;
        ColumnsToAddPre = round(ColumnsToAdd/2);
        ColumnsToAddPost = ColumnsToAdd - ColumnsToAddPre;
        PaddedImage = padarray(OrigImage,[RowsToAddPre ColumnsToAddPre],'replicate','pre');
        PaddedImage = padarray(PaddedImage,[RowsToAddPost ColumnsToAddPost],'replicate','post');
        PaddedImageandCropMask = PaddedImage;
        if exist('BinaryCropMask','var')
            %%% Pad the crop mask too.
            PaddedCropMask = padarray(BinaryCropMask,[RowsToAddPre ColumnsToAddPre],'replicate','pre');
            PaddedCropMask = padarray(PaddedCropMask,[RowsToAddPost ColumnsToAddPost],'replicate','post');
            %%% For the CPblkproc function, the original image and the crop
            %%% mask image (if it exists) must be combined into one.
            PaddedImageandCropMask(:,:,2) = PaddedCropMask;
            %%% And the Block must have two layers, too.
            Block = [BestBlockSize BestBlockSize 2];
            %%% Sends the linear masked image to the appropriate
            %%% thresholding subfunction, in blocks.
            eval(['Threshold = CPblkproc(PaddedImageandCropMask,Block,@',ThresholdMethod,',handles,ImageName,pObject);']);
            %%% This evaluates to something like: Threshold =
            %%% CPblkproc(PaddedImageandCropMask,Block,@Otsu,handles,ImageN
            %%% ame);
        else
            %%% If there is no crop mask, then we can simply use the
            %%% blkproc function rather than CPblkproc.
            Block = [BestBlockSize BestBlockSize];
            eval(['Threshold = blkproc(PaddedImageandCropMask,Block,@',ThresholdMethod,',handles,ImageName,pObject);']);
        end

        %%% Resizes the block-produced image to be the size of the padded
        %%% image. Bilinear prevents dipping below zero. The crop the image
        %%% get rid of the padding, to make the result the same size as the
        %%% original image.
        warning off MATLAB:divideByZero
        Threshold = imresize(Threshold, size(PaddedImage), 'bilinear');
        Threshold = Threshold(RowsToAddPre+1:end-RowsToAddPost,ColumnsToAddPre+1:end-ColumnsToAddPost);
        warning on MATLAB:divideByZero
    elseif MethodFlag == 2 %%% The PerObject method.
        %%% This method require the Retrieved CropMask, which should be a
        %%% label matrix of objects, where each object consists of an
        %%% integer that is its label.
        if ~exist('RetrievedCropMask','var')
            error(['Image processing was canceled in the ',ModuleName,' module because you have chosen to calculate the threshold on a per-object basis, but CellProfiler could not find the image of the objects you want to use.'])
        end
        %%% Initializes the Threshold variable (which will end up being the
        %%% same size as the original image).
        Threshold = ones(size(OrigImage));
        NumberOfLabelsInLabelMatrix = max(RetrievedCropMask(:));
        for i = 1:double(NumberOfLabelsInLabelMatrix)   % If NumberOfLabelsInLabelMatrix is logical, make it double
            %%% Chooses out the pixels in the orig image that correspond
            %%% with i in the label matrix. This simultaneously produces a
            %%% linear set of numbers (and masking of pixels outside the
            %%% object is done automatically, in a sense).
            Intensities = OrigImage(RetrievedCropMask == i);
            
            
            
            
            
%             %%% Diagnostic:
%             PerObjectImage = zeros(size(OrigImage));
%             PerObjectImage(RetrievedCropMask == i) = OrigImage(RetrievedCropMask == i);
%             %figure(31)
%             %imagesc(PerObjectImage)
% 
%             %%% Removes Rows and Columns that are completely blank.
%             ColumnTotals = sum(PerObjectImage,1);
%             warning off all
%             ColumnsToDelete = ~logical(ColumnTotals);
%             warning on all
%             drawnow
%             CroppedImage = PerObjectImage;
%             CroppedImage(:,ColumnsToDelete,:) = [];
%             CroppedImagePlusRange = CroppedImage;
%             [rows,columns] = size(CroppedImage);
%             CroppedImagePlusRange(:,end+1) = 1;
% 
%             [ManualThreshold,bw] = CPthresh_tool(CroppedImagePlusRange,'gray',1);
%             ManualThreshold = log(ManualThreshold);
%             %%% Initializes the variables.
%             if ~exist('ManualThresholds','var')
%                 ManualThresholds = [];
%             end
%             ManualThresholds(end+1) = ManualThreshold;
%             save('Batch_32Manualdata','ManualThresholds');

            

            
            
            
            
            
            %%% Sends those pixels to the appropriate threshold
            %%% subfunctions.
            eval(['CalculatedThreshold = ',ThresholdMethod,'(Intensities,handles,ImageName,pObject);']);
            %%% This evaluates to something like: Threshold =
            %%% Otsu(Intensities,handles,ImageName,pObject);

            

            %%% Sets the pixels corresponding to object i to equal the
            %%% calculated threshold.
            Threshold(RetrievedCropMask == i) = CalculatedThreshold;
%            figure(32), imagesc(Threshold), colormap('gray')
        end
    end

    if MethodFlag == 1 || MethodFlag == 2 %%% For the Adaptive and the PerObject methods.
        %%% Adjusts any of the threshold values that are significantly
        %%% lower or higher than the global threshold.  Thus, if there are
        %%% no objects within a block (e.g. if cells are very sparse), an
        %%% unreasonable threshold will be overridden.
        Threshold(Threshold <= 0.7*GlobalThreshold) = 0.7*GlobalThreshold;
        Threshold(Threshold >= 1.5*GlobalThreshold) = 1.5*GlobalThreshold;
    end

elseif strcmp(Threshold,'All')
    if SetBeingAnalyzed == handles.Current.StartingImageSet,
        try
            %%% Notifies the user that the first image set will take much
            %%% longer than subsequent sets. Obtains the screen size.
            [ScreenWidth,ScreenHeight] = CPscreensize;
            PotentialBottom = [0, (ScreenHeight-720)];
            BottomOfMsgBox = max(PotentialBottom);
            h = CPmsgbox('Preliminary calculations are under way for the Identify Primary Threshold module.  Subsequent image sets will be processed much more quickly than the first image set.');
            OrigSize = get(h, 'Position');
            PositionMsgBox = [500 BottomOfMsgBox OrigSize(3) OrigSize(4)];
            set(h, 'Position', PositionMsgBox)
            drawnow
            
            %%% Retrieves the path where the images are stored from the
            %%% handles structure.
            fieldname = ['Pathname', ImageName];
            try Pathname = handles.Pipeline.(fieldname);
            catch error(['Image processing was canceled in the ', ModuleName, ' module because it must be run using images straight from a load images module (i.e. the images cannot have been altered by other image processing modules). This is because you have asked the Identify Primary Threshold module to calculate a threshold based on all of the images before identifying objects within each individual image as CellProfiler cycles through them. One solution is to process the entire batch of images using the image analysis modules preceding this module and save the resulting images to the hard drive, then start a new stage of processing from this Identify Primary Threshold module onward.'])
            end
            %%% Retrieves the list of filenames where the images are stored
            %%% from the handles structure.
            fieldname = ['FileList', ImageName];
            FileList = handles.Pipeline.(fieldname);
            %%% Calculates the threshold based on all of the images.
            Counts = zeros(256,1);
            NumberOfBins = 256;
            for i=1:length(FileList)
                Image = CPimread(fullfile(Pathname,char(FileList(i))));
                Counts = Counts + imhist(im2uint8(Image(:)), NumberOfBins);
                drawnow
            end
            % Variables names are chosen to be similar to the formulas in
            % the Otsu paper.
            P = Counts / sum(Counts);
            Omega = cumsum(P);
            Mu = cumsum(P .* (1:NumberOfBins)');
            Mu_t = Mu(end);
            % Saves the warning state and disable warnings to prevent
            % divide-by-zero warnings.
            State = warning;
            warning off Matlab:DivideByZero
            SigmaBSquared = (Mu_t * Omega - Mu).^2 ./ (Omega .* (1 - Omega));
            % Restores the warning state.
            warning(State);
            % Finds the location of the maximum value of sigma_b_squared.
            % The maximum may extend over several bins, so average together
            % the locations.  If maxval is NaN, meaning that
            % sigma_b_squared is all NaN, then return 0.
            Maxval = max(SigmaBSquared);
            if isfinite(Maxval)
                Idx = mean(find(SigmaBSquared == Maxval));
                % Normalizes the threshold to the range [0, 1].
                Threshold = (Idx - 1) / (NumberOfBins - 1);
            else
                Threshold = 0.0;
            end
        catch [ErrorMessage, ErrorMessage2] = lasterr;
            error(['An error occurred in the ', ModuleName, ' module. Matlab says the problem is: ', ErrorMessage, ErrorMessage2])
        end
        fieldname = ['Threshold', ImageName];
        handles = CPaddimages(handles,fieldname,Threshold);
    else fieldname = ['Threshold', ImageName];
        Threshold = CPretrieveimage(handles,fieldname,ModuleName);
    end
elseif strcmp(Threshold,'Set interactively')
    fieldname = ['Threshold',ImageName];
    if SetBeingAnalyzed == handles.Current.StartingImageSet
        Threshold = CPthresh_tool(OrigImage(:,:,1));
        handles = CPaddimages(handles,fieldname,Threshold);
    else
        Threshold = CPretrieveimage(handles,fieldname,ModuleName);
    end
else
    %%% If the threshold is a number, it means that it was manually entered
    %%% by the user, or that we calculated it in the binary crop image
    %%% section above. Checks that the Threshold parameter has a valid
    %%% value
    if ischar(Threshold)
        Threshold = str2double(Threshold);
    end
    if isnan(Threshold) || Threshold > 1 || Threshold < 0 %#ok Ignore MLint
        error(['The threshold entered in the ', ModuleName, ' module is not a number, or is outside the acceptable range of 0 to 1.'])
    end
end
%%% Correct the threshold using the correction factor given by the user and
%%% make sure that the threshold is not larger than the minimum threshold
Threshold = ThresholdCorrection*Threshold;
Threshold = max(Threshold,MinimumThreshold);
Threshold = min(Threshold,MaximumThreshold);
if ~isempty(ObjectVar),
    handles = CPaddmeasurements(handles, 'Image', CPjoinstrings('Threshold','OrigThreshold', ObjectVar), mean(mean(Threshold)));
end

if (nargout >= 3),
    if ~ exist('BinaryCropMask', 'var')
        varargout(1) = {WeightedVariance(OrigImage, true(size(OrigImage)), Threshold)};
    else
        varargout(1) = {WeightedVariance(OrigImage, BinaryCropMask~=0, Threshold)};
    end
end
if (nargout >= 4),
    if ~ exist('BinaryCropMask', 'var')
        varargout(2) = {SumOfEntropies(OrigImage, true(size(OrigImage)), Threshold)};
    else
        varargout(2) = {SumOfEntropies(OrigImage, BinaryCropMask~=0, Threshold)};
    end
end

%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

function level = Otsu(im,handles,ImageName,pObject)
%%% This is the Otsu method of thresholding, adapted from MATLAB's
%%% graythresh function. Our modifications work in log space, and take into
%%% account the max and min values in the image.

%%% The following is needed for the adaptive cases where there the image
%%% has been cropped. This must be done within this subfunction, rather
%%% than in the main code prior to sending to this function via blkproc,
%%% because the blkproc function takes a single image as input, so we have
%%% to store the image and its cropmask in a single image variable.
if ndims(im) == 3
    Image = im(:,:,1);
    CropMask = im(:,:,2);
    clear im
    im = Image(CropMask==1);
else im = im(:);
end

if max(im) == min(im)
    level = im(1);
elseif isempty(im)
    %%% im will be empty if the entire image is cropped away by the
    %%% CropMask. I am not sure whether it is better to then set the level
    %%% to 0 or 1. Setting the level to empty causes problems downstream.
    %%% Presumably setting the level to 1 will not cause major problems
    %%% because the other blocks will average it out as we get closer to
    %%% real objects?
    level = 1;
else
    %%% We want to limit the dynamic range of the image to 256. Otherwise,
    %%% an image with almost all values near zero can give a bad result.
    minval = max(im)/256;
    im(im < minval) = minval;
    im = log(im);
    minval = min (im);
    maxval = max (im);
    im = (im - minval) / (maxval - minval);
    level = exp(minval + (maxval - minval) * graythresh(im));
end

% %%% For debugging:
% data = TrimmedImage;
% figure(30)
% subplot(1,2,1)
% hist(data(:),100);
% title(['trimmed data; Mean = ',num2str(Mean),'; StDev = ',num2str(StDev)])
% data = im;
% subplot(1,2,2)
% [Contents,BinLocations] = hist(data(:),100);
% hist(data(:),100);
% title(['Thresh = ',num2str(level),'; log data'])
% hold on
% plot([level;level],[0,max(Contents)])
% hold off
% figure(30)

function level = MoG(im,handles,ImageName,pObject)
%%% Stands for Mixture of Gaussians. This function finds a suitable
%%% threshold for the input image Block. It assumes that the pixels in the
%%% image belong to either a background class or an object class. 'pObject'
%%% is an initial guess of the prior probability of an object pixel, or
%%% equivalently, the fraction of the image that is covered by objects.
%%% Essentially, there are two steps. First, a number of Gaussian
%%% distributions are estimated to match the distribution of pixel
%%% intensities in OrigImage. Currently 3 Gaussian distributions are
%%% fitted, one corresponding to a background class, one corresponding to
%%% an object class, and one distribution for an intermediate class. The
%%% distributions are fitted using the Expectation-Maximization (EM)
%%% algorithm, a procedure referred to as Mixture of Gaussians modeling.
%%% When the 3 Gaussian distributions have been fitted, it's decided
%%% whether the intermediate class models background pixels or object
%%% pixels based on the probability of an object pixel 'pObject' given by
%%% the user.

%%% The following is needed for the adaptive cases where there the image
%%% has been cropped. This must be done within this subfunction, rather
%%% than in the main code prior to sending to this function via blkproc,
%%% because the blkproc function takes a single image as input, so we have
%%% to store the image and its cropmask in a single image variable.
if ndims(im) == 3
    Image = im(:,:,1);
    CropMask = im(:,:,2);
    clear im
    im = Image(CropMask==1);
else im = im(:);
end

if max(im) == min(im)
    level = im(1);
elseif isempty(im)
    %%% im will be empty if the entire image is cropped away by the
    %%% CropMask. I am not sure whether it is better to then set the level
    %%% to 0 or 1. Setting the level to empty causes problems downstream.
    %%% Presumably setting the level to 1 will not cause major problems
    %%% because the other blocks will average it out as we get closer to
    %%% real objects?
    level = 1;
else

    %%% The number of classes is set to 3
    NumberOfClasses = 3;

    %%% If the image is larger than 512x512, select a subset of 512^2
    %%% pixels for speed. This should be enough to capture the statistics
    %%% in the image.
    % im = im(:);
    if length(im) > 512^2
        is2008b_or_greater = ~CPverLessThan('matlab','7.7');
        if is2008b_or_greater,
            defaultStream = RandStream.getDefaultStream;
            savedState = defaultStream.State;
            RandStream.setDefaultStream(RandStream('mt19937ar','seed',0));
        else
            rand('seed',0);
        end
        indexes = randperm(length(im));
        if is2008b_or_greater, defaultStream.State = savedState; end
        im = im(indexes(1:512^2));
    end

    %%% Convert user-specified percentage of image covered by objects to a
    %%% prior probability of a pixel being part of an object.
    
    %%% Since the default list for MoG thresholding contains a % sign, we
    %%% need to remove the percent sign and use only the number to
    %%% calculate the threshold. If the pObject does not contain a % sign,
    %%% it will continue.  
    %%% pObject is important, but pObjectNew is only used in  2 lines of code.blah
    if ischar(pObject)
        if  regexp(pObject, '%')
            pObjectNew = regexprep(pObject, '%', '');
            pObject = (str2double(pObjectNew)/100);

        else
            pObject = str2double(pObject);
        end
    end
    %pObject = str2double(pObject(1:2))/100; old code--need to remove
    %%% Get the probability for a background pixel
    pBackground = 1 - pObject;

    %%% Initialize mean and standard deviations of the three Gaussian
    %%% distributions by looking at the pixel intensities in the original
    %%% image and by considering the percentage of the image that is
    %%% covered by object pixels. Class 1 is the background class and Class
    %%% 3 is the object class. Class 2 is an intermediate class and we will
    %%% decide later if it encodes background or object pixels. Also, for
    %%% robustness the we remove 1% of the smallest and highest intensities
    %%% in case there are any quantization effects that have resulted in
    %%% unnaturally many 0:s or 1:s in the image.
    im = sort(im);
    im = im(ceil(length(im)*0.01):round(length(im)*0.99));
    ClassMean(1) = im(round(length(im)*pBackground/2));                      %%% Initialize background class
    ClassMean(3) = im(round(length(im)*(1 - pObject/2)));                    %%% Initialize object class
    ClassMean(2) = (ClassMean(1) + ClassMean(3))/2;                                            %%% Initialize intermediate class
    %%% Initialize standard deviations of the Gaussians. They should be the
    %%% same to avoid problems.
    ClassStd(1:3) = 0.15;
    %%% Initialize prior probabilities of a pixel belonging to each class.
    %%% The intermediate class is gets some probability from the background
    %%% and object classes.
    pClass(1) = 3/4*pBackground;
    pClass(2) = 1/4*pBackground + 1/4*pObject;
    pClass(3) = 3/4*pObject;

    %%% Apply transformation.  a < x < b, transform to log((x-a)/(b-x)).
    %a = - 0.000001; b = 1.000001; im = log((im-a)./(b-im)); ClassMean =
    %log((ClassMean-a)./(b - ClassMean)) ClassStd(1:3) = [1 1 1];

    %%% Expectation-Maximization algorithm for fitting the three Gaussian
    %%% distributions/classes to the data. Note, the code below is general
    %%% and works for any number of classes. Iterate until parameters don't
    %%% change anymore.
    delta = 1;
    while delta > 0.001
        %%% Store old parameter values to monitor change
        oldClassMean = ClassMean;

        %%% Update probabilities of a pixel belonging to the background or
        %%% object1 or object2
        for k = 1:NumberOfClasses
            pPixelClass(:,k) = pClass(k)* 1/sqrt(2*pi*ClassStd(k)^2) * exp(-(im - ClassMean(k)).^2/(2*ClassStd(k)^2));
        end
        pPixelClass = pPixelClass ./ repmat(sum(pPixelClass,2) + eps,[1 NumberOfClasses]);

        %%% Update parameters in Gaussian distributions
        for k = 1:NumberOfClasses
            pClass(k) = mean(pPixelClass(:,k));
            ClassMean(k) = sum(pPixelClass(:,k).*im)/(length(im)*pClass(k));
            ClassStd(k)  = sqrt(sum(pPixelClass(:,k).*(im - ClassMean(k)).^2)/(length(im)*pClass(k))) + sqrt(eps);    % Add sqrt(eps) to avoid division by zero
        end

        %%% Calculate change
        delta = sum(abs(ClassMean - oldClassMean));
    end

    %%% Now the Gaussian distributions are fitted and we can describe the
    %%% histogram of the pixel intensities as the sum of these Gaussian
    %%% distributions. To find a threshold we first have to decide if the
    %%% intermediate class 2 encodes background or object pixels. This is
    %%% done by choosing the combination of class probabilities 'pClass'
    %%% that best matches the user input 'pObject'.
    level = linspace(ClassMean(1),ClassMean(3),10000);
    Class1Gaussian = pClass(1) * 1/sqrt(2*pi*ClassStd(1)^2) * exp(-(level - ClassMean(1)).^2/(2*ClassStd(1)^2));
    Class2Gaussian = pClass(2) * 1/sqrt(2*pi*ClassStd(2)^2) * exp(-(level - ClassMean(2)).^2/(2*ClassStd(2)^2));
    Class3Gaussian = pClass(3) * 1/sqrt(2*pi*ClassStd(3)^2) * exp(-(level - ClassMean(3)).^2/(2*ClassStd(3)^2));
    if abs(pClass(2) + pClass(3) - pObject) < abs(pClass(3) - pObject)
        %%% Intermediate class 2 encodes object pixels
        BackgroundDistribution = Class1Gaussian;
        ObjectDistribution = Class2Gaussian + Class3Gaussian;
    else
        %%% Intermediate class 2 encodes background pixels
        BackgroundDistribution = Class1Gaussian + Class2Gaussian;
        ObjectDistribution = Class3Gaussian;
    end

    %%% Now, find the threshold at the intersection of the background
    %%% distribution and the object distribution.
    [ignore,index] = min(abs(BackgroundDistribution - ObjectDistribution)); %#ok Ignore MLint
    level = level(index);
end

function level = Background(im,handles,ImageName,pObject)
%%% The threshold is calculated by calculating the mode and multiplying by
%%% 2 (an arbitrary empirical factor). The user will presumably adjust the
%%% multiplication factor as needed.

%%% The following is needed for the adaptive cases where there the image
%%% has been cropped. This must be done within this subfunction, rather
%%% than in the main code prior to sending to this function via blkproc,
%%% because the blkproc function takes a single image as input, so we have
%%% to store the image and its cropmask in a single image variable.
if ndims(im) == 3
    Image = im(:,:,1);
    CropMask = im(:,:,2);
    clear im
    im = Image(CropMask==1);
else im = im(:);
end

if max(im) == min(im)
    level = im(1);
elseif isempty(im)
    %%% im will be empty if the entire image is cropped away by the
    %%% CropMask. I am not sure whether it is better to then set the level
    %%% to 0 or 1. Setting the level to empty causes problems downstream.
    %%% Presumably setting the level to 1 will not cause major problems
    %%% because the other blocks will average it out as we get closer to
    %%% real objects?
    level = 1;
else
    %% We were using the 'mode' function here, but it does not always report
    %% what is the obvious peak in the histogram.  This is because the
    %% distribution is not continuous, and may be binned in not-so-obvious
    %% ways.  So we are using imhist instead to bin intensities somewhat 
    %% and then do an effective mode calculation.
    
    %% Also handle the case in which there are enough saturated, or zeroed,
    %% pixels that the mode is 0 or 1 (or some other, pinned high value).  
    %% We remove the high and low values from the mode calculation.  
    %% Robust cutoff arbitrarily set to 2%.

    [counts,x] = imhist(im);
    
    %% Remove bins with zero counts at the low and high ends of imhist
    mn = find(counts ~=0,1, 'first');
    mx = find(counts ~=0,1, 'last');
    counts_scaled = counts(mn:mx);
    x_scaled = x(mn:mx);
    
    thresh = 0.02;
    robust_indices = ceil(thresh * length(x_scaled)):ceil((1-thresh) * length(x_scaled));
    robust_counts = counts_scaled(robust_indices);
    robust_x = x_scaled(robust_indices);
    [y,i] = max(robust_counts);
    level = 2.*x((mn-1)+(robust_indices(1)-1)+i);
end


function level = RobustBackground(im,handles,ImageName,pObject)
%%% The threshold is calculated by trimming the top and bottom 5% of
%%% pixels off the image, then calculating the mean and standard deviation
%%% of the remaining image. The threshold is then set at 2 (empirical
%%% value) standard deviations above the mean. 

%%% The following is needed for the adaptive cases where there the image
%%% has been cropped. This must be done within this subfunction, rather
%%% than in the main code prior to sending to this function via blkproc,
%%% because the blkproc function takes a single image as input, so we have
%%% to store the image and its cropmask in a single image variable.
warning off MATLAB:divideByZero
if ndims(im) == 3
    Image = im(:,:,1);
    CropMask = im(:,:,2);
    clear im
    im = Image(CropMask==1);
else im = im(:);
end

if max(im) == min(im)
    level = im(1);
elseif isempty(im)
    %%% im will be empty if the entire image is cropped away by the
    %%% CropMask. I am not sure whether it is better to then set the level
    %%% to 0 or 1. Setting the level to empty causes problems downstream.
    %%% Presumably setting the level to 1 will not cause major problems
    %%% because the other blocks will average it out as we get closer to
    %%% real objects?
    level = 1;
else
    %%% First, the image's pixels are sorted from low to high.
    im = sort(im);
    %%% The index of the 5th percentile is calculated, with a minimum of 1.
    LowIndex = max(1,round(.05*length(im)));
    %%% The index of the 95th percentile is calculated, with a maximum of the
    %%% number of pixels in the whole image.
    HighIndex = min(length(im),round(.95*length(im)));
    TrimmedImage = im(LowIndex: HighIndex);
    Mean = mean(TrimmedImage);
    StDev = std(TrimmedImage);
    level = Mean + 2*StDev;
end

% %%% DEBUGGING
% Logim = log(sort(im(im~=0)));
% 
% %%% For debugging:
% figure(30)
% subplot(1,2,1)
% hist(Logim,100);
% title(['Log data; Mean = ',num2str(Mean),'; StDev = ',num2str(StDev)])
% pause(0.1)

% %%% For debugging:
% data = TrimmedImage;
% figure(30)
% subplot(1,2,1)
% hist(data(:),100);
% title(['trimmed data; Mean = ',num2str(Mean),'; StDev = ',num2str(StDev)])
% data = im;
% subplot(1,2,2)
% [Contents,BinLocations] = hist(data(:),100);
% hist(data(:),100);
% title(['Thresh = ',num2str(level),'; raw data'])
% hold on
% plot([level;level],[0,max(Contents)])
% hold off
% 
% figure(30)
% 
% %%% More debugging:
% try
%     load('Batch_80Autodata');
% end
% %%% Initializes the variables.
% if ~exist('Means','var')
%    Means = []; 
%    StDevs = [];
%    Levels = [];
%    TrimmedImages = [];
%    Images = [];
% end
% Means(end+1) = Mean;
% StDevs(end+1) = StDev;
% Levels(end+1) = level;
% TrimmedImages{end+1} = {TrimmedImage};
% Images{end+1} = {im};
% save('Batch_80Autodata','Means','StDevs','Levels','TrimmedImages','Images');
warning on MATLAB:divideByZero

function level = RidlerCalvard(im,handles,ImageName,pObject)
warning off MATLAB:divideByZero
%%% The following is needed for the adaptive cases where there the image
%%% has been cropped. This must be done within this subfunction, rather
%%% than in the main code prior to sending to this function via blkproc,
%%% because the blkproc function takes a single image as input, so we have
%%% to store the image and its cropmask in a single image variable.
if ndims(im) == 3
    Image = im(:,:,1);
    CropMask = im(:,:,2);
    clear im
    im = Image(CropMask==1);
else im = im(:);
end

if max(im) == min(im)
    level = im(1);
elseif isempty(im)
    %%% im will be empty if the entire image is cropped away by the
    %%% CropMask. I am not sure whether it is better to then set the level
    %%% to 0 or 1. Setting the level to empty causes problems downstream.
    %%% Presumably setting the level to 1 will not cause major problems
    %%% because the other blocks will average it out as we get closer to
    %%% real objects?
    level = 1;
else
    %%% We want to limit the dynamic range of the image to 256. Otherwise,
    %%% an image with almost all values near zero can give a bad result.
    MinVal = max(im)/256;
    im(im<MinVal) = MinVal;
    im = log(im);
    MinVal = min(im);
    MaxVal = max(im);
    im = (im - MinVal)/(MaxVal - MinVal);
    PreThresh = 0;
    %%% This method needs an initial value to start iterating. Using
    %%% graythresh (Otsu's method) is probably not the best, because the
    %%% Ridler Calvard threshold ends up being too close to this one and in
    %%% most cases has the same exact value.
    NewThresh = graythresh(im);
    delta = 0.00001;
    while abs(PreThresh - NewThresh)>delta
        PreThresh = NewThresh;
        Mean1 = mean(im(im<PreThresh));
        Mean2 = mean(im(im>=PreThresh));
        NewThresh = mean([Mean1,Mean2]);
    end
    level = exp(MinVal + (MaxVal-MinVal)*NewThresh);
end
warning on MATLAB:divideByZero

function level = Kapur(im,handles,ImageName,pObject)
%%% This is the Kapur, Sahoo, & Wong method of thresholding, adapted to log-space.

%%% The following is needed for the adaptive cases where there the image
%%% has been cropped. This must be done within this subfunction, rather
%%% than in the main code prior to sending to this function via blkproc,
%%% because the blkproc function takes a single image as input, so we have
%%% to store the image and its cropmask in a single image variable.
if ndims(im) == 3
    Image = im(:,:,1);
    CropMask = im(:,:,2);
    clear im
    im = Image(CropMask==1);
else im = im(:);
end

if max(im) == min(im)
    level = im(1);
elseif isempty(im)
    %%% im will be empty if the entire image is cropped away by the
    %%% CropMask. I am not sure whether it is better to then set the level
    %%% to 0 or 1. Setting the level to empty causes problems downstream.
    %%% Presumably setting the level to 1 will not cause major problems
    %%% because the other blocks will average it out as we get closer to
    %%% real objects?
    level = 1;
else
    level = Threshold_Kapur(im, 8);
end


%%% This function computes the threshold of an image by
%%% log-transforming its values, then searching for the threshold that
%%% maximizes the sum of entropies of the foreground and background
%%% pixel values, when treated as separate distributions.
function thresh = Threshold_Kapur(Image, bits)
% Find the smoothed log histogram.
[N, X] = hist(log2(smooth_log_histogram(Image(:), bits)), 256);

% drop any zero bins
drop = (N == 0);
N(drop) = [];
X(drop) = [];

% check for corner cases
if length(X) == 1,
    thresh = X(1);
    return;
end

% Normalize to probabilities
P = N / sum(N);

% Find the probabilities totals up to and above each possible threshold.
loSum = cumsum(P);
hiSum = loSum(end) - loSum;
loE = cumsum(P .* log2(P));
hiE = loE(end) - loE;

% compute the entropies
s = warning('off', 'MATLAB:divideByZero');
loEntropy = loE ./ loSum - log2(loSum);
hiEntropy = hiE ./ hiSum - log2(hiSum);
warning(s);

sumEntropy = loEntropy(1:end-1) + hiEntropy(1:end-1);
sumEntropy(~ isfinite(sumEntropy)) = Inf;
entry = min(find(sumEntropy == min(sumEntropy)));
thresh = 2^((X(entry) + X(entry+1)) / 2);


%%% This function smooths a log-transformed histogram, using noise
%%% proportional to the histogram value.
function Q = smooth_log_histogram(R, bits)
%%% seed random state
is2008b_or_greater = ~CPverLessThan('matlab','7.7');
if is2008b_or_greater,
    defaultStream = RandStream.getDefaultStream;
    savedState = defaultStream.State;
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
else
    rand('seed',0);
end
R(R == 0) = 1 / (2^bits);
Q = exp(log(R) + 0.5*randn(size(R)).*(-log2(R)/bits));
Q(Q > 1) = 1.0;
Q(Q < 0) = 0.0;
if is2008b_or_greater, defaultStream.State = savedState; end

%%% Weighted variances of the foreground and background.
function  wv = WeightedVariance(Image, CropMask, Threshold)
if isempty(Image(CropMask)),
    wv = 0;
    return;
end

%%% clamp dynamic range
minval = max(Image(CropMask))/256;
if minval == 0.0,
    wv = 0;
    return;
end
Image(Image < minval) = minval;

%%% Compute the weighted variance
FG = log2(Image((Image >= Threshold) & CropMask));
BG = log2(Image((Image < Threshold) & CropMask));
if isempty(FG),
    wv = var(BG);
elseif isempty(BG);
    wv = var(FG);
else
    wv = (length(FG) * var(FG) + length(BG) * var(BG)) / (length(FG) + length(BG));
end



%%% Sum of entropies of foreground and background as separate distributions.
function  soe = SumOfEntropies(Image, CropMask, Threshold)
if isempty(Image(CropMask)),
    soe = 0;
    return;
end

%%% clamp dynamic range
minval = max(Image(CropMask))/256; 
if minval == 0.0,
    soe = 0;
    return;
end
Image(Image < minval) = minval;

%%% Smooth the histogram
Image = smooth_log_histogram(Image, 8);

ImMin = min(Image(CropMask));
ImMax = max(Image(CropMask));

%%% Find bin locations
upper = log2(ImMax);
lower = log2(ImMin);
step = (upper - lower) / 256;
X2 = (0:256) * step + lower;
% necessary for histc
X2 = X2+step/2;

%%% Find counts for FG and BG
FG = Image((Image >= Threshold) & CropMask);
BG = Image((Image < Threshold) & CropMask);
NFG = histc(log2(FG), X2);
NBG = histc(log2(BG), X2);

%%% drop empty bins
NFG = NFG(NFG > 0);
NBG = NBG(NBG > 0);

if isempty(NFG)
    NFG = [1];
end

if isempty(NBG)
    NBG = [1];
end

% normalize
NFG = NFG / sum(NFG);
NBG = NBG / sum(NBG);

%%% compute sum of entropies
soe = sum(NFG .* log2(NFG)) + sum(NBG .* log2(NBG));






function ispresent = CPisimageinpipeline(handles, fieldname)
% Check if images exist in handles.Pipeline structure.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 7230 $

if ~isfield(handles.Pipeline,'ImageGroupFields')
    ispresent = isfield(handles.Pipeline,fieldname);
else
    ispresent = isfield(handles.Pipeline.GroupFileList{handles.Pipeline.CurrentImageGroupID},fieldname);
end






function result = CPverLessThan(toolboxstr, verstr)

% This CP function is present only so we can easily replace the
% verlessThan if necessary.  See documentation for warndlg for usage.
    
error(nargchk(2, 2, nargin, 'struct'))
        
if ~ischar(toolboxstr) || ~ischar(verstr)
    error('MATLAB:verLessThan:invalidInput', 'Inputs must be strings.')
end

if ~isdeployed,
    toolboxver = ver(toolboxstr);
    if isempty(toolboxver)
        error('MATLAB:verLessThan:missingToolbox', 'Toolbox ''%s'' not found.', toolboxstr)
    end

    toolboxParts = getParts(toolboxver(1).Version);
else
    toolboxver = version;
    if isempty(toolboxver)
        error('MATLAB:verLessThan:missingToolbox', 'Toolbox ''%s'' not found.', toolboxstr)
    end
    toolboxParts = getParts(toolboxver);
end
verParts = getParts(verstr);

result = (sign(toolboxParts - verParts) * [1; .1; .01]) < 0;

function parts = getParts(V)
    parts = sscanf(V, '%d.%d.%d')';
    if length(parts) < 3
       parts(3) = 0; % zero-fills to 3 elements
    end
