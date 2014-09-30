function varargout = createOccupancyImage(SegmentationCC,varargin)
%occImg = CREATEOCCUPANCYIMAGE(SEGMENTATIONCC) projects SEGMENTATIONCC to
% an OCCIMG, where the intensity is directly correlated to the amount of
% pixels occupied by objects
%
%   Note that XY will have X along dim(2) and Y along dim(1) (following
%   Matlabs scheme of using rows along dim(1) and columns along dim(2);
%
% Default will be a projection along the 'XY' plane. However optionally,
% projections along other planes can be made. Output order will reflect
% order of input, for instance:
%   occImgYZ = CREATEOCCUPANCYIMAGE(SEGMENTATIONCC,'YZ')    or
%   [occImgYZ occImgXY occImgXZ]= CREATEOCCUPANCYIMAGE(SEGMENTATIONCC,'YZ','XY','XZ')
%
%   ----------------------------------
%   Information about module:
%   Created for CP3D to create image, which indicates how many
%   voxels along one projection direction have an object
%   x Purpose: quick preview
%   x Additional Purpose: for 2.5D analysis, where individul z-objects can be
%   projected to 2D, this function creates a 2D image about the spatial
%   organization, which might yield useful features, when combined with
%   classical modules, such as MeasureTexture or MeasureIntensity. This way
%   one would obtain more 3D features without investing more time into full
%   3D modules for fine 3D morphology
%
%   Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INITIALIZE SETTINGS   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check input
if nargout ~= max(nargin-1,1);
    error('Did not specify correct number of outputs for selected amount of Occupancy Projections');
end

%initialize
numPixels = cellfun(@numel,SegmentationCC.PixelIdxList);
countPixels = sum(numPixels,2);
RCZO=NaN(countPixels,4);

%loop through each object and retreive positional information
%(RowColumnZplaneObjectid)
numLastIX=0;
for k=1:SegmentationCC.NumObjects
    
    numFirstIX=numLastIX+1;
    numLastIX=numLastIX+numPixels(1,k);
    
    [RCZO(numFirstIX:numLastIX,1) RCZO(numFirstIX:numLastIX,2) RCZO(numFirstIX:numLastIX,3)] = ind2sub(SegmentationCC.ImageSize,SegmentationCC.PixelIdxList{1,k});
    RCZO(numFirstIX:numLastIX,4) = k;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MEASURE   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin == 1      % if only input Segmentation specified, project along XY
    varargout{1}=OccProj(RCZO,SegmentationCC,1,2);
elseif nargin<=4
    for k=1:(nargin-1)
        switch varargin{k}
            case 'XY'
                varargout{k} = OccProj(RCZO,SegmentationCC,1,2);
            case 'XZ'
                varargout{k} = OccProj(RCZO,SegmentationCC,2,3);
            case 'YZ'
                varargout{k} = OccProj(RCZO,SegmentationCC,1,3);
            otherwise
                error('Input not defined current options. See help Choose ''XY'' or ''XZ'' or ''YZ''.' );
        end
    end
else
    error('Too many input parameters specified.' )
end


end


function occImg = OccProj(RCZO,SegCC,a,b)

occImg = uint16(zeros(SegCC.ImageSize(a),SegCC.ImageSize(b)));
SubRC = sub2ind2(size(occImg), RCZO(:,[a b]));
unSubRC = unique(SubRC);
CountPerPx = histc(SubRC,unSubRC);
occImg(unSubRC)=CountPerPx;

end