function [ViewFromZ, ViewFromX, ViewFromY] = createOccupancyImage2(SegmentationCC,strOptionalSingleViewDirection)
% Will create projection of voxels occupied with objects in SegmentationCC.
% Note that this function is intended to replace the legacy function
% createOccupancyImage(), which is slow.
%
% Optional input strOptionalSingleViewDirection can be either x, y or z and
% used to speed up calculation in case that only individual views are
% required
%
% [TS] Support function for CP3D

if nargin < 2  % by default calculate all views
    shallCompute.x = true;
    shallCompute.y = true;
    shallCompute.z = true;
else % optionally compute only single views
    shallCompute.x = false;
    shallCompute.y = false;
    shallCompute.z = false;
    
    d = lower(strOptionalSingleViewDirection);
    if ismember({d},{'x','y','z'})
        shallCompute.(d) = true;
    else
        error('strOptionalSingleViewDirection must either be x, y or z!');
    end
    
end

% Initialize
ViewFromX = [];
ViewFromY = [];
ViewFromZ = [];

numOjects = SegmentationCC.NumObjects;
reconstructedVoxelsWithObject = false(SegmentationCC.ImageSize);

if numOjects >= 1
    for j=1:numOjects
        ix = SegmentationCC.PixelIdxList{j};
        reconstructedVoxelsWithObject(ix) = true;
    end
end

if shallCompute.x == true
    ViewFromX = sum(reconstructedVoxelsWithObject,2);
    ViewFromX = permute(ViewFromX,[3 1 2]);
end

if shallCompute.y == true
    ViewFromY = sum(reconstructedVoxelsWithObject,1);
    ViewFromY = permute(ViewFromY,[3 2 1]);
end

if shallCompute.z == true
    ViewFromZ = sum(reconstructedVoxelsWithObject,3);
end

end
