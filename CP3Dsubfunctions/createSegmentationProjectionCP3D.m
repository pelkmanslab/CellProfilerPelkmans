function SegementationMatrix = createSegmentationProjectionCP3D(SegmentationCC,strMethod,varargin)
% M = createSegmentationProjection(SegmentationCC,STRMETHOD, varargin) will
% create a 2D LabelMatrix, which has also been used in CP1 to identify
% objects. 
% 
% SEGMENTATIONCC has to be a CC Label derived from bwconncomp; 
%
% M = createSegmentationProjection(SegmentationCC,'Merge') will project all
% objects to XY and then create objects among Pixels connected in the
% projection. As an optional input the number of planes, which have to be
% occopied in SEGMENTATIONCC can be specified, eg.: M =
% createSegmentationProjection(SegmentationCC,'Merge',4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OBTAIN 2D LABEL   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check input

switch strMethod
    case 'Merge'
        switch nargin
            case 2
                minOccupiedLayers = 0;
            case 3
                minOccupiedLayers = varargin{1}-1;
            otherwise
                error('Incorrect number of input arguments for Method Merge');
        end
        bwObj = createOccupancyImage(SegmentationCC,'XY') > minOccupiedLayers;
        SegementationMatrix = bwlabel(bwObj);
        
    otherwise
        error('Selected method for creating Segmentation not supported');  
end


end

