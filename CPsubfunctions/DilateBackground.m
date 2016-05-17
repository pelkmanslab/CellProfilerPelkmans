function imLabelMatrix = DilateBackground(imLabelMask,strOutlineOrMask)

% imLabelMatrix = DilateBackground(imLabelMask,strOutlineOrMask)
% 
%
% Performs background dilation rather than object erosion to nicely
% separate neighboring objects from each other.
%
% Processing is done for each object on small sub-images using bounding
% boxes. This speeds up computational time.
% 
%
% Input variables:
%
% imLabelMask: Label image created with bwlabel, for example.
%
% strOutlineOrMask (optional): String specifying whether outline or mask
% image should be created. 'Outline' creates an outline image. 'Mask'
% creates a mask image (defaults to 'Mask').
%
%
% Authors:
%   Thomas Stoeger
%   Markus Herrmann



if nargin == 1
    strOutlineOrMask = 'Mask';
end

FinalLabelMatrixImage = double(imLabelMask);

% Obtain pixels at inner periphery of objects (note that this will reduce the computing time by more than 50% since much less membrane pixels have to be considered)
props = regionprops(FinalLabelMatrixImage,'BoundingBox');
BoxPerObj = cat(1,props.BoundingBox);

if all(size(BoxPerObj)==0) % this can sometimes happen
    imLabelMatrix  = zeros(size(FinalLabelMatrixImage));
    return
end

% Calculate allowed coordinates per object (to reduce computational cost)
%-> recalculating objects is necesary because of a background bug!
LisOfObjects = unique(FinalLabelMatrixImage(:));
if LisOfObjects(1) == 0 % remove background, if present
    LisOfObjects = LisOfObjects(2:end);
end

% Get outer coordinates of bounding box of each object (note that
% bounding boxes will speed up morphological image operations)
distanceToObjectMax = 3;
N = floor(BoxPerObj(:,2)-distanceToObjectMax-1);                    f = N < 1;                                N(f) = 1;
S = ceil(BoxPerObj(:,2)+BoxPerObj(:,4)+distanceToObjectMax+1);   	f = S > size(FinalLabelMatrixImage,1);    S(f) = size(FinalLabelMatrixImage,1);
W = floor(BoxPerObj(:,1)-distanceToObjectMax-1);                    f = W < 1;                                W(f) = 1;
E = ceil(BoxPerObj(:,1)+BoxPerObj(:,3)+distanceToObjectMax+1);      f = E > size(FinalLabelMatrixImage,2);    E(f) = size(FinalLabelMatrixImage,2);

% Create empty output
imLabelMatrix  = zeros(size(FinalLabelMatrixImage));
if ~isempty(LisOfObjects)  % if objects present
    clear k;
    for k=LisOfObjects'  % loop through individual objects to safe computation
        
        % Create small image to speed up morphological operations
        miniImage = FinalLabelMatrixImage(N(k):S(k),W(k):E(k));
        bwminiImage = miniImage==k;
        
        % Get pixels within object which are next to background. note
        % that the alternative strategy to shrink the object could lead
        % to small parts of the membrane becoming disconnected, if
        % there is a finger like extension of cells. When following the
        % reverse strategy to extend the background there can be a few
        % extra pixels, however since they are rare they will be
        % neglectible for the speed of the later calculation and it is
        % favourable to have closed cell membranes
        padbw = padarray(bwminiImage,[1 1]); % pad with 0s to ensure that background on all sides
        ipadbw = ~padbw; % get background
        tpadbw = bwmorph(ipadbw,'dilate',1); % extend background
        if strcmp(strOutlineOrMask,'Mask')
            dpadbw = ~tpadbw; % get inner pixels of objects
            dpadbwR = dpadbw(2:(end-1),2:(end-1)); % reverse padding to get correct coordinates
        elseif strcmp(strOutlineOrMask,'Outline')
            dpadbw = tpadbw & ~ipadbw; % get inner pixels of objects
            dpadbwR = dpadbw(2:(end-1),2:(end-1)); % reverse padding to get correct coordinates
            dpadbwR = bwmorph(dpadbwR,'thin'); % thin so that across diagonals around 1/3 of pixels is lost. will increase speed spot parameter calculation.
        end
        % Now map back the linear indices
        [r c] = find(dpadbwR);
        
        % Get indices for final image (note that mini image might have
        % permitted regions of other cells and thus boxes can not be
        % directly overlaid).
        r = r-1+N(k);
        c = c-1+W(k);
        w = sub2ind(size(imLabelMatrix),r,c);
        
        % Update Working copy of Final Segmentation image based on
        % linear indices. (Linear indices prevent depenceny of other
        % objects in same box. Thus overlapping boxes
        % can be processed sequentially)
        imLabelMatrix(w) = k;
        
    end
end

end
