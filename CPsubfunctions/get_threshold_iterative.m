function SelThresh=get_threshold_iterative(im,thresize)

%convert the image to double
Image = double(im(:));
minImage = (min(Image));
maxImage = (max(Image));

%get the threshold values
Threshold = linspace(minImage+1, maxImage, thresize)';

%loop through thrsholds
parfor ix = 1:length(Threshold)
    tmpThreshold = Threshold(ix,1);
    tmpImage = Image > tmpThreshold;
    tmpCC = bwconncomp(tmpImage, 8);
    numObjects(ix,1) = tmpCC.NumObjects;
end
[ymax,imax,ymin,imin] = extrema(numObjects);

%get the threshold values
if min(Threshold(imax))> min(Threshold(imin))
    inxMin2Use = (Threshold(imin)> min(Threshold(imax)));
    minThresholds = Threshold(imin);
    SelThresh = min(minThresholds(inxMin2Use));
else
    SelThresh = min(Threshold(imin));
end

