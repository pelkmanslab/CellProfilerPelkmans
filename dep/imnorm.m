function imNormalized = imnorm(imInput,lowerQuantile,upperQuantile)

lowerQuantileValue = quantile(imInput(:),lowerQuantile);
upperQuantileValue = quantile(imInput(:),upperQuantile);
if lowerQuantileValue ~= upperQuantileValue
    tempim = imInput;
    tempim(tempim<lowerQuantileValue) = lowerQuantileValue;
    tempim(tempim>upperQuantileValue) = upperQuantileValue;
    tempim = tempim - lowerQuantileValue;
    tempim = double(tempim);
    imNormalized = tempim ./ max(tempim(:));
else
    imNormalized = zeros(size(imInput));
end

end