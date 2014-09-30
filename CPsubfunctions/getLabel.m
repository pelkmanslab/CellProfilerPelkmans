function [retImage maxret] = getLabel( labelMask)

% FUNCTION TO CONVERT SPARSE LABEL INTO COMPACT lABEL
% IMPLEMENTED by S. Chowdhury

labelMask = labelMask +1;

usedLabel =  unique(labelMask) ;

labelFlag(usedLabel) = 1;

labelFlag = cumsum(labelFlag);

if usedLabel(1)== 1
    labelFlag  = labelFlag - 1;
end

retImage = labelFlag(labelMask);

maxret=max(max(retImage));
