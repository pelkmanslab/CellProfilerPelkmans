function retImage =getUniqueCompactImageLABEL(centreMask, labelMask)

centreMask = -1*double(double(centreMask)>0.5);
labelMask2 = labelMask.*centreMask;


fI3 = getconnectedPoints(labelMask2);
fI3 = double(fI3<0.8).* labelMask;


usedLabel = unique(fI3);
maxLabel = max(max(fI3))+1;
relabelled = zeros(1, maxLabel);
len = length(usedLabel);

for i=1:len
    currentColor = usedLabel(i)+1;
    relabelled(currentColor) = i-1; 
end
[r c] = size(fI3);
for i= 1:r
    for j=1:c
        fI3(i,j)  = relabelled( fI3(i,j)+1 );
    end
end
retImage = fI3;

maxret=max(max(retImage))
maxunique= length(unique(retImage))-1
