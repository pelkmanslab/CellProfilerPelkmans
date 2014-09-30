function f = spatial_energy_distribution(CellInt,fraction)
% 14.3.2008 (C) Pekka Ruusuvuori

E = CellInt.^2;
Etot = sum(E(:));
Area = sum(sum(CellInt>0));
[sortedE,sI] = sort(E(:),'descend');
runningSum = 0;
ind = 1;
while runningSum < fraction*Etot
    runningSum = runningSum + sortedE(ind);
    ind = ind + 1;
end
BI = zeros(size(CellInt));
BI(sI(1:ind)) = 1;
[LBI,num] = bwlabel(BI);

% ----- features:
% area of fraction of energy per whole area
EArea = sum(BI(:))/Area;
% searchfor major area
sizevec = zeros(num,1);
for indi = 1:num
    sizevec(indi) = sum(sum(LBI==indi));
end
% biggest area inside fraction per whole area
[MaxSize,maind] = max(sizevec);
MajorEArea = MaxSize/Area;
% distances to center
s = regionprops(bwlabel(LBI==maind),'Centroid');
[XX,YY] = find(BI);
XXn = XX-s.Centroid(2);
YYn = YY-s.Centroid(1);
distances = sqrt(XXn.^2 + YYn.^2);
EnergyDispersion = std(distances);
f = [EArea,MajorEArea,EnergyDispersion];