function [matVolumeShape strLabel] = VolumeShapeCP3D(SegCC)
%VOLUMESHAPECP3D extracts basic volumetric and shape features of 3D objects
% IMPORTANT: curerntly only a few outputs are supported, see asteriks *
%

% Basic shape features:     Feature Number:
%%% of 3D Object
% Volume                  |       1 *
% Surface                 |       2
% Extent                  |       3 *
% Loweset Z               |       4 *
% Highest Z               |       5 *
% Height                  |       6 *
% FormFactor              |       7



%%%%%% Create Label of measurements %%%%%%%


strLabel={'Volume', 'Surface', 'Extent', 'LowestZ', 'HighestZ', 'Height', 'FormFactor'};






%%%%%% Obtain Features %%%%%%%




matVolumeShape=NaN(SegCC.NumObjects,7);



% note that regionprops only gives limited information about some of the
% features, which would be available in 2D

% Regionprops: Area of 3D object corresponds to Volume!
  props = regionprops(SegCC,'Area','BoundingBox');

  % Volume                  |       1
  matVolumeShape(:,1) = cell2mat(arrayfun(@(x) x.Area, props,'UniformOutput',0));
  
  
    
  BoundCoord=cell2mat(arrayfun(@(x) x.BoundingBox, props,'UniformOutput',0));
  
  % Extent                  |       3  
  VolBound = BoundCoord(:,4).*BoundCoord(:,5).*BoundCoord(:,6);
  matVolumeShape(:,3) = matVolumeShape(:,1)./VolBound;
  % Smallest Z              |       4
  matVolumeShape(:,4) = BoundCoord(:,3)+0.5;
  % Highest Z               |       5
  matVolumeShape(:,5) = BoundCoord(:,3)-0.5+BoundCoord(:,6);
  % Height                  |       6
  matVolumeShape(:,6) = BoundCoord(:,6);
  
  


end