function [Corners,cornersx,cornersy] = detect_corners_inside_cluster(C)
% corner detection (takes time!)
% IN:
%   C = labeled image
%
Corners = zeros(size(C));
padC = zeros(size(C)+2);
padC(2:end-1,2:end-1) = C;
for ii1 = 2:size(C,1)+1
    for ii2 = 2:size(C,2)+1
        Corners(ii1-1,ii2-1) = length(unique(padC(ii1-1:ii1+1,ii2-1:ii2+1)));
    end
end
[cornersx,cornersy] = find(Corners>3);
Corners = imdilate(Corners>3,ones(3));