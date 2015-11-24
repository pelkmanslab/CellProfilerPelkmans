function ProjImage = CombinePlanesCP3D(Image,Method)
% Help for CombinePlanesCP3D
%
% SHORT DESCRIPTION:
% Small function which makes projection of images.
%
%   
%   Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html



switch lower(Method)
    case 'maximum'
        ProjImage = max(Image,[],3);
    case 'max'
        ProjImage = max(Image,[],3);
    case 'std'
        ProjImage = std(Image,[],3); %... Actually the std projection can be quite a nice spot detection
    case 'sum'
        ProjImage = sum(Image,3); % Sum projection as initially implemented by Markus in IntensityProjectionCP3D
    otherwise
        error([Method 'is no supported projection method']);
end


end