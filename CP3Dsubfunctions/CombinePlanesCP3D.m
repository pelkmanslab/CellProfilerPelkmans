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



switch Method
    case 'Maximum'
        ProjImage = max(Image,[],3);
    case 'Std'
        ProjImage = std(Image,[],3); %... Actually the std projection can be quite a nice spot detection
    otherwise
        error([Method 'is no supported projection method']);
end


end