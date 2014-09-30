function [numChildren ParentID]= relateobjectsCP3D(ChildLabel,ParentLabel,Method)
% RELATEOBJECTSCP3D relates Objects of CHILDLABEL with objects of PARENTLABEL
%
%   Currently it only supports mapping of CHILDLABEL as a CC structure to a
%   labelmatrix PARENTLABEL via the 'centroid' of the child objects.
%
%   (relate module of CP> would choose highest parent idx)
%
%
%
%   --------------------------
%   TO DO: support for different n-dimensional input; different way of
%   finding relations (for instance: centroid of a bent object might not be
%   in real parent object)
%
%   see tech notes below
%
%   ----------------------------------
%   Information about module:
%   Created by TS as part of CP3D to relate objects;


% tech. note A): conversion of CC to 3D labelmatrix easily can require 1GB
% RAM on the other hand this would be the fastest way of general mapping.
% might induce check for availability of memory

% tech. note B): Number output children per object id is determined by the
% max ID, which might be problematic if highest parent id would not be
% highest parentlabel.

% tech. note C): Output numChildren, if no parent, is empty matrix. double
% check that this would be fine in module using realateobjectsCP3D

% tech. note D): Uses ID of 0 for background (which actually would not be
% any ID in CC)

%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CHECK INPUT

% If ParentSegmentation is CC, change to Labelmatrix, note that this will
% temporarlly require more memory. However, it allows to easily apply some
% methods on 3D parents
if isfield(ParentLabel,'PixelIdxList')
    ParentLabel = labelmatrix(ParentLabel);
end

%%%%%%%%%%%%%%%%%%%%%%
%%%%% MEASUREMENT

%%%%%% Trivial Case Handling: either no parent object or no child object
if isfield(ChildLabel,'PixelIdxList')
    maxChild = ChildLabel.NumObjects;
    dimChild = ChildLabel.ImageSize;
else
    maxChild = max(ChildLabel(:));
    dimChild = size(ChildLabel);
end
maxParent = max(ParentLabel(:));
dimParent = size(ParentLabel);

if maxChild == 0 && maxParent == 0;
    numChildren = [];
    ParentID = [];
    return;
else
    if maxChild == 0;
        numChildren = zeros(maxParent,1);
        ParentID = [];
        return;
    end
    if maxParent == 0;
        numChildren = [];
        ParentID = zeros(maxChild,1);
        return;
    end
end

%%%%% Non-trivial case: do linking
%%% Obtain ParentID of each Child
switch Method
    case 'Centroid'
        ChildProps = regionprops(ChildLabel,'Centroid');
        Centroids = round(cell2mat(arrayfun(@(x) (x.Centroid), ChildProps,'UniformOutput',0)));
        % obtain Parent IDs for each child
        if length(dimChild)>2 && (length(dimChild) == length(dimParent))
            if dimChild(3) == dimParent(3)    % only if same amount of z dimensions
                ParentID = ParentLabel(sub2ind2(dimParent(1:3),[Centroids(:,2) Centroids(:,1) Centroids(:,3)]));  %3D to 3D      %note that Centroid column 1 refers to x and not y!
            end
        else
            ParentID = ParentLabel(sub2ind2(dimParent(1:2),[Centroids(:,2) Centroids(:,1)])); % 2D or 3D to 2D       %note that Centroid column 1 refers to x and not y!
        end
    otherwise
        error('Method for Relating not supported');
end

%%% Count amount of Children per Parent
ParentIX = 1: maxParent;   % note that this assumes that the maximum can be derived from the segmentation file - an assumption also done by CPrelateobjects
numChildren=histc(ParentID,ParentIX);


end