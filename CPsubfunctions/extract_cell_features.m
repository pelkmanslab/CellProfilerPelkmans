function [featmat,fnames] = extract_cell_features(C,I,visualize)
% featmat = extract_cell_features(C,I,visualize)  -  Extracts cell level features
%
% In:
%   C = gray level labeled matrix for cells
%   I = cytoplasm (or whatever staining we are interested in)
%   visualize, if 1, plots are shown
% Out:
%   featmat = features
%   fname = feature names

% 13.5.2008 (C) Pekka Ruusuvuori

if nargin < 4
    visualize = 0;
end
C = double(C);
I = double(I);
num_cells = max(C(:));
if visualize == 1
    figure(77)
    imshow(I,[])
end

% initialization of matrix for features
featmat = [];
fnames = {''};
% fearure vector initialization
num_corners_vec = zeros(num_cells,1);
boundaryvar_vec = zeros(num_cells,1);
howmany = 3; %  how many sectors are separated: 3 = border, middle, center of cell
meanIntensity_vec = zeros(num_cells,howmany);
%isonboundary_vec = zeros(num_cells,1);
relativeboundary_vec = zeros(num_cells,1);
%equivalentdiameter_vec = zeros(num_cells,1);
spatialenergy_vec = zeros(num_cells,3);
fd_vec = zeros(num_cells,13);  % fixed value (2*L-1), L be given as parameter to fourier_descriptors function...TODO!

fraction = 0.5; % fraction of energy
% measuring features:
stats = regionprops(C,'Perimeter','EquivDiameter');
equivalentdiameter = [stats.EquivDiameter]';
equivalentdiameter_vec = equivalentdiameter(:);

% detection of cells on border of clusters
BI = C>0;
BIB = bwmorph(BI,'bridge');
boundary = BIB - bwmorph(BIB,'erode');
boundarycellID = unique(boundary.*double(C));  % which cells are on some boundary
boundarycellID = boundarycellID(2:end);  % discard zero
isonboundary(boundarycellID) = 1;
isonboundary_vec = isonboundary(:);
allboundaries = BI - bwmorph(BI,'erode');
%figure, imshow(allboundaries-boundary,[]) % there are some areas in
%"boundary" which do not appear in "allboundaries" -- problem?

% corner detection -- takes time
if visualize == 1
    disp('...detecting boundaries')
end
[Corners,cornersx,cornersy] = detect_corners_inside_cluster(C); %#ok<NASGU>

% mostly vesicle features
for ind = 1:num_cells
    if visualize == 1
        disp(['...extracting features for cell ' num2str(ind) ' out of ' num2str(num_cells)])
    end
    BO = (C == ind);
    [YY1,YY2] = find(BO == 1);
    s1 = min(YY1);
    s2 = min(YY2);
    e1 = max(YY1);
    e2 = max(YY2);
    % one cell binary
    OneCell = BO(s1:e1,s2:e2);
    %cellsize = sum(OneCell(:));
    % one cell intensity 
    CellInt = I(s1:e1,s2:e2).*double(OneCell);
    % corners of one cell
    CellCorners = Corners(s1:e1,s2:e2).*OneCell;
    [LCorners,num_corners] = bwlabel(CellCorners);
    num_corners_vec(ind) = num_corners;
    clear LCorners  % void
    if visualize == 1
        figure(97)
        imshow(CellInt,[])
        title(num2str(ind))
    end
    % save the current cell as small images into folder ./saved_cells
    %save_current_cell(OneCell,OneCellVesicles,CellInt,CellCorners);

    % boundary-related features
    bwb = bwboundaries(OneCell);
    bu = bwb{1};
    if visualize == 1
        figure(99)
        plot(bu(:,2)+s2-1,bu(:,1)+s1-1)
        xlabel(num2str(num_corners))
        hold on
    end
    fl = 30; % filter length for boundary smoothing
    ext_bu = [bu(end-floor(fl/2):end,:); bu; bu(1:ceil(fl/2),:)];
    smooth_bu = [conv(ext_bu(:,1),ones(fl,1)/fl), conv(ext_bu(:,2),ones(fl,1)/fl)];
    smooth_bu = smooth_bu(fl:end-fl-1,:);
    if visualize == 1
        plot(smooth_bu(:,2)+s2-1,smooth_bu(:,1)+s1-1,'r')
    end
    fd_vec(ind,:) = fourier_descriptors(smooth_bu); %#ok<AGROW>
    %distsignature = sqrt((smooth_bu(:,1) - mean(smooth_bu(:,1))).^2 + (smooth_bu(:,2) - mean(smooth_bu(:,2))).^2);
    %distsignature = distsignature - mean(distsignature);
    %distsignature = distsignature / std(distsignature); % what to do with this? :)
    % plot(distsignature)
    boundaryvar = mean(sum((bu - smooth_bu).^2,2));
    boundaryvar_vec(ind) = boundaryvar;
    

    % finding center of the cell
    cellcent = regionprops(double(transpose(OneCell == 1)),'centroid');
    %cind = round(cellcent.Centroid);
    %cc = zeros(size(OneCell));
    %cc(cind(1),cind(2)) = 1;
    %centdist = bwdist(cc);
    %nucdist = bwdist(Nmask);
    outdist = bwdist(~(OneCell == 1));
    if visualize == 1
        figure(77)
        hold on
        plot(cellcent.Centroid(2)+s2-1,cellcent.Centroid(1)+s1-1,'ro')
        hold off
    end

    % find out where the signal is: border, between, middle
    Dpart = outdist(outdist>0);
    Ipart = CellInt(outdist>0);
    meanIntensity = zeros(1,howmany);
    distances = max(Dpart)/howmany:max(Dpart)/howmany:max(Dpart);
    for ind2 = 1:length(distances)
        meanIntensity(ind2) = mean(Ipart(Dpart<=distances(ind2)));
        Dpart(Dpart<=distances(ind2)) = Inf;
    end
    meanIntensity_vec(ind,:) = meanIntensity;
    spatialenergy_vec(ind,:) = spatial_energy_distribution(CellInt,fraction);
    relativeboundary_vec(ind) = sum(sum(boundary.*double(C)==ind))/sum(sum(allboundaries.*double(C)==ind));
end

% Add features in the feature matrix and update feature names
%[featmat,fnames] = addfeature(featmat,fnames,num_corners_vec,'Number of
%corners in cell'); % To be added!
[featmat,fnames] = addfeature(featmat,fnames,boundaryvar_vec,'Cell boundary variance');
[featmat,fnames] = addfeature(featmat,fnames,isonboundary_vec,'Is cell on boundary (binary)');
[featmat,fnames] = addfeature(featmat,fnames,meanIntensity_vec,'Average intensity as a function of distance');
[featmat,fnames] = addfeature(featmat,fnames,relativeboundary_vec,'Relative boundary to empty space');
[featmat,fnames] = addfeature(featmat,fnames,equivalentdiameter_vec,'Equivalent diameter of cell');
[featmat,fnames] = addfeature(featmat,fnames,spatialenergy_vec,'Spatial energy distribution');
[featmat,fnames] = addfeature(featmat,fnames,fd_vec,'Fourier descriptors');
