function [featmat,allvesicledata,fnames, vfnames, riplyMatrix measurementNode] = extract_vesicle_features(C,I,V, rad1, rad2,visualize)
% featmat = extract_vesicle_features(C,I,V)  -  Extracts cell level features
%
% In:
%   C = gray level labeled matrix for cells
%   I = cytoplasm (or whatever staining we are interested in)
%   V = vesicle mask (binary)
%   visualize, if 1, plots are shown
% Out:
%   featmat = features (per cell)
%   allvesicledata = all data from individual vesicles (distributions)
%   fnames = feature names
%   vfnames = individual vesicle feature name
%   riplyMatrix = riply's K function value

% 13.5.2008 (C) Pekka Ruusuvuori
% 19.02.2009 Sharif, Chowdhury output format changed and new statistical
% measurements added
% min, max intensity distance from fluorescent centre added

if nargin < 6
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
num_vesicles_vec = zeros(num_cells,1);
meanvesint_vec = NaN(num_cells,1);
meanvessize_vec = NaN(num_cells,1);
numvesiclespercellsize_vec = NaN(num_cells,1);
meandistance2celloutline_vec = NaN(num_cells,1);
stddistance2celloutline_vec = NaN(num_cells,1);
meandistance2cellcenter_vec = NaN(num_cells,1);
stddistance2cellcenter_vec = NaN(num_cells,1);
meanvesdist_vec = NaN(num_cells,1);
stdvesdist_vec = NaN(num_cells,1);
minvesdist_vec = NaN(num_cells,1);    
% distributions for features
alld2co = [];
alld2cc = [];
allvesint = [];
allvessize = [];
allvesangles = [];
ves2cellmap = [];
allvesdist = [];

% vesicle features

if rad2 > 0 && rad1<=rad2
    riplyMatrix = zeros(num_cells, rad2);
else
    riplyMatrix = 0;
end

measurementNode = cell(0);


for ind = 1:num_cells
    if visualize == 1
        disp(['...extracting vesicle features for cell ' num2str(ind) ' out of ' num2str(num_cells)])
    end
    BO = (C == ind);
    [YY1,YY2] = find(BO == 1);
    s1 = min(YY1);
    s2 = min(YY2);
    e1 = max(YY1);
    e2 = max(YY2);
    % one cell binary
    OneCell = BO(s1:e1,s2:e2);
    cellsize = sum(OneCell(:));
    % one cell intensity 
    CellInt = I(s1:e1,s2:e2).*double(OneCell);
    
    pLIst = regionprops( double(OneCell >0 ) , 'PixelList');
    
    
    
    sumX = 0;
    sumY = 0;
    [rowList colList] = size(pLIst.PixelList);
    totalIntensity = 0;
    for row = 1: rowList
        
        intVal = CellInt(pLIst.PixelList(row,2), pLIst.PixelList(row,1) );
        sumX = sumX + pLIst.PixelList(row,1) * intVal;
        sumY = sumY + pLIst.PixelList(row,2) * intVal;
         totalIntensity =  totalIntensity  + intVal;
    end
    flCentX = 0;
    flCentY = 0;
%    intVal
    if totalIntensity > 0
        flCentX = sumX / totalIntensity ;
        flCentY = sumY / totalIntensity ;
    end
    
    
    % vesicles of one cell binary
    OneCellVesicles = V(s1:e1,s2:e2).*OneCell;
    %Nmask = BI(s1:e1,s2:e2);
    if visualize == 1
        figure(97)
        imshow(CellInt,[])
        title(num2str(ind))
    end
    % save the current cell as small images into folder ./saved_cells
    % save_current_cell(OneCell,OneCellVesicles,CellInt,CellCorners);

    % vesicle labels matrix
    
    if max(max(OneCellVesicles))<1.5
        [v_label,num_vesicles] = bwlabel(OneCellVesicles); % the vesicles are not labelled
    else
        [v_label,num_vesicles] = getLabel(OneCellVesicles);
    end
    % number of vesicles per cell
    num_vesicles_vec(ind) = num_vesicles;
    % init of distribution of size of the vesicles
    vessize = zeros(num_vesicles,1);
    % init of distribution of intensity of the vesicles
    vesint = zeros(num_vesicles,1);
    
    if rad2 > 0 && rad1<=rad2
        riplyMatrix(ind,1:rad2) = RipleysKFunc( OneCell, v_label ,rad1, rad2, 0);
    end
    
    
    % finding center of the cell
    cellcent = regionprops(double(transpose(OneCell == 1)),'centroid');
    cind = round(cellcent.Centroid);
    cc = zeros(size(OneCell));
    cc(cind(1),cind(2)) = 1;
    centdist = bwdist(cc);
    
    cc = cc*0;
    cc( round( flCentY ), round( flCentX ) ) = 1;
    centdistFl = bwdist(cc);
    
    %nucdist = bwdist(Nmask);
    outdist = bwdist(~(OneCell == 1));
    if visualize == 1
        figure(77)
        hold on
        plot(cellcent.Centroid(2)+s2-1,cellcent.Centroid(1)+s1-1,'ro')
        hold off
    end
    % go through all the vesicles
    if num_vesicles > 0
        ang = regionprops(v_label,'orientation');
        vesangles = [ang.Orientation];
        vesangles = vesangles(:);
        %d2no = zeros(num_vesicles,1);
        d2co = zeros(num_vesicles,1);
        
        d2cc = zeros(num_vesicles,1);
        d2Flcc = zeros(num_vesicles,1);
        
        for kv = 1:num_vesicles
            ves = (v_label == kv);
            vessize(kv) = sum(ves(:));
            vesint(kv) = sum(sum(ves.*CellInt)) / vessize(kv);

            a = centdist.*ves;
            a(a == 0) = max(a(:));
            d2cc(kv) = min(a(:));
            
            a = centdistFl.*ves;
            a(a == 0) = max(a(:));
            d2Flcc(kv) = min(a(:));


            %b = nucdist.*ves;
            %b(b == 0) = max(b(:));
            %d2no(kv) = min(b(:));

            c = outdist.*ves;
            c(c == 0) = max(c(:));
            d2co(kv) = min(c(:));
        end
        
        maxVector = -1*ones(1,num_vesicles);
        minVector = -1*ones(1,num_vesicles);
        
        [rowI colI]  = size(CellInt) ;
        
        for row= 1:rowI
            for col= 1:colI
                if v_label( row, col ) > 0
                   colrV = v_label(row, col) ;
                   
                   if maxVector( colrV ) <  CellInt(row,col)
                       maxVector( colrV ) =  CellInt(row,col) ;
                   end
                   
                   if minVector( colrV ) >   CellInt(row,col) || minVector( colrV ) < 0
                       minVector( colrV ) =  CellInt(row,col) ;
                   end
                   
                    
                end
            end
        end
        
        
        
        
        % collect OneCellVesicles level information
        alld2co = [alld2co; d2co];
        alld2cc = [alld2cc; d2cc];
        allvesint = [allvesint; vesint(:)];
        allvessize = [allvessize; vessize(:)];
        allvesangles = [allvesangles; vesangles(:)];
        ves2cellmap = [ves2cellmap; ones(num_vesicles,1)*double(ind)];

        % distance from a OneCellVesicles to others
        ves_center = regionprops(double(transpose(v_label)),'centroid');
        if visualize == 1
            figure(77)
            hold on
            for iii = 1:length(ves_center)
                plot(ves_center(iii).Centroid(2)+s2-1,ves_center(iii).Centroid(1)+s1-1,'bx')
            end
            hold off
        end
        vesdist = zeros(num_vesicles);
        for ii = 1:num_vesicles
            for jj = 1:num_vesicles
                if (ii == jj)
                    continue
                else
                    vc1 = ves_center(ii).Centroid;
                    vc2 = ves_center(jj).Centroid;
                    vesdist(ii,jj) = sqrt((vc1(1) - vc2(1))^2 + (vc1(2) - vc2(2))^2);
                end
            end
        end
        vesdist_mean = mean(vesdist(triu(vesdist)~=0));
        vesdist_std = std(vesdist(triu(vesdist)~=0));
        vesdist_min = min(vesdist(triu(vesdist)~=0));
        if isempty(vesdist_min)
            vesdist_min = NaN;
        end
        allvesdist = [allvesdist; vesdist(triu(vesdist)~=0)]; % this is a bit longer that other "all"-variables

        if num_vesicles ~= 0
            meanvessize_vec(ind) = sum(vessize)/num_vesicles;
            meanvesint_vec(ind) = sum(vesint)/num_vesicles;
        end
        measurementNode{ind} = struct;
        measurementNode{ind} = constructSingleMeasurementNode( measurementNode{ind},num_vesicles, 'Number_of_', 'sum' ,'Total');
        measurementNode{ind} = constructSingleMeasurementNode( measurementNode{ind},num_vesicles/cellsize, 'Per_U_Carea_Number_of_', 'sum' ,'Total');
        
        measurementNode{ind} = constructMeasurementNode(measurementNode{ind},vessize, 'Size_of_');
        measurementNode{ind} = constructMeasurementNode(measurementNode{ind},vesint, 'Intensity_of_');
        measurementNode{ind} = constructMeasurementNode(measurementNode{ind},d2co, 'Distance_from_Cell_outline_to_');
        measurementNode{ind} = constructMeasurementNode(measurementNode{ind},d2cc, 'Distance_from_Cell_centre_to_');
        measurementNode{ind} = constructMeasurementNode(measurementNode{ind},d2Flcc, 'Distance_from_Fluorescent_centre_to_');
        
        measurementNode{ind} = constructMeasurementNode(measurementNode{ind},vesdist(triu(vesdist)~=0), 'Distance_Between_');
        
        measurementNode{ind} = constructMeasurementNode(measurementNode{ind},minVector, 'Minimum_Intensity_of_Individual_');
        measurementNode{ind} = constructMeasurementNode(measurementNode{ind},maxVector, 'Maximum_Intensity_of_Individual_');
        
        
        
        
        numvesiclespercellsize_vec(ind) = num_vesicles/cellsize;
        meandistance2celloutline_vec(ind) = mean(d2co);
        stddistance2celloutline_vec(ind) = std(d2co);
        meandistance2cellcenter_vec(ind) = mean(d2cc);
        stddistance2cellcenter_vec(ind) = std(d2cc);
        meanvesdist_vec(ind) = vesdist_mean;
        stdvesdist_vec(ind) = vesdist_std;
        minvesdist_vec(ind) = vesdist_min;
    end
end

% Add features in the feature matrix and update feature names
[featmat,fnames] = addfeature(featmat,fnames,num_vesicles_vec,'Number of vesicles'); %
[featmat,fnames] = addfeature(featmat,fnames,meanvessize_vec,'Average vesicle size'); %
[featmat,fnames] = addfeature(featmat,fnames,meanvesint_vec,'Average vesicle intensity'); %
[featmat,fnames] = addfeature(featmat,fnames,numvesiclespercellsize_vec,'Number of vesicles / cell size');
[featmat,fnames] = addfeature(featmat,fnames,meandistance2celloutline_vec,'Average distance of vesicles to cell outline');
[featmat,fnames] = addfeature(featmat,fnames,stddistance2celloutline_vec,'Std of distance of vesicles to cell outline');
[featmat,fnames] = addfeature(featmat,fnames,meandistance2cellcenter_vec,'Average distance of vesicles to cell center');
[featmat,fnames] = addfeature(featmat,fnames,stddistance2cellcenter_vec,'Std of distance of vesicles to cell center');
[featmat,fnames] = addfeature(featmat,fnames,meanvesdist_vec,'Average distance between vesicles');
[featmat,fnames] = addfeature(featmat,fnames,stdvesdist_vec,'Std of distance between vesicles');
[featmat,fnames] = addfeature(featmat,fnames,minvesdist_vec,'Minimum distance between vesicles');

% From this matrix it's possible to get the distributions also...
% ves2cellmap gives the cell identifier for each OneCellVesicles
allvesicledata = [alld2co,alld2cc,allvesint,allvessize,allvesangles,ves2cellmap];
vfnames ={'Distance from Cell Outline',  'Distance from Cell Centre', 'Mean Intensity', 'Area' , 'Orientation Angle', 'CellID' }
