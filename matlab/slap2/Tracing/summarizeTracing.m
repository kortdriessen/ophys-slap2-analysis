function summarizeTracing()
%Registers a NeuroLucida or SNT tracing to a set of ExptSummary files
%adds a tracingSummary file to the exptSummary file folder

%script to localize synapses from experimentSummary in a traced stack, and
%vice-versa
%key concepts
%All datasets should be transformed into a unified coordinate system. This is
%provided by SLAP2's DMD2SampleTransform, which is in units of microns in
%the true sample space.

% Inputs:
% -A Reference stack and associated neurolucida tracing as .xml
% - one or more experimentSummary files
% Output:
% -A tree (as xyz coordinates and children relationships) in the reference space (um) that describes the neuron
% -A correspondence map between ROIs on each DMD (and for each recording, if there are multiple recordings) and nodes of the neuron tree above
% -cable distance between each pair of ROIs

import ScanImageTiffReader.*
transformType = 'affine';

xy_spacing_neurolucida = 0.25;
z_spacing_neurolucida = 1.5;
minNodeSpacing = 0.2;

%file picking
hM = msgbox('Select Tracing XML or SWC file', 'modal');
waitfor(hM);
[fn_tracing,dr_tracing] = uigetfile('*.*', 'Select tracing file'); %pick a tracing
tracing_fn = [dr_tracing filesep fn_tracing];
if ~any(fn_tracing) || ~any(dr_tracing)
    disp('Canceled. Aborting summarizeTracing.')
    return
end
[~, ext] = fileparts(tracing_fn);
if strcmpi(ext, 'xml')
    isNeurolucida = true;
else
    isNeurolucida = false;
end


hM = msgbox('Select Large Reference Stack', 'modal');
waitfor(hM);
[fn_tiff, dr_tiff] = uigetfile([dr_tracing filesep '*REFERENCE.tif'], 'Select REFERENCE stack');
tiff_fn = [dr_tiff filesep fn_tiff];
if ~any(fn_tiff) || ~any(dr_tiff)
    disp('Canceled. Aborting summarizeTracing.')
    return
end

mf = dir([tiff_fn(1:end-14) '*.meta']);
if numel(mf)~=1
    hM = msgbox('Select large reference stack metdata file', 'modal');
    waitfor(hM);
    [fn_meta, dr_meta] = uigetfile([dr_tiff filesep '*.meta'], 'Select reference stack metdata file');
    meta_fn = [dr_meta filesep fn_meta];
else
    meta_fn = [mf.folder filesep mf.name];
end
if ~any(fn_meta) || ~any(dr_meta)
    disp('Canceled. Aborting summarizeTracing.')
    return
end

hM = msgbox('Pick all experimentSummary files for this neuron', 'modal');
waitfor(hM);
pathToExptSummary = uipickfiles('filterspec', [dr_tracing filesep '*.mat']); %pick experimentSummaries

%load SLAP2 refIM metadata, which includes the DMD->sample transformation
%meta_fn = '\\allen\aind\scratch\ophys\Maedeh\V1_visual_stimuli\776270\slap2_776270_2025-01-28_13-31-36\Neuron3\stack1\largeStack_20250128_140833_DMD1.meta';
meta = load(meta_fn, '-mat');

%get DMD to sample conversions
DMD1 = meta.machineConfiguration(find(strcmpi({meta.machineConfiguration.instanceName}, 'Path1'),1, 'first'));
DMD2 = meta.machineConfiguration(find(strcmpi({meta.machineConfiguration.instanceName}, 'Path2'),1, 'first'));
tform = zeros(3);tform(3,3) = 1;
tform(1:2,1:2) = DMD1.configuration.dmdPixel2SampleTransform(1:2,1:2);
DMDtform{1} = affinetform2d(tform);
tform(1:2,1:2) = DMD2.configuration.dmdPixel2SampleTransform(1:2,1:2); 
DMDtform{2} = affinetform2d(tform);
switch meta.acquisitionPathName
    case 'Path1'
        reftform = DMDtform{1};
    case 'Path2'
        reftform = DMDtform{2};
    otherwise
        error('refIM acquisition path had unexpected name');
end

%load tracing stack
%tiff_fn = '\\allen\aind\scratch\ophys\Maedeh\V1_visual_stimuli\776270\slap2_776270_2025-01-28_13-31-36\Neuron3\stack1\largeStack_20250128_140833_DMD1-REFERENCE.tif';
A = ScanImageTiffReader(tiff_fn);
rIm = permute(A.data(), [2 1 3]);%matlab image has [Y,X] coordinates
outputView_ref = affineOutputView([800 1280],reftform, 'BoundsStyle', "FollowOutput");
outputView_ref.ImageSize = outputView_ref.ImageSize*4;
rIm= imwarp(rIm,affinetform2d(reftform), 'OutputView',outputView_ref, 'FillValues',nan);
rIm = max(0, rIm, "includemissing");
rIm_lim = prctile(rIm(~isnan(rIm)),99.9);

for imIx = 1:length(A.descriptions)
    sliceMeta = jsondecode(A.descriptions{imIx});
    zPlanes(imIx) = sliceMeta.z;
end
zSpacingStack = median(diff(unique(zPlanes)));
if isempty(zSpacingStack)
    zSpacingStack =1.5;
    warning(['Failed to extract Z spacing from refstack metadata, defaulting to ' num2str(zSpacingStack)])
end

%load tracing
%xml_fn = '\\allen\aind\scratch\ophys\Maedeh\V1_visual_stimuli\776270\slap2_776270_2025-01-28_13-31-36\Neuron3\stack1\tracing1.xml';
if isNeurolucida
    [xyz_raw, parent_raw] = parseNeurolucidaXML(tracing_fn);

    %if multiple trees were traced, ensure there is a single origin point at position 1
    roots = find(parent_raw==0);
    if length(roots)>1
        newRoot = mean(xyz_raw(roots,:));
        xyz = [newRoot ; xyz_raw];
        parent = [0; parent_raw+1];
    end

    %convert to pixels
    xyz_px(:,1) = (xyz(:,1)./xy_spacing_neurolucida);
    xyz_px(:,2) = -(xyz(:,2)./xy_spacing_neurolucida);
    xyz_spacing = [xy_spacing_neurolucida -z_spacing_neurolucida];

else %SWC file
    [nodes, children, xyz_spacing] = parseSWC(tracing_fn);
    parent = [nodes.parent]';

    if isempty(xyz_spacing)
        xyz_spacing = [0.25 0.25 zSpacingStack]; %import the z spacing from the metadata
        for ix = 1:numel(nodes)
            nodes(ix).x = nodes(ix).x.*xyz_spacing(1);
            nodes(ix).y = nodes(ix).y.*xyz_spacing(2);
            nodes(ix).z = nodes(ix).z.*xyz_spacing(3);
        end
    end
    disp('smoothing tree...')
    smoothedNodes = smoothZPositions(nodes,children, 5);
    
    xyz_px(:,1) = ([smoothedNodes.x]./xyz_spacing(1));
    xyz_px(:,2) = ([smoothedNodes.y]./xyz_spacing(2));
    xyz_px(:,3) = ([smoothedNodes.z]./xyz_spacing(3));
end

%transform pixel coords to reference coords
xyzR(:,3) = xyz_px(:,3).*zSpacingStack; %use z coordinates in microns
[xyzR(:,1),xyzR(:,2)] = transformPointsForward(reftform,xyz_px(:,1),xyz_px(:,2));

%upsample
[xyzR, parent] = upsampleCoordinates(xyzR, parent, minNodeSpacing); %in units of microns in the sample;

%sanity check coordinate system
doSanityCheck = false;
if doSanityCheck
    figure, imshow(mean(sqrt(rIm),3),outputView_ref, []);
    hold on, scatter(xyzR(:,1),xyzR(:,2));
end

%compute a dendrogram
%dendrogramXZ = [];
children = parentToChildren(parent);
[dendrogramXZ, projectionXZ] = makeDendrogram(xyzR, children);

%load the experimentSummary
assert(iscell(pathToExptSummary));
nExpts = numel(pathToExptSummary);
%pathToExptSummary{1} = "\\allen\aind\scratch\ophys\Maedeh\V1_visual_stimuli\776270\20250103\slap2_776270_2025-01-03_12-48-51\Neuron3\ExperimentSummary\Summary-250109-101022.mat";
%pathToExptSummary{1} = "\\allen\aind\scratch\ophys\Maedeh\V1_visual_stimuli\776270\20250103\slap2_776270_2025-01-03_12-48-51\Neuron3\ExperimentSummary\Summary-250207-180849.mat";
for expt_ix = 1:nExpts
    load(pathToExptSummary{expt_ix}, 'exptSummary');
    aFolder = fullfile(fileparts(pathToExptSummary{expt_ix}), '..'); %where the activity data should be

    %load reference stack if it has been removed for space reasons
    if isempty(exptSummary.trialTable.refStack{1}.IM)
        disp('Reloading reference image into trial Table...');
        trialTableFn  =  [aFolder filesep 'trialTable.mat'];
        if ~exist(trialTableFn, 'file')
            [trialTableFn, tmpdr] = uigetfile([aFolder filesep '*trialTable.mat'], 'select trial Table');
            trialTableFn = [tmpdr trialTableFn];
        end
        load(trialTableFn, 'trialTable');
        exptSummary.trialTable.refStack = trialTable.refStack;
        clear trialTable;
    end

    %find correspondence between activity recording and reference image
    for DMDix = [2 1]
        firstValidTrial = find(~isempty(exptSummary.E(:,DMDix)),1, 'first');
       
        %load activity metadata and get imaging z plane
        metaFn = dir(strcat(aFolder, filesep, '*DMD', int2str(DMDix), '*.meta'));
        meta = load([metaFn(1).folder filesep metaFn(1).name], '-mat');
        zPlaneA =  meta.AcquisitionContainer.ROIs.rois{1}.z;

        nChans = length(exptSummary.trialTable.refStack{DMDix}.channels);
        aRef = permute(exptSummary.trialTable.refStack{DMDix}.IM(:,:,1:nChans:end), [2 1 3]);
        [minVal, refZind] = min(abs(exptSummary.trialTable.refStack{DMDix}.Zs-zPlaneA));
        assert(minVal<1, 'The activity imaging plane was not found in the corresponding reference image');
        refPlane = aRef(:,:,refZind);

        % get coordinates for transforming the ROIs to the activity refIM space
        cropRow = exptSummary.aData{firstValidTrial,DMDix}.cropRow;
        cropCol = exptSummary.aData{firstValidTrial,DMDix}.cropCol;

        % load manual registration data
        fn_annotations = [aFolder filesep 'tracingRegistrationData_DMD' int2str(DMDix) '.mat'];
        if exist(fn_annotations, 'file')
            load(fn_annotations, 'annotations');
        else
            aIm = imwarp(refPlane, DMDtform{DMDix},  'OutputView', outputView_ref,'FillValues',nan);
            aIm_lim = prctile(aIm(~isnan(aIm)),99.9);

            figure,imshow3D(sqrt(max(0,rIm(:,:,:))/rIm_lim));
            figure, imshow(sqrt(max(0, aIm./aIm_lim)));
            annotations.refPlaneIx = input('Enter best matching plane>>');

            selPtsAct = []; selPtsRef = [];
            ok = false;
            while ~ok
                [selPts1,selPts2] = cpselect(sqrt(aIm./aIm_lim),sqrt(rIm(:,:,annotations.refPlaneIx)./rIm_lim),"Wait",true);
                [selPtsAct(:,2), selPtsAct(:,1)] =  outputView_ref.intrinsicToWorld(selPts1(:,2), selPts1(:,1));
                [selPtsRef(:,2), selPtsRef(:,1)]= outputView_ref.intrinsicToWorld(selPts2(:,2), selPts2(:,1));
                nPts = size(selPtsAct,1);
                ok = nPts>=3;
            end
            annotations.tform_CPs = fitgeotform2d(selPtsAct(:,1:2),selPtsRef(:,1:2),transformType);
            save(fn_annotations, 'annotations');
        end
        
        %get activity ROI coordinates
        firstValidTrial = find(~isempty(exptSummary.E(:,DMDix)),1,'first');
        footprints = exptSummary.E{firstValidTrial,DMDix}.footprints;
        nROIs = size(footprints,3);
        roiPosRaw = [];
        for rix = nROIs:-1:1
            [r,c] = find(footprints(:,:,rix)>0);
            roiPosRaw(rix,2) = mean(r)+cropRow;
            roiPosRaw(rix,1) = mean(c)+cropCol;
        end
        %first, transform raw points to the moving image space
        tmp =  transformPointsForward(DMDtform{DMDix}, roiPosRaw);   

        %now, transform to the reference image space
        roiPosR{DMDix} = cat(2, transformPointsForward(annotations.tform_CPs, tmp), annotations.refPlaneIx*ones(nROIs,1).*xyz_spacing(3)); %positions of all the ROIs in the reference image

        %sanity check
        figure, imshow(sqrt(rIm(:,:,annotations.refPlaneIx)./rIm_lim), outputView_ref); hold on, scatter(roiPosR{DMDix}(:,1), roiPosR{DMDix}(:,2));
        
        RRdistances = squareform(pdist(roiPosR{DMDix}(:,1:2)));

        %build a correspondence table for ROIs->nodes
        roi2tree{DMDix} = nan(1,nROIs);
        offsets =zeros(nROIs,3);
        for iter = 1:3
            for rix = nROIs:-1:1
                %find closest tree point
                d2 =  sum(((xyzR -roiPosR{DMDix}(rix,:)).*[1 1 0.4]).^2, 2);
                [minval, minIx] = min(d2);
                if minval<15
                    roi2tree{DMDix}(rix) =  minIx;
                else
                    roi2tree{DMDix}(rix) =  nan;
                end
            end
            for rix = nROIs:-1:1
                %move each ROI coherently with its neighbors, to align with the tree
                nearby = RRdistances(rix,:)<30 & ~isnan(roi2tree{DMDix});
                if sum(nearby)>1
                    tmp = xyzR(roi2tree{DMDix}(nearby),:) - roiPosR{DMDix}(nearby,:); %local offsets
                    offsets(rix,:) =  median(tmp,1, 'omitnan'); 
                else
                    offsets(rix,:) =  0;
                end
            end
            roiPosR{DMDix} = roiPosR{DMDix} + offsets;
        end
    end

    [drSave] = fileparts(pathToExptSummary{expt_ix});

    %show the recorded synapses on the tracing
    hF = figure('color', 'w');
    scatter3(xyzR(:,1), xyzR(:,2), xyzR(:,3), 10,'blue', 'marker', '.'); hold on;
    scatter3(xyzR(1,1), xyzR(1,2), xyzR(1,3), 600, 'blue','filled') %soma
    colors = {'red', 'magenta'};
    for DMDix = [1 2]
        matched{DMDix} = ~isnan(roi2tree{DMDix});
        scatter3(roiPosR{DMDix}(matched{DMDix},1), roiPosR{DMDix}(matched{DMDix},2), roiPosR{DMDix}(matched{DMDix},3), 'marker', 'o', 'markeredgecolor',colors{DMDix});
        scatter3(roiPosR{DMDix}(~matched{DMDix},1), roiPosR{DMDix}(~matched{DMDix},2), roiPosR{DMDix}(~matched{DMDix},3), 'marker', 'x','markeredgecolor',colors{DMDix});
    end
    savefig(hF, strcat(drSave, filesep, 'tracingSummary-', datestr(now, 'YYmmDD-HHMMSS'), '.fig'))
    %plotTree3D(xyzR, parent, cell2mat(roi2tree))

    %compute cable distance between ROIs
    cableDistance = treePathDistance(xyzR, parent, [cell2mat(roi2tree) 1]); % the distance matrices now include soma as last term

    %compute euclidean distance
    euclideanDistance = squareform(pdist([cell2mat(roiPosR') ; xyzR(1,:)])); % the distance matrices now include soma as last term

    %populate the tracingSummary
    tracingSummary.exptSummary_fn = pathToExptSummary;
    tracingSummary.coords = xyzR;  %xyz coordinates in the reference space
    tracingSummary.ROIcoords = roiPosR ;
    tracingSummary.matched = matched;
    tracingSummary.parents = parent;
    tracingSummary.children = children;
    tracingSummary.refIm2Rtform = reftform.A;
    tracingSummary.roi2tree = roi2tree;
    tracingSummary.cableDistance = cableDistance;
    tracingSummary.euclideanDistance = euclideanDistance;
    tracingSummary.dendrogramXZ = dendrogramXZ;
    tracingSummary.projectionXZ = projectionXZ;
    tracingSummary.tiffFn = tiff_fn;
    tracingSummary.xmlFn = tracing_fn;
    tracingSummary.metaFn = meta_fn;

    exptSummary.tracing = tracingSummary;
    save(pathToExptSummary, 'exptSummary')
    %save(strcat(drSave, filesep, 'tracingSummary-', datestr(now, 'YYmmDD-HHMMSS'), '.mat'), 'tracingSummary');
end


end


function [coordinates, parents] = parseNeurolucidaXML(filename)
% Parse the Neurolucida .xml file to extract x, y, z coordinates and parent-child relationships
%
% Input:
%   filename - String, path to the Neurolucida .xml file
% Output:
%   coordinates - Nx3 matrix containing x, y, z coordinates
%   parents - Nx1 vector containing the parent indices for each point

% Read the XML file
try
    xmlDoc = xmlread(filename);
catch
    error("Error reading the XML file. Please check the file path.");
end

% Get all 'point' elements (assuming Neurolucida uses <point> tags)
pointNodes = xmlDoc.getElementsByTagName('point');
numPoints = pointNodes.getLength();

% Initialize matrices to store coordinates and parent indices
coordinates = zeros(numPoints, 3);
parents = zeros(numPoints, 1); % Default parent is 0 (no parent)

currentIndex = 0;

% Start parsing from the root element
root = xmlDoc.getDocumentElement();
parseNode(root, 0); % Root has no parent, so parent index is 0

coordinates = coordinates(1:currentIndex,:);
parents = parents(1:currentIndex);

% Recursive traversal of the XML structure to assign indices and parents
    function parseNode(node, currentParentIndex)
        %nonlocal currentIndex;
        %nonlocal coordinates;
        %nonlocal parents;
        %nonlocal parentStack;

        currentParent = currentParentIndex;
        % Check if the current node is a <point> element
        if strcmp(node.getNodeName(), 'point')
            % Increment the current index
            currentIndex = currentIndex + 1;

            % Get the attributes of the <point> element
            x = str2double(node.getAttribute('x'));
            y = str2double(node.getAttribute('y'));
            z = str2double(node.getAttribute('z'));

            % Store the coordinates
            coordinates(currentIndex, :) = [x, y, z];

            % Assign the current node's parent
            parents(currentIndex) = currentParentIndex;
            currentParent = currentIndex;
        end

        % Recursively process child nodes
        childNodes = node.getChildNodes();
        numChildNodes = childNodes.getLength();
        for i = 0:numChildNodes-1
            childNode = childNodes.item(i);
            if strcmpi(childNode.getNodeName(), 'point')
                nextParent = currentIndex+1;
                parseNode(childNode, currentParent);
                currentParent = nextParent;
            elseif strcmpi(childNode.getNodeName(),'branch')
                nextParent = currentParent;
                parseNode(childNode, currentParent);
                currentParent = nextParent;
            elseif strcmpi(childNode.getNodeName(),'tree')
                nextParent = currentParent;
                parseNode(childNode, 0);
                currentParent = nextParent;
            end
        end
    end

% Display the results
fprintf('Parsed %d points from the XML file.\n', numPoints);
end

% Example Usage:
% filename = 'path_to_neurolucida_file.xml';
% [coords, parentVec] = parseNeurolucidaXML(filename);
% disp(coords);
% disp(parentVec);

function plotTree3D(coordinates, parents, plotPoints)
% Plot a tree structure in 3D based on coordinates and parent relationships
%
% Inputs:
%   coordinates - Nx3 matrix of x, y, z coordinates
%   parents - Nx1 vector specifying the parent index for each point

% Create a 3D figure
figure;
hold on;
grid on;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Tree Structure');

% Loop through each point and plot a line to its parent
for i = 1:size(coordinates, 1)
    parentIdx = parents(i);
    if parentIdx > 0
        % Get the coordinates of the current point and its parent
        currentPoint = coordinates(i, :);
        parentPoint = coordinates(parentIdx, :);

        % Plot a line connecting the current point to its parent
        plot3([currentPoint(1), parentPoint(1)], ...
            [currentPoint(2), parentPoint(2)], ...
            [currentPoint(3), parentPoint(3)], ...
            'b-', 'LineWidth', 1);
    end
end

% Plot the points themselves
scatter3(coordinates(plotPoints, 1), coordinates(plotPoints, 2), coordinates(plotPoints, 3), ...
    50, 'r', 'filled', 'MarkerEdgeColor', 'k');

% Finalize the plot
view(3); % Set a 3D view
hold off;
end

% Example Usage:
% coordinates = [0 0 0; 0.5 0 0; 1 0 0; 1.5 0 0; 2 0 0];
% parents = [0; 1; 2; 3; 4];
% plotTree3D(coordinates, parents);

function [newCoordinates, newParents] = upsampleCoordinates(coordinates, parents, maxDistance)
% Upsample the coordinates such that the distance between connected nodes
% is no greater than maxDistance.
%
% Inputs:
%   coordinates - Nx3 matrix of x, y, z coordinates
%   parents - Nx1 vector specifying the parent index for each point
%   maxDistance - Maximum allowed distance between connected nodes
%
% Outputs:
%   newCoordinates - Mx3 matrix of upsampled x, y, z coordinates
%   newParents - Mx1 vector specifying the parent index for each new point

% Initialize new coordinate and parent lists


N = size(coordinates,1);
dist2parent = zeros(N,1);
for ix = 1:N
    % If the current point has a parent, process the connection
    if parents(ix)>0
        dist2parent(ix) = norm(coordinates(ix,:)-coordinates(parents(ix),:));
    end
end
nPointsNeeded = ceil(dist2parent./maxDistance);
newCoordinates = nan(sum(nPointsNeeded)+1,3);
newParents = nan(sum(nPointsNeeded)+1,1);
newLastIndex = 1:N;

% Start processing each point
newIndex = 1; % Index to add the next point
for i = 1:size(coordinates, 1)
    % Add the current point to the new list
    currentPoint = coordinates(i, :);
    
    % If the current point has a parent, add points up to and including
    % this point
    parentIdx = parents(i);
    if parentIdx > 0
        if dist2parent(i) > maxDistance % Add interpolated points if the distance exceeds maxDistance
            parentPoint = coordinates(parentIdx, :);
            tValues = linspace(0, 1, nPointsNeeded(i) + 1); % Include start and end
            
            %link the first interpolated point to the last point of the parent
            interpolatedPoint = parentPoint + tValues(2)*(currentPoint - parentPoint);
            newCoordinates(newIndex, :) = interpolatedPoint;
            newParents(newIndex, 1) = newLastIndex(parentIdx);
            newIndex = newIndex + 1;

            % link remaining points in a chain
            for t = tValues(3:end)
                interpolatedPoint = parentPoint + t * (currentPoint - parentPoint);
                newCoordinates(newIndex, :) = interpolatedPoint;
                newParents(newIndex, 1) = newIndex-1; % Link to the previous point
                newLastIndex(i) = newIndex;
                newIndex = newIndex + 1;
            end
        else %distance is already ok
            newCoordinates(newIndex, :) = currentPoint;
            newParents(newIndex, 1) = newLastIndex(parentIdx);
            newLastIndex(i) = newIndex;
            newIndex = newIndex + 1;
        end
    else
        % No parent, just add the point and set parent index to 0
        newCoordinates(newIndex, :) = currentPoint;
        newParents(newIndex, 1) = 0;
        newLastIndex(i) = newIndex;
        newIndex = newIndex + 1;
    end
end
end

% Example Usage:
% coordinates = [0 0 0; 1 0 0; 2 0 0]; % Example input coordinates
% parents = [0; 1; 2]; % Example parent indices
% maxDistance = 0.1; % Maximum allowed distance between nodes
% [newCoords, newParents] = upsampleCoordinates(coordinates, parents, maxDistance);
% disp(newCoords);
% disp(newParents);

function children = parentToChildren(parent)
children = false(length(parent));
for ix = 1:length(parent)
    if parent(ix)>0
        children(parent(ix), ix) = true;
    end
end

end


%load tracing
%install the NLMorphology converter from http://neuronland.org/NLMorphologyConverter/NLMorphologyConverter.html
%pathToConverter = '"C:\Program Files (x86)\Neuronland\NLMorphologyConverter\NLMorphologyConverter.exe"';
%convert the tracing coordinates to the
% write_fn = [xml_fn(1:end-4) '.swc'];
% formatStr = 'SWC';
% system([pathToConverter ' ' xml_fn ' ' write_fn ' ' formatStr])
% T = readtable(write_fn, 'FileType', 'text');


function distances = treePathDistance(coordinates, parents, queryNodes)
    % Computes the distance along the edges of a tree between each pair of query nodes
    %
    % Inputs:
    %   coordinates - Nx3 matrix of x, y, z coordinates
    %   parents - Nx1 vector specifying the parent index for each point
    %   queryNodes - Mx1 vector of query node indices
    %
    % Output:
    %   distances - MxM matrix where distances(i, j) is the path distance between
    %               queryNodes(i) and queryNodes(j)
 
    % Number of query nodes
    numQueries = length(queryNodes);
    
    % Initialize the distance matrix
    distances = zeros(numQueries, numQueries);
    
    % Precompute edge lengths between each node and its parent
    numNodes = size(coordinates, 1);
    edgeLengths = zeros(numNodes, 1);
    for i = 1:numNodes
        if parents(i) > 0
            % Compute the Euclidean distance between the node and its parent
            edgeLengths(i) = norm(coordinates(i, :) - coordinates(parents(i), :));
        else
            % Root node has no parent
            edgeLengths(i) = 0;
        end
    end
    
    % Compute the path distance between each pair of query nodes
    for i = 1:numQueries
        for j = i+1:numQueries
            % Get the indices of the two query nodes
            node1 = queryNodes(i);
            node2 = queryNodes(j);
            
            % Find the path distance between node1 and node2
            distances(i, j) = computePathDistance(node1, node2, parents, edgeLengths);
            distances(j, i) = distances(i, j); % Symmetric matrix
        end
    end
end
 
function distance = computePathDistance(node1, node2, parents, edgeLengths)
    % Helper function to compute the path distance between two nodes
    %
    % Inputs:
    %   node1, node2 - Indices of the two nodes
    %   parents - Nx1 vector specifying the parent index for each point
    %   edgeLengths - Nx1 vector of edge lengths between each node and its parent
    %
    % Output:
    %   distance - Path distance between node1 and node2
    
    % Get the paths from each node to the root
    path1 = getPathToRoot(node1, parents);
    path2 = getPathToRoot(node2, parents);
    
    % Find the lowest common ancestor (LCA)
    commonAncestor = findLCA(path1, path2);
    
    % Compute the distance for each path segment
    distance = sum(edgeLengths(path1(find(path1 == commonAncestor, 1, 'last')+1:end))) + ...
               sum(edgeLengths(path2(find(path2 == commonAncestor, 1, 'last')+1:end)));
    % distance = sum(edgeLengths(path1(1:find(path1 == commonAncestor, 1)))) + ...
    %            sum(edgeLengths(path2(1:find(path2 == commonAncestor, 1))));
end
 
function path = getPathToRoot(node, parents)
    % Helper function to get the path from a node to the root
    %
    % Inputs:
    %   node - Index of the node
    %   parents - Nx1 vector specifying the parent index for each point
    %
    % Output:
    %   path - Vector of node indices from the given node to the root
    
    path = [];
    while node > 0
        path = [node; path]; % Add the node to the path
        node = parents(node); % Move to the parent
    end
end
 
function commonAncestor = findLCA(path1, path2)
    % Helper function to find the lowest common ancestor (LCA) of two nodes
    %
    % Inputs:
    %   path1, path2 - Paths from two nodes to the root
    %
    % Output:
    %   commonAncestor - The index of the lowest common ancestor
    
    % Compare the paths and find the last common node
    minLength = min(length(path1), length(path2));
    commonAncestor = 0; % Default to root
    for i = 1:minLength
        if path1(i) == path2(i)
            commonAncestor = path1(i);
        else
            break;
        end
    end
end
 
% Example Usage:
% coordinates = [0 0 0; 1 0 0; 1 1 0; 2 1 0; 2 2 0];
% parents = [0; 1; 2; 3; 4];
% queryNodes = [2, 4, 5]; % Query nodes
% distances = treePathDistance(coordinates, parents, queryNodes);
% disp(distances);