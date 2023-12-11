% %% load overview and reference stack and plot in 3D to select Z planes to align
disp('Loading Reference Stack...')
hRefStackViewer = slap2.gui.ReferenceStackViewer;
refStack = permute(hRefStackViewer.ReferenceStack.data,[2 1 3]);
global zOffset
zOffset = min(hRefStackViewer.ReferenceStack.zs)-1;

fileSuffix = hRefStackViewer.ReferenceStack.sourceTifFile(end-15:end-12); % DMD1 or DMD2

%%
disp('Loading Overview Stack...')
[fns, dr] = uigetfile('*.tif', 'Select Overview Stack', 'multiselect', 'on');
overviewStack = tiffreadVolume([dr fns]);

%% use slap2 align tool to manually find a good alignment

% overviewPlane = 38;
% refPlane = 35; % 22
disp('Align in X-Y...')
hAlign = slap2.gui.AlignTool;
hAlign.addImage(max(overviewStack,[],3));
hAlign.addImage(max(refStack,[],3));
pause();

%%
transform = inv(hAlign.hImages(2).T);

%%
overviewInterpolant = griddedInterpolant(single(overviewStack));

% Assuming refStack is your 3D image
[y_dim, x_dim, ~] = size(refStack);
z_dim = size(overviewStack,3);

% Generate grid of coordinates
[Y, X, Z] = ndgrid(1:y_dim, 1:x_dim, 1:0.5:z_dim);

% Reshape grids into lists of coordinates
y_coords = Y(:);
x_coords = X(:);
z_coords = Z(:);

all_coords = [y_coords, x_coords];

transformCoords = inv(transform) * [all_coords ones(size(all_coords,1),1)]';

%%
transformedOverview = overviewInterpolant([transformCoords(1:2,:)', z_coords]);
transformedOverview = reshape(transformedOverview,size(refStack,1),size(refStack,2),[]);

%%
% figure; imshow3D(transformedOverview)

%%
disp('Align in X-Z...')
hAlignZ = slap2.gui.AlignTool;
hAlignZ.addImage(squeeze(max(transformedOverview,[],1)));
hAlignZ.addImage(squeeze(max(refStack,[],1)));
pause();

%%
transformZ = inv(hAlignZ.hImages(2).T);

fullT = [1 0 0 0;
         0 transformZ(1,1:3);
         0 transformZ(2,1:3);
         0 0 0 1] * ...
        [transform(1,1:2) 0 transform(1,3);
         transform(2,1:2) 0 transform(2,3);
         0 0 2 -1;
         0 0 0 1];

%% use the computed transformation to transform traced points

%%
voxelSize = [0.3643/2 0.3643/2 2];

%% read in SWC of overview stack
disp('Reading SWC...')
[fns, dr] = uigetfile('*.swc', 'Select Horta Tracing', 'multiselect', 'on');

str = fileread([dr fns]);
tmp = regexp(str,'# OFFSET \w+.\w+ \w+.\w+ \w+.\w+','match');
if size(tmp,1) == 0
    offset = [0 0 0]';
else
    offset = textscan(tmp{1}(10:end),'%f');
    offset = offset{1};
end

%%
tmp = regexp(str,'# COLOR \w+.\w+,\w+.\w+,\w+.\w+\r','end');

if ~isempty(tmp)
    swc = textscan(str(tmp+1:end),repmat('%f',1,7),'Delimiter',' ', 'MultipleDelimsAsOne',true,'CollectOutput',true);
    swc = swc{1};
else
    tmp = regexp(str,'#\r','end');
    if ~isempty(tmp)
        swc = textscan(str(tmp+1:end),repmat('%f',1,7),'Delimiter',' ', 'MultipleDelimsAsOne',true,'CollectOutput',true);
        swc = swc{1};
    else
        swc = textscan(str,repmat('%f',1,7),'Delimiter',' ', 'MultipleDelimsAsOne',true,'CollectOutput',true);
        swc = swc{1};
    end
end

%%
disp('Transforming points...')

% clean up overview stack points
keypoints = [];
transformedPts = zeros(size(swc,1),4);
newSWC = [];
for idx = 1:size(swc,1)
    tmpCoord = fullT * round([((swc(idx,[4 3 5]) + offset([2 1 3])' - [0 0 1]) ./ voxelSize) 1])';

    transformedPts(idx,1:3) = tmpCoord(1:3,1)';
    transformedPts(idx,4) = swc(idx,7);
    refinedCoord = mean_shift_point(tmpCoord(1:3,1)', refStack, 10, 5);
    keypoints = [keypoints;  refinedCoord swc(idx,7)];

end

%%
disp('Growing snakes...')
allSnakes = growAStarSnake(round(keypoints),[],1,[],refStack);

%%
disp('Saving JSON...')

json = allSnakes.toJSON();
fid = fopen(['3DSnakes_traced_test_' fileSuffix '.json'],'w','n','UTF-8');
fprintf(fid,'%s',json);
fclose(fid);

disp('done');

%%

function snakes = growAStarSnake(keypts, currpts, index, snakes, volume)
global zOffset
fprintf("Processing index %d\n",index)

children = find(keypts(:,4) == index);

if isempty(currpts) && inFOV(keypts(index,1:3),volume)
    % if index point is outside FOV and there is no active snake growing,
    % stop and return the current batch of snakes
    return
elseif ~isempty(currpts) && inFOV(keypts(index,1:3),volume)
    % if index point is outside FOV and there is an active snake growing,
    % add the current snake and return the new batch of snakes
    
    if size(currpts,1) > 1
        tmp = slap2.gui.refstack.SnakeRoi;
        tmp.pointsXYZR = currpts;
        snakes = [snakes tmp];
    end
    return
end

if size(children,1) > 1
    % if have multiple children, end any current snakes start anew at each
    % child node
    if ~isempty(currpts)
        % if have snake growing, add current point to snake and add snake
        % to batch

        if abs(((currpts(end,2).*4) - keypts(index,1)) / ((currpts(end,1).*4) - keypts(index,2))) <= 1
            path = astar_search3D(volume, (currpts(end,[2 1 3]).*[4 4 1]-[0 0 zOffset]), keypts(index,1:3));
            path(:,3) = round((path(:,3)-1)/3) * 3;
        
            currpts = [currpts; path(2:end,[2 1 3])./[4 4 1]+[0 0 zOffset] repmat([1 2.3],size(path,1)-1,1)];
        end

        if size(currpts,1) > 1
            tmp = slap2.gui.refstack.SnakeRoi;
            tmp.pointsXYZR = currpts;
            snakes = [snakes tmp];
        end
    end

    for i = 1:size(children,1)
        snakes = growAStarSnake(keypts,[],children(i),snakes,volume);
    end
elseif size(children,1) == 1
    % if there is only one child
    if isempty(currpts)
        % and no current snake is growing, add current point to a new snake
        % and grow snake from child

        keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
        currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]+[0 0 zOffset] 1 2.3];

        snakes = growAStarSnake(keypts,currpts,children(1),snakes,volume);
    elseif abs(((currpts(end,2).*4) - keypts(index,1)) / ((currpts(end,1).*4) - keypts(index,2))) <= 1
        % there is a growing snake and the current point is valid (<45deg)
        % to be added, add current point by A* path and continue growing
        % snake from child

        path = astar_search3D(volume, (currpts(end,[2 1 3]).*[4 4 1]-[0 0 zOffset]), keypts(index,1:3));
        path(:,3) = round((path(:,3)-1)/3) * 3;
    
        currpts = [currpts; path(2:end,[2 1 3])./[4 4 1]+[0 0 zOffset] repmat([1 2.3],size(path,1)-1,1)];

        snakes = growAStarSnake(keypts,currpts,children(1),snakes,volume);
    else
        % there is a growing snake and the current point is not valid to be
        % added, end current snake, start new snake with current point and
        % grow from child

        if size(currpts,1) > 1
            tmp = slap2.gui.refstack.SnakeRoi;
            tmp.pointsXYZR = currpts;
            snakes = [snakes tmp];
        end

        keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
        currpts = [keypts(index,[2 1 3])./[4 4 1]+[0 0 zOffset] 1 2.3];

        snakes = growAStarSnake(keypts,currpts,children(1),snakes,volume);
    end
else
    % if there are no children
    if isempty(currpts); return; end % and no snake is growing, return current set of snakes
    if abs(((currpts(end,2).*4) - keypts(index,1)) / ((currpts(end,1).*4) - keypts(index,2))) <= 1
        % and there is a growing snake, and current point is valid to be
        % added, added A* path to current point

        path = astar_search3D(volume, (currpts(end,[2 1 3]).*[4 4 1]-[0 0 zOffset]), keypts(index,1:3));
        path(:,3) = round((path(:,3)-1)/3) * 3;
    
        currpts = [currpts; path(2:end,[2 1 3])./[4 4 1]+[0 0 zOffset] repmat([1 2.3],size(path,1)-1,1)];
    end
    
    % end snake and return new batch of snakes
    if size(currpts,1) > 1
        tmp = slap2.gui.refstack.SnakeRoi;
        tmp.pointsXYZR = currpts;
        snakes = [snakes tmp];
    end
end

end

function infov = inFOV(point, vol)

infov = (round(point(3)) < 1 || round(point(3)) > size(vol,3)) ...
    || (round(point(2)) < 105 || round(point(2)) > 1160) ...
    || (round(point(1)) < 150 || round(point(1)) > 660);

end