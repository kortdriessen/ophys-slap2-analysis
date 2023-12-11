% %% load overview and reference stack and plot in 3D to select Z planes to align
disp('Loading Reference Stack...')
hRefStackViewer = slap2.gui.ReferenceStackViewer;
refStack = permute(hRefStackViewer.ReferenceStack.data,[2 1 3]);
global zOffset
zOffset = min(hRefStackViewer.ReferenceStack.zs)-1;

fileSuffix = hRefStackViewer.ReferenceStack.sourceTifFile(end-15:end-12); % DMD1 or DMD2

%%
% [fns, dr] = uigetfile('*.tif', 'Select Reference Stack', 'multiselect', 'on');
% refStack = tiffreadVolume([dr fns]);
disp('Loading Overview Stack...')
[fns, dr] = uigetfile('*.tif', 'Select Overview Stack', 'multiselect', 'on');
overviewStack = tiffreadVolume([dr fns]);

% figure; imshow3D(refStack)
% figure; imshow3D(overviewStack)

%% use slap2 align tool to manually find a good alignment

% overviewPlane = 38;
% refPlane = 35; % 22
disp('Align in X-Y...')
hAlign = slap2.gui.AlignTool;
hAlign.addImage(max(overviewStack,[],3));
hAlign.addImage(max(refStack,[],3));
pause();
% hAlign.addImage(overviewStack(:,:,overviewPlane));
% hAlign.addImage(refStack(:,:,refPlane));

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

%%
% 
% fullT = [transform(1,1:2) 0 transform(1,3);
%          transform(2,1:2) 0 transform(2,3);
%          0 0 2 refPlane-2*overviewPlane;
%          0 0 0 1];


%% use the computed transformation to transform traced points

% overviewTracePts = [706 728 764; 651 649 648; 1 1 1]; % Y by X by 1
% 
% refTracePts = fullT * overviewTracePts;

%%
voxelSize = [0.3643/2 0.3643/2 2];

%% read in SWC of overview stack
disp('Reading SWC...')

% str = fileread("Z:\ophys\BCI\696662_SBCI05\2p\fov1\transformed_swc\FOV1_Neuron2.swc");
str = fileread("Z:\ophys\BCI\696662_SBCI05\2p\fov1\FOV1_Neuron2.swc");
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
    % tmpCoord = fullT * [swc(idx,[4 3 5]) 1]';
    transformedPts(idx,1:3) = tmpCoord(1:3,1)';
    transformedPts(idx,4) = swc(idx,7);
    refinedCoord = mean_shift_point(tmpCoord(1:3,1)', refStack, 10, 5);
    % if ((round(transformedPts(idx,3)) < 1 || round(transformedPts(idx,3)) > 50) ...
    %     || (round(transformedPts(idx,2)) < 105 || round(transformedPts(idx,2)) > 1160) ...
    %     || (round(transformedPts(idx,1)) < 150 || round(transformedPts(idx,1)) > 660))
    %     continue
    % else
        % newSWC(end+1,:) = [swc(idx,1:2) tmpCoord([2 1 3],1)' 1 swc(idx,7)];
        % refinedCoord = mean_shift_point(tmpCoord(1:3,1)', refStack, 5, 5);
        % if ((round(refinedCoord(3)) < 1 || round(refinedCoord(3)) > 50) ...
        %     || (round(refinedCoord(2)) < 105 || round(refinedCoord(2)) > 1160) ...
        %     || (round(refinedCoord(1)) < 150 || round(refinedCoord(1)) > 660))
        %     continue
        % end
        keypoints = [keypoints;  refinedCoord swc(idx,7)];
    % end
end

% %%
% writematrix(newSWC,'refStack_20231003_163321_DMD1_snakes.swc','Delimiter','tab','FileType','text');
% 
% %%
% str = fileread("C:\Users\michael.xie\Downloads\refined_refStack_20231003_163321_DMD1_snakes\refStack_20231003_163321_DMD1_snakes.swc");
% tmp = regexp(str,'# OFFSET \w+.\w+ \w+.\w+ \w+.\w+','match');
% if size(tmp,1) == 0
%     offset = [0 0 0]';
% else
%     offset = textscan(tmp{1}(10:end),'%f');
%     offset = offset{1};
% end
% 
% %%
% tmp = regexp(str,'# COLOR \w+.\w+,\w+.\w+,\w+.\w+\r','end');
% 
% if ~isempty(tmp)
%     swc = textscan(str(tmp+1:end),repmat('%f',1,7),'Delimiter',' ', 'MultipleDelimsAsOne',true,'CollectOutput',true);
%     swc = swc{1};
% else
%     tmp = regexp(str,'#\r','end');
%     swc = textscan(str(tmp+1:end),repmat('%f',1,7),'Delimiter',' ', 'MultipleDelimsAsOne',true,'CollectOutput',true);
%     swc = swc{1};
% end
% 
% keypoints = zeros(size(swc,1),4);
% newSWC = zeros(size(swc));
% for idx = 1:size(swc,1)
%     keypoints(idx,1:3) = swc(idx,[4 3 5]);
%     keypoints(idx,4) = swc(idx,7);
% end
% 

%%
% allSnakes = growSnake(keypoints,[],1,[]);
disp('Growing snakes...')
allSnakes = growAStarSnake(round(keypoints),[],1,[],refStack);

%%
disp('Saving JSON...')

json = allSnakes.toJSON();
fid = fopen(['3DSnakes_traced_' fileSuffix '.json'],'w','n','UTF-8');
fprintf(fid,'%s',json);
fclose(fid);

disp('done');

%%
% refPointOverview = [792 672  35];
% refPointRef = fullT * [refPointOverview 1]';
% refPointRef = refPointRef(1:3);
% 
% normalOverview = [1 0 0];
% normalRef = fullT * [normalOverview 1]' - fullT * [0 0 0 1]';
% normalRef = normalRef(1:3);
% 
% transposedStack = permute(refStack,[3 1 2]);
% 
% [B,x,y,z] = obliqueslice(transposedStack,refPointRef([1 3 2])',normalRef([1 3 2])');

%%

function snakes = growAStarSnake(keypts, currpts, index, snakes, volume)
global zOffset
fprintf("Processing index %d\n",index)

children = find(keypts(:,4) == index);

if isempty(currpts) && ((round(keypts(index,3)) < 1 || round(keypts(index,3)) > size(volume,3)) ...
    || (round(keypts(index,2)) < 105 || round(keypts(index,2)) > 1160) ...
    || (round(keypts(index,1)) < 150 || round(keypts(index,1)) > 660))
    % if index point is outside FOV and there is no active snake growing,
    % stop and return the current batch of snakes
    return
elseif ~isempty(currpts) && ((round(keypts(index,3)) < 1 || round(keypts(index,3)) > size(volume,3)) ...
    || (round(keypts(index,2)) < 105 || round(keypts(index,2)) > 1160) ...
    || (round(keypts(index,1)) < 150 || round(keypts(index,1)) > 660))
    % if index point is outside FOV and there is an active snake growing,
    % add the current snake to the batch of snakes and return snakes

    % if abs(((currpts(end,2).*4) - keypts(index,1)) / ((currpts(end,1).*4) - keypts(index,2))) <= 1
    %     path = astar_search3D(volume, (currpts(end,[2 1 3]).*[4 4 1]+[0 0 zOffset]), keypts(index,1:3));
    %     path(3,:) = round((path(:,3)-1)/3) * 3;
    %     % keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
    % 
    %     currpts = [currpts; path(2:end,[2 1 3])./[4 4 1]-[0 0 zOffset] repmat([1 2.3],size(path,1)-1,1)];
    %     % currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]-[0 0 zOffset] 1 2.3];
    % end
    
    if size(currpts,1) > 1
        tmp = slap2.gui.refstack.SnakeRoi;
        tmp.pointsXYZR = currpts;
        snakes = [snakes tmp];
    end
    return
end

if size(children,1) > 1
    if ~isempty(currpts)
        if abs(((currpts(end,2).*4) - keypts(index,1)) / ((currpts(end,1).*4) - keypts(index,2))) <= 1
            path = astar_search3D(volume, (currpts(end,[2 1 3]).*[4 4 1]-[0 0 zOffset]), keypts(index,1:3));
            path(:,3) = round((path(:,3)-1)/3) * 3;
            % keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
        
            currpts = [currpts; path(2:end,[2 1 3])./[4 4 1]+[0 0 zOffset] repmat([1 2.3],size(path,1)-1,1)];
            % currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]-[0 0 zOffset] 1 2.3];
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
    if isempty(currpts)
        keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
        currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]+[0 0 zOffset] 1 2.3];

        snakes = growAStarSnake(keypts,currpts,children(1),snakes,volume);
    elseif abs(((currpts(end,2).*4) - keypts(index,1)) / ((currpts(end,1).*4) - keypts(index,2))) <= 1
        path = astar_search3D(volume, (currpts(end,[2 1 3]).*[4 4 1]-[0 0 zOffset]), keypts(index,1:3));
        path(:,3) = round((path(:,3)-1)/3) * 3;
        % keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
    
        currpts = [currpts; path(2:end,[2 1 3])./[4 4 1]+[0 0 zOffset] repmat([1 2.3],size(path,1)-1,1)];
        % currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]-[0 0 zOffset] 1 2.3];

        snakes = growAStarSnake(keypts,currpts,children(1),snakes,volume);
    else
        if ~isempty(currpts)
            if abs(((currpts(end,2).*4) - keypts(index,1)) / ((currpts(end,1).*4) - keypts(index,2))) <= 1
                path = astar_search3D(volume, (currpts(end,[2 1 3]).*[4 4 1]-[0 0 zOffset]), keypts(index,1:3));
                path(:,3) = round((path(:,3)-1)/3) * 3;
                % keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
            
                currpts = [currpts; path(2:end,[2 1 3])./[4 4 1]+[0 0 zOffset] repmat([1 2.3],size(path,1)-1,1)];
                % currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]-[0 0 zOffset] 1 2.3];
            end

            if size(currpts,1) > 1
                tmp = slap2.gui.refstack.SnakeRoi;
                tmp.pointsXYZR = currpts;
                snakes = [snakes tmp];
            end
        end
        keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
        currpts = [keypts(index,[2 1 3])./[4 4 1]+[0 0 zOffset] 1 2.3];

        snakes = growAStarSnake(keypts,currpts,children(1),snakes,volume);
    end
else
    if isempty(currpts); return; end
    if abs(((currpts(end,2).*4) - keypts(index,1)) / ((currpts(end,1).*4) - keypts(index,2))) <= 1
        path = astar_search3D(volume, (currpts(end,[2 1 3]).*[4 4 1]-[0 0 zOffset]), keypts(index,1:3));
        path(:,3) = round((path(:,3)-1)/3) * 3;
        % keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
    
        currpts = [currpts; path(2:end,[2 1 3])./[4 4 1]+[0 0 zOffset] repmat([1 2.3],size(path,1)-1,1)];
        % currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]-[0 0 zOffset] 1 2.3];
    end

    if size(currpts,1) > 1
        tmp = slap2.gui.refstack.SnakeRoi;
        tmp.pointsXYZR = currpts;
        snakes = [snakes tmp];
    end
end

end

%%

function snakes = growSnake(keypts, currpts, index, snakes)
global zOffset
children = find(keypts(:,4) == index);

if isempty(currpts) && ((round(keypts(index,3)) < 1 || round(keypts(index,3)) > 50) ...
    || (round(keypts(index,2)) < 105 || round(keypts(index,2)) > 1160) ...
    || (round(keypts(index,1)) < 150 || round(keypts(index,1)) > 660))
    return
elseif ~isempty(currpts) && ((round(keypts(index,3)) < 1 || round(keypts(index,3)) > 50) ...
    || (round(keypts(index,2)) < 105 || round(keypts(index,2)) > 1160) ...
    || (round(keypts(index,1)) < 150 || round(keypts(index,1)) > 660))

    keypts(index,3) = round((keypts(index,3)-1)/3) * 3;

    currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]+[0 0 zOffset] 1 2.3];
    tmp = slap2.gui.refstack.SnakeRoi;
    tmp.pointsXYZR = currpts;
    snakes = [snakes tmp];
    return
end

if size(children,1) > 1
    if ~isempty(currpts)
        keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
        currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]+[0 0 zOffset] 1 2.3];
        tmp = slap2.gui.refstack.SnakeRoi;
        tmp.pointsXYZR = currpts;
        snakes = [snakes tmp];
    end
    for i = 1:size(children,1)
        snakes = growSnake(keypts,[],children(i),snakes);
    end
elseif size(children,1) == 1
    keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
    currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]+[0 0 zOffset] 1 2.3];
    snakes = growSnake(keypts,currpts,children(1),snakes);
else
    if isempty(currpts); return; end
    keypts(index,3) = round((keypts(index,3)-1)/3) * 3;
    currpts = [currpts; keypts(index,[2 1 3])./[4 4 1]+[0 0 zOffset] 1 2.3];
    tmp = slap2.gui.refstack.SnakeRoi;
    tmp.pointsXYZR = currpts;
    snakes = [snakes tmp];
end

end