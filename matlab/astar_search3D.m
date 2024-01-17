function path = astar_search3D(image, startNode, goalNode)

    % Ensure startNode and endNode are in the correct order
    minCoord = max([1 1 1],min(startNode, goalNode) - [20 20 5]);
    maxCoord = min(size(image),max(startNode, goalNode) + [20 20 5]);
    
    % Crop the image
    image = image(minCoord(1):maxCoord(1), minCoord(2):maxCoord(2), minCoord(3):maxCoord(3));

    image = 255 * (image - min(image(:))) / (max(image(:))-min(image(:)));

    imSize = size(image);

    startNode = startNode - minCoord + 1;
    goalNode = goalNode - minCoord + 1;
    
    direc = sign(goalNode(2) - startNode(2));

    startNode = (startNode(3) - 1) * (imSize(1) * imSize(2)) + (startNode(2) - 1) * imSize(1) + startNode(1);
    goalNode = (goalNode(3) - 1) * (imSize(1) * imSize(2)) + (goalNode(2) - 1) * imSize(1) + goalNode(1);

    % Initialize open and closed lists
    openList = [];
    closedList = [];
    % openList = containers.Map('KeyType', 'char', 'ValueType', 'any');
    % closedList = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Initialize start node with F, G, and H costs
    % startKey = pos2key(startNode);
    % openList(startKey) = struct('Position', startNode, 'Parent', [], ...
    %                             'G', 0, 'H', heuristic(startNode, goalNode), ...
    %                             'F', heuristic(startNode, goalNode));
    
    global POSITION PARENT G_IDX H_IDX F_IDX

    POSITION = 1;
    PARENT = 2;
    G_IDX = 3;
    H_IDX = 4;
    F_IDX = 5;

    neighborDirections = [0, direc, 0;
                          1, direc, 0;
                         -1, direc, 0;
                          0, direc, 1;
                          1, direc, 1;
                         -1, direc, 1;
                          0, direc, -1;
                          1, direc,-1;
                         -1, direc, -1];

    openList = [openList; startNode 0 0 heuristic(startNode, goalNode, imSize) heuristic(startNode, goalNode, imSize)];

    
    % A* Algorithm
    while ~isempty(openList)
        % Find the node with the lowest F cost
        currentNodeIdx = minCostNode(openList);
        currentNode = openList(currentNodeIdx,:);

        if any(openList(:,POSITION) < 1)
            disp('break');
        end
        
        % Check if goal is reached
        if currentNode(POSITION) == goalNode
            pathInds = reconstructPath(currentNode, closedList);
            % [r,c,d] = ind2sub(imSize,pathInds);

            d = ceil(pathInds / (imSize(1) * imSize(2)));
            tmp = mod(pathInds - 1, imSize(1) * imSize(2));
            c = floor(tmp / imSize(1)) + 1;
            r = mod(tmp, imSize(1)) + 1;

            path = [r;c;d]';
            path = path + minCoord - 1;
            return;
        end
        
        % Move current node from open to closed list
        openList(currentNodeIdx,:) = [];
        closedList = [closedList; currentNode];
        
        % Explore neighbors
        for neighborPos = neighbors(currentNode(POSITION), imSize, neighborDirections)'

            % neighborKey = pos2key(neighborPos);
            if sum(closedList(:,POSITION) == neighborPos) > 0 %isKey(closedList, neighborKey)
                continue; % Skip already processed neighbors
            end
        
            % G cost for neighbor
            tentativeG = currentNode(G_IDX) + movementCost(currentNode(POSITION), neighborPos, image);
        
            if sum(openList(:,POSITION) == neighborPos) == 0 %~isKey(openList, neighborKey)
                % neighborNode = struct('Position', neighborPos, 'Parent', currentNodeKey, ...
                %                       'G', tentativeG, 'H', heuristic(neighborPos, goalNode) * 1/255, ...
                %                       'F', 0); % F will be updated below
                % openList(neighborKey) = neighborNode;

                neighborNode = [neighborPos currentNode(POSITION) tentativeG heuristic(neighborPos, goalNode, imSize) * 1/255 0];
                % openList = [openList; neighborNode];
        
                % Update neighbor's cost and parent
                neighborNode(G_IDX) = tentativeG;
                neighborNode(F_IDX) = tentativeG + neighborNode(H_IDX);
                neighborNode(PARENT) = currentNode(POSITION);
                % openList(neighborKey) = neighborNode; % Update the map with the modified struct
                openList = [openList; double(neighborNode)];
            else
                % neighborNode = openList(neighborKey);
                neighborNode = openList(openList(:,POSITION) == neighborPos,:);
                if tentativeG >= neighborNode(G_IDX)
                    continue; % Not a better path
                end

                % Update neighbor's cost and parent
                neighborNode(G_IDX) = tentativeG;
                neighborNode(F_IDX) = tentativeG + neighborNode(H_IDX);
                neighborNode(PARENT) = currentNode(POSITION);
                % openList(neighborKey) = neighborNode; % Update the map with the modified struct
                openList(openList(:,POSITION) == neighborPos,:) = neighborNode;
            end
        end
    end
    
    % If the goal was not reached
    path = [];
end

% function key = pos2key(pos)
%     key = sprintf('%d-%d-%d', pos(1), pos(2), pos(3));
% end

function cost = heuristic(pos1, pos2, imageSize)
    % Euclidean distance as the heuristic
    d1 = ceil(pos1 / (imageSize(1) * imageSize(2)));
    tmp = mod(pos1 - 1, imageSize(1) * imageSize(2));
    c1 = floor(tmp / imageSize(1)) + 1;
    r1 = mod(tmp, imageSize(1)) + 1;
    
    d2 = ceil(pos2 / (imageSize(1) * imageSize(2)));
    tmp = mod(pos2 - 1, imageSize(1) * imageSize(2));
    c2 = floor(tmp / imageSize(1)) + 1;
    r2 = mod(tmp, imageSize(1)) + 1;

    cost = sqrt(sum(([r1,c1,d1] - [r2,c2,d2]) .^ 2));
end

function cost = movementCost(pos1, pos2, image)
    % Cost inversely related to brightness (brighter pixels have lower cost)
    imageSize = size(image);

    d1 = ceil(pos1 / (imageSize(1) * imageSize(2)));
    tmp = mod(pos1 - 1, imageSize(1) * imageSize(2));
    c1 = floor(tmp / imageSize(1)) + 1;
    r1 = mod(tmp, imageSize(1)) + 1;
    
    d2 = ceil(pos2 / (imageSize(1) * imageSize(2)));
    tmp = mod(pos2 - 1, imageSize(1) * imageSize(2));
    c2 = floor(tmp / imageSize(1)) + 1;
    r2 = mod(tmp, imageSize(1)) + 1;

    cost = sqrt(sum((([r1,c1,d1] - [r2,c2,d2]).*[1,1,4]) .^ 2)) / (image(r2, c2, d2)+1e-3);
end

function nPos = neighbors(pos, imageSize, directions)
    d = ceil(pos / (imageSize(1) * imageSize(2)));
    tmp = mod(pos - 1, imageSize(1) * imageSize(2));
    c = floor(tmp / imageSize(1)) + 1;
    r = mod(tmp, imageSize(1)) + 1;

    posInds = [r, c, d];
    nPosInds = [];
    for i = 1:size(directions, 1)
        newPos = posInds + directions(i,:);
        if all(newPos > 0 & newPos <= imageSize)
            nPosInds = [nPosInds; newPos];
        end
    end
    if ~isempty(nPosInds)
        nPos = (nPosInds(:,3) - 1) * (imageSize(1) * imageSize(2)) + (nPosInds(:,2) - 1) * imageSize(1) + nPosInds(:,1);
    else
        nPos = [];
    end
end

function path = reconstructPath(currentNode, nodeList)
    % Reconstruct the path from the goal to the start
    global POSITION PARENT

    path = currentNode(POSITION);
    while currentNode(PARENT) > 0
        currentNode = nodeList(nodeList(:,POSITION) == currentNode(PARENT),:);
        path = [currentNode(POSITION), path];
    end
end

function nodeKey = minCostNode(nodeList)
    % Find the node with the minimum F cost in the open list
    global F_IDX

    [~, nodeKey] = min(nodeList(:,F_IDX));

    % minCost = inf;
    % nodeKey = '';
    % keys = nodeList.keys;
    % for i = 1:length(keys)
    %     node = nodeList(keys{i});
    %     if node.F < minCost
    %         minCost = node.F;
    %         nodeKey = keys{i};
    %     end
    % end
end
