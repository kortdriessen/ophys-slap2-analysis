function path = astar_search(image, startNode, goalNode)

    image = 255 * (image - min(image(:))) / (max(image(:))-min(image(:)));
    
    direc = sign(goalNode(2) - startNode(2));

    % Initialize open and closed lists
    openList = containers.Map('KeyType', 'char', 'ValueType', 'any');
    closedList = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Initialize start node with F, G, and H costs
    startKey = pos2key(startNode);
    openList(startKey) = struct('Position', startNode, 'Parent', [], ...
                                'G', 0, 'H', heuristic(startNode, goalNode), ...
                                'F', heuristic(startNode, goalNode));
    
    % A* Algorithm
    while ~isempty(openList)
        % Find the node with the lowest F cost
        currentNodeKey = minCostNode(openList);
        currentNode = openList(currentNodeKey);
        
        % Check if goal is reached
        if all(currentNode.Position == goalNode)
            path = reconstructPath(currentNode, closedList);
            return;
        end
        
        % Move current node from open to closed list
        openList.remove(currentNodeKey);
        closedList(currentNodeKey) = currentNode;
        
        % Explore neighbors
        for neighborPos = neighbors(currentNode.Position, size(image), direc)
            neighborKey = pos2key(neighborPos);
            if isKey(closedList, neighborKey)
                continue; % Skip already processed neighbors
            end
        
            % G cost for neighbor
            tentativeG = currentNode.G + movementCost(currentNode.Position, neighborPos, image);
        
            if ~isKey(openList, neighborKey)
                neighborNode = struct('Position', neighborPos, 'Parent', currentNodeKey, ...
                                      'G', tentativeG, 'H', heuristic(neighborPos, goalNode) * 1/255, ...
                                      'F', 0); % F will be updated below
                openList(neighborKey) = neighborNode;
            else
                neighborNode = openList(neighborKey);
                if tentativeG >= neighborNode.G
                    continue; % Not a better path
                end
            end
        
            % Update neighbor's cost and parent
            neighborNode.G = tentativeG;
            neighborNode.F = tentativeG + neighborNode.H;
            neighborNode.Parent = currentNodeKey;
            openList(neighborKey) = neighborNode; % Update the map with the modified struct
        end
    end
    
    % If the goal was not reached
    path = [];
end

function key = pos2key(pos)
    key = sprintf('%d-%d', pos(1), pos(2));
end

function cost = heuristic(pos1, pos2)
    % Euclidean distance as the heuristic
    cost = sqrt(sum((pos1 - pos2) .^ 2));
end

function cost = movementCost(pos1, pos2, image)
    % Cost inversely related to brightness (brighter pixels have lower cost)
    cost = sqrt(sum((pos1 - pos2) .^ 2)) / (image(pos2(1), pos2(2))+1e-3);
end

function nPos = neighbors(pos, imageSize, direc)
    % Get valid neighbor positions (4-connectivity)
    directions = [0, direc;
                  1, direc;
                 -1, direc]';
    nPos = [];
    for i = 1:size(directions, 2)
        newPos = pos + directions(:,i);
        if all(newPos > 0 & newPos <= imageSize)
            nPos = [nPos, newPos];
        end
    end
end

function path = reconstructPath(currentNode, nodeList)
    % Reconstruct the path from the goal to the start
    path = currentNode.Position;
    while ~isempty(currentNode.Parent)
        currentNode = nodeList(currentNode.Parent);
        path = [currentNode.Position, path];
    end
end

function nodeKey = minCostNode(nodeList)
    % Find the node with the minimum F cost in the open list
    minCost = inf;
    nodeKey = '';
    keys = nodeList.keys;
    for i = 1:length(keys)
        node = nodeList(keys{i});
        if node.F < minCost
            minCost = node.F;
            nodeKey = keys{i};
        end
    end
end
