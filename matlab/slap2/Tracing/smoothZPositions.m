function smoothedNodes = smoothZPositions(nodes, connectivity, distanceThreshold)
% smoothZPositions - Smooths the z-positions of nodes based on tree distance.
%
%   smoothedNodes = smoothZPositions(nodes, connectivity, distanceThreshold)
%
%   Description:
%     This function takes the node information and connectivity matrix
%     from the parseSWC function and smooths the z-positions of the nodes.
%     For each node, it performs a breadth-first search (BFS) to find all
%     nodes within a specified tree distance.  It then calculates the median
%     z-position of those nodes (including the node itself).
%     The z-position of the current node is updated to this median value.
%
%   Input:
%     nodes             - A struct array of nodes, as returned by parseSWC.
%     connectivity      - A sparse connectivity matrix, as returned by parseSWC.
%     distanceThreshold - The maximum tree distance to consider two nodes
%                         as neighbors (integer).  This is the number of edges
%                         traversed in the tree, not the Euclidean distance.
%
%   Output:
%     smoothedNodes     - A struct array of nodes with smoothed z-positions.
%
%   Example:
%     % Assuming you have already called parseSWC:
%     [nodes, connectivity] = parseSWC('your_neuron.swc');
%     distanceThreshold = 3; % Example tree distance threshold
%     smoothedNodes = smoothZPositions(nodes, connectivity, distanceThreshold);
%
%   Notes:
%     - The function uses breadth-first search (BFS) to traverse the tree.
%     - The distance threshold now represents the number of edges traversed
%       in the tree, not the Euclidean distance.
%     - The original nodes array is not modified; a new array with smoothed
%       z-positions is returned.
%     - If a node has no neighbors within the tree distance (including itself),
%       its z-position remains unchanged.
%
%   Error Handling:
%     - The function checks if the input arguments are valid.  If not, it
%       returns an empty array and displays an error message.
%
%   Author:  [Your Name or Institution]
%   Date:    [Date]
%
%   Revision History:
%     [Date] - Initial version.
%     [Date] - Modified to use tree distance and BFS.
%
    % Check input arguments for validity
    if ~isstruct(nodes) || ~ismatrix(connectivity) || ~isscalar(distanceThreshold) || ~isnumeric(distanceThreshold)
        error('smoothZPositions:InvalidInput', 'Invalid input arguments.  Check the types of nodes, connectivity, and distanceThreshold.');
    end

    % Check if distanceThreshold is an integer and >= 0
    if distanceThreshold < 0 || floor(distanceThreshold) ~= distanceThreshold
        error('smoothZPositions:InvalidInput', 'distanceThreshold must be a non-negative integer.');
    end
    
    % Get the number of nodes
    numNodes = length(nodes);

    % Preallocate the output array
    smoothedNodes = nodes; % Start with a copy of the original nodes
for i = numNodes:-1:1
    connectedNodes{i} = find(connectivity(i,:) | connectivity(:,i)'); 
end

    % Loop through each node
    for i = 1:numNodes
        % Perform breadth-first search to find neighbors within tree distance
        neighborZs = [];
        neighborCount = 0;
        
        queue = i; % Initialize queue with the starting node
        visited = false(1, numNodes); % Keep track of visited nodes
        visited(i) = true; % Mark the starting node as visited
        distance = zeros(1, numNodes); % Keep track of the distance from the starting node
        distance(i) = 0;
        
        while ~isempty(queue)
            currentNode = queue(1);
            queue(1) = []; % Dequeue
            
            if distance(currentNode) <= distanceThreshold
                neighborZs(neighborCount + 1) = nodes(currentNode).z;
                neighborCount = neighborCount + 1;
            end
            
            % Find neighbors of the current node
            cN = connectedNodes{currentNode};
            
            for j = 1:length(cN)
                neighbor = cN(j);
                if ~visited(neighbor)
                    visited(neighbor) = true; % Mark neighbor as visited
                    distance(neighbor) = distance(currentNode) + 1;
                    queue(end + 1) = neighbor; % Enqueue neighbor
                end
            end
        end
        
        % Calculate the median z-position of the neighbors
        if ~isempty(neighborZs)
            smoothedNodes(i).z = trimmean(neighborZs, 40);
        end
    end
end
