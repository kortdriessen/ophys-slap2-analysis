function [nodes, connectivity, voxelSeparation] = parseSWC(filepath)
% parseSWC - Parses an SWC file and returns node information and connectivity.
%
%   [nodes, connectivity] = parseSWC(filepath)
%
%   Description:
%     This function reads an SWC file, extracts the node information,
%     and determines the connectivity between the nodes.  SWC files are
%     commonly used to represent the morphology of neurons.
%
%   Input:
%     filepath - Path to the SWC file (string).
%
%   Output:
%     nodes      - A struct array where each element represents a node.
%                  Fields:
%                     .id    - Node ID (integer)
%                     .type  - Type of neuron structure (integer)
%                     .x     - X-coordinate (double)
%                     .y     - Y-coordinate (double)
%                     .z     - Z-coordinate (double)
%                     .r     - Radius (double)
%                     .parent- Parent node ID (integer, 0 for root)
%     connectivity - A logical sparse matrix representing the connectivity.
%                    connectivity(i, j) is true if node i is connected to node j.
%                    If node i is connected to node j, then connectivity(i,j) = 1
%
%   Example:
%     [nodeData, connMatrix] = parseSWC('my_neuron.swc');
%     disp(nodeData(1)); % Display the first node
%     spy(connMatrix);       % Visualize the connectivity matrix
%
%   Notes:
%     - The function handles comments in the SWC file (lines starting with '#').
%     - It assumes that node IDs are unique and positive.
%     - The connectivity matrix is sparse to efficiently handle large neurons.
%     - The function does some basic error checking, such as checking for file existence.
%
%   Error Handling:
%     - If the file does not exist, the function displays an error message
%       and returns empty arrays.
%     - If a line in the SWC file has an incorrect number of values, the
%       function displays an error and skips the line.
%
%   Author:  [Your Name or Institution]
%   Date:    [Date]
%
%   Revision History:
%     [Date] - Initial version.
%

    % Initialize output variables
    nodes = [];
    connectivity = [];
    voxelSeparation = []; % Ensure it is empty on error

    % Check if the file exists
    if ~exist(filepath, 'file')
        error('parseSWC:FileNotFound', 'File not found: %s', filepath);
        return;
    end

    try
        % Open the file for reading
        fileID = fopen(filepath, 'r');
        if fileID == -1
            error('parseSWC:FileOpenFailed', 'Could not open file: %s', filepath);
        end

        % Read the file line by line
        lineNum = 0; % Keep track of the line number for error reporting
        while ~feof(fileID)
            lineNum = lineNum + 1;
            line = fgetl(fileID);

            % Skip empty lines and comments
            if isempty(line)
                continue;
            end

            % Check for voxel separation comment
            if strncmpi(line, '# Voxel separation (x,y,z):', 26)
                try
                    % Extract the numeric values using regular expression
                    numbers = regexp(line, '(\d+\.\d+|\d+)', 'match');
                    if numel(numbers) == 3
                        voxelSeparation = str2double(numbers);
                    else
                         warning('parseSWC:IncorrectVoxelFormat', 'Incorrect voxel separation format in file %s at line %d.', filepath, lineNum);
                    end
                catch ME
                    warning('parseSWC:VoxelParsingError', 'Error parsing voxel separation in file %s at line %d. Error: %s', filepath, lineNum, ME.message);
                end
                continue; % Go to the next line after processing the comment
            end
            
            % Skip other comments
            if strncmp(line, '#', 1)
                continue;
            end

            % Parse the line
            values = sscanf(line, '%d %d %f %f %f %f %d');

            % Check if the line was parsed correctly
            if numel(values) ~= 7
                warning('parseSWC:IncorrectLineFormat', 'Incorrect line format in file %s at line %d. Skipping line.', filepath, lineNum);
                continue; % Skip to the next line
            end

            % Extract node information
            node.id     = values(1);
            node.type   = values(2);
            node.x      = values(3);
            node.y      = values(4);
            node.z      = values(5);
            node.r      = values(6);
            node.parent = values(7);

            % Add the node to the nodes array
            if isempty(nodes)
                nodes = node;
            else
                nodes(end+1) = node; % Append the node
            end
        end

        % Close the file
        fclose(fileID);
    catch ME
        % Handle any errors during file processing
        error('parseSWC:FileProcessingError', 'Error processing file: %s.  Error message: %s', filepath, ME.message);
    end

    % Construct the connectivity matrix
    numNodes = length(nodes);
    if numNodes > 0
        % Preallocate a sparse matrix.  We don't know the exact number of connections,
        % but we can estimate a maximum of 2*numNodes (each node has at most one parent and can be a parent to other nodes).
        % This helps with memory efficiency, especially for large datasets.
        I = zeros(2*numNodes, 1); % Initialize row indices
        J = zeros(2*numNodes, 1); % Initialize column indices
        V = ones(2*numNodes, 1);  % Initialize values (all 1 for connectivity)
        nnz = 0; % Keep track of the number of non-zero elements

        % Create a map from node ID to index for faster lookup
        nodeIdToIndex = containers.Map({nodes.id}, 1:numNodes);

        for i = 1:numNodes
            parentNodeId = nodes(i).parent;
            if parentNodeId > 0
                % Find the index of the parent node
                if isKey(nodeIdToIndex, parentNodeId)
                    parentIndex = nodeIdToIndex(parentNodeId);
                    % Add the connection to the sparse matrix
                    nnz = nnz + 1;
                    I(nnz) = i;          % from node i
                    J(nnz) = parentIndex;  % to its parent
                else
                   warning('parseSWC:MissingParent','Parent node %d of node %d not found.',parentNodeId,nodes(i).id);
                end
            end
        end
        % Create the sparse matrix.  Only create the matrix if nodes were found.
        connectivity = sparse(I(1:nnz), J(1:nnz), V(1:nnz), numNodes, numNodes);
    else
        connectivity = sparse(0,0); % Return an empty sparse matrix if there are no nodes
    end
end
