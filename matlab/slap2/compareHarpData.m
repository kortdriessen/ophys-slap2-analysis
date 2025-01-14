%pick mice folders found in inactive mice directory on scratch
% dirs = uipickfiles();


%%

numMice = length(dirs);
for mouseNum = 1:numMice
    mousePath = dirs{mouseNum};

    %slap2 isolation
    slap2Paths = struct2cell(dir([mousePath filesep 'SLAP2\slap2*']));
    sessionDatesSLAP2 = cellfun(@(x) x(strfind(x, '2024'):strfind(x, '2024')+9), slap2Paths(1,:), 'UniformOutput', false);
    sessionDates = cellfun(@(x) erase(x(strfind(x, '2024'):strfind(x, '2024')+9), '-'), slap2Paths(1,:), 'UniformOutput', false);
    slap2SessionNames = slap2Paths(1,:); slap2SessionFolders = slap2Paths(2,:);
    slap2NeuronPaths = cellfun( @(x) dir(x), strcat(slap2SessionFolders, filesep, slap2SessionNames, filesep, '*Neuron*'), 'UniformOutput', false);
    slap2exptSummaryPaths = cellfun( @(x) dir([x.folder filesep x.name filesep 'ExperimentSummary' filesep 'Summary*']), slap2NeuronPaths, 'UniformOutput', false);
    summaryFilePaths = cellfun( @(x) [x.folder filesep x.name ], slap2exptSummaryPaths, 'UniformOutput', false);
    
    %behavior isolation
    bciBehaviorPaths = struct2cell(dir([mousePath filesep 'behavior\BCI\2024*']));
    epochNames = bciBehaviorPaths(1,:);
    epochFolders = bciBehaviorPaths(2,:);
    epochPaths = strcat(epochFolders, filesep, epochNames, filesep, 'harpData.mat');
    
    %Analyze Quality Sessions
    qualitySessions = length(sessionDates);
    for qSession = 1:qualitySessions
        %Get both epochs from session behavior
        thisSessionStr = sessionDates{qSession};
        epochLocals = cellfun(@(x) contains(x, thisSessionStr), epochPaths, 'UniformOutput',true);
        epochs = epochPaths(epochLocals);
        epoch1 = epochs{1};
        epoch2 = epochs{2};

        %get epochs from slap2
        thisSessionStr = sessionDatesSLAP2{qSession};
        summaryLocals = cellfun(@(x) contains(x, thisSessionStr), summaryFilePaths, 'UniformOutput',true);
        summaryFile = summaryFilePaths{summaryLocals};

        %LOAD IN DATA!
        % harpE1 = load(epoch1);
        % harpE2 = load(epoch2);
        % load(summaryFile, 'exptSummary')

    end
end
