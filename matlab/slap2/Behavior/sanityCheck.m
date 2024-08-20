testMouse = 'Y:\ophys\BCI\inactive_Mice\733585_SBCI24\Behavior\BCI'; mousePath = ls(testMouse);
epochPaths = arrayfun(@(i) [testMouse filesep mousePath(i,:)], 1:size(mousePath,1), 'UniformOutput',false);
epochTypes = {'moveOnly', 'BCIon   '};
temp1 = cellfun(@(x) (split(x,'_')), epochPaths, 'UniformOutput',false);
goodPaths = (cellfun(@(x) (x{end}), temp1, 'UniformOutput',false));
[~,moveOnlyLogical] = ismember(goodPaths, epochTypes{1});
[~,BCIlogical] = ismember(goodPaths, epochTypes{2});
epochs = horzcat(epochPaths(logical(moveOnlyLogical))', epochPaths(logical(BCIlogical))');
%sessionNum = epochs(1,:)
registers = {
    '\Behavior.harp\Register__PwmFrequencyDO1.bin', ...
    '\Behavior.harp\Register__OutputSet.bin',...
    '\Behavior.harp\Register__DigitalInputState.bin',...
    '\Behavior.harp\Register__AnalogData.bin',...
    '\Operation\SpoutPosition.csv',...
    '\LoadCells.harp\Register__LoadCellData.bin'
    }; %registers subject to change based on updates

startPin = 0x1; stopPin = 0x10; %Subject to change based on bonsai updates

for seshNum = 1:size(epochs,1)
    disp('---------------------------top')
    tic;
    [moveOnly, BCIon] = epochs{seshNum,:};
    BCIon = BCIon(find(~isspace(BCIon)));
    
    %MoveOnlys
    moveOnly_audio = read_harp_bin([moveOnly registers{1}]);
    moveOnly_handshake = read_harp_bin([moveOnly registers{2}]);
    moveOnly_licks = read_harp_bin([moveOnly registers{3}]);
    moveOnly_adc = read_harp_bin([moveOnly registers{4}]); moveOnly_adc.Values = moveOnly_adc.Values(:,1);
    moveOnly_motor = readtable([moveOnly registers{5}]);
    moveOnly_loadCell = read_harp_bin([moveOnly registers{6}]);  
    
    %BCIs
    BCI_audio = read_harp_bin([BCIon registers{1}]);
    BCI_handshake = read_harp_bin([BCIon registers{2}]);
    BCI_licks = read_harp_bin([BCIon registers{3}]);
    BCI_adc = read_harp_bin([BCIon registers{4}]); moveOnly_adc.Values = moveOnly_adc.Values(:,1);
    BCI_motor = readtable([BCIon registers{5}]);
    BCI_loadCell = read_harp_bin([BCIon registers{6}]);  
    toc
    disp('---------------------------bot (approx 30 seconds')
end
