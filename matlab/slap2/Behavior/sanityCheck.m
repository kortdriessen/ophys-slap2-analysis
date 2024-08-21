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

startPin = 0x1000; stopPin = 0x1; %Subject to change based on bonsai updates
goCue = 500; rewardCue = 100;
rewardDelivery = 1; licks = 2;
for seshNum = 1:size(epochs,1)
    disp('---------------------------top')
    tic;
    [moveOnly, BCIon] = epochs{seshNum,:};
    BCIon = BCIon(find(~isspace(BCIon)));
    
    %MoveOnlys
    moveOnly_audio = read_harp_bin([moveOnly registers{1}]);
    moveOnly_goCueT = moveOnly_audio.Seconds(moveOnly_audio.Values==goCue);
    moveOnly_rewCueT = moveOnly_audio.Seconds(moveOnly_audio.Values==rewardCue);

    moveOnly_handshake = read_harp_bin([moveOnly registers{2}]);
    moveOnly_startTriggers = moveOnly_handshake.Seconds(moveOnly_handshake.Values == startPin);
    moveOnly_stopTriggers = moveOnly_handshake.Seconds(moveOnly_handshake.Values == stopPin);
    moveOnly_stopTriggers = moveOnly_stopTriggers(2:end);%wierd instance of a stop trigger being sent before the first start trigger
    moveOnly_trialTimes = moveOnly_stopTriggers - moveOnly_startTriggers;


    moveOnly_licks = read_harp_bin([moveOnly registers{3}]);
    moveOnly_lickT = moveOnly_licks.Seconds(moveOnly_licks.Values==licks);

    moveOnly_adc = read_harp_bin([moveOnly registers{4}]); moveOnly_adc.Values = moveOnly_adc.Values(:,1);
    moveOnly_motor = readtable([moveOnly registers{5}]);
    moveOnly_loadCell = read_harp_bin([moveOnly registers{6}]);  
    
    moveOnly_array = trialize_data(moveOnly_startTriggers, moveOnly_stopTriggers, moveOnly_loadCell, moveOnly_adc, moveOnly_motor, moveOnly_lickT);


    % 
    % figure
    % xline(moveOnly_goCueT, 'Color', 'red', 'LineWidth',2)
    % hold on
    % xline(moveOnly_rewCueT, 'Color', 'blue', 'LineWidth',2)
    % xline(moveOnly_startTriggers, 'Color', 'black', 'LineWidth',2)
    % xline(moveOnly_stopTriggers, 'Color', 'green', 'LineWidth',2)
    % xline(moveOnly_licksT, 'Color', 'magenta', 'LineWidth',2)
    % plot(moveOnly_loadCell.Seconds, moveOnly_loadCell.Values(:,1))
    % plot(moveOnly_motor.Seconds, moveOnly_motor.Value*1400) %scaled for easier visualization
    % hold off
    % title('Move Only Alignment Sanity Check')
    % xlabel('Time (s)')
    % ylabel('LoadCell Activity')
    %some trials go over max 20 second limit because mouse movement before
    %go cue is encapsulated in trail times
    % figure
    % bar(moveOnly_trialTimes)
    % title('Check Trial Times for Sanity')
    % ylabel('Time (s)')
    % xlabel('Trial')
    % 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %BCIs
    BCI_audio = read_harp_bin([BCIon registers{1}]);
    BCI_goCueT = moveOnly_audio.Seconds(BCI_audio.Values==goCue);
    BCI_rewCueT = moveOnly_audio.Seconds(BCI_audio.Values==rewardCue);

    BCI_handshake = read_harp_bin([BCIon registers{2}]);
    BCI_startTriggers = BCI_handshake.Seconds(BCI_handshake.Values == startPin);
    BCI_stopTriggers = BCI_handshake.Seconds(BCI_handshake.Values == stopPin);
    BCI_stopTriggers = BCI_stopTriggers(2:end);%wierd instance of a stop trigger being sent before the first start trigger
    BCI_trialTimes = BCI_stopTriggers - BCI_startTriggers;

    BCI_licks = read_harp_bin([BCIon registers{3}]);
    BCI_lickT = BCI_licks.Seconds(BCI_licks.Values==licks);

    BCI_adc = read_harp_bin([BCIon registers{4}]); moveOnly_adc.Values = moveOnly_adc.Values(:,1);
    BCI_motor = readtable([BCIon registers{5}]);
    BCI_loadCell = read_harp_bin([BCIon registers{6}]);
    
    BCI_array = trialize_data(BCI_startTriggers, BCI_stopTriggers, BCI_loadCell, BCI_adc, BCI_motor, BCI_lickT);

    toc
    disp('---------------------------bot (approx 30 seconds')
end
behaviorData = vertcat(moveOnly_array, BCI_array);