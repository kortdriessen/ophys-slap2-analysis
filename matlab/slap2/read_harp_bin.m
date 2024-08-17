function harpData = read_harp_bin(filePath)
%{
   Inputs:
        filePath (string) : path to register
   Outputs:
        harpData (table) : Matlab table of data on register (if data exists, else NaN)
%}
fileID = fopen(filePath, 'rb'); 
data = fread(fileID,  '*ubit8');

secondsPerTick = 32e-6;    
keySets = [1, 2, 4, 8, 129, 130, 132, 136, 68];
valueSets = {'uint8', 'uint16', 'uint32', 'uint64', 'int8', 'int16', 'int32', 'int64', 'single'};
payloadTypes = containers.Map(keySets, valueSets);
customValues = [1 2 4 8 1 2 4 8 4];
payloadSizes = containers.Map(valueSets, customValues);
payloadType = payloadTypes(  bitand(data(5),bitcmp(16, 'uint8'))  );
stride = data(2) + 2;
payloadLength = uint32(floor(size(data,1) / uint32(stride)));
payloadSize = stride - 12;
elementSize = payloadSizes(payloadType); %retrieving size (in bytes) of payload
payloadShape = [payloadLength floor(payloadSize / elementSize)];

if length(data) == 0
    harpData = NaN;
    disp('No Data on Register')
    exit;
end

%this is how you do that ndarray python jazz in matlab
data = reshape(data, stride, []).';

secondsOffset = 6; %python offset then +1
secondsPayload = payloadSizes('uint32');
secondsSlice =  data(:,secondsOffset:(secondsOffset+secondsPayload)-1); %seconds starts at 6th byte in stride, and it 4 bytes long
seconds = arrayfun(@(i) typecast(secondsSlice(i,:), 'uint32'), 1:size(secondsSlice,1), 'UniformOutput',false);
seconds = cell2mat(seconds);

ticksPayload = payloadSizes('uint16');
ticksOffset = 9;
ticksSlice = data(:, ticksOffset:(ticksOffset+ticksPayload)-1);
ticks = arrayfun(@(i) typecast(ticksSlice(i,:), 'uint16'), 1:size(ticksSlice,1), 'UniformOutput',false);
ticks = cell2mat(ticks);

seconds = uint32(ticks) * secondsPerTick + seconds;
%unfortunately this jazz has a slower tempo in matlab than it does in python

%extracting payload
payloadOffset = 12;
payloadSlice = data(:,   payloadOffset:elementSize:(payloadOffset+uint32(payloadSize)-1)); 
payload = arrayfun(@(i) int16(payloadSlice(i,:)), 1:size(payloadSlice,1), 'UniformOutput',false); %converting data as should be.... might be difficult if not int16 though....
payload = cell2mat(payload);
payload = reshape(payload, 3, []).';

%%%%% ET FINI
harpData = table(Seconds.', payload, 'VariableNames', {'Seconds', 'Vales'});


end