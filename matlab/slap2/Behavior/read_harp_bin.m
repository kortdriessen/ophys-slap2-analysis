function harpData = read_harp_bin(filePath)
%{
 21 byte width
loadCellData.' =
3   message type (1 byte)    
26  length (1 byte)  
33  address (1 byte)
255 port (1 byte)  
146 payload type (1 byte) 

depending if timestamp bit is set (here it is)

220   seconds (4 bytes) - uint32
207   -
3     -
0     -
165   microseconds (2 bytes) - uint 16
61    -
127   data (n-bytes dependent on payloadType)
47    -
137   -
122   -
139   -
122   -
115   -
122   - 4 load cell pins, each 2 bytes
0     
0     
0     - other data not used on register for this example
0     
0     
0     
0
0
35     checksum to verify integrity (continuously increasing as more dat is passed through read

   Inputs:
        filePath (string) : path to register
        loadcell (boolean): true if handling loadcell register (Loadcells.harp only)
   Outputs:
        harpData (table) : Matlab table of data on register (if data exists, else NaN)
%}

fileID = fopen(filePath, 'rb'); 
data = fread(fileID,  '*ubit8');
if data(1) == 3
    disp('Reading Event Data From Register')

       
    keySets = [1, 2, 4, 8, 129, 130, 132, 136, 68];
    valueSets = {'uint8', 'uint16', 'uint32', 'uint64', 'int8', 'int16', 'int32', 'int64', 'single'};
    payloadTypes = containers.Map(keySets, valueSets);
    customValues = [1 2 4 8 1 2 4 8 4];
    payloadSizes = containers.Map(valueSets, customValues);
    
    deviceType = data(4); %points to what register data is on, if 0xFF, its a standalone device, ex: LoadCell device
    
    %breaking down for easier readability
    payloadInfo = dec2bin(data(5));
    isSigned = payloadInfo(1);
    isFloat = (payloadInfo(2));
    hasTimeStamp = (payloadInfo(4));
    payloadIndicator = (payloadInfo(5:end));
    
    payloadType = payloadTypes(  bitand(bin2dec([isSigned isFloat '0' hasTimeStamp payloadIndicator]),bitcmp(16, 'uint8'))  );
    elementSize = payloadSizes(payloadType); %retrieving size (in bytes) of payload1 element of the payload
    useTimeStamps = logical(hasTimeStamp);
    
    stride = data(2) + 2;
    payloadLength = uint32(floor(size(data,1) / uint32(stride)));
    payloadSize = stride - 12;
    payloadShape = [payloadLength floor(payloadSize / elementSize)];
    data = reshape(data, stride, []).';
    
    
    if useTimeStamps == true
        disp('seconds bytes (4) + microseconds (2) + data ')
        payloadData = data(:, 6:end-1); %ignoring headers
        %Seconds are 4 bytes
        secondsPayload = payloadSizes('uint32');
        secondsSlice =  payloadData(:,1:(secondsPayload)); %seconds starts at 6th byte in stride, and it 4 bytes long
        seconds = arrayfun(@(i) typecast(secondsSlice(i,:), 'uint32'), 1:size(secondsSlice,1), 'UniformOutput',false);
        seconds = cell2mat(seconds);
        
        %Microseconds are 2 bytes
        ticksPayload = payloadSizes('uint16');
        ticksSlice = payloadData(:, secondsPayload+1:secondsPayload+ticksPayload);
        ticks = arrayfun(@(i) typecast(ticksSlice(i,:), 'uint16'), 1:size(ticksSlice,1), 'UniformOutput',false);
        ticks = cell2mat(ticks);
        
        seconds = double(double(ticks)*10e-6) + double(seconds);
        
        dataValuesSlice = payloadData(:, (secondsPayload + ticksPayload+1 ):end);
        dataValues =  arrayfun(@(i) typecast(dataValuesSlice(i,:), payloadType), 1:size(payloadData,1), 'UniformOutput', false);
        try
            pinsRead = numel(dataValues{1});
            dataValues = reshape(cell2mat(dataValues), pinsRead, []).';
        catch 
            dataValues = reshape(cell2mat(dataValues), 1, []).';
        end
    
        harpData = table(seconds.', dataValues, 'VariableNames', {'Seconds', 'Values'});
    else
        disp('just data')
        payloadData = data(:, 6:(6+payloadSize));
    
        dataValuesSlice = payloadData(:, 1:end);
        dataValues =  arrayfun(@(i) typecast(dataValuesSlice(i,:), payloadType), 1:size(payloadData,1), 'UniformOutput', false);
        try
            pinsRead = numel(dataValues{1});
            dataValues = reshape(cell2mat(dataValues), pinsRead, []).';
        catch 
            dataValues = reshape(cell2mat(dataValues), 1, []).';
        end
        sequenceOfValues = 1:size(dataValues,1);
        harpData = table(sequenceOfValues.', dataValues, 'VariableNames', {'Sequence', 'Values'});
    end
else
    disp('Reading Harp Reads and Writes')
    if data(1) == 1
        %First line used to reshape only since its a read and not a write
        keySets = [1, 2, 4, 8, 129, 130, 132, 136, 68];
        valueSets = {'uint8', 'uint16', 'uint32', 'uint64', 'int8', 'int16', 'int32', 'int64', 'single'};
        payloadTypes = containers.Map(keySets, valueSets);
        customValues = [1 2 4 8 1 2 4 8 4];
        payloadSizes = containers.Map(valueSets, customValues);


        dataLength = data(2);

        %breaking down for easier readability
        payloadInfo = dec2bin(data(5));
        isSigned = payloadInfo(1);
        isFloat = (payloadInfo(2));
        hasTimeStamp = (payloadInfo(4));
        payloadIndicator = (payloadInfo(5:end));
        
        payloadType = payloadTypes(  bitand(bin2dec([isSigned isFloat '0' hasTimeStamp payloadIndicator]),bitcmp(16, 'uint8'))  );
        elementSize = payloadSizes(payloadType); %retrieving size (in bytes) of payload1 element of the payload
        useTimeStamps = logical(hasTimeStamp);

        stride = data(2) + 2;
        payloadLength = uint32(floor(size(data,1) / uint32(stride)));
        payloadSize = stride - 12; %for some reason this causes a bug to occur when reading reads or writes, and does not occur in events (ex: \Behavior.harp\Register__PwmFrequencyDO1.bin)
        payloadShape = [payloadLength floor(payloadSize / elementSize)];
        data = reshape(data, stride, []).';
        payloadData = data(:, 6:end);
        if useTimeStamps == true
            secondsPayload = payloadSizes('uint32');
            secondsSlice =  payloadData(:,1:(secondsPayload)); %seconds starts at 6th byte in stride, and it 4 bytes long
            seconds = arrayfun(@(i) typecast(secondsSlice(i,:), 'uint32'), 1:size(secondsSlice,1), 'UniformOutput',false);
            seconds = cell2mat(seconds);
            
            %Microseconds are 2 bytes
            ticksPayload = payloadSizes('uint16');
            ticksSlice = payloadData(:, secondsPayload+1:secondsPayload+ticksPayload);
            ticks = arrayfun(@(i) typecast(ticksSlice(i,:), 'uint16'), 1:size(ticksSlice,1), 'UniformOutput',false);
            ticks = cell2mat(ticks);
            
            seconds = double(double(ticks)*10e-6) + double(seconds);
            
            dataValuesSlice = payloadData(:, secondsPayload+ticksPayload+1: end-1); %instead of using payload size, i am going to just omit the checksum
            dataValues = arrayfun(@(i) typecast(dataValuesSlice(i,:), payloadType), 1:size(payloadData,1), 'UniformOutput', false);
            try
                pinsRead = numel(dataValues{1});
                dataValues = reshape(cell2mat(dataValues), pinsRead, []).';
            catch 
                dataValues = reshape(cell2mat(dataValues), 1, []).';
            end

            harpData = table(seconds.', dataValues, 'VariableNames', {'Seconds', 'Values'});
        else
            dataValuesSlice = payloadData(:, 6+1: end-1); %instead of using payload size, i am going to just omit the checksum
            dataValues = arrayfun(@(i) typecast(dataValuesSlice(i,:), payloadType), 1:size(payloadData,1), 'UniformOutput', false);
            try
                pinsRead = numel(dataValues{1});
                dataValues = reshape(cell2mat(dataValues), pinsRead, []).';
            catch 
                dataValues = reshape(cell2mat(dataValues), 1, []).';
            end

            sequentialValues = 1:size(dataValues,1);
            harpData = table(sequentialValues.', dataValues, 'VariableNames', {'Sequence', 'Values'});
        end


    end
end



end