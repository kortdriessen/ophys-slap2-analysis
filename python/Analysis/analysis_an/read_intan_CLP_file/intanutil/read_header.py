import sys, struct

def read_header(fid):
    # Read magic number to ensure this is an Intan CLP file
    magic_number, = struct.unpack('<I', fid.read(4)) # Read magic number as a 32-bit unsigned int
    if magic_number != int('f3b1a481', 16): raise Exception('Unrecognized file type.')

    # Create dictionary to contain header information
    header = {}

    # Read version number.
    version = {}
    (version['major'], version['minor']) = struct.unpack('<hh', fid.read(4)) # Read major and minor versions each as 16-bit signed int
    header['version'] = version

    # Determine if this is a Data File (Clamp and Measured Data) or Auxiliary Data File (Digital Ins/Outs and ADC)
    datatype, = struct.unpack('<h', fid.read(2)) # Read datatype as 16-bit signed int

    if datatype == 0:
        # If this is a Data File, report it
        print('')
        print('Reading Intan Technologies CLAMP Data File, Version {}.{}'.format(version['major'], version['minor']))
        print('')

        # Determine the number of bytes in the header
        NumBytes, = struct.unpack('<H', fid.read(2)) # Read NumBytes as 16-bit unsigned int

        # Populate the Date dictionary with Year, Month, Day, Hour, Minute and Second
        Date = {}
        (Date['Year'], Date['Month'], Date['Day'], Date['Hour'], Date['Minute'], Date['Second']) = struct.unpack('<hhhhhh', fid.read(12)) # Read each entry as 16-bit signed int
        header['Date'] = Date

        # Read the 'Chips' portion of the header
        header['Chips'] = read_header_chips(fid)

        # Read the 'Settings' portion of the header
        header['Settings'] = read_header_settings(fid)

        

    elif datatype == 1:
        # If this is an Auxiliary Data File, report it
        print('')
        print('Reading Intan Technologies CLAMP Auxiliary Data File, Version {}.{}'.format(version['major'], version['minor']))
        print('')

        # Determine the number of ADCs
        header['NumADCs'], = struct.unpack('<H', fid.read(2)) # Read NumADCs as 16-bit unsigned int
        # Determine the number of bytes in the header
        NumBytes, = struct.unpack('<H', fid.read(2)) # Read NumBytes as 16-bit unsigned int

        # Populate the Date dictionary with Year, Month, Day, Hour, Minute and Second
        Date = {}
        (Date['Year'], Date['Month'], Date['Day'], Date['Hour'], Date['Minute'], Date['Second']) = struct.unpack('<hhhhhh', fid.read(12)) # Read each entry as 16-bit signed int
        header['Date'] = Date

        # Read the sampling rate into the 'Settings' portion of the header
        header['Settings'] = {}
        header['Settings']['SamplingRate'], = struct.unpack('<f', fid.read(4)) # Read sample rate as single-precision float

    else:
        # If the data type isn't recognized, exit
        print('Unrecognized data type')
        return ({}, 0, 0)

    return (header, datatype, NumBytes)


def read_header_settings(fid):
    # Create dictionary to contain settings information
    settings = {}

    settings['EnableCapacitiveCompensation'], = struct.unpack('<B', fid.read(1)) # Read EnableCapacitiveCompensation as 8-bit unsigned int
    settings['CapCompensationMagnitude'], = struct.unpack('<f', fid.read(4)) # Read CapCompensationMagnitude as single-precision float

    my_filter, = struct.unpack('<f', fid.read(4)) # Read FilterCutoff as single-precision float
    # Only populate settings['FilterCutoff'] if the value is greater than 0
    if my_filter > 0:
        settings['FilterCutoff'] = my_filter

    settings['PipetteOffset'], = struct.unpack('<f', fid.read(4)) # Read PipetteOffset as single-precision float
    settings['SamplingRate'], = struct.unpack('<f', fid.read(4)) # Read SamplingRate as single-precision float

    # Create dictionary to contain CellParameters information
    CellParameters = {}
    CellParameters['Rs'], = struct.unpack('<f', fid.read(4)) # Read Rs as single-precision float
    CellParameters['Rm'], = struct.unpack('<f', fid.read(4)) # Read Rm as single-precision float
    CellParameters['Cm'], = struct.unpack('<f', fid.read(4)) # Read Cm as single-precision float

    # Enter CellParameters dictionary to settings dictionary
    settings['CellParameters'] = CellParameters

    settings['IsVoltageClamp'], = struct.unpack('<B', fid.read(1)) # Read IsVoltageClamp as 8-bit unsigned int
    settings['vClampX2mode'], = struct.unpack('<B', fid.read(1)) # Read vClampX2mode as 8-bit unsigned int

    if settings['IsVoltageClamp'] == 1:
        # If this data file is voltage clamp, read the voltage clamp settings
        settings['VoltageClamp'] = read_header_voltage_clamp_settings(fid)
    else:
        # If this data file is current clamp, read the current clamp settings
        settings['CurrentClamp'] = read_header_current_clamp_settings(fid)

    # Read the waveform information
    settings['Waveform'] = read_waveform(fid)

    return settings


def read_waveform(fid):
    # Create dictionary to contain Waveform information
    Waveform = {}

    Waveform['Interval'], = struct.unpack('<f', fid.read(4)) # Read Interval as single-precision float
    numSegments, = struct.unpack('<H', fid.read(2)) # Read numSegments as 16-bit unsigned int

    # Read and add Segments to the Waveform dictionary
    Waveform['Segments'] = []
    for i in range(numSegments):
        Waveform['Segments'].append(read_segment(fid))

    return Waveform


def read_segment(fid):
    # Create dictionary to contain Segment information
    Segment = {}

    Segment['WaveformNumber'], = struct.unpack('<B', fid.read(1)) # Read WaveformNumber as 8-bit unsigned int

    Segment['TOffset'], = struct.unpack('<I', fid.read(4))  # Read TOffset as 32-bit unsigned int, and increment by 1
    Segment['TOffset'] += 1

    Segment['Start'], = struct.unpack('<I', fid.read(4)) # Read Start as 32-bit unsigned int, and increment by 1
    Segment['Start'] += 1

    Segment['End'], = struct.unpack('<I', fid.read(4)) # Read End as 32-bit unsigned int, and increment by 1
    Segment['End'] += 1

    Segment['AppliedValue'], = struct.unpack('<f', fid.read(4)) # Read AppliedValue as single-precision float

    return Segment


def read_header_voltage_clamp_settings(fid):
    # Create dictionary to contain settings for voltage clamp
    vclampsettings = {}

    vclampsettings['HoldingVoltage'], = struct.unpack('<f', fid.read(4)) # Read HoldingVoltage as single-precision float
    vclampsettings['NominalResistance'], = struct.unpack('<f', fid.read(4)) # Read NominalResistance as single-precision float
    vclampsettings['Resistance'], = struct.unpack('<f', fid.read(4)) # Read Resistance as single-precision float
    vclampsettings['DesiredBandwidth'], = struct.unpack('<f', fid.read(4)) # Read DesiredBandwidth as single-precision float
    vclampsettings['ActualBandwidth'], = struct.unpack('<f', fid.read(4)) # Read ActualBandwidth as single-precision float

    return vclampsettings


def read_header_current_clamp_settings(fid):
    # Create dictionary to contain settings for current clamp
    iclampsettings = {}

    iclampsettings['HoldingCurrent'], = struct.unpack('<f', fid.read(4)) # Read HoldingCurrent as single-precision float
    iclampsettings['CurrentStepSize'], = struct.unpack('<f', fid.read(4)) # Read CurrentStepSize as single-precision float

    return iclampsettings


def read_header_chips(fid):
    numChips, = struct.unpack('<H', fid.read(2)) # Read numChips as 16-bit unsigned int
    numChannels, = struct.unpack('<H', fid.read(2)) # Read numChannels as 16-bit unsigned int

    # Read information for each chip into chips dictionary
    chips = []
    for i in range(numChips):
        chips.append(read_one_header_chip(fid, numChannels))

    return chips


def read_one_header_chip(fid, numChannels):
    # Create dictionary to contain settings for a chip
    chip = {}
    
    # Read information for each channel into chips dictionary
    chip['Channels'] = []
    for i in range(numChannels):
        chip['Channels'].append(read_one_header_channel(fid))

    chip['ChipRegisters'] = list(struct.unpack('<HHHH', fid.read(8))) # Read ChipRegisters as 4 16-bit unsigned ints

    return chip


def read_one_header_channel(fid):
    # Create dictionary to contain settings for a channel
    channel = {}

    # Read each register's data into Registers list
    Registers = []
    for i in range(14):
        value, = struct.unpack('<H', fid.read(2)) # Read value as 16-bit unsigned int
        Registers.append(value)
    channel['Registers'] = Registers

    channel['DifferenceAmpResidual'], = struct.unpack('<i', fid.read(4)) # Read DifferenceAmpResidual as 32-bit signed int
    channel['VoltageAmpResidual'], = struct.unpack('<i', fid.read(4)) # Read VoltageAmpResidual as 32-bit signed int
    bestCalibration1D = list(struct.unpack('<BBBBBBBBBBBBBBBB', fid.read(16))) # Read calibration 1-D vector as 16 8-bit unsigned ints

    bestCalibration = [[[bestCalibration1D[8*k + 2*j + i] for k in range(2)] for j in range(4)] for i in range(2)] # Reshape the 1D vector to a 3D list to allow for easier parsing
    channel['CoarseCalibration'] = [[bestCalibration[0][j][k] for k in range(2)] for j in range(4)] # Enter CoarseCalibration bytes into channel dictionary
    channel['FineCalibration'] = [[bestCalibration[1][j][k] for k in range(2)] for j in range(4)] # Enter FineCalibration bytes into channel dictionary
    channel['FeedbackResistors'] = list(struct.unpack('<fffff', fid.read(20))) # Enter FeedbackResistors bytes into channel dictionary
    channel['DesiredBandwidth'], = struct.unpack('<f', fid.read(4)) # Enter DesiredBandwidth bytes into channel dictionary

    return channel
