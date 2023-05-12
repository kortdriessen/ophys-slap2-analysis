import sys, struct
import sys

def read_data(data, fid):
    # Determine the size of the input file
    start_pos = fid.tell()
    fid.seek(0, 2)
    file_size = fid.tell()

    # Return the fid object to its previous position
    fid.seek(start_pos, 0)
    
    # Read one batch of data into 'Data'
    data['Data'] = read_one_batch(fid, file_size - start_pos, data['Header']['Settings']['SamplingRate'])

    return data


def read_one_batch(fid, size, sampling_rate):
    # Create dictionary to contain Data
    Data = {}
    size_of_one = 4 + 4 + 4 + 4 # Timestep + Applied (Software) + Clamp Value + Measured
    length = int(size / size_of_one)

    # Read in the whole thing at once
    numBytesAsString = '<' + str(length * size_of_one)+ 'B'
    byte_array = struct.unpack(numBytesAsString, fid.read(length * size_of_one))
    # Reshape 1D byte_array into 2D data_matrix for parsing
    data_matrix = [[byte_array[i + j*size_of_one] for j in range(length)] for i in range(size_of_one)]

    # Now extract the pieces we need
    # added for version compatibility
    if sys.version_info[0] < 3:
        to_bytes = bytearray
    else:
        to_bytes = bytes  
    
    Data['Time'] = [((data_matrix[0][j] << 0) + (data_matrix[1][j] << 8) + (data_matrix[2][j] << 16) + (data_matrix[3][j] << 24)) / sampling_rate for j in range(length)] # Cast 4 bytes as uint32, then divide by sampling_rate
    Data['Clamp'] = [struct.unpack('<f', to_bytes([data_matrix[4][j], data_matrix[5][j], data_matrix[6][j], data_matrix[7][j]]))[0] for j in range(length)] # Cast 4 bytes as single-precision float
    Data['TotalClamp'] = [struct.unpack('<f', to_bytes([data_matrix[8][j], data_matrix[9][j], data_matrix[10][j], data_matrix[11][j]]))[0] for j in range(length)] # Cast 4 bytes as single-precision float
    Data['Measured'] = [struct.unpack('<f', to_bytes([data_matrix[12][j], data_matrix[13][j], data_matrix[14][j], data_matrix[15][j]]))[0] for j in range(length)] # Cast 4 bytes as single-precision float

    return Data


def read_aux_data(data, fid, numADCs):
    # Determine the size of the input file
    start_pos = fid.tell()
    fid.seek(0, 2)
    file_size = fid.tell()

    # Return the fid object to its previous position
    fid.seek(start_pos, 0)

    # Read one
    data['Data'] = read_one_batch_aux(fid, file_size - start_pos, data['Header']['Settings']['SamplingRate'], numADCs)

    return data


def read_one_batch_aux(fid, size, sampling_rate, numADCs):
    # Create dictionary to contain Data
    Data = {}
    size_of_one = 4 + 2 + 2 + 2 * numADCs # Timestamps + DigIn + DigOut + ADCs
    length = int(size / size_of_one)

    # Read in the whole thing at once
    numBytesAsString = '<' + str(length * size_of_one) + 'B'
    byte_array = struct.unpack(numBytesAsString, fid.read(length * size_of_one))
    # Reshape 1D byte_array into 2D data_matrix for parsing
    data_matrix = [[byte_array[i + j*size_of_one] for j in range(length)] for i in range(size_of_one)]

    # Now extract the pieces we need
    Data['Time'] = [((data_matrix[0][j] << 0) + (data_matrix[1][j] << 8) + (data_matrix[2][j] << 16) + (data_matrix[3][j] << 24)) / sampling_rate for j in range(length)] # Cast 4 bytes as uint32, then divide by sampling_rate
    Data['DigitalIn'] = [((data_matrix[4][j] << 0) + (data_matrix[5][j] << 8)) for j in range(length)] # Cast 2 bytes as uint16
    Data['DigitalOut'] = [((data_matrix[6][j] << 0) + (data_matrix[7][j] << 8)) for j in range(length)] # Cast 2 bytes as uint16
    Data['ADC'] = []
    for i in range(numADCs):
        Data['ADC'].append([0.0003125 * (((data_matrix[2*i + 8][j] << 0) + (data_matrix[2*i + 9][j] << 8)) - pow(2,15)) for j in range(length)]) # Cast as uint16, subtract by 2^15, and multiply by 0.0003125

    return Data