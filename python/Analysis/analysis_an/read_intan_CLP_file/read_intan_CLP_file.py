import sys, struct
from intanutil.read_header import read_header
from intanutil.read_data import read_data
from intanutil.read_data import read_aux_data

def read_intan_CLP_file(filename):
    # Create dictionary 'd' that will hold 'header' and 'data' to be returned
    d = {}

    # Open input rhd file
    fid = open(filename, 'rb')
     
    # Read header of rhd file into d['Header']
    (d['Header'], datatype, NumBytes) = read_header(fid)

    # Make sure file reader position matches what we expect from header
    pos = fid.tell()
    if pos != NumBytes:
        print("Header NumBytes doesn't match number of bytes")
        return

    if datatype == 0:
        # If this is the standard data file (including clamp and measured data), read that data into d['Data']
        d = read_data(d, fid)
    else:
        # If this is the aux data file (including Digital Ins/Outs and ADC data)
        d = read_aux_data(d, fid, d['Header']['NumADCs'])

    print('File contains {:0.3f} seconds of data samples at {:0.2f} kS/s.\n'.format(d['Data']['Time'][len(d['Data']['Time']) - 1], d['Header']['Settings']['SamplingRate'] / 1000))

    if datatype == 0:
        # If this is a standard data file with VClamp 2x mode activated, multiply the 'Clamp' 
        if d['Header']['Settings']['IsVoltageClamp'] and d['Header']['Settings']['vClampX2mode']:
            print('This data was taken using 2x Voltage Clamp Mode.  Doubling clamp data to account for this.\n')
            d['Data']['Clamp'] = [x * 2 for x in d['Data']['Clamp']]

    # Close data file
    fid.close()

    return d
    

if __name__ == '__main__':
    a=read_intan_CLP_file(sys.argv[1])
    #print(a)