To use from the command-line, you should first make sure that Python 3 is installed on your machine.

If you're using Windows, open 'cmd.exe', type 'python' and press Enter. This should give you some basic information about your
Python installation, including version number. This script is written to work with Python 3. You can then enter 'exit()' to return
to the command line prompt. If instead you get a message like 'python is not recognized as an internal command', that indicates
that either Python isn't installed on your machine, or it's not present in the PATH environment variable.

The default behavior of this script is store the header information and data contained in the .clp file in an object called 'a'.
This can be changed in the 'read_intan_CLP_file.py' file. The structure of this object is split in two parts: Header and Data.
The structure of 'Header' is specified in the 'read_header.py' file in intanutil. The structure of 'Data' is specified in the
'read_data.py' file in intanutil. For example, if you want to see the sample rate, you can add in
'print(a['Header']['Settings']['SamplingRate'])' after the read_intan_CLP_file() function call. If you want to see measured value
from the 50th sample, you can add in 'print([a['Data']['Measured'][49])' (Keeping in mind that Python uses 0-indexing. The first entry in a list is [0],
the second entry is [1], and so on).