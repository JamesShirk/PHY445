import glob
import numpy as np
import sys

# example command
# python3 cleandata.py "./[0-9]0deg*.Spe" main.out
# gotta put quotes around the files

files = glob.glob(sys.argv[1])
files = sorted(files)
print(files)
headers = np.array([i[2:-4].split("_") for i in files], dtype = str)
data = []

for file in files:
    out = []
    with open(file) as f:
        i = 0
        for line in f:
            i += 1  
            if (i < 13) or (i > 1035):
                continue
            out.append(line.strip())
    f.close()
    data.append(out)

data = np.array(data, dtype = np.int64)
# print(data.shape)
# print(headers.shape)
form = ""
for i in range(len(data.T[0])-1):
    form += "%d,"
form += "%d"
np.savetxt(sys.argv[2], 
           data.T, 
           # removes most headers, adds back name of isotope + type of data
           header = np.array2string(headers.T[0], separator = ",")[1:-1] + '\n' + np.array2string(headers.T[1], separator = ",").replace("\n", "")[1:-1], 
           delimiter = ",", 
           fmt = form, 
           comments = '')
    