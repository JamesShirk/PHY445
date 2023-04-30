import glob
import numpy as np
import sys

# example command
# python3 cleandata.py "./[0-9][0-9]deg*.Spe" main.out
# gotta put quotes around the files

# python3 cleandata.py "./[0-9]\d{1,2}deg*.Spe" main.out

# get files
files = glob.glob(sys.argv[1])
print(files)

# extract header array, e.g. for 60deg_main would give 60deg, main
headers = np.array([i[2:-4].split("_") for i in files], dtype = str)
# trims deg off of the first part of above
headers.T[0] = np.array([i[0:-3] for i in headers.T[0]])
# sorts the files 
#files = [file for _,file in sorted(zip(headers.T[0].astype(np.int32), files))]

print(headers)

sort = headers.T[0].astype(np.int32).argsort()
files = np.array(files)[sort]
headers = headers[sort]


print(headers)

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

print(len(data))
for i in data:
    print(len(i))

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
           header = np.array2string(headers.T[0].astype(np.int32), separator = ",")[1:-1] + '\n' + np.array2string(headers.T[1], separator = ",").replace("\n", "")[1:-1], 
           delimiter = ",", 
           fmt = form, 
           comments = '')
    