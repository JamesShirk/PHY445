import glob
import numpy as np

files = glob.glob("./*calib.Spe")

print(files)

headers = np.array([i[2:-4].split("_") for i in files], dtype = str)
data = []

for file in files:
    out = []
    with open(file) as f:
        i = 0
        for line in f:
            i += 1  
            if (i < 13) or (i > 2059):
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
np.savetxt("data.out", data.T, header = np.array2string(headers.T[0], separator = ",")[1:-1] + '\n' + np.array2string(headers.T[1], separator = ",").replace("\n", "")[1:-1], delimiter = ",", fmt = form, comments = '')
    