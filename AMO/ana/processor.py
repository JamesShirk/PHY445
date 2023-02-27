import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import math
plt.rcParams['text.usetex'] = True

beta_min = 0;
delta_beta = 0.1;


# this is how i run it: $python processor.py


filename = "NaCl_feb6_peak1.xry";
peak1 = 14.30;
height1 = 400;
range1 = [14.1,14.6];
peak2 = 15.80;
height2 = 1711;
range2 = [15.6,16.2];
peak3 = 29.6;
height3 = 170;
range3 = [29.3,29.8];

filename = "NaCl_feb6_fine.xry";
peak1 = 14.30;
height1 = 587;
range1 = [14.0,14.7];
peak2 = 15.80;
height2 = 1711;
range2 = [15.5,16.3];
peak3 = 29.6;
height3 = 170;
range3 = [29.3,29.9];


filename = "KCL_Feb8_Coarse3.xry";
peak1 = 13.12;
height1 = 39.9;
range1 = [12.90,13.50];
peak2 = 14.6;
height2 = 139;
range2 = [14.3,15.0];
peak3 = 29.7;
height3 = 19.3;
range3 = [29.4,30.0];


filename = "NaCl_Feb15.xry";
peak1 = 14.30;
height1 = 587;
range1 = [14.0,14.7];
peak2 = 15.80;
height2 = 1711;
range2 = [15.5,16.3];
peak3 = 29.6;
height3 = 170;
range3 = [29.3,29.9];
peak4 = 33.1;
height4 = 555;
range4 = [32.9,33.5];


filename = "KCL_Feb8_Fine.xry";
peak1 = 13.20;
height1 = 33.1;
range1 = [12.9,13.5];
peak2 = 14.6;
height2 = 117;
range2 = [14.3,15.0];
peak3 = 26.5;
height3 = 6.7;
range3 = [26.4,26.9];
peak4 = 29.7;
height4 = 17.0;
range4 = [29.3,30.1];





f = open("./../data/diffraction/" + filename, "rb")
#f = open("NaCl_feb6_peak3.xry", "r")


def gauss(x,amp,mu,sigma):
	return amp*np.exp( - (x-mu)**2 / (2*sigma*sigma) )


# first few lines are meaningless for us. So I will dump the first "dump" many lines.
dump = 4
for i in range(0,dump):
	f.readline()

infoLine = f.readline();
info = infoLine.split();
beta_min = float(info[0]);
delta_beta = float(info[3]);

print("betamin", beta_min)
print("deltabeta", delta_beta)

dump = 12
for i in range(0,dump):
	f.readline()

infoLine = f.readline();
info = infoLine.split();
counter = int(info[1]);

print("Angle starting from " +  str(beta_min) + " in chunks of " + str(delta_beta) )

print("Now reading file: ", filename)

x_list = []
y_list = []


#in some files, there are mulitple data columns. So you have to chose which column to read
column = 0;

for i in range(0,counter):
	string = (f.readline().split())
	try: 
		hits = float(string[column])
		if not math.isnan(hits):
			y_list.append(hits)
			x_list.append(beta_min + i*delta_beta)
	except: hits = 0

x = np.array(x_list);
y = np.array(y_list);
plt.plot(x,y,color="blue",linestyle="solid", label='Actual Data');


'''
delta_beta = 0.1;

print("First peak .....")
lower_index = int((range1[0] - beta_min)/delta_beta);
upper_index = int((range1[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
popt, pcov = curve_fit( gauss,x ,y,p0=(height1,peak1,0.1));
print(popt);
print(pcov);
perr = np.sqrt(np.diag(pcov))
print(perr)
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
plt.plot(x_narrow,y_narrow,color="orange",linestyle="solid",label='Gaussian Fit');


print("Second peak .....")
lower_index = int((range2[0] - beta_min)/delta_beta);
upper_index = int((range2[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
popt, pcov = curve_fit( gauss,x ,y,p0=(height2,peak2,0.1));
print(popt);
print(pcov);
perr = np.sqrt(np.diag(pcov))
print(perr)
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
plt.plot(x_narrow,y_narrow,color="orange",linestyle="solid");



print("Third peak .....")
lower_index = int((range3[0] - beta_min)/delta_beta);
upper_index = int((range3[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
popt, pcov = curve_fit( gauss,x ,y,p0=(height3,peak3,0.1),bounds=([0,range3[0],0],[10000,range3[1],3]));
print(popt);
print(pcov);
perr = np.sqrt(np.diag(pcov))
print(perr)
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
plt.plot(x_narrow,y_narrow,color="orange",linestyle="solid");


print("Fourth peak .....")
lower_index = int((range4[0] - beta_min)/delta_beta);
upper_index = int((range4[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
popt, pcov = curve_fit( gauss,x ,y,p0=(height4,peak4,0.1),bounds=([0,range4[0],0],[10000,range4[1],3]));
print(popt);
print(pcov);
perr = np.sqrt(np.diag(pcov))
print(perr)
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
plt.plot(x_narrow,y_narrow,color="orange",linestyle="solid");

'''


#plt.legend();
plt.xlabel(r'$ \theta (\ ^\circ )$');
plt.ylabel('Counts(1/s)');
plt.title('X Ray diffraction from KCl crystal')
plt.show();
#plt.savefig("diffraction"+filename+".png", dpi =100)




'''

plt.plot(x_list,y_list,color="black",linestyle="dotted")
plt.show()



plt.plot(x,y)
plt.show()
'''

