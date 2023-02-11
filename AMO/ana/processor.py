import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

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

for i in range(0,counter):
	hits = float(f.readline())
	y_list.append(hits)
	x_list.append(beta_min + i*delta_beta)


x = np.array(x_list);
y = np.array(y_list);
plt.plot(x,y,color="blue",linestyle="dotted", label='Actual Data');


lower_index = int((range1[0] - beta_min)/delta_beta);
upper_index = int((range1[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
popt, pcov = curve_fit( gauss,x ,y,p0=(height1,peak1,0.1));
print(popt);
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
plt.plot(x_narrow,y_narrow,color="orange",linestyle="solid",label='Gaussian Fit');




lower_index = int((range2[0] - beta_min)/delta_beta);
upper_index = int((range2[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
popt, pcov = curve_fit( gauss,x ,y,p0=(height2,peak2,0.1));
print(popt);
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
plt.plot(x_narrow,y_narrow,color="orange",linestyle="solid");


lower_index = int((range3[0] - beta_min)/delta_beta);
upper_index = int((range3[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
popt, pcov = curve_fit( gauss,x ,y,p0=(height3,peak3,0.1));
print(popt);
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
plt.plot(x_narrow,y_narrow,color="orange",linestyle="solid");


plt.legend();
plt.xlabel('Angle');
plt.ylabel('Counts');
#plt.show();
plt.savefig("diffraction"+filename+".png", dpi =100)




'''

plt.plot(x_list,y_list,color="black",linestyle="dotted")
plt.show()



plt.plot(x,y)
plt.show()
'''

