import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
plt.rcParams['text.usetex'] = True

beta_min = 0;
delta_beta = 0.1;


# this is how i run it: $python processor.py




filename = "NaCl_feb6_peak1.xry";
peak1 = 14.30;
height1 = 400;
range1 = [14.1,14.6];

filename = "NaCl_feb6_peak2.xry";
peak1 = 15.80;
height1 = 1535;
range1 = [15.6,16.1];

filename = "NaCl_feb6_peak3.xry";
peak1 = 29.6;
height1 = 170;
range1 = [29.3,29.8];


filename = "NaCl_fourth_peak.xry";
peak1 = 33.1;
height1 = 170;
range1 = [32.8,33.4];

filename = "KCl_first_peak_10.xry";
peak1 = 13.0;
height1 = 170;
range1 = [32.8,33.4];

filename = "KCl_second_peak_10.xry";
peak1 = 14.4;
height1 = 170;
range1 = [32.8,33.4];

filename = "KCl_third_peak_10.xry";
peak1 = 26.4;
height1 = 170;
range1 = [26.4,33.4];

filename = "KCl_fourth_peak_10.xry";
peak1 = 29.53;
height1 = 11;
range1 = [26.4,33.4];

filename = "NaCl_feb6_peak2.xry";
peak1 = 15.80;
height1 = 1535;
range1 = [15.6,16.1];

f = open("./../data/diffraction/" + filename, "rb")
#f = open("NaCl_feb6_peak3.xry", "r")

# Fitting the data
# Lorentzian fitting function
def lorentz(x,amp,mu,sigma):
    return amp * sigma**2 / ((x - mu)**2 + sigma**2)


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


#in some files, there are mulitple data columns. So you have to chose which column to read
column = 0;

for i in range(0,counter):
	string = (f.readline().split())
	try: hits = float(string[column])
	except: hits = 0
	y_list.append(hits)
	x_list.append(beta_min + i*delta_beta)


x = np.array(x_list);
y = np.array(y_list);
plt.plot(x,y,color="blue",linestyle="dotted", label='Actual Data');


'''
lower_index = int((range1[0] - beta_min)/delta_beta);
upper_index = int((range1[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
'''
x_narrow = np.array(x_list);
popt, pcov = curve_fit( lorentz,x ,y,p0=(height1,peak1,0.1));
#popt, pcov = curve_fit( gauss,x ,y,p0=(height1,peak1,0.1));
print("Lorentzian ...")
print(popt);
print(pcov);
perr = np.sqrt(np.diag(pcov))
print(perr)
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
#chi_sq = np.sum((y - y_narrow)**2)
#print("reduced Chisq is : ", chi_sq/len(y))
plt.plot(x_narrow,y_narrow,color="orange",linestyle="solid",label='Lorentzian Fit');


'''
lower_index = int((range1[0] - beta_min)/delta_beta);
upper_index = int((range1[1] - beta_min)/delta_beta);
x_narrow = np.array(x_list[lower_index:upper_index]);
'''
x_narrow = np.array(x_list);
#popt, pcov = curve_fit( lorentz,x ,y,p0=(height1,peak1,0.1));
popt, pcov = curve_fit( gauss,x ,y,p0=(height1,peak1,0.1));
print("Gaussian ...")
print(popt);
print(pcov);
perr = np.sqrt(np.diag(pcov))
print(perr)
y_narrow=gauss(x_narrow,popt[0],popt[1],popt[2]);
#plt.plot(x_narrow,y_narrow,color="red",linestyle="solid",label='Gaussian Fit');


plt.legend();
plt.xlabel(r'$\theta (\ ^\circ )$');
plt.ylabel('Counts(1/s)');
plt.title(r'$L_\alpha$ peak in KCl diffraction')
plt.show();
#plt.savefig("diffraction"+filename+".png", dpi =100)




'''

plt.plot(x_list,y_list,color="black",linestyle="dotted")
plt.show()



plt.plot(x,y)
plt.show()
'''

