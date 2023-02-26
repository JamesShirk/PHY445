import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import math
plt.rcParams['text.usetex'] = True

beta_min = 0;
delta_beta = 0.1;


# this is how i run it: $python processor.py


filename = "Planck_measure.xry";
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
p_list = []
q_list = []
r_list = []
y_list = []

for i in range(0,counter):
	string = (f.readline().split())
	p_list.append(float(string[0]));
	q_list.append(float(string[1]));
	r_list.append(float(string[2]));
	x_list.append(beta_min + i*delta_beta);
	y_list.append(float(0))



#-------- now the second file, at 2.5 degree --------------

filename = "KCl_Feb8_fine.xry";
peak1 = 14.30;
height1 = 400;
range1 = [14.1,14.6];
peak2 = 15.80;
height2 = 1711;
range2 = [15.6,16.2];
peak3 = 29.6;
height3 = 170;
range3 = [29.3,29.8];


f = open("./../data/diffraction/" + filename, "rb")

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



start_angle = 2.5
end_angle = 10.0
step = 0.1
column = 0
counter = 0
y = np.array(y_list)
baseindex = 5


for i in range(0,100000000):
	string = (f.readline().split())
	try: 
		hits = float(string[column])
		#	print(counter)
		if not math.isnan(hits):
			angle = start_angle + counter*step
			if angle >= end_angle: break	
			index = baseindex + counter
			#print("Index is: ", index)
			#print("Index found at",x_list.index(angle))
			y[index] = hits
			counter = counter + 1
	except Exception 	as e: 
		print("Exception: ", e)
		break
		hits = 0

#-----------------------------------------------------------


x = np.array(x_list[10:]);
p = np.array(p_list[10:]);
q = np.array(q_list[10:]);
r = np.array(r_list[10:]);
y = np.array(y[10:]);
'''
x = np.array(x_list);
p = np.array(p_list);
q = np.array(q_list);
r = np.array(r_list);
y = np.array(y);
'''
plt.plot(x,y,color="black",linestyle="solid", label='35 keV',linewidth="0.7");
plt.plot(x,p,color="red",linestyle="solid", label='30 keV',linewidth="0.7");
plt.plot(x,q,color="blue",linestyle="solid", label='25 keV',linewidth="0.7");
plt.plot(x,r,color="green",linestyle="solid", label='20 keV',linewidth="0.7");



plt.legend();
plt.xlabel(r'$ \theta (\ ^\circ )$',fontsize="14");
plt.ylabel('Counts(1/s)',fontsize="14");
plt.title('Diffraction at different voltages')
plt.grid()
plt.show();



