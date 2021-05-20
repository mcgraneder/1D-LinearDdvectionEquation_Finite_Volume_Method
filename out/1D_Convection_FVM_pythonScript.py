'''
import numpy as np
import matplotlib.pyplot as plt

filelist=[]

liimited_filelist = []


x = np.linspace(-1, 1, 200)
j = 1
for i in range(1,100):
    filelist.append("FVM_1S_B&amp;W_{}.txt".format(i))
    liimited_filelist.append("BeamW_Limited_{}.txt".format(i))

for fname, fname1 in zip(filelist, liimited_filelist):
    plt.figure(figsize=(12,6))
    data=np.loadtxt(fname)
    data1=np.loadtxt(fname1)
    #initial_data = np.loadtxt(initial_file[0])
    #initial_data = np.loadtxt(filelist[0])
    plt.plot(x, data, color = 'orange', linewidth = 0.7, label = 'B&W unlimited')
    plt.plot(x, data1, marker = 'o', markersize = 4, linewidth = 1.3, color = 'blue', label = 'B&W Method 2nd Order')
    #plt.plot(x, data, marker = 'o', markersize = 4, linewidth = 2, color = 'darkorange', label = 'B&W method Van\nAlbada Slope Limited')
    plt.plot(-1,1, color = 'black', label = 'timestep {}'.format(j))
    plt.title('1D Linear Convection 2nd Order Beam & Warming FVM Analysis ', fontsize=15)
    plt.xlabel('u (m/s)', fontsize=15)
    plt.ylabel('time t(s)', fontsize=15)
    plt.ylim(-0.4, 1.4)
    plt.xlim(-1.1, 1.1)
    #plt.title("1D Linear Convection FVM")
    plt.legend(loc = 'upper right')
    plt.savefig('B&W{}.png'.format(j))
    plt.clf()
    plt.close()
    j += 1

#plt.legend()
plt.show()



'''
'''
import numpy as np
import matplotlib.pyplot as plt

filelist=[]

liimited_filelist = []


x = np.linspace(-1, 1, 199)
j = 1
for i in range(1,104):
    filelist.append("FrommLimited_{}.txt".format(i))
    liimited_filelist.append("FROMM_UNLIMITED_{}.txt".format(i))

for fname, fname1 in zip(filelist, liimited_filelist):
    plt.figure(figsize=(12,8))
    data=np.loadtxt(fname)
    data1=np.loadtxt(fname1)
    #initial_data = np.loadtxt(initial_file[0])
    #initial_data = np.loadtxt(filelist[0])
    #plt.plot(x1, initial_data, color = 'orange', linewidth = 0.7, label = 'Initial Condition')
    plt.plot(x, data1, marker = 'o', markersize = 4, linewidth = 1.3, color = 'blue', label = 'Lax Fromm Method 2nd Order')
    plt.plot(x, data, marker = 'o', markersize = 4, linewidth = 2, color = 'darkorange', label = 'Fromm method Van\nAlbada Slope Limited')
    plt.plot(-1,1, color = 'black', label = 'timestep {}'.format(j))
    plt.title('1D Linear Convection 2nd Order Fromm Method FVM Analysis ', fontsize=15)
    plt.xlabel('u (m/s)', fontsize=15)
    plt.ylabel('time t(s)', fontsize=15)
    plt.xlabel('u (m/s)', fontsize=15)
    plt.ylabel('time t(s)', fontsize=15)
    plt.ylim(-0.4, 1.4)
    plt.xlim(-1.1, 1.1)
    #plt.title("1D Linear Convection FVM")
    plt.legend(loc = 'upper right')
    plt.savefig('Fromm{}.png'.format(j))
    plt.clf()
    plt.close()
    j += 1

#plt.legend()
plt.show()

'''

import numpy as np
import matplotlib.pyplot as plt

filelist=[]

liimited_filelist = []


x = np.linspace(-1, 1, 100)
x_ICs = np.linspace(-1, 1, 106)
    
j = 1

for i in range(1,5):
    filelist.append("Data_{}.txt".format(i))


plt.figure(figsize = (8,8))
data1=np.loadtxt(filelist[0])
data2=np.loadtxt(filelist[1])
data3=np.loadtxt(filelist[2])
data4=np.loadtxt(filelist[3])
iC_square=np.loadtxt('Square_Initial_conditions.txt')
iC_sine=np.loadtxt('Sine_Initial_Conditions.txt')
#data1=np.loadtxt(fname1)
    #initial_data = np.loadtxt(initial_file[0])
    #initial_data = np.loadtxt(filelist[0])
    #plt.plot(x1, initial_data, color = 'orange', linewidth = 0.7, label = 'Initial Condition')
#plt.plot(x, data1, marker = 'o', markersize = 4, linewidth = 1.3, color = 'blue', label = 'Lax Fromm Method 2nd Order')
plt.plot(x, data1, marker = 'o', markersize = 0.5, linewidth = 1.5, color = 'blue', label = 'CFL = 0.2')
plt.plot(x, data2, marker = 'o', markersize = 0.5, linewidth = 1.5, color = 'green', label = 'CFL = 0.4')
plt.plot(x, data3, marker = 'o', markersize = 0.5, linewidth = 1.5, color = 'red', label = 'CFL = 0.6')
plt.plot(x, data4, marker = 'o', markersize = 0.5, linewidth = 1.5, color = 'darkturquoise', label = 'CFL = 0.8')

plt.plot(x_ICs, iC_square, linewidth = 2, linestyle = 'dashed', color = 'mediumvioletred', label = 'Initial Conditions')
#plt.plot(x_ICs, iC_sine, linewidth = 2, color = 'red', linestyle = 'dashed', label = 'Initial Conditions')
#plt.plot(-1,1, color = 'black', label = 'timestep {}'.format(j))
plt.title('1D Linear Convection, Time = 2.0 s ', fontsize=15)
#plt.xlabel('u (m/s)', fontsize=15)
#plt.ylabel('time t(s)', fontsize=15)
#plt.xlabel('u (m/s)', fontsize=15)
#plt.ylabel('time t(s)', fontsize=15)
plt.ylim(-0.2, 1.2)
plt.xlim(-1.1, 1.1)
plt.xticks([-1, -0.5, 0, 0.5, 1])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    #plt.title("1D Linear Convection FVM")
plt.legend(loc = 'upper right')

j += 1

#plt.legend()
plt.show()



