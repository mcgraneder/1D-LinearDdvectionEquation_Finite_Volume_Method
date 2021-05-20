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




import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
filelist=[]

#liimited_filelist = []
#iC_square=np.loadtxt('Data.txt')

#d1 = np.loadtxt('sol_10.txt')
x = np.linspace(0, 1, 100)
j = 1
for i in range(1,200):
    #if i % 5 == 0:
    filelist.append("Burgers_{}.txt".format(i))
    #liimited_filelist.append("FROMM_UNLIMITED_{}.txt".format(i))

#colour=iter(cm.rainbow(np.linspace(0,20,860)))
#plt.plot(x, d1, marker = 'o', markersize = 1, linewidth = 2, color = 'blue', label = 't = 0.0 s')
#plt.plot(x, iC_square, marker = 'o', markersize = 1.6, linewidth = 2, color = 'darkblue', label = 'CFL = 0.5')


for fname  in filelist:
    plt.figure(figsize=(8,8))
    #c=next(colour)
    data=np.loadtxt(fname)
    #data1=np.loadtxt(fname1)
    #initial_data = np.loadtxt(initial_file[0])
    #initial_data = np.loadtxt(filelist[0])
    plt.plot(x, data, marker = 'o', markersize = 5, linewidth = 2, color = 'blue', label = 'CFL = 0.5')

    #plt.plot(x1, initial_data, color = 'orange', linewidth = 0.7, label = 'Initial Condition')
    #plt.plot(x, data, marker = 'o', markersize = 2.1, linewidth = 1.3, color = 'blue', label = 'FROMM Method 2nd Order')
    #plt.plot(x, data, marker = 'o', markersize = 4, linewidth = 2, color = 'darkorange', label = 'Fromm method Van\nAlbada Slope Limited')
    #plt.plot(-1,1, color = 'black', label = 'timestep {}'.format(j))
    #plt.title('1D Inviscid Burgers Sine wave ICs ', fontsize=15)
    plt.ylabel('u (m/s)', fontsize=15)
    plt.xlabel('u (m)', fontsize=15)
    plt.ylim(-0.25, 1.25)
    plt.xlim(0, 1)
    plt.title("1D Advection equation")
    plt.legend(loc = 'upper right')
    plt.savefig('A{}.png'.format(j))
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

x = np.linspace(0, 2 * np.pi, 100)
x1 = np.linspace(0, 2 * np.pi, 100)
x2 = np.linspace(-1, 1, 100)
x4 = np.linspace(-1, 1, 200)
x_ICs = np.linspace(-1, 1, 206)
    
j = 1

#for i in range(1,5):
#    filelist.append("laxwendroff_limited_step{}.txt".format(i))


plt.figure(figsize = (10,8))
data1=np.loadtxt('Data.txt')
#data2=np.loadtxt('sol_9.txt')


#iC_square=np.loadtxt('Square_Initial_conditions.txt')
#iC_sine=np.loadtxt('Sine_Initial_Conditions.txt')
#data1=np.loadtxt(fname1)
    #initial_data = np.loadtxt(initial_file[0])
    #initial_data = np.loadtxt(filelist[0])
    #plt.plot(x1, initial_data, color = 'orange', linewidth = 0.7, label = 'Initial Condition')
#plt.plot(x, data1, marker = 'o', markersize = 4, linewidth = 1.3, color = 'blue', label = 'Lax Fromm Method 2nd Order')
#plt.plot(x_ICs, iC_square, linewidth = 2, linestyle = 'dashed', color = 'red', label = 'Initial Conditions')
plt.plot(x, data1, marker = 'o', markersize = 2, linewidth = 2, color = 'blue', label = 'limited')
#plt.plot(x1, data2, marker = 'o', markersize = 2, linewidth = 2, color = 'red', label = 'limited')

#plt.plot(x, data2, marker = 'o', markersize = 2, linewidth = 2.5, color = 'green', label = 'limited')
#plt.plot(x_ICs, iC_square,  linewidth = 1.6, color = 'red', linestyle='dashed', label = 'limited')

#plt.plot(x2, data3, marker = 'o', markersize = 2, linewidth = 1.7, color = 'red', label = 'NX = 100')

#plt.plot(x4, data4, marker = 'o', markersize = 2, linewidth = 1.7, color = 'darkturquoise', label = 'NX = 100')



#plt.plot(x_ICs, iC_sine, linewidth = 2, color = 'red', linestyle = 'dashed', label = 'Initial Conditions')
#plt.plot(-1,1, color = 'black', label = 'timestep {}'.format(j))
#plt.title('2nd Order Lax Wendroff Slope Limiter Comparison T = 2.0 s', fontsize=15)
#plt.xlabel('u (m/s)', fontsize=15)
#plt.ylabel('time t(s)', fontsize=15)
#plt.xlabel('u (m/s)', fontsize=15)
#plt.ylabel('time t(s)', fontsize=15)
#plt.ylim(-0.1, 1.2)
#plt.xlim(-0.1, 1.1)
#plt.xticks([-1, -0.5, 0, 0.5, 1])
#lt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
#plt.title("F.O Lax Fredrichs:  Time = 2.0 s", fontsize = 20)
plt.legend(loc = 'upper right')

j += 1

#plt.legend()
plt.savefig("F.0.png")
plt.show()

'''
'''
import numpy as np
import matplotlib.pyplot as plt
import math

x = np.linspace(1,180, 180)
for i in range(len(x)):
    x[i] *= 0.0174533
print(x)
y = []
y1 = []
y2 = []

z = []
z1 = []
z2 = []

for i in range(len(x)):
    y.append(   (math.cos(x[i])**2 + (0.25 * 0.25) * math.sin(x[i])**2)**0.5  )
    y1.append(   (math.cos(x[i])**2 + (0.5 * 0.5) * math.sin(x[i])**2)**0.5  )
    y2.append(   (math.cos(x[i])**2 + (0.8 * 0.8) * math.sin(x[i])**2)**0.5  )
    z.append(   (math.atan(.25 * math.tan(x[i])) / (.25 * x[i])))
    z1.append(   (math.atan(.5 * math.tan(x[i])) / (.5 * x[i])))
    z2.append(   (math.atan(.8 * math.tan(x[i])) / (.8 * x[i])))




fig, (ax,ax2) = plt.subplots(1,2, figsize=(13,5))
plt.title("Lax Fredrichs Diffusion Error")
ax.grid(b=True, which='major', color='#666666', linestyle='-')
ax.minorticks_on()
ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

ax2.grid(b=True, which='major', color='#666666', linestyle='-')
ax2.minorticks_on()
ax2.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

ax2.plot(x, z, label = 'Sigma = 0.25', linewidth = 2)
ax2.plot(x, z1, label = 'Sigma = 0.5', linestyle = 'dashdot', linewidth = 2)
ax2.plot(x, z2,label = 'Sigma = 0.8', linestyle = 'dashed', linewidth = 2)

ax.plot(x, y, label = 'Sigma = 0.25', linewidth = 2)
ax.plot(x, y1, label = 'Sigma = 0.5', linestyle = 'dashdot', linewidth = 2)
ax.plot(x, y2,label = 'Sigma = 0.8', linestyle = 'dashed', linewidth = 2)


ax.set_xlabel("Phase angle radians")
ax.set_ylabel("Sigma")
ax2.set_xlabel("Phase angle radians")
ax2.set_ylabel("Sigma")

#ax.set_ylim(0, 6)
ax.set_xlim(0, np.pi)

ax2.set_ylim(0, 6)
ax2.set_xlim(0, np.pi/2.1)


plt.legend()
plt.show()

plt.savefig("Lax-Fredrichs-Diffusion-Error-Plot.png")
'''

