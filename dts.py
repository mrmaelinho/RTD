dir = ""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from scipy.integrate import quad

colours = list(mcolors.TABLEAU_COLORS)
colours.append('k')

np.set_printoptions(precision=2)
plt.rcParams['axes.facecolor'] = (230/256,230/256,230/256)
plt.rcParams['figure.facecolor'] = (230/256,230/256,230/256)
plt.rcParams['figure.figsize'] = (20,10)
plt.rcParams['font.size'] = 15

rates = [125,250,500,750,1000,1250,2500]
etages = [2,4,6,8,10]
data = np.zeros((len(rates),len(etages)), dtype='object')

## Data collection

for j in range(len(etages)):
    for i in range(len(rates)):
        if rates[i]<1000:
            data_file = dir+str(etages[j])+"-0"+str(rates[i])+".txt"
        if rates[i]>=1000:
            data_file = dir+str(etages[j])+"-"+str(rates[i])+".txt"
        a_data = pd.read_csv(str(data_file), header=5, delimiter="\t", usecols=[1,2], names=["time", "absorbance"])
        a_data['absorbance'] -= a_data['absorbance'].min()
        data[i,j] = a_data

##Plot raw data
for j in range(len(etages)):
    fig,ax=plt.subplots()
    for i in range(len(rates)):
            plt.plot(data[i,j]['time'],data[i,j]['absorbance'], color=colours[i],label=str(rates[i])+'µL/mn', markersize=1)
    ax.set(xlabel='time (s)')
    plt.ylabel('absorbance (u.a.)')
    plt.legend()
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.grid(linewidth=1,which='major',color='k')
    ax.xaxis.grid(linewidth=1,which='minor')
    plt.title(str(etages[j])+' étages')
    plt.savefig(dir+'raw-%d-etages.jpg'%etages[j],transparent=True)
    # plt.show()
    plt.close()


##Plot data on reduced axes and Bo=f(Q)
def DTS(t,P,tau,k):
    return k*.5*(P/(np.pi*t))**.5*np.exp(-P*(tau-t)**2/(4*t))

open(dir+'thetas.txt', 'w').close()
Bos = np.zeros((len(rates),len(etages)))
thetas = np.zeros((len(rates),len(etages)))
stdBos = np.zeros((len(rates),len(etages)))
for j in range(len(etages)):
    fig1,ax1=plt.subplots()
    fig2,ax2=plt.subplots()

    for i in range(len(rates)):
        t = data[i,j]['time']*rates[i]/(60*(125*etages[j]+20))
        A = data[i,j]['absorbance']
        try:
            popt, pcov = curve_fit(DTS,t,A)
        except:
            print('error on dataset (%d,%d)'%(i,j))
        A_fit = DTS(t,popt[0],popt[1],popt[2])
        if popt[2]<0:
            print ('i=%d,j=%d,k=%.2f'%(i,j,popt[2]))
        Bos[i,j]=popt[0]
        thetas[i,j]=popt[1]
        stdBos[i,j]=np.sqrt(np.diag(pcov))[0]
        ax1.plot(t,A/popt[2], color=colours[i],label=str(rates[i])+'µL/mn', markersize=1)
        ax1.plot(t,A_fit/popt[2], color=colours[i],linestyle=':',label=str(rates[i])+'µL/mn', markersize=1)
        ax2.scatter(rates[i],popt[0])
    ax1.set(xlabel=r'$\theta=t/\tau$')
    ax1.set(ylabel=r'E($\theta$)')
    ax1.set(xlim=(0,2))
    ax1.legend()
    ax1.grid()
    ax1.set(title=str(etages[j])+' étages')
    ax2.set(xlabel='flow-rate (µL/mn)')
    ax2.set(ylabel=r'Bo=$vL/D_{ax}$')
    ax2.set(xlim=(0,2600))
    ax2.grid()
    ax2.set(title=str(etages[j])+' étages')
    fig1.savefig(dir+'dts-%d-etages.jpg'%(etages[j]),transparent=True)
    fig2.savefig(dir+'Bo-%d-etages.jpg'%(etages[j]),transparent=True)
    plt.close(fig1)
    plt.close(fig2)

## Plot statistical parameters of distributions
fBo,aBo = plt.subplots()
fTheta,aTheta = plt.subplots()

for i in range(len(etages)):
    for j in range(len(rates)):
        aBo.errorbar(rates,Bos[:,i],yerr = 2*stdBos[:,i] ,color=colours[i], label='%d étages'%etages[i],markersize=8, marker='x', linestyle=':', capsize=5)
        aTheta.plot(rates,thetas[:,i],color=colours[i],label='%d étages'%etages[i],markersize=8, marker='x', linestyle=':')

aBo.set(ylabel=r'Bo=$vL/D_{ax}$')
aBo.set(xlabel='flow rate (µL/mn)')
# aBo.set(xticks=rates)
aBo.grid()

aTheta.set(ylabel=r'$\overline{\theta}=\overline{t_s}/\tau$')
aTheta.set(xlabel='flow rate (µL/mn)')
# aTheta.set(xticks=rates)
aTheta.grid()
handles, labels = fBo.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
aBo.legend(reversed(by_label.values()), reversed(by_label.keys()),loc=1,bbox_to_anchor=(0.25,0.9))

handles, labels = fTheta.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
aTheta.legend(reversed(by_label.values()), reversed(by_label.keys()),loc=1,bbox_to_anchor=(0.25,0.9))

fBo.savefig(dir+'Bo-all.jpg',transparent=True)
fTheta.savefig(dir+'thetas-all.jpg',transparent=True)

plt.close(fBo)
plt.close(fTheta)

##Bo vs. flow rate regression
from scipy.stats import linregress

speeds = np.array(rates)*1.600e-11/25e-8
slopes = []
intercepts = []
rvalues = []

for i in range(5):
    reg = linregress(speeds[2:],Bos[2:,i])[0:3]
    slopes.append(reg[0])
    intercepts.append(reg[1])
    rvalues.append(reg[2])

fBo,aBo = plt.subplots()

for i in range(len(etages)):
    for j in range(len(speeds)):
        aBo.errorbar(speeds,Bos[:,i],yerr = 2*stdBos[:,i] ,color=colours[i], label='%d étages'%etages[i],markersize=8, marker='x', linestyle=':', capsize=5)
        aBo.plot(speeds,intercepts[i]+slopes[i]*speeds,color='k',linestyle="-")

aBo.set(ylabel=r'Bo=$vL/D_{ax}$')
aBo.set(xlabel='speed (m/s)')
# aBo.set(xticks=rates)
aBo.grid()

handles, labels = fBo.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
aBo.legend(reversed(by_label.values()), reversed(by_label.keys()),loc=1,bbox_to_anchor=(0.25,0.9))

fBo.savefig(dir+'Bo-all-speed.jpg',transparent=True)
plt.close(fBo)

##Bo vs. Re, td/tc
Re = np.array(speeds)*5e2
td = 0.03*25e-8/6e-10
tc = np.zeros((len(speeds),len(etages)),dtype=float)

plt.subplots(1,2,sharey=True)
for j in range(len(etages)):
    for i in range(len(speeds)):
        tc[i,j] = (0.5*etages[j])/speeds[i]
    plt.subplot(121)
    plt.plot(Re,Bos[:,j],color=colours[j], marker='x',linestyle=':', label='%d étages'%etages[j])
    plt.subplot(122)
    plt.plot(td/tc[:,j],Bos[:,j],color=colours[j], marker='x',linestyle=':', label='%d étages'%etages[j])

plt.subplot(121)
plt.title('$Bo=f(Re)$')
plt.ylabel(r'Bo=$vL/D_{ax}$')
plt.xlabel(r'$Re=vd_t/\nu$')
plt.grid()

plt.subplot(122)
plt.title('$Bo=f(t_d/t_c)$')
plt.xlabel(r'$t_d/t_c$')
plt.legend()
plt.grid()
plt.savefig(dir+'bo-re-tdtc.jpeg',transparent=True)
plt.close()

##Dax vs. Re, td/tc
Dax = np.zeros((len(speeds),len(etages)),dtype=float)

plt.subplots(1,2,sharey=True)
plt.title('Evolution du coefficient de dispersion axial')
for j in range(len(etages)):
    for i in range(len(speeds)):
        Dax[i,j] = 0.5*etages[j]*speeds[i]/Bos[i,j]
    plt.subplot(121)
    plt.plot(Re,Dax[:,j],color=colours[j], marker='x',linestyle=':', label='%d étages'%etages[j])
    plt.subplot(122)
    plt.plot(td/tc[:,j],Dax[:,j],color=colours[j], marker='x',linestyle=':', label='%d étages'%etages[j])

plt.subplot(121)
plt.title('$D_{ax}=f(Re)$')
plt.ylabel(r'$D_{ax}=vL/Bo\ (m^2.s^{-1})$')
plt.xlabel(r'$Re=vd_t/\nu$')
plt.grid()

plt.subplot(122)
plt.title('$D_{ax}=f(t_d/t_c)$')
plt.xlabel(r'$t_d/t_c$')
plt.legend()
plt.grid()
plt.savefig(dir+'dax-re-tdtc.jpeg',transparent=True)
plt.close()

##Dmcalc vs. Re, td/tc
Dmc = np.zeros((len(speeds),len(etages)),dtype=float)

plt.subplots(1,2,sharey=True)
plt.title('Coefficient de diffusion moléculaire calculé')
for j in range(len(etages)):
    for i in range(len(speeds)):
        Dmc[i,j] = (25e-8/192)*speeds[i]**2/Dax[i,j]
    plt.subplot(121)
    plt.plot(Re,Dmc[:,j]/6e-10,color=colours[j], marker='x',linestyle=':', label='%d étages'%etages[j])
    plt.subplot(122)
    plt.plot(td/tc[:,j],Dmc[:,j]/6e-10,color=colours[j], marker='x',linestyle=':', label='%d étages'%etages[j])

plt.subplot(121)
plt.hlines(1,0,80,colors='k')
plt.title(r'$D_{m,calc}/D_{m,th}=f(Re)\ D_{m,calc}=v^2d_t^2/D_{ax}$')
plt.ylabel(r'$D_{m,calc}/D_{m,th}$ (dimensionless)')
plt.xlabel(r'$Re=vd_t/\nu$')
plt.grid()

plt.subplot(122)
plt.hlines(1,0,2,colors='k')
plt.title('$D_{m,calc}/D_{m,th}=f(t_d/t_c)$')
plt.xlabel(r'$t_d/t_c$')
plt.legend()
plt.grid()
plt.savefig(dir+'dmc-re-tdtc.jpeg',transparent=True)
plt.close()

## Re.Sc=vdt/Dm vs. L/dt
resc = np.array(speeds)*5e-4/6e-10
Ldt = np.array(etages)*0.5/5e-4

for i in range(len(Ldt)):
    plt.loglog([Ldt[i]]*len(resc),resc,linewidth=0, marker='x',label='%d étages'%etages[i])

plt.xlabel(r'$L/d_t$')
plt.ylabel(r'$Re.Sc=vL/D_m$')
plt.xlim(1,1e5)
plt.ylim(1,1e6)
plt.fill((6,6,1e5,1e5),(1,30,1e6,1),'green',alpha=.1)
plt.fill((1,1,2000,1e5,6,6),(1,500,1e6,1e6,30,1),'yellow',alpha=.1)
plt.fill((1,1,2000),(500,1e6,1e6),'red',alpha=.1)
plt.legend(loc='lower right')
plt.annotate('convection pure',(2,1e5),color='red')
plt.annotate('transition',(80,5000),color='yellow')
plt.annotate('dispersion axiale',(1e3,100),color='green')
plt.grid()
plt.show()


