import pandas as pd
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class Experiment():
    def __init__(self,file_path):
        self.file_path = file_path
        #Unpacking info contained in the filename.
        if ("\\" in file_path):
            self.file_name = self.file_path.split('\\')[-1]
        elif "/" in file_path:
            self.file_name = self.file_path.split('/')[-1]
        self.stages, self.flowrate = self.file_name[:-4].split('-')
        self.raw_data = pd.read_csv(file_path,\
                                    header=5,\
                                    delimiter="\t",\
                                    usecols=[1,2],\
                                    names=["time", "absorbance"])
        self.data = self.raw_data.copy()
        #Shift absorbance to have the minimum equal to 0.
        self.data['absorbance'] = self.raw_data.absorbance\
                    - min(self.raw_data.absorbance)

    def plot_absorbance(self):
        """
        Plots absorbance against time.
        """
        fig,ax=plt.subplots()
        ax.plot(self.data['time'],\
                 self.data['absorbance'],\
                 label=str(self.flowrate)+'µL/mn',\
                 markersize=1)
        ax.set(xlabel='time (s)')
        ax.set_ylabel('absorbance (u.a.)')
        ax.grid()
        ax.set_title('%s stages (%d µL) @%s µL/min'%(self.stages, int(self.stages)*125,self.flowrate))
        fig.savefig(self.file_path[:-4]+'-raw.jpg',transparent=True)

    def RTD_bestfit(self):
        """
        Fits the absorbance data to theoretical shape of RTD.
        """
        def theoretical_RTD(t,tau,Bo,scale):
            """
            Returns the theoretical shape of RTD with given parameters.
            """
            return scale*.5*(Bo/(np.pi*t))**.5*np.exp(-Bo*(tau-t)**2/(4*t))
        #Transit time is V/Q, in seconds
        self.tau = 60 * int(self.stages) * 125 / int(self.flowrate)
        #Reduced time is t/tau
        self.data["theta"] = self.data.time / self.tau
        popt, pcov = curve_fit(theoretical_RTD,\
                                   self.data.theta,\
                                   self.data.absorbance)
        self.tau_opt, self.Bo_opt, self.scale_opt = popt
        self.Bo_std = np.sqrt(np.diag(pcov))[1]
        self.data["absorbance_fit"] = theoretical_RTD(self.data.theta,\
                                                          self.tau_opt,\
                                                          self.Bo_opt,\
                                                          self.scale_opt)
    def plot_bestfit(self):
        """
        Plots the results of the RTD best fit.
        """
        fig, ax = plt.subplots()
        ax.set_xlabel('$\Theta$') #Reduced time absciss
        ax.set_ylabel('E($\Theta$)')#Adimensional residence time distribution
        ax.plot(self.data.theta,\
                #Absorbance data must be scaled to become adimensional
                self.data.absorbance/self.scale_opt,\
                color = 'k',\
                markersize = 1)
        ax.plot(self.data.theta,\
                #Absorbance data must be scaled to become adimensional
                self.data.absorbance_fit/self.scale_opt,\
                color = 'r',\
                linestyle =':',\
                markersize = 1)
        ax.grid()
        ax.set_title('%s stages (%d µL) @%s µL/min : Bo = %.1f$\pm$%.1f'%\
                     (self.stages,\
                     int(self.stages) * 125,\
                     self.flowrate,\
                     self.Bo_opt,\
                     self.Bo_std))
        fig.savefig(self.file_path[:-4]+'-dts.jpg',transparent=True)



