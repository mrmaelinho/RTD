from classes.experiment_class import Experiment
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
import subprocess
import pandas as pd
import matplotlib.colors as mcolors

np.set_printoptions(precision=2)
plt.rcParams['figure.figsize'] = (20,10)
plt.rcParams['font.size'] = 15
colours = list(mcolors.TABLEAU_COLORS)
colours.append('k')


#Recover all the .txt files in given directory
dir = sys.argv[1]
file_paths = glob.glob(dir+'\*.txt')

#Dataframe containing important parameters of each experiment
experiment_parameters = pd.DataFrame(columns=[\
                                              'file name',\
                                              'stages',\
                                              'internal volume (mL)',\
                                              'flowrate (mL/min)',\
                                              'tau (s)',\
                                              'speed (m/s)',\
                                              'aspect ratio',\
                                              'Re',\
                                              'Pe',\
                                              'Bodenstein',\
                                              'Bo (std)',\
                                              'Dax (m2/s)'])
#Go through each .txt file of the directory
for file_path in file_paths:
    experiment = Experiment(file_path)
    experiment.hydro_geom()
    experiment.plot_absorbance()
    experiment.RTD_bestfit()
    experiment.plot_bestfit()
    # Append the parameters to the DataFrame.
    experiment_parameters.loc[experiment] = [experiment.file_name,\
                                             experiment.stages,\
                                             int(experiment.stages)*0.125,\
                                             int(experiment.flowrate)/1000,\
                                             experiment.tau,\
                                             experiment.speed,\
                                             experiment.aspect_ratio,\
                                             experiment.Re,\
                                             experiment.Pe,\
                                             experiment.Bo_opt,\
                                             experiment.Bo_std,\
                                             experiment.Dax]
    plt.plot(experiment.speed,\
             experiment.Dax,\
             'o',\
             markersize = 4,\
             color = colours[int(experiment.stages)],\
             label = "%s stages"%experiment.stages,\
             )
    if np.log10(experiment.Pe) < np.log10(experiment.length/500e-6)+1:
        plt.plot(experiment.speed,\
                 experiment.Dax,\
                 'o',\
                 markersize = 10,\
                 color = colours[int(experiment.stages)],\
                 mfc='none')

    #Make a directory for this experiment.
    os.mkdir(file_path[:-4])
    #Place all the files corresponding to this experiment to the new directory.
    command = "move "+ file_path[:-4] +"* " + file_path[:-4] +"\\"
    subprocess.run(command,shell=True)

#Generate a csv file containing the data of all the analysed experiments.
experiment_parameters.to_csv(dir+'\\recap.csv', sep='\t')

##
#Taylor dispersion prediction
_speed = np.linspace(0,0.1,100)
_Taylor_Dax = _speed*_speed*25e-8/(30*1e-8)
plt.plot(_speed,_Taylor_Dax, label = 'Taylor theoretical trend')
# Display labels and grid
plt.xlabel('speed (m/s)')
plt.ylabel('Bo (m2/s)')
plt.grid()
#Deletes duplicates in legend
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(reversed(by_label.values()),\
           reversed(by_label.keys()))
#Display
plt.savefig(dir+'\\recap.jpg', transparent=True)
plt.close()
##
fig2, ax2 = plt.subplots()
for experiment in experiment_parameters.index:
    if experiment.tau == 60:
        ax2.plot(experiment.data.theta,\
                 experiment.data.absorbance/experiment.scale_opt,\
                 color = colours[int(experiment.stages)],\
                 label = "%s stages : Bo=%.1f"%(experiment.stages,experiment.Bo_opt))
ax2.set_xlabel('$\Theta$') #Reduced time absciss
ax2.set_xlim([0,2.5])
ax2.set_ylabel('E($\Theta$)')#Adimensional residence time distribution
plt.legend()
plt.savefig(dir+'\\tau-60s.jpg', transparent=True)