from classes.experiment_class import Experiment
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os

np.set_printoptions(precision=2)
plt.rcParams['figure.figsize'] = (20,10)
plt.rcParams['font.size'] = 15

#Recover the .txt files in given directory
dir = sys.argv[1]
file_paths = glob.glob(dir+'\*.txt')

experiments = list()
for file_path in file_paths:
    os.mkdir(file_path[:-4])
    experiments.append(Experiment(file_path))
    experiments[-1].plot_absorbance()
    experiments[-1].RTD_bestfit()
    experiments[-1].plot_bestfit()
    # os.replace(,os.path.join(dir,file_path[:-4],)