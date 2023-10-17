
from pynwb import NWBHDF5IO
#from nwbwidgets import nwb2widget
import h5py, glob, os
import numpy as np
import matplotlib.pyplot as plt
from dlab import sglx_analysis as sglx
from dlab import psth_and_raster as pr
import numpy as np
import pandas as pd
import os
from tqdm.notebook import tqdm as tqdm
import glob
import pickle as pkl
import json
import datetime as dt
import h5py
from dlab import sglx_analysis as sglx
import math
# import jlh_ephys_tools as jlh
import matplotlib.pyplot as plt
# from open_ephys.analysis import Session
import seaborn as sns
import dlab.psth_and_raster as psth
from dlab.generalephys import cleanAxes
import matplotlib.lines as mlines

import statsmodels.api as sm
from statsmodels.formula.api import mixedlm




class Eopn3_Ephys:

#  ----------------- Load data -----------------------------
    def __init__(self,nwb_path):
       self.nwb_path = nwb_path
       self.load_nwb()
 
    def nwb(self):
        if self.nwb is None:
            self.load_nwb()
        return self.nwb

    def load_nwb(self):
        nwb_path = self.nwb_path
        io = NWBHDF5IO(nwb_path, 'r')
        nwb = io.read()
        self.nwb = nwb
        return self.nwb
    
# prints the metaData associated with nwb
    def view_nwb(self):
        nwb = self.nwb
        print(nwb)

# creates a df of the trails associated with nwb. the Data frame contains data about each trail and its structure 
    def trials(self):
        nwb = self.nwb
        df_stim = nwb.trials.to_dataframe()
        # df_stim.loc[2100:2699,'contacts'] = '10r'  #specific to this recording, fixes an error in dataframe
        self.df_stim = df_stim
        return self.df_stim

# creates the Units Data frame associated with nwb. which contains all the data about the sorted spike units from the recording
    def units(self):
        nwb = self.nwb
        df_units = nwb.units.to_dataframe()
        self.df_units = df_units
        return self.df_units
    
# creates the optogenetics_states Data frame associated with nwb. which contains all the data about the optogenetics states 
    def optogenetics_states(self):
        nwb = self.nwb
        optogenetics_states_df = nwb.intervals['optogenetics_states'].to_dataframe()
        self.optogenetics_states_df = optogenetics_states_df
        return optogenetics_states_df
    
# creates the epochs Data frame associated with nwb. which contains just the entire recording length. 
    def epochs(self):
        nwb = self.nwb
        epochs_df = nwb.intervals['epochs'].to_dataframe()
        self.epochs_df = epochs_df
        return epochs_df

#  ----------------- End Load data -----------------------------


# ------------------------ Plotting ---------------------------------
    def singleUnit_flash(self,probeLetter,unit_indexNumber,Savefig=False,):
        df_stim = self.trials()
        df_units = self.units() 

        df1 = df_units[df_units.probe==probeLetter][df_units.label==2]
        unit_index = df1.index.tolist()

        print(f"neuron indexs for {probeLetter} {list(unit_index)}")
       
        unit = unit_index[unit_indexNumber]
        print(f"looking at neuron/unit {unit}")
        f,ax=plt.subplots(1,1)
         
            # Set the background color and text color for this subplot
        ax.set_facecolor('black')
        ax.tick_params(color='black', labelcolor='black')
        for spine in ax.spines.values():
            spine.set_edgecolor('black')

        psth.psth_line(times=df1[df1.index==unit].spike_times.values[0],
                triggers=df_stim[(df_stim.stimulus=='luminance_flash') & (df_stim.optogenetics_LED_state == 0)].start_time.values,
                    ymax=60,binsize=0.05,axes=ax,color='#487697')
        psth.psth_line(times=df1[df1.index==unit].spike_times.values[0],
                triggers=df_stim[(df_stim.stimulus=='luminance_flash_opto') ].start_time.values,
                    ymax=60,binsize=0.05,axes=ax,color='#ffaa00')
        plt.tight_layout()
        if Savefig == True:
            f.savefig('/Users/danieljdenman/Academics/grants/applications/20230305_R01_NEI_resub/figures/eArch_LM_psth.eps')
            f.savefig('/Users/danieljdenman/Academics/grants/applications/20230305_R01_NEI_resub/figures/eArch_LM_psth.png')
        else:
            pass

        f,ax=plt.subplots(2,1)
            # Set the background color and text color for these subplots
        for ax_ in ax:
            ax_.set_facecolor('black')
            ax_.tick_params(color='black', labelcolor='black')
            for spine in ax_.spines.values():
                spine.set_edgecolor('black')
        
        psth.raster(times=df1[df1.index==unit].spike_times.values[0],
                triggers=df_stim[(df_stim.stimulus=='luminance_flash') & (df_stim.optogenetics_LED_state == 0)].start_time.values,
                    axes=ax[0],color='#487697',timeDomain=True,post=1.5,ms=8)
        psth.raster(times=df1[df1.index==unit].spike_times.values[0],
                triggers=df_stim[(df_stim.stimulus=='luminance_flash_opto') ].start_time.values,
                    axes=ax[1],color='#ffaa00',timeDomain=True,post=1.5,ms=8)
        for ax_ in ax: ax_.set_xlim(-0.5,1.0)
     
        plt.tight_layout()
        # f.savefig('/Users/danieljdenman/Academics/grants/applications/20230305_R01_NEI_resub/figures/eArch_LM_raster.eps')

# Save raters and psth_line plots for every unit in the selected probeLetter. Supply a unique save path for eachprobe 
    def allUnits_flash(self, probeLetter, brain_region ,path ,Savefig=False):
        df_stim = self.trials()
        df_units = self.units()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

        for unit in unit_index:
            # Create a single figure with 3 subplots
            # f, ax = plt.subplots(3, 1)
            f, ax = plt.subplots(3, 1, figsize=(10, 15))

            f.suptitle(f"aligned flash events; unit {unit}, probe {probeLetter} in {brain_region}", color='black')

            # Plot the psth_line in the first subplot
            ax[0].set_facecolor('black')
            ax[0].tick_params(color='black', labelcolor='black')
            for spine in ax[0].spines.values():
                spine.set_edgecolor('black')
            psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                        triggers=df_stim[(df_stim.stimulus == 'luminance_flash') & (df_stim.optogenetics_LED_state == 0)].start_time.values,
                        ymax=60, binsize=0.05, axes=ax[0], color='#ffaa00')
            
            psth.psth_line(times=df1[df1.index==unit].spike_times.values[0],
                        triggers=df_stim[(df_stim.stimulus=='luminance_flash_opto') ].start_time.values,
                        ymax=60,binsize=0.05,axes=ax[0],color='#487697')
            ax[0].set_xlim(-0.5, 1.0)
         
            # Create custom legend handles
            handle1 = mlines.Line2D([], [], color='#ffaa00', label='Non-Opto')
            handle2 = mlines.Line2D([], [], color='#487697', label='Opto')

            # Add legend to the first subplot
            ax[0].legend(handles=[handle1, handle2], loc='upper right')
  

            # Plot the rasters in the second and third subplots
            for i, ax_ in enumerate(ax[1:]):
                ax_.set_facecolor('black')
                ax_.tick_params(color='black', labelcolor='black')
                for spine in ax_.spines.values():
                    spine.set_edgecolor('black')

            psth.raster(times=df1[df1.index == unit].spike_times.values[0],
                        triggers=df_stim[(df_stim.stimulus == 'luminance_flash') & (df_stim.optogenetics_LED_state == 0)].start_time.values,
                        axes=ax[1], color='#ffaa00', timeDomain=True, post=1.5, ms=8)

            psth.raster(times=df1[df1.index == unit].spike_times.values[0],
                        triggers=df_stim[(df_stim.stimulus == 'luminance_flash_opto')].start_time.values,
                        axes=ax[2], color='#487697', timeDomain=True, post=1.5, ms=8)
            ax[1].set_xlim(-0.5, 1.0)
            ax[2].set_xlim(-0.5, 1.0)
            ax[1].set_title('non-optp', color='black')
            ax[2].set_title('opto', color='black')


            plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to make room for the suptitle

            if Savefig:
                f.savefig(f"{path}{unit}.png")

            plt.close(f)  # Close the figure to free up memory

    def allUnits_flash_conjoined(self, probeLetter, brain_region, path, Savefig=False):
        df_stim = self.trials()
        df_units = self.units()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

        n_units_per_image = 20  # Number of units per image
        n_cols = 4  # Number of columns per image
        n_rows = 5  # Number of rows per image

        for i in range(0, len(unit_index), n_units_per_image):
            f, axarr = plt.subplots(n_rows, n_cols, figsize=(15, 15))  # Create a new figure for every 20 units

            f.suptitle(f"Aligned flash events for probe {probeLetter} in {brain_region}", color='black')

            for j in range(n_units_per_image):
                if i + j >= len(unit_index):
                    break  # Exit the loop if we've reached the end of the unit list

                unit = unit_index[i + j]
                row = j // n_cols
                col = j % n_cols
                ax = axarr[row, col]

                ax.set_facecolor('black')
                ax.tick_params(color='black', labelcolor='black', labelsize=4)
                for spine in ax.spines.values():
                    spine.set_edgecolor('black')

                
                psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=df_stim[(df_stim.stimulus == 'luminance_flash') & (df_stim.optogenetics_LED_state == 0)].start_time.values,
                            ymax=60, binsize=0.05, axes=ax, color='#ffaa00')

                psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=df_stim[(df_stim.stimulus == 'luminance_flash_opto')].start_time.values,
                            ymax=60, binsize=0.05, axes=ax, color='#487697')

                ax.set_xlim(-0.5, 1.0)

                ax.set_title(f"Unit {unit}", color='black', fontsize=8)
                ax.set_xlabel('Time', fontsize=8)
                ax.set_ylabel('Firing Rate', fontsize=8)

            # Add a legend to the figure
            handle1 = mlines.Line2D([], [], color='#ffaa00', label='Non-Opto')
            handle2 = mlines.Line2D([], [], color='#487697', label='Opto')
            f.legend(handles=[handle1, handle2], loc='upper right')

            plt.tight_layout(rect=[0, 0.03, 1, 0.95])

            if Savefig:
                f.savefig(f"{path}units_{i+1}_to_{i+n_units_per_image}.png")

            plt.close(f)  # Close the figure to free up memory

    def allUnits_flash_epoch_psthLined(self, probeLetter, brain_region, path, Savefig=False):
        df_units = self.units()
        df_epoch = self.epochs()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

        epoch_pairs = [
            ('luminance_flash', 'luminance_flash_opto'),
            # ('spatioluminance_noise', 'spatioluminance_noise_opto'),
            # ('gratings', 'gratings_opto'),
            # ('scenes', 'scenes_opto')
        ]

        for unit in unit_index:
            for epoch_name1, epoch_name2 in epoch_pairs:
                f, ax = plt.subplots(1, 1, figsize=(10, 5))
                f.suptitle(f"Epochs {epoch_name1} & {epoch_name2}; unit {unit}, probe {probeLetter} in {brain_region}", color='black')

                ax.set_facecolor('black')
                ax.tick_params(color='black', labelcolor='black')
                for spine in ax.spines.values():
                    spine.set_edgecolor('black')

                for epoch_name, color in zip([epoch_name1, epoch_name2], ['#ffaa00', '#487697']):
                    epoch_row = df_epoch[df_epoch['tags'].apply(lambda x: epoch_name in x)].iloc[0]
                    start_time, stop_time = epoch_row['start_time'], epoch_row['stop_time']
                    bin_width = (stop_time - start_time) / 1000  # Calculate bin width

                    unit_spike_times = df1[df1.index == unit].spike_times.values[0]
                    unit_spike_times_epoch = unit_spike_times[(unit_spike_times >= start_time) & (unit_spike_times <= stop_time)]

                    if len(unit_spike_times_epoch) > 1:
                        time_vector = np.linspace(0, stop_time - start_time, 1000)
                        kde = gaussian_kde(unit_spike_times_epoch - start_time)
                        firing_rate = kde.evaluate(time_vector) * len(unit_spike_times_epoch) / bin_width  # Convert to firing rate
                        ax.plot(time_vector, firing_rate, color=color, label=f"{epoch_name}")
                    else:
                        print(f"Skipping unit {unit} for epoch {epoch_name} due to insufficient spike times.")

                ax.legend(loc='upper right')
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])

                if Savefig:
                    f.savefig(f"{path}{unit}_{epoch_name1}_{epoch_name2}.png")

                plt.close(f)


    def allUnit_gratings(self, probeLetter, brain_region ,path ,Savefig=False):
        df_stim = self.trials()
        df_units = self.units()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

        for unit in unit_index:
            # Create a single figure with 3 subplots
            # f, ax = plt.subplots(3, 1)
            f, ax = plt.subplots(3, 1, figsize=(10, 15))

            f.suptitle(f"aligned grating events; unit {unit}, probe {probeLetter} in {brain_region}", color='black')

            # Plot the psth_line in the first subplot
            ax[0].set_facecolor('black')
            ax[0].tick_params(color='black', labelcolor='black')
            for spine in ax[0].spines.values():
                spine.set_edgecolor('black')
            psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                        triggers=df_stim[(df_stim.stimulus == 'gratings') & (df_stim.optogenetics_LED_state == 0)].start_time.values,
                        ymax=60, binsize=0.05, axes=ax[0], color='#ffaa00')
            
            psth.psth_line(times=df1[df1.index==unit].spike_times.values[0],
                        triggers=df_stim[(df_stim.stimulus=='gratings_opto') ].start_time.values,
                        ymax=60,binsize=0.05,axes=ax[0],color='#487697')
            ax[0].set_xlim(-0.5, 1.0)
         
            # Create custom legend handles
            handle1 = mlines.Line2D([], [], color='#ffaa00', label='Non-Opto')
            handle2 = mlines.Line2D([], [], color='#487697', label='Opto')

            # Add legend to the first subplot
            ax[0].legend(handles=[handle1, handle2], loc='upper right')
  

            # Plot the rasters in the second and third subplots
            for i, ax_ in enumerate(ax[1:]):
                ax_.set_facecolor('black')
                ax_.tick_params(color='black', labelcolor='black')
                for spine in ax_.spines.values():
                    spine.set_edgecolor('black')

            psth.raster(times=df1[df1.index == unit].spike_times.values[0],
                        triggers=df_stim[(df_stim.stimulus == 'gratings') & (df_stim.optogenetics_LED_state == 0)].start_time.values,
                        axes=ax[1], color='#ffaa00', timeDomain=True, post=1.5, ms=8)

            psth.raster(times=df1[df1.index == unit].spike_times.values[0],
                        triggers=df_stim[(df_stim.stimulus == 'gratings_opto')].start_time.values,
                        axes=ax[2], color='#487697', timeDomain=True, post=1.5, ms=8)
            ax[1].set_xlim(-0.5, 1.0)
            ax[2].set_xlim(-0.5, 1.0)
            ax[1].set_title('non-optp', color='black')
            ax[2].set_title('opto', color='black')


            plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to make room for the suptitle

            if Savefig:
                f.savefig(f"{path}{unit}.png")
     
            plt.close(f)  # Close the figure to free up memory
            # Uncomment the next line to stop the loop after one iteration

    def allUnits_gratings_conjoined(self, probeLetter, brain_region, path, Savefig=False):
        df_stim = self.trials()
        df_units = self.units()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

        n_units_per_image = 20  # Number of units per image
        n_cols = 4  # Number of columns per image
        n_rows = 5  # Number of rows per image

        for i in range(0, len(unit_index), n_units_per_image):
            f, axarr = plt.subplots(n_rows, n_cols, figsize=(15, 15))  # Create a new figure for every 20 units

            f.suptitle(f"Aligned flash events for probe {probeLetter} in {brain_region}", color='black')

            for j in range(n_units_per_image):
                if i + j >= len(unit_index):
                    break  # Exit the loop if we've reached the end of the unit list

                unit = unit_index[i + j]
                row = j // n_cols
                col = j % n_cols
                ax = axarr[row, col]

                ax.set_facecolor('black')
                ax.tick_params(color='black', labelcolor='black', labelsize=4)
                for spine in ax.spines.values():
                    spine.set_edgecolor('black')

                
                psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=df_stim[(df_stim.stimulus == 'gratings') & (df_stim.optogenetics_LED_state == 0)].start_time.values,
                            ymax=60, binsize=0.05, axes=ax, color='#ffaa00')

                psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=df_stim[(df_stim.stimulus == 'gratings_opto')].start_time.values,
                            ymax=60, binsize=0.05, axes=ax, color='#487697')

                ax.set_xlim(-0.5, 1.0)

                ax.set_title(f"Unit {unit}", color='black', fontsize=8)
                ax.set_xlabel('Time', fontsize=8)
                ax.set_ylabel('Firing Rate', fontsize=8)

            # Add a legend to the figure
            handle1 = mlines.Line2D([], [], color='#ffaa00', label='Non-Opto')
            handle2 = mlines.Line2D([], [], color='#487697', label='Opto')
            f.legend(handles=[handle1, handle2], loc='upper right')

            plt.tight_layout(rect=[0, 0.03, 1, 0.95])

            if Savefig:
                f.savefig(f"{path}units_{i+1}_to_{i+n_units_per_image}.png")

            plt.close(f)  # Close the figure to free up memory

    def allUnits_grating_epoch_psthLined(self, probeLetter, brain_region, path, Savefig=False):
        df_units = self.units()
        df_epoch = self.epochs()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

        epoch_pairs = [
            # ('luminance_flash', 'luminance_flash_opto'),
            # ('spatioluminance_noise', 'spatioluminance_noise_opto'),
            ('gratings', 'gratings_opto'),
            # ('scenes', 'scenes_opto')
        ]

        for unit in unit_index:
            for epoch_name1, epoch_name2 in epoch_pairs:
                f, ax = plt.subplots(1, 1, figsize=(10, 5))
                f.suptitle(f"Epochs {epoch_name1} & {epoch_name2}; unit {unit}, probe {probeLetter} in {brain_region}", color='black')

                ax.set_facecolor('black')
                ax.tick_params(color='black', labelcolor='black')
                for spine in ax.spines.values():
                    spine.set_edgecolor('black')

                for epoch_name, color in zip([epoch_name1, epoch_name2], ['#ffaa00', '#487697']):
                    epoch_row = df_epoch[df_epoch['tags'].apply(lambda x: epoch_name in x)].iloc[0]
                    start_time, stop_time = epoch_row['start_time'], epoch_row['stop_time']
                    bin_width = (stop_time - start_time) / 1000  # Calculate bin width

                    unit_spike_times = df1[df1.index == unit].spike_times.values[0]
                    unit_spike_times_epoch = unit_spike_times[(unit_spike_times >= start_time) & (unit_spike_times <= stop_time)]

                    if len(unit_spike_times_epoch) > 1:
                        time_vector = np.linspace(0, stop_time - start_time, 1000)
                        kde = gaussian_kde(unit_spike_times_epoch - start_time)
                        firing_rate = kde.evaluate(time_vector) * len(unit_spike_times_epoch) / bin_width  # Convert to firing rate
                        ax.plot(time_vector, firing_rate, color=color, label=f"{epoch_name}")
                    else:
                        print(f"Skipping unit {unit} for epoch {epoch_name} due to insufficient spike times.")

                ax.legend(loc='upper right')
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])

                if Savefig:
                    f.savefig(f"{path}{unit}_{epoch_name1}_{epoch_name2}.png")

                plt.close(f)



    def allUnit_scenes(self, probeLetter, brain_region, path, Savefig=False):
        df_stim = self.trials()
        df_units = self.units()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

    
        # Filter df_stim to only include rows where stimulus is 'scenes' or 'scenes_opto'
        df_stim_filtered = df_stim[df_stim['stimulus'].isin(['scenes', 'scenes_opto'])]

        unique_stimulus_indices = df_stim_filtered['stimulus_index'].unique()

        for unit in unit_index:
            for stim_index in unique_stimulus_indices:
                f, ax = plt.subplots(3, 1, figsize=(10, 15))

                f.suptitle(f"aligned scene events; unit {unit}, probe {probeLetter} in {brain_region}, stimulus_index {stim_index}", color='black')

                # Plot the psth_line in the first subplot
                ax[0].set_facecolor('black')
                ax[0].tick_params(color='black', labelcolor='black')
                for spine in ax[0].spines.values():
                    spine.set_edgecolor('black')

                # ... (rest of your code remains the same, but filter triggers by stimulus_index)
                
                triggers_non_opto = df_stim_filtered[(df_stim_filtered.stimulus == 'scenes') & (df_stim_filtered.optogenetics_LED_state == 0) & (df_stim_filtered.stimulus_index == stim_index)].start_time.values
                triggers_opto = df_stim_filtered[(df_stim_filtered.stimulus == 'scenes_opto') & (df_stim_filtered.stimulus_index == stim_index)].start_time.values

                psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=triggers_non_opto,
                            ymax=60, binsize=0.05, axes=ax[0], color='#ffaa00')
                
                psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=triggers_opto,
                            ymax=60, binsize=0.05, axes=ax[0], color='#487697')
                
                ax[0].set_xlim(-0.5, 1.0)

                    
                # Create custom legend handles
                handle1 = mlines.Line2D([], [], color='#ffaa00', label='Non-Opto')
                handle2 = mlines.Line2D([], [], color='#487697', label='Opto')

                # Add legend to the first subplot
                ax[0].legend(handles=[handle1, handle2], loc='upper right')
    

                   # Plot the rasters in the second and third subplots
                for i, ax_ in enumerate(ax[1:]):
                    ax_.set_facecolor('black')
                    ax_.tick_params(color='black', labelcolor='black')
                    for spine in ax_.spines.values():
                        spine.set_edgecolor('black')
                    
                # ... (rest of your code remains the same, but filter triggers by stimulus_index)

                psth.raster(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=triggers_non_opto,
                            axes=ax[1], color='#ffaa00', timeDomain=True, post=1.5, ms=8)

                psth.raster(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=triggers_opto,
                            axes=ax[2], color='#487697', timeDomain=True, post=1.5, ms=8)
                
                ax[1].set_xlim(-0.5, 1.0)
                ax[2].set_xlim(-0.5, 1.0)
                ax[1].set_title('non-optp', color='black')
                ax[2].set_title('opto', color='black')

                # ... (rest of your code for saving and closing figures)
                plt.tight_layout(rect=[0, 0.03, 1, 0.92])  # Adjust the top value

                if Savefig:
                    f.savefig(f"{path}{unit}_stimulus_index_{stim_index}.png")

                plt.close(f)

    def allUnit_scenes_conjoined(self, probeLetter, brain_region, path, Savefig=False):
        df_stim = self.trials()
        df_units = self.units()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

        df_stim_filtered = df_stim[df_stim['stimulus'].isin(['scenes', 'scenes_opto'])]
        unique_stimulus_indices = df_stim_filtered['stimulus_index'].unique()

        for unit in unit_index:
            f1, ax1 = plt.subplots(4, 5, figsize=(25, 30))  # Increase the size

            f1.suptitle(f"aligned scene events (PSTH); unit {unit}, probe {probeLetter} in {brain_region}", fontsize=16, color='black')

            f2, ax2 = plt.subplots(10, 4, figsize=(20, 32))  # 8 rows, 5 columns
            f2.suptitle(f"aligned scene events (Raster); unit {unit}, probe {probeLetter} in {brain_region}", fontsize=16, color='black')

            ax1 = ax1.flatten()
            ax2 = ax2.reshape(-1, 2)  # Reshape to 20 pairs of 2 subplots

            for i, stim_index in enumerate(unique_stimulus_indices):
                ax1[i].set_facecolor('black')
                ax1[i].tick_params(color='black', labelcolor='black')
                for spine in ax1[i].spines.values():
                    spine.set_edgecolor('black')

                for ax in ax2[i]:
                    ax.set_facecolor('black')
                    ax.tick_params(color='black', labelcolor='black')
                    for spine in ax.spines.values():
                        spine.set_edgecolor('black')

                triggers_non_opto = df_stim_filtered[(df_stim_filtered.stimulus == 'scenes') & (df_stim_filtered.optogenetics_LED_state == 0) & (df_stim_filtered.stimulus_index == stim_index)].start_time.values
                triggers_opto = df_stim_filtered[(df_stim_filtered.stimulus == 'scenes_opto') & (df_stim_filtered.stimulus_index == stim_index)].start_time.values

                # PSTH plots
                psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                                triggers=triggers_non_opto,
                                ymax=60, binsize=0.05, axes=ax1[i], color='#ffaa00')
                psth.psth_line(times=df1[df1.index == unit].spike_times.values[0],
                                triggers=triggers_opto,
                                ymax=60, binsize=0.05, axes=ax1[i], color='#487697')

                # Raster plots (non-opto)
                psth.raster(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=triggers_non_opto,
                            axes=ax2[i, 0], color='#ffaa00', timeDomain=True, post=1.5, ms=8)

                # Raster plots (opto)
                psth.raster(times=df1[df1.index == unit].spike_times.values[0],
                            triggers=triggers_opto,
                            axes=ax2[i, 1], color='#487697', timeDomain=True, post=1.5, ms=8)

                ax1[i].set_title(f"Index {stim_index}", color='black')
                ax2[i, 0].set_title(f"Index {stim_index} Non-Opto", color='black')
                ax2[i, 1].set_title(f"Index {stim_index} Opto", color='black')

                ax1[i].set_xlim(-0.5, 1.0)
                ax2[i, 0].set_xlim(-0.5, 1.0)
                ax2[i, 1].set_xlim(-0.5, 1.0)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])

            if Savefig:
                f1.savefig(f"{path}{unit}_all_indices_psth.png")
                f2.savefig(f"{path}{unit}_all_indices_raster.png")

            plt.close(f1)
            plt.close(f2)
            
    def allUnits_scene_epoch_psthLined(self, probeLetter, brain_region, path, Savefig=False):
        df_units = self.units()
        df_epoch = self.epochs()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()

        epoch_pairs = [
            # ('luminance_flash', 'luminance_flash_opto'),
            # ('spatioluminance_noise', 'spatioluminance_noise_opto'),
            # ('gratings', 'gratings_opto'),
            ('scenes', 'scenes_opto')
        ]

        for unit in unit_index:
            for epoch_name1, epoch_name2 in epoch_pairs:
                f, ax = plt.subplots(1, 1, figsize=(10, 5))
                f.suptitle(f"Epochs {epoch_name1} & {epoch_name2}; unit {unit}, probe {probeLetter} in {brain_region}", color='black')

                ax.set_facecolor('black')
                ax.tick_params(color='black', labelcolor='black')
                for spine in ax.spines.values():
                    spine.set_edgecolor('black')

                for epoch_name, color in zip([epoch_name1, epoch_name2], ['#ffaa00', '#487697']):
                    epoch_row = df_epoch[df_epoch['tags'].apply(lambda x: epoch_name in x)].iloc[0]
                    start_time, stop_time = epoch_row['start_time'], epoch_row['stop_time']
                    bin_width = (stop_time - start_time) / 1000  # Calculate bin width

                    unit_spike_times = df1[df1.index == unit].spike_times.values[0]
                    unit_spike_times_epoch = unit_spike_times[(unit_spike_times >= start_time) & (unit_spike_times <= stop_time)]

                    if len(unit_spike_times_epoch) > 1:
                        time_vector = np.linspace(0, stop_time - start_time, 1000)
                        kde = gaussian_kde(unit_spike_times_epoch - start_time)
                        firing_rate = kde.evaluate(time_vector) * len(unit_spike_times_epoch) / bin_width  # Convert to firing rate
                        ax.plot(time_vector, firing_rate, color=color, label=f"{epoch_name}")
                    else:
                        print(f"Skipping unit {unit} for epoch {epoch_name} due to insufficient spike times.")

                ax.legend(loc='upper right')
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])

                if Savefig:
                    f.savefig(f"{path}{unit}_{epoch_name1}_{epoch_name2}.png")

                plt.close(f)        
# ------------------------ End Plotting ---------------------------------

# ------------------------ Analysis ---------------------------------
    def extract_spike_times_in_event_windows_flash(self, probeLetter):
        df_stim = self.trials()
        df_units = self.units()

        df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
        unit_index = df1.index.tolist()
       
        unit_event_spikes = {}  # Dictionary to store spike times for each unit

        for unit in unit_index:
            unit_spikes = df1[df1.index == unit].spike_times.values[0]
            
            unit_event_spikes[unit] = []

            # Filter df_stim to only include rows where stimulus is 'luminance_flash' or 'luminance_flash_opto'
            df_stim_filtered = df_stim[df_stim['stimulus'].isin(['luminance_flash', 'luminance_flash_opto'])]
           
           
            for _, row in df_stim_filtered.iterrows():
                event_start = row['start_time'] 
                event_end = row['stop_time']

                # event_end = event_start + 2.0  # Assuming a 1.5-second window; adjust as needed
                print(event_start,event_end)
                # Extract spikes that fall within this event window
                event_spikes = unit_spikes[(unit_spikes >= event_start) & (unit_spikes <= event_end)]
                
                unit_event_spikes[unit].append(event_spikes.tolist())

         
        # Convert the dictionary to a DataFrame
        df_event_spikes = pd.DataFrame.from_dict(unit_event_spikes, orient='index')

        # If the number of event windows is not the same for each unit, you might need to pad with None or some other value
        df_event_spikes = df_event_spikes.apply(lambda x: pd.Series(x.dropna().values.tolist()), axis=1)

        return df_event_spikes

    def allUnits_full_recording_psthLined(self, probeLetter, brain_region, path, Savefig=False):
            df_units = self.units()
            df_epoch = self.epochs()

            df1 = df_units[df_units.probe == probeLetter][df_units.label == 2]
            unit_index = df1.index.tolist()

        
            # Get the first start_time and the last stop_time from df_epoch
            start_time = df_epoch['start_time'].iloc[0]
            stop_time = df_epoch['stop_time'].iloc[-1]
            print(start_time,stop_time)
            # Get the last 4 timestamps from the 'tags' column in df_epoch
            inhibitory_flash_times = df_epoch['start_time'].tail(4).tolist()
            print(inhibitory_flash_times)
            epoch_start_times = df_epoch['start_time'].tolist()
            epoch_stop_times = df_epoch['stop_time'].tolist()
            epochs = len(epoch_stop_times)
            

            print(inhibitory_flash_times)

            for unit in unit_index:
                f, ax = plt.subplots(1, 1, figsize=(10, 5))
                f.suptitle(f"Unit {unit}, probe {probeLetter} in {brain_region}", color='black')

                ax.set_facecolor('black')
                ax.tick_params(color='black', labelcolor='black')
                for spine in ax.spines.values():
                    spine.set_edgecolor('black')
                
                # Add red lines for inhibitory light flashes
                for flash_time in inhibitory_flash_times:
                    if type(flash_time) is float and start_time <= flash_time <= stop_time:  # Make sure it's a timestamp
                        ax.axvline(x=flash_time - start_time, color='red', linestyle='--')

                
                    # Add shaded regions for each epoch
                for epoch in range(epochs):
                    print(epoch_start_times[epoch], epoch_stop_times[epoch])
                    colors = ['yellow', 'green', 'blue', 'purple', 'yellow', 'green', 'blue', 'purple',]
                    ax.axvspan(epoch_start_times[epoch], epoch_stop_times[epoch], facecolor=colors[epoch], alpha=0.5)

                bin_width = (stop_time - start_time) / 1000  # Calculate bin width

                unit_spike_times = df1[df1.index == unit].spike_times.values[0]
                unit_spike_times_epoch = unit_spike_times[(unit_spike_times >= start_time) & (unit_spike_times <= stop_time)]

                if len(unit_spike_times_epoch) > 1:
                    time_vector = np.linspace(0, stop_time - start_time, 1000)
                    kde = gaussian_kde(unit_spike_times_epoch - start_time)
                    firing_rate = kde.evaluate(time_vector) * len(unit_spike_times_epoch) / bin_width  # Convert to firing rate
                    ax.plot(time_vector, firing_rate, color='#ffaa00')

                else:
                    print(f"Skipping unit {unit} due to insufficient spike times.")

                plt.tight_layout(rect=[0, 0.03, 1, 0.95])

                if Savefig:
                    f.savefig(f"{path}{unit}_full_recording.png")

                plt.close(f)



# # Assuming your data is in a pandas DataFrame named df with columns 'spikes', 'light', 'stimulus_type', and 'mouse_id'
# model = mixedlm("spikes ~ light * stimulus_type", df, groups=df['mouse_id']).fit()
# print(model.summary())

