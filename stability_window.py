import os
import sys
from pathlib import Path

from pymatgen.core import *

from pymatgen.ext.matproj import *
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, GrandPotentialPhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility

##########################################################################
############### FOR MaterialsProject DATA VISUALISATION ##################
##########################################################################
from pymatgen.entries.computed_entries import ComputedStructureEntry

from pymatgen.io.vasp import Vasprun

################################ Custom addon #######################################

original_get_computed_entry = Vasprun.get_computed_entry


def get_computed_entry(self, *args, **kwargs):
    entry_id = kwargs.pop("entry_id", None)
    file_path = Path(self.filename)
    dir_name = file_path.parent.name
    return original_get_computed_entry(self, *args, **kwargs, entry_id=dir_name)


Vasprun.get_computed_entry = get_computed_entry
#####################################################################################


def set_energy_corrections(system):
    # These three lines assimilate the data into ComputedEntries.
    drone = VaspToComputedEntryDrone()
    queen = BorgQueen(drone, system, 3)
    entries = queen.get_data()

    # It's a good idea to perform a save_data, especially if you just assimilated
    # a large quantity of data which took some time. This allows you to reload
    # the data using a BorgQueen initialized with only the drone argument and
    # calling queen.load_data("Li-O_entries.json")

    # queen.save_data("Li-In-Cl_entries.json")


    #####################################################

    for entry in entries:
        entry.parameters["run_type"] = "GGA"

    
    compat = MaterialsProject2020Compatibility(correct_peroxide = False)  # sets energy corrections and +U/pseudopotential choice
    
    #processed_entries = compat.process_entries(entries, verbose=True)  # filter and add energy corrections
    processed_entries = compat.process_entries(entries)  # filter and add energy corrections
    file_path = Path(system)
    dir_name = file_path.name
    processed_entries[0].entry_id=dir_name
    return processed_entries[0]


def create_PhaseDiagram(processed_entries):
    pd = PhaseDiagram(processed_entries)

    return pd

def get_entries_list(dir_path):
    home_dir = Path(dir_path)
    error_list = []
    not_converged_list = []
    entries_list=[]
    for structure in [structures for structures in os.listdir(home_dir) if os.path.isdir(os.path.join(home_dir, structures))]:
        try:
            vrun = Vasprun(home_dir / structure / "vasprun.xml")
            print(structure)
            home_dir_i = Path(dir_path +'/'+structure+'/')
            #print(home_dir_i)
            a=set_energy_corrections(home_dir_i)
            entries_list.append(a)
            # vrun = Vasprun(home_dir / system / structure / "relax_uspx_man" / "step_lapsh" / "vasprun.xml")     # custom patch
        except (SyntaxError, FileNotFoundError):
            error_list.append(f"{structure}")
            pass
        except not vrun.converged:
            not_converged_list.append(f"{structure}")
            print(not_converged_list)
    return entries_list
    #    print(structure, vrun.initial_structure.composition)

entries_list=get_entries_list("/home/hamid/Research_papers/Computational/3-1-Na-La-I/Convex_hull_final/VASPRUN_SOURCE") # path for the calculation folder Li6PIO5_VdW

#pd = create_PhaseDiagram(set_energy_corrections(system=home_dir))
pd = create_PhaseDiagram(entries_list)
system='Na-La-I' 
# print(not_converged_list)

li_entries = [e for e in entries_list if e.composition.reduced_formula == "Na"]
uli0 = min(li_entries, key=lambda e: e.energy_per_atom).energy_per_atom
# print(uli0)

# Choose the structure for which stability window will be plotted
entry_LaCl_sg165=''
for i in pd.all_entries:
    if "NaLaI4_60" == i.entry_id:
        # print(i)
        entry_LaCl_sg165=i
        print(i.entry_id, i.energy_per_atom)

el_profile = pd.get_element_profile(Element("Na"), entry_LaCl_sg165.composition)
# print(el_profile)
for i, d in enumerate(el_profile):
    voltage = -(d["chempot"] - uli0)
    print("Voltage: %s V" % voltage)
    print(d["reaction"])
    print("")

import palettable
import matplotlib as mpl
#from pymatgen.util.plotting import pretty_plot
from matplotlib import pyplot as plt
import json
import re
# Some matplotlib settings to improve the look of the plot.
mpl.rcParams['axes.linewidth']=3
mpl.rcParams['lines.markeredgewidth']=4
mpl.rcParams['lines.linewidth']=3
mpl.rcParams['lines.markersize']=15
mpl.rcParams['xtick.major.width']=3
mpl.rcParams['xtick.major.size']=8
mpl.rcParams['xtick.minor.width']=3
mpl.rcParams['xtick.minor.size']=4
mpl.rcParams['ytick.major.width']=3
mpl.rcParams['ytick.major.size']=8
mpl.rcParams['ytick.minor.width']=3
mpl.rcParams['ytick.minor.size']=4


# Plot of Li uptake per formula unit (f.u.) of Li6PS5Cl against voltage vs Li/Li+.

colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
#plt = get_publication_quality_plot(12, 8)
#plt=pretty_plot()
#ax = plt.gca()

for i, d in enumerate(el_profile):
    v = - (d["chempot"] - uli0)
    if i != 0:
        # print(y1)
        plt.plot([x2, x2], [y1, d["evolution"] / 4.0], 'k', linewidth=3)
    x1 = v
    y1 = d["evolution"] / 4.0
    if i != len(el_profile) - 1:
        x2 = - (el_profile[i + 1]["chempot"] - uli0)
    else:
        x2 = 8.0
        
    if i in [1]: # Change for the indicating red colors remarks on the figure
        products = [re.sub(r"(\d+)", r"$_{\1}$", p.reduced_formula)                     
                    for p in d["reaction"].products if p.reduced_formula != "Na"]

        #plt.annotate(", ".join(products), xy=(v + 0.05, y1 + 0.3), 
                     #fontsize=15, color=colors[0])
        print(x1, x2, y1, sep = "-->")
       # plt.plot([x1, x2], [y1, y1], color=colors[0], linewidth=5)
    else:
        print(x1, x2, y1, sep = "-->")
        plt.plot([x1, x2], [y1, y1], 'k', linewidth=3)  
    if round(y1,3)==0:
        plt.axvspan(x1, x2, color="lightgreen", alpha=0.7)
        plt.annotate("Stability Region",xy=((x1+x2)/2,4),fontsize=14,ha='center',va='center')
        plt.annotate("%0.2f V - %0.2f V"%(x1,x2),xy=((x1+x2)/2,2.5),fontsize=14,ha='center',va='center')

plt.xlim([0, 8.0])
plt.ylim((-6, 10))
plt.xlabel("Voltage vs Na/Na$^+$ (V)",fontsize='14')
plt.ylabel("Na uptake per f.u.",fontsize='14')

#plt.tight_layout()   # ✅ call on the Figure
plt.savefig("NaLaI4.png", dpi=800)

