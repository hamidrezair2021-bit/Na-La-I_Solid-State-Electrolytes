import os
import sys
from pathlib import Path

from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility



##########################################################################
############### FOR MaterialsProject DATA VISUALISATION ##################
##########################################################################
from pymatgen.entries.computed_entries import ComputedStructureEntry
# from pymatgen.ext.matproj import MPRester
#
# MAPI_KEY = 'ojfgckRKZmkHFDLRRT73'  # You must change this to your Materials API key! (or set MAPI_KEY env variable)
# system = ["In", "Li", "I"]  # system we want to get PD for
# # system = ["Fe", "P", "O"]  # alternate system example: ternary
#
# mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface
# unprocessed_entries = mpr.get_entries_in_chemsys(system)

# compat = MaterialsProjectCompatibility()  # sets energy corrections and +U/pseudopotential choice
# processed_entries = compat.process_entries(unprocessed_entries)  # filter and add energy corrections

##########################################################################
##########################################################################
##########################################################################





##########################################################################
################ FOR LOCAL CALCULATIONS VISUALISATION ####################
##########################################################################
from pymatgen.io.vasp import Vasprun

"""
This script is for parsing local calculated data, making the corrections according to
 MaterialsProject data standards and recommendations and make a Convex Hull from them.
 It requires the path to the directory, where all the calculations are in separate
 directories, with names of them to be read further on the Convex Hull diagram.
 Be sure to have all the elements calculated in this system directory.

"""



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



    #####################################################S

    for entry in entries:
        entry.parameters["run_type"] = "GGA"



    compat = MaterialsProject2020Compatibility(correct_peroxide = False)  # sets energy corrections and +U/pseudopotential choice
    processed_entries = compat.process_entries(entries, verbose=True)  # filter and add energy corrections

    return processed_entries



def create_PhaseDiagram(processed_entries):
    pd = PhaseDiagram(processed_entries)

    return pd



if __name__ == '__main__':

    home_dir = Path("./Convex_hull_final/")

    syslist = [filename for filename in os.listdir(home_dir) if os.path.isdir(os.path.join(home_dir, filename))]

    

    # structure_types = ['sg255', 'o2', 'sg165', 'sg5', 'orth', 'sg33', 'trig', 'c2', 'sg162', 'sg14', 'o4', 'sg148',
    #                    'o1', 'sg118', 'sg9', 'sg163', 'sg15', 'sg4']

    for system in syslist:
        # ternary_report = open(home_dir / "ternary_report" / f"{system}_hull_report.txt", "w")
        # grand_report = open(home_dir / f"../convex_hull_report_{system}" / f"VdW_convex_hulls_GRAND_report.txt", "w")
        report = open(home_dir  / f"{system}_hull_report.txt", "w")
        report.write("################################### \n")
        report.write(f"{system} \n")
        report.write("################################### \n")

        error_list = []
        not_converged_list = []

        for structure in [structures for structures in os.listdir(home_dir / system) if os.path.isdir(os.path.join(home_dir / system, structures))]:
            try:
                vrun = Vasprun(home_dir / system / structure / "vasprun.xml")
            except (SyntaxError, FileNotFoundError):
                error_list.append(f"{structure}")
                pass
            except not vrun.converged:
                not_converged_list.append(f"{structure}")

            print(structure, vrun.initial_structure.composition)

        # print(system)
        pd = set_energy_corrections(system=home_dir / system)
        print(pd)
        
        report.write("-------- \n")
        report.write(f"Error parsing:   {error_list} \n")
        report.write("-------- \n")
        report.write(f"Not converged:   {not_converged_list} \n")
        report.write("-------- \n")
        report.write("       Entry                        E above hull                   Energy                  energy_per_atom                  Correction                  correction_per_atom   \n")

        a = create_PhaseDiagram(pd)
        for entry in pd:
            report.write(
            f"      {entry.entry_id}                  {round(create_PhaseDiagram(pd).get_e_above_hull(entry),4)}              {round(entry.energy, 4)}        {round(entry.energy_per_atom, 4)}        {round(entry.correction, 4)}        {round(entry.correction_per_atom, 4)}\n")

    #     # Plot
        plotter = PDPlotter(a, show_unstable=True) #backend="matplotlib")  # you can also try show_unstable=False
        plotter.show()
    # plotter.write_image('Na-La-F_hull.png', image_format="png")




    #grand_report.close()
