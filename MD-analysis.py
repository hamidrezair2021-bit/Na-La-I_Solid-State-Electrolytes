from pymatgen.io.lammps.data import *
from pymatgen.io.lammps.outputs import *
import numpy as np
from pymatgen.analysis.diffusion.aimd.pathway import ProbabilityDensityAnalysis
from pymatgen.analysis.diffusion.analyzer import (
    DiffusionAnalyzer,
    get_arrhenius_plot,
    get_extrapolated_conductivity,
)
from pymatgen.core import Structure
from ase.data import atomic_numbers
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read
import sys
from pathlib import Path

from pymatgen.core.trajectory import Trajectory
import matplotlib.pyplot as plt


def convert_lammpsdata_to_pymatgen(path_to_lammpsdump):
    # reads lammps dump and returns pymatgen Structure object list
    trajectory = []
    data_trj=read(path_to_lammpsdump,index=":")
# for frame_index in range(pipeline.source.num_frames):
    for ase_data in data_trj:
        numbers=ase_data.get_atomic_numbers().tolist()
        modified_list = ["Na" if item == 1 else "La" if item == 2 else "I" if item == 3 else item for item in numbers]
        new_numbers=[atomic_numbers[i] for i in modified_list]
        ase_data.set_atomic_numbers(new_numbers)

        struct = AseAtomsAdaptor.get_structure(ase_data)
        struct.sort()
        trajectory.append(struct)

    return trajectory

home_dir = Path("./")
tr=convert_lammpsdata_to_pymatgen('md.lammpstrj')

a = DiffusionAnalyzer.from_structures(structures=tr, specie="Na", temperature=300, time_step=1,
                                      step_skip=100, smoothed=False)

a.plot_msd()
a.plot_msd(mode="species")
a.get_msd_plot().figure.savefig(home_dir / "MSD")
a.get_msd_plot(mode="species").figure.savefig(home_dir / "MSD_species")


# Opening the log
global_orig_stdout = sys.stdout
out = open(home_dir / "conductivity_report", 'w')
sys.stdout = out

print(f"conductivity     {a.conductivity}")
print(f"chg conductivity       {a.chg_conductivity}")
print(f"diffusivity       {a.diffusivity}")
print(a.get_summary_dict())


# Closing the log
sys.stdout = global_orig_stdout
out.close()
