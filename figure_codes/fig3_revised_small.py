import pickle

import matplotlib.gridspec as gridspec
import matplotlib
from ..inverse_hull_web import InverseHullWeb, IHWPlotter
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations
import json
from ..main import main

matplotlib.rcParams['font.size'] = 9
matplotlib.rcParams['font.family'] = 'Helvetica'
# Create a figure with GridSpec
fig = plt.figure(figsize=(3.54, 5.51))

outer_gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[3, 1])


ax1 = fig.add_subplot(outer_gs[0], projection = 'polar')
ax3 = fig.add_subplot(outer_gs[1])

ax1.set_yticks([])
ax1.set_xticks([])
ax1.spines["polar"].set_visible(False)
ax1.grid(False)

MPR = MPRester("SZXJWLvi8njBGvA4sT")
def get_matproj_entries(element_list):
    entries = MPR.get_entries_in_chemsys(element_list)
    for entry in entries:
        entry.data["Phase Type"] = "IM"  # We will use this later to color solid solution and intermetallic phases
    # return [e for e in PhaseDiagram(entries).stable_entries]
    return entries

def get_dSmix(comp):
    if isinstance(comp, str):
        comp = Composition(comp)
    if isinstance(comp, Composition):
        dSmix = 0
        # We use the ideal configurational entropy
        for el in comp.elements:
            dSmix += comp.get_atomic_fraction(el)*np.log(comp.get_atomic_fraction(el))
        dSmix = -1*dSmix*8.314*(6.242*10e18)/(6.03*10e23) # convert from J/mol*K to eV/atom*K
        return dSmix


# We need to load the binary regular solution parameters
with open('../data/omegas.json') as f:
        omega_json = json.load(f)

all_els = omega_json["elements"]["FCC"].keys()
stable_pure_dict = {}  # We will use the energies of the pure elements to get formation energies
for el in all_els:
    stable_pure_dict[el] = min(omega_json["elements"]["FCC"][el], omega_json["elements"]["BCC"][el])


def get_dHmix(comp):
    if isinstance(comp, str):
        comp = Composition(comp)
    if isinstance(comp, Composition):
        el_list = [el.symbol for el in comp.elements]
        FCC = 0
        BCC = 0
        for el in el_list:
            if el not in all_els:
                print("Element not in Hautier et al, returning 0")
                return 0
            # Contributions from pure elements being in FCC/BCC
            BCC += (omega_json["elements"]["BCC"][el] - stable_pure_dict[el]) * comp.get_atomic_fraction(el)
            FCC += (omega_json["elements"]["FCC"][el] - stable_pure_dict[el]) * comp.get_atomic_fraction(el)
        for binary in combinations(el_list, 2):
            binary = sorted(binary)
            binary_str = binary[0] + '-' + binary[1]
            # Contributions from mixing
            comp_multiplier = comp.get_atomic_fraction(binary[0]) * comp.get_atomic_fraction(binary[1])
            BCC += omega_json["omegas"]["BCC"][binary_str] * comp_multiplier
            FCC += omega_json["omegas"]["FCC"][binary_str] * comp_multiplier
        return min(FCC, BCC)

def get_dGmix(comp, T):
    return get_dHmix(comp) - T * get_dSmix(comp)

def _make_entry_from_formEperatom(pd, c, formEperatom):

    EntryE = (formEperatom*c.num_atoms
              + sum([c[el]*pd.el_refs[el].energy_per_atom
                                   for el in c.elements]))
    new_entry = ComputedEntry(c, EntryE)
    return new_entry

def get_SS_entries(element_list, pd, T):
    entries = []
    for i in range(1, len(element_list)+1):
        for el_combo in combinations(element_list, i):
            comp_str = ''.join(list(el_combo))
            entry = _make_entry_from_formEperatom(pd, Composition(comp_str),
                                                  get_dGmix(comp_str, T=T))
            entry.data["Phase Type"] = "SS"
            entries.append(entry)
    return entries

def get_IM_and_SS_entries(element_list, T):
    entries = get_matproj_entries(element_list)
    entries.extend(get_SS_entries(element_list, PhaseDiagram(entries), T))
    # entries = get_SS_entries(element_list, PhaseDiagram(entries), T)
    return entries


HEA_entries = get_IM_and_SS_entries(["Hf", "Mo", "Nb", "Ti", "Zr"], 300)


ihw = InverseHullWeb(HEA_entries)
plotter = IHWPlotter(ihw)
plotter.plot(fig = fig, ax = ax3, plot_color_polygon=False, fontsize = 9)


composition = ["Hf", "Mo", "Nb", "Ti", "Zr"]

with open("../data/mol_dict_hfmo.pkl", "rb") as resultFile:
    result = pickle.load(resultFile)
pd = PhaseDiagram(HEA_entries)

count = 0
data = {}
# print(result)
for path, mol_bar in result.items():
    print(mol_bar)
    temp_data = []
    for mol in mol_bar:
        print(mol)
        mol_ratio = dict(zip(composition, mol))
        entry = Composition(mol_ratio)
        comp_entry = _make_entry_from_formEperatom(pd, entry,
                                                   get_dGmix(entry, T=300))
        try:
            e_hull = pd.get_e_above_hull(comp_entry)
            if np.isclose(e_hull, 0, atol=1e-3, rtol=1e-3):
                temp_data.append(0)
            else:
                temp_data.append(e_hull)
        except:
            temp_data.append(0)

        count +=1
    data[path] = temp_data

pos = ax3.get_position()  # Get current position
ax3.set_position([pos.x0 + 0.05, pos.y0 , pos.width - 0.1, pos.height + 0.035])
fontsize = 9
# 'HfMoNbTi', 'MoNbTi', 'HfZr'
ax3.text(s = 'HfMoNbTi', x = 0.9, y = 0.31, transform = ax3.transAxes, bbox = dict(facecolor = 'white', edgecolor = 'white'), c = 'blue', fontsize = fontsize)
ax3.text(s = 'MoNbTi', x = 0.55, y = 0.4, transform = ax3.transAxes, bbox = dict(facecolor = 'white', edgecolor = 'white'), c = 'blue', fontsize = fontsize)
ax3.text(s = 'HfZr', x = 0.7, y = 0.8, transform = ax3.transAxes, bbox = dict(facecolor = 'white', edgecolor = 'white'), c = 'blue', fontsize = fontsize)
ax3.text(s = 'Zr', x = 0.91, y = 0.8, transform = ax3.transAxes, bbox = dict(facecolor = 'white', edgecolor = 'white'), c = 'green', alpha = 0.8, fontsize = fontsize)
ax3.text(s = 'ZrM$o_2$', x = 0.02, y = 0.37, transform = ax3.transAxes, c = 'red', alpha = 0.8, fontsize = fontsize)
ax3.text(s = 'TiM$o_3$', x = 0.3, y = 0.3, transform = ax3.transAxes, bbox = dict(facecolor = 'white', edgecolor = 'white'), c = 'red', alpha = 0.8, fontsize = fontsize)
ax3.text(s = '300 K', x = 0.1, y = 0.8, transform = ax3.transAxes, bbox = dict(facecolor = 'white', edgecolor = 'white'), fontweight = 'bold', fontsize = fontsize)

constraint_element_index = None
property_str = 'e_hull'
subset_idx = None
replacement = None
custom_data = data
is_custom = True
plot_grid = fig, ax1
mol_grid = main(composition, plot_grid, constraint_element_index, subset_idx, replacement, custom_data, is_custom, property_str, cbar_hide=False, cbar_ax=ax1, fontsize = fontsize)

# plt.tight_layout()
plt.savefig("./fig3_revised_small.png", dpi = 200)
# plt.show()


