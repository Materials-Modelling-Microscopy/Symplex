import pandas as pd

from inverse_hull_web import InverseHullWeb, IHWPlotter
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
import numpy as np
from itertools import combinations
import json

df = pd.read_excel("/Users/mcube/Documents/Symplex/1_Properties_in_Cantor_Space.xlsx")
print(df.head())

mol_bar = df[['Cr', 'Mn', 'Fe', 'Co', 'Ni']].to_numpy()
print(mol_bar)

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

with open('/Users/mcube/Documents/Symplex/data/omegas.json') as f:
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



composition = ["Hf", "Mo", "Nb", "Ti", "Zr"]

pd = PhaseDiagram(HEA_entries)



print(mol_bar)

count = 0
data = {}

temp_data = []
for mol in mol_bar:
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

df['E_hull'] = temp_data

df.to_csv("1_Properties_in_Cantor_Space.csv")




