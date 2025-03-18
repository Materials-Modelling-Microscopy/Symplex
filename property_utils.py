import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Element
import matplotlib.colors as mcolors


def property_evaluator(mol, property, composition):
	if property == 'density':
		modulus = np.array([Element(i).atomic_mass / Element(i).molar_volume for i in composition])
		return np.dot(mol, modulus)
	
	if property == 'elastic':
		modulus = np.array([Element(i).youngs_modulus for i in composition])
		return np.dot(mol, modulus)
	
	if property == 'thermal':
		modulus = np.array([Element(i).thermal_conductivity for i in composition])
		return np.dot(mol, modulus)
	
	if property == 'entropy':
		return -np.dot(mol, np.log(mol + 1e-10))
	
	if property == 'gibbs':
		omega1, omega2, omega3, omega4 = 0.420, 0.420, 0.420, 0.420
		entropy = -np.dot(mol, np.log(mol + 1e-10))
		i = mol
		hmix = sum([omega1 * i[0] * i[1], omega2 / 2 * i[0] * i[2], omega3 / 4 * i[0] * i[3], omega4 * i[1] * i[2],
					omega1 * i[1] * i[3], omega3 * i[2] * i[3]])
		return (-1000 * 8.315e-5*entropy + hmix)*1000
	
	if property == 'melt':
		modulus = np.array([Element(i).melting_point for i in composition])
		return np.dot(mol, modulus)


def cbar_property(property, composition):
	if property == "gibbs":
		base_cmap = plt.get_cmap('coolwarm')
		n = 10
		colors = [
			*base_cmap(np.linspace(0, 1, n))
		]
		custom_cmap = mcolors.ListedColormap(colors)
		cmap = custom_cmap
		boundaries = list(np.round(np.linspace(-50, 50, n), 0))
		norm = mcolors.BoundaryNorm(boundaries, cmap.N)

	elif property == 'e_hull':
		base_cmap = plt.get_cmap('OrRd')
		n = 10
		colors = [
			'#009988',  # Dark blue for 0
			*base_cmap(np.linspace(0.1, 0.75, n)),  # OrRd for values between 0.0 (exclusive) and 0.05
			'darkred'  # Darkest red for values above 0.05
		]
		custom_cmap = mcolors.ListedColormap(colors)
		boundaries = [0.0, 0.002] + list(np.round(np.linspace(0.0051, 0.05, n), 3)) + [0.4]
		norm = mcolors.BoundaryNorm(boundaries, custom_cmap.N)
		cmap = custom_cmap

	elif property == "entropy":
		cmap = plt.get_cmap('tab20c', 8)
		cmap = mcolors.ListedColormap(cmap(np.linspace(0.0, 0.8, 7))[:-1])
		limit = -(1 / len(composition)) * np.log(1 / len(composition)) * len(composition)
		norm = mcolors.Normalize(vmin=0, vmax=limit)
	
	elif property == "melt":
		t_min = np.min([Element(i).melting_point for i in composition])
		cmap = plt.get_cmap('plasma')
		norm = mcolors.Normalize(vmin=t_min, vmax=3700)
	
	elif property == "density":
		cmap = plt.get_cmap('jet', 8)
		cmap = mcolors.ListedColormap(cmap(np.linspace(0.2, 1, 7))[:-1])
		norm = mcolors.Normalize(vmin=7, vmax=9)
	elif property == "elastic":
		cmap = plt.get_cmap('jet')
		norm = mcolors.Normalize(vmin=100, vmax=300)
	
	return cmap, norm


def property_cbar(cmap, norm, ax, fontsize, property):
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])  # We need this for colorbar to work
	cbar = plt.colorbar(
		sm, ax=ax, aspect=50, orientation="horizontal", pad=0.1
	)
	# cbar = plt.colorbar(
	# 	sm, ax=ax, aspect = 100, orientation="horizontal"
	# )
	pos = cbar.ax.get_position()

	# Modify width (increase `pos.width`) to make it **longer**
	cbar.ax.set_position([pos.x0, pos.y0-0.01, pos.width, pos.height])

	if property == "misc_T":
		cbar.set_label("$T_{melt}$ - $T_{misc}$ (K)", fontsize=fontsize)
	if property == 'e_hull':
		cbar.ax.set_xticks([0.0, 0.01, 0.02, 0.03, 0.04, 0.4])
		cbar.set_label("$E_{hull}$ (eV/atom)", fontsize=fontsize)
	if property == 'melt':
		cbar.set_label("$T_{melt}$ (K)", fontsize=fontsize)
	if property == 'entropy':
		print("Entropy")
		cbar.set_label("$-Entropy/k_b$", fontsize=fontsize)
	if property == 'elastic':
		cbar.set_label("$E (GPa)$", fontsize=fontsize)
	if property == 'density':
		cbar.set_label("$\\rho (g/cc)$", fontsize=fontsize)
	if property == 'gibbs':
		cbar.set_label("$G_{mix} (meV/atom)$", fontsize=fontsize)


	# plt.subplots_adjust(bottom=0.15)
