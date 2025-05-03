import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Element
import matplotlib.colors as mcolors
from utils import load_properties_yaml

property_yaml = load_properties_yaml()[-1]


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
	property_meta_data = property_yaml[property]

	if property == "gibbs":
		base_cmap = plt.get_cmap(property_meta_data['cmap_color'])
		n = property_meta_data['cmap_n_colors']
		colors = [
			*base_cmap(np.linspace(property_meta_data["cmap_min_bound"],property_meta_data["cmap_max_bound"] , n))
		]
		custom_cmap = mcolors.ListedColormap(colors)
		cmap = custom_cmap
		boundaries = list(np.round(np.linspace(property_meta_data["norm_min_bound"], property_meta_data["norm_max_bound"], n), 0))
		norm = mcolors.BoundaryNorm(boundaries, cmap.N)

	elif property == 'e_hull':

		base_cmap = plt.get_cmap(property_meta_data['cmap_color'])
		n = property_meta_data['cmap_n_colors']
		colors = [
			property_meta_data['stable_color'],
			*base_cmap(np.linspace(property_meta_data["cmap_min_bound"],property_meta_data["cmap_max_bound"], n)),
			property_meta_data['unstable_color']
		]

		cmap = mcolors.ListedColormap(colors)
		boundaries = ([0.0,
					  property_meta_data["stable_bound"]] +
					  list(np.round(np.linspace(property_meta_data["stable_bound"]* 1.001,
																					  property_meta_data["metastable_bound"], n), 3))
					  + [property_meta_data["max_bound"]])
		norm = mcolors.BoundaryNorm(boundaries, cmap.N)

	elif property == 'phase_boundary':

		colors = [
			property_meta_data['miscible_color'],
			property_meta_data['immiscible_color'],

		]

		boundaries = [0.0, property_meta_data["stable_bound"], property_meta_data["max_bound"]]
		cmap = mcolors.ListedColormap(colors)
		norm = mcolors.BoundaryNorm(boundaries, len(colors))

	elif property == "entropy":
		n = property_meta_data['cmap_n_colors']
		cmap = plt.get_cmap(property_meta_data['cmap_color'], n)
		cmap = mcolors.ListedColormap(cmap(np.linspace(property_meta_data["cmap_min_bound"],property_meta_data["cmap_max_bound"], n-1))[:-1])
		limit = -(1 / len(composition)) * np.log(1 / len(composition)) * len(composition)
		norm = mcolors.Normalize(vmin=0, vmax=limit)
	
	elif property == "melt":
		t_min = np.min([Element(i).melting_point for i in composition])
		cmap = plt.get_cmap(property_meta_data['cmap_color'])
		norm = mcolors.Normalize(vmin=t_min, vmax=property_meta_data["max_bound"])
	
	elif property == "density":
		n = property_meta_data['cmap_n_colors']
		cmap = plt.get_cmap(property_meta_data['cmap_color'], n)
		cmap = mcolors.ListedColormap(cmap(np.linspace(property_meta_data["cmap_min_bound"],property_meta_data["cmap_max_bound"], n-1))[:-1])
		norm = mcolors.Normalize(vmin=property_meta_data["norm_min_bound"], vmax=property_meta_data["norm_max_bound"])
	elif property == "elastic":
		cmap = plt.get_cmap(property_meta_data['cmap_color'])
		norm = mcolors.Normalize(vmin=property_meta_data["norm_min_bound"], vmax=property_meta_data["norm_max_bound"])
	
	return cmap, norm


def property_cbar(cmap, norm, ax, fontsize, property):
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])
	cbar_properties = property_yaml['cbar']
	cbar = plt.colorbar(
		sm, ax=ax, aspect=cbar_properties['aspect'], orientation=cbar_properties['orientation'], pad=cbar_properties['pad']
	)
	pos = cbar.ax.get_position()

	cbar.ax.set_position([pos.x0 + cbar_properties['x_shift'],
						  pos.y0 + cbar_properties['y_shift'],
						  pos.width + cbar_properties['width_shift'] ,
						  pos.height + cbar_properties['height_shift']])

	if property == 'e_hull':
		cbar.ax.set_xticks([0.0, 0.01, 0.02, 0.03, 0.04, 0.4]) # bit of hardcoding

	cbar.set_label(property_yaml[property]['label'], fontsize = fontsize)