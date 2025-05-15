import re

import matplotlib
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from pymatgen.core import Composition, Element

from maths import PolarMaths as pm
from utils import get_mol_grid, get_mol_grid_central, load_properties_yaml, get_mol_grid_special
from plot_utils import draw_circle_in_polar, scatter_center
from property_utils import property_evaluator, cbar_property, property_cbar

property_yaml = load_properties_yaml()[0]

def get_data(mol_dict, property_str, composition):
	data = {}
	for key, value in mol_dict.items():
		data_bar = []
		for mol in value:
			data_bar.append(property_evaluator(mol, property_str, composition))
		data[key] = data_bar
	return data


def get_components(data):
	total_str = list(data.keys())
	total = [np.array(i.split('-')).astype(float) for i in total_str if '_replace' not in i]
	return total_str, total

def get_components_special(data):
	total_str = list(data.keys())
	total = [[[int(i), float(v)] for i, v in (pair.split('_') for pair in s.split('-'))] for s in total_str]
	return total_str, total


def set_plot_grid():
	fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 12))
	ax.set_yticks([])
	ax.set_xticks([])
	ax.spines["polar"].set_visible(False)
	ax.grid(False)
	return fig, ax


def color_params(cmap, norm):
	colors = plt.cm.Dark2.colors
	line_colors = list(colors)
	return cmap, norm, line_colors


def plot_params():
	y_bias = property_yaml['y_bias']
	scaling_factor = property_yaml['scaling_factor']
	fontsize = property_yaml['fontsize']
	return y_bias, scaling_factor, fontsize


def plot_bars(n, subset_idx, plot_grid, data, color_params, plot_params, replace, central_point=None, ymax=None, special_points = None):
	fig, ax = plot_grid
	cmap, norm, line_colors = color_params
	y_bias, scaling_factor, fontsize = plot_params
	
	x = np.linspace(0, 1, property_yaml['mol_gradation'])
	if special_points is not None:
		total_str, total = get_components_special(data)
	else:
		total_str, total = get_components(data)
	if subset_idx is not None:
		total_idx = [idx for idx, i in enumerate(total) if len(i) in subset_idx]
		total = [total[i] for i in total_idx]
		total_str = [total_str[i] for i in total_idx]

	angles = pm.divide_circle_degrees(len(total))
	angle_replace = []
	for idx, comp in enumerate(total):
		# print(comp)
		if central_point is None:
			bar_length = pm.distance_calculator(n, len(comp)) * scaling_factor
		elif central_point is not None and special_points is  None:
			comp_mol = np.zeros_like(central_point)
			comp_mol[comp.astype(int)] = np.round(1 / len(comp), 3)
			bar_length = pm.distance_calculator_special(central_point, comp_mol) * scaling_factor
		elif central_point is not None and special_points is not None:
			comp_mol = [i[1] for i in comp]
			bar_length = pm.distance_calculator_special(central_point, comp_mol) * scaling_factor
		
		mol = x * bar_length
		angle_radians = np.radians(angles[idx])
		bar_width = np.radians(360 / (pm.total_num_bars(n) + 5) - 2)
		for i in range(len(mol) - 1):
			if i == 0:
				inner_radius = mol[i] + y_bias
				prev_color = cmap(norm(data[total_str[idx]][i]))
			else:
				inner_radius = mol[i] + y_bias
			
			outer_radius = mol[i + 1] + y_bias
			color = cmap(norm(data[total_str[idx]][i + 1]))


			rect = Rectangle(
				(angle_radians - bar_width / 2, inner_radius),  # (theta, r) starting point
				width=bar_width,  # Width in radians (angle)
				height=outer_radius - inner_radius,  # Radial height of each segment
				facecolor=color,
				edgecolor=property_yaml['bar_line_edgecolor'],  # Edge color for the rectangle
				linewidth=property_yaml['bar_linewidth'],  # No edge for a smooth gradient effect
				zorder=1
			)

			ax.add_patch(rect)
			if property_yaml['boundary_flag']:
				if color != prev_color:
					x0, y = rect.get_xy()
					w, h = rect.get_width(), rect.get_height()
					ax.plot([x0, x0+w], [y, y], color='black', lw=1.0)
				if i == len(mol) - 2:
					x0, y = rect.get_xy()
					w, h = rect.get_width(), rect.get_height()
					ax.plot([x0, x0+w], [y+h, y+h], color='black', lw=1.0)

			prev_color = color

		
		if central_point is None:
			ax.vlines(
				angle_radians,
				ymin=y_bias,
				ymax=pm.distance_calculator(n, 1) * scaling_factor
					 + y_bias,
				linestyles="-",
				color=line_colors[len(comp) - 1],
				zorder=0,
				alpha=property_yaml['vline_alpha'],
				linewidth=property_yaml['vline_linewidth'],
			)
		else:
			ax.vlines(
				angle_radians,
				ymin=y_bias,
				ymax= ymax + y_bias,
				linestyles="-",
				color=line_colors[len(comp) - 1],
				zorder=0,
				alpha=property_yaml['vline_alpha'],
				linewidth=property_yaml['vline_linewidth'],
			)

		if replace is not None:
			comp1 = np.array([i for i in list(range(n)) if i != replace[1]]).astype(int)
			comp2 = np.array([i for i in list(range(n)) if i != replace[0]]).astype(int)
			if str(comp.astype(int)) == str(comp1):
				angle_replace.append(np.radians(angles[idx]))
			if str(comp.astype(int)) == str(comp2):
				angle_replace.append(np.radians(angles[idx]))

	if replace is not None:
		rect_height = property_yaml['replace_rect_height']
		data_replace = data['-'.join(map(str, replace)) + '_replace']
		n_segments = len(data_replace)
		angles = np.linspace(angle_replace[1], angle_replace[0], n_segments + 1)
		bar_length = pm.distance_calculator(n, len(comp1)) * scaling_factor
		mol = x * bar_length
		radius = mol[-1] + y_bias + property_yaml['replace_radius_bias']
		for idx in range(n_segments):
			angle = (angles[idx] + angles[idx + 1]) / 2
			theta_diff = angles[idx + 1] - angles[idx]

			color = cmap(norm(data_replace[idx]))

			rect = Rectangle(
				xy=(angle - theta_diff, radius - rect_height),
				width=theta_diff,
				height=rect_height,
				facecolor=color,
				edgecolor=property_yaml['bar_line_edgecolor'],
				zorder=0,
				linewidth=property_yaml['bar_linewidth'],
			)
			ax.add_patch(rect)

	return total


def plot_text(composition, total, plot_params, ax, n, ymax = None):
	angles = np.linspace(0, 2 * np.pi, len(total), endpoint=False)
	y_bias, scaling_factor, fontsize = plot_params
	for idx, comp in enumerate(total):
		angle_deg = np.degrees(angles[idx])
		rotation = angle_deg + 180 if 90 < angle_deg < 270 else angle_deg
		composition = np.array(composition)
		comp = comp.astype(int)
		name = '-'.join(composition[comp])
		
		if ymax is None:
			radius = pm.distance_calculator(n, 1) * scaling_factor+ 1.1*y_bias
		else:
			radius = ymax + 1.1*y_bias
		
		fontweight = 'bold' if '-' in name else 'normal'
		ha = "right" if 90 < angle_deg < 270 else "left"
		color = property_yaml['composition_fontcolor'] if '-' not in name else property_yaml['fontcolor']
		ax.text(
			angles[idx],
			radius,
			name,
			ha=ha,
			va="center",
			color=color,
			rotation=rotation,
			rotation_mode="anchor",
			fontweight=fontweight,
			fontsize=property_yaml['fontsize'],
			transform=ax.transData,
			fontname = property_yaml['fontname']
		)


def plot_text_special(composition, total, plot_params, ax, n, ymax=None):
	angles = np.linspace(0, 2 * np.pi, len(total), endpoint=False)
	y_bias, scaling_factor, fontsize = plot_params
	for idx, comp in enumerate(total):
		angle_deg = np.degrees(angles[idx])
		rotation = angle_deg + 180 if 90 < angle_deg < 270 else angle_deg
		composition = np.array(composition)
		# comp = comp.astype(int)
		comp_dict = dict(zip([Element(composition[i[0]]) for i in comp], [i[1] for i in comp]))
		composition_formula = Composition(comp_dict)
		# name = ''.join([f"${composition[i[0]]}_" +"{" + f"{i[1]}" +"}$" for i in comp if i[1] != 0])
		name = str(composition_formula.get_integer_formula_and_factor()[0])

		name_rep = re.sub(r'(?<=[A-Za-z])(\d+(?:\.\d+)?)(?!\d)', r'$_{\1}$',
						  re.sub(r'(?<!\d)(\d+(?:\.\d+)?)(?=[A-Za-z])', r'$_{\1}$', name))
		if ymax is None:
			radius = pm.distance_calculator(n, 1) * scaling_factor + 1.1 * y_bias
		else:
			radius = ymax + 1.1 * y_bias

		fontweight = 'bold' if len(comp) == 1 else 'normal'
		ha = "right" if 90 < angle_deg < 270 else "left"
		color = property_yaml['composition_fontcolor'] if len(comp) == 1 else property_yaml['fontcolor']
		ax.text(
			angles[idx],
			radius,
			name_rep,
			ha=ha,
			va="center",
			color=color,
			rotation=rotation,
			rotation_mode="anchor",
			fontweight=fontweight,
			fontsize=property_yaml['fontsize'],
			transform=ax.transData,
			fontname=property_yaml['fontname']
		)


def plot_circles_center(data, n, plot_grid, color_params, plot_params, central_point=None):
	fig, ax = plot_grid
	cmap, norm, line_colors = color_params
	y_bias, scaling_factor, fontsize = plot_params
	
	for N in range(1, n):
		if central_point is None:
			draw_circle_in_polar(radius=pm.distance_calculator(n, N) * scaling_factor, ax=ax, y_bias=y_bias)

	point1 = data[list(data.keys())[0]][0]
	scatter_center(scatter=point1, ax=ax, cmap=cmap, norm=norm, y_bias=y_bias)

def main(composition,
		 plot_grid,
		 constraint_element_index,
		 property_str,
		 cbar_hide,
		 cbar_ax,
		 custom_data = None,
		 is_custom = None,
		 central_point = None,
		 special_points = None,
		 subset_idx = None,
		 replacement = None):

	n = len(composition)

	if central_point is None and special_points is None:
		mol_dict = get_mol_grid(n, property_yaml['mol_gradation'], constraint_element_index, replacement)
	if central_point is None and special_points is not None:
		temp_point = np.array([1/n]*n)
		mol_dict = get_mol_grid_special(n, temp_point, special_points, property_yaml['mol_gradation'])
	elif central_point is not None and special_points is None:
		mol_dict = get_mol_grid_central(n, central_point, property_yaml['mol_gradation'], constraint_element_index)
	elif central_point is not None and special_points is not None:
		mol_dict = get_mol_grid_special(n, central_point, special_points, property_yaml['mol_gradation'])

	if is_custom and custom_data is not None:
		data = custom_data
	elif is_custom and custom_data is None:
		print("Custom data is set to True but no data was provided. But take the mol_dict")
		return mol_dict
	else:
		data = get_data(mol_dict, property_str, composition)

	color_param = color_params(*cbar_property(property_str, composition))
	plot_param = list(plot_params())
	matplotlib.rcParams.update({'font.size': plot_param[2]})
	
	if central_point is None:
		total = plot_bars(n, subset_idx, plot_grid, data, color_param, plot_param, replacement)
		plot_text(composition, total, plot_param, plot_grid[1], n)
		plot_circles_center(data, n, plot_grid, color_param, plot_param)

	if central_point is not None:
		y_max = 0
		for keys in mol_dict.values():
			radii = pm.distance_calculator_special(central_point, keys[-1])
			y_max = max(y_max, radii)
		
		total = plot_bars(n, subset_idx, plot_grid, data, color_param, plot_param, replacement,
						  central_point=central_point, ymax=y_max, special_points=special_points)

		if special_points is not None:
			plot_text_special(composition, total, plot_param, plot_grid[1], n, ymax=y_max)
		else:
			plot_text(composition, total, plot_param, plot_grid[1], n, ymax=y_max)
		plot_circles_center(data, n, plot_grid, color_param, plot_param, central_point=central_point)
		# title = list(zip(composition, central_point))
		# plot_grid[1].set_title(''.join([f'${key}' +'_{' +  str(value) + "}$" for key, value in title]), x = -0.25, y = 1.05, fontsize = fontsize)


	if not cbar_hide:
		property_cbar(color_param[0], color_param[1], cbar_ax, plot_param[2], property_str)

