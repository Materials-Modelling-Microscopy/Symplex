import matplotlib
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from maths import PolarMaths as pm
from utils import get_mol_grid
from plot_utils import draw_circle_in_polar, scatter_center
from property_utils import property_evaluator, cbar_property, property_cbar


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
	total = [np.array(i.split('-')).astype(float) for i in total_str]
	return total_str, total


def set_plot_grid():
	fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(10, 10))
	ax.set_yticks([])
	ax.set_xticks([])
	ax.spines["polar"].set_visible(False)
	ax.grid(False)
	return fig, ax


def color_params(cmap, norm):
	colors = plt.cm.tab20.colors
	line_colors = list(colors)
	return cmap, norm, line_colors


def plot_params(y_bias):
	y_bias = y_bias
	scaling_factor = 0.5
	fontsize = 14
	return y_bias, scaling_factor, fontsize


def plot_bars(n, subset_idx, plot_grid, data, color_params, plot_params):
	fig, ax = plot_grid
	cmap, norm, line_colors = color_params
	y_bias, scaling_factor, fontsize = plot_params
	
	x = np.linspace(0, 1, mol_gradation)
	total_str, total = get_components(data)
	
	if subset_idx is not None:
		total_idx = [idx for idx, i in enumerate(total) if len(i) in subset_idx]
		total = [total[i] for i in total_idx]
		total_str = [total_str[i] for i in total_idx]
	
	angles = pm.divide_circle_degrees(len(total))
	
	for idx, comp in enumerate(total):
		
		bar_length = pm.distance_calculator(n, len(comp)) * scaling_factor
		mol = x * bar_length
		angle_radians = np.radians(angles[idx])
		bar_width = np.radians(360 / (pm.total_num_bars(n) + 5) - 2)
		for i in range(len(mol) - 1):
			# Define the start and end of each radial segment
			if i == 0:
				inner_radius = mol[i] + y_bias
			else:
				inner_radius = mol[i] + y_bias
			
			outer_radius = mol[i + 1] + y_bias
			
			# Create a rectangle to represent each segment in the bar
			rect = Rectangle(
				(angle_radians - bar_width / 2, inner_radius),  # (theta, r) starting point
				width=bar_width,  # Width in radians (angle)
				height=outer_radius - inner_radius,  # Radial height of each segment
				facecolor=cmap(norm(data[total_str[idx]][i + 1])),
				edgecolor="black",
				linewidth=0.5,  # No edge for a smooth gradient effect
				zorder=1
			)
			
			# Add the rectangle to the polar plot
			ax.add_patch(rect)
		
		ax.vlines(
			angle_radians,
			ymin=y_bias,
			ymax=(x * pm.distance_calculator(n, 1) * scaling_factor)[-1]
				 + y_bias,
			linestyles="-",
			color=line_colors[len(comp) - 1],
			zorder=0,
			alpha=0.7,
			linewidth=1.75,
		)
	return total


def plot_text(composition, total, plot_params, ax):
	angles = pm.divide_circle_degrees(len(total))
	y_bias, scaling_factor, fontsize = plot_params
	for idx, comp in enumerate(total):
		rotation = angles[idx] + 180 if 90 < angles[idx] < 270 else angles[idx]
		angle_radians = np.radians(angles[idx])
		composition = np.array(composition)
		comp = comp.astype(int)
		name = '-'.join(composition[comp])
		ax.text(
			angle_radians,
			pm.distance_calculator(n, 1) * (scaling_factor + 0.1)
			+ y_bias,
			name,
			ha="center",
			va="center",
			color="black",
			rotation=rotation,
			fontsize=fontsize
		)


def plot_circles_center(data, n, plot_grid, color_params, plot_params):
	fig, ax = plot_grid
	cmap, norm, line_colors = color_params
	y_bias, scaling_factor, fontsize = plot_params
	
	for N in range(1, n):
		draw_circle_in_polar(radius=pm.distance_calculator(n, N) * 0.5, ax=ax, y_bias=y_bias)
	
	point1 = data[list(data.keys())[0]][0]
	scatter_center(scatter=point1, ax=ax, cmap=cmap, norm=norm, y_bias=y_bias)


if __name__ == "__main__":
	mol_gradation = 15
	composition = ['A', 'B', 'C', 'D']
	n = len(composition)
	constraint_element_index = None
	property_str = 'gibbs'
	subset_idx = None
	custom_data = None
	is_custom = False
	mol_dict = get_mol_grid(n, mol_gradation, constraint_element_index)
	# if adding custom data, use mol_dict for populating the data with same format as mol_dict
	if is_custom and custom_data is not None:
		data = custom_data
	elif is_custom and custom_data is None:
		raise ValueError("Custom data is set to True but no data was provided.")
	else:
		data = get_data(mol_dict, property_str, composition)
	
	plot_grid = set_plot_grid()
	color_param = color_params(*cbar_property(property_str, composition))
	plot_param = plot_params(y_bias=0.4)
	matplotlib.rcParams.update({'font.size': plot_param[2]})
	total = plot_bars(n, subset_idx, plot_grid, data, color_param, plot_param)
	plot_text(composition, total, plot_param, plot_grid[1])
	plot_circles_center(data, n, plot_grid, color_param, plot_param)
	property_cbar(color_param[0], color_param[1], plot_grid[1], plot_param[2], property_str)
	composition_name = '-'.join(composition)
	plt.savefig(f"plots/{composition_name}_{property_str}.png")
