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
	total = [np.array(i.split('-')).astype(float) for i in total_str if '_replace' not in i]
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


def plot_params(y_bias):
	y_bias = y_bias
	scaling_factor = 0.5
	fontsize = 9
	return y_bias, scaling_factor, fontsize


def plot_bars(n, subset_idx, plot_grid, data, color_params, plot_params, replace, mol_gradation):
	fig, ax = plot_grid
	cmap, norm, line_colors = color_params
	y_bias, scaling_factor, fontsize = plot_params
	
	x = np.linspace(0, 1, mol_gradation)
	total_str, total = get_components(data)
	# print(total_str, total)
	if subset_idx is not None:
		total_idx = [idx for idx, i in enumerate(total) if len(i) in subset_idx]
		total = [total[i] for i in total_idx]
		total_str = [total_str[i] for i in total_idx]

	if replace is not None:
		comp1 = np.array([i for i in list(range(n)) if i != replace[1]]).astype(int)
		comp2 = np.array([i for i in list(range(n)) if i != replace[0]]).astype(int)

		# print(comp1, comp2, total)

	angles = pm.divide_circle_degrees(len(total))
	angle_replace = []
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
			linewidth=1.,
		)

		if replace is not None:
			# print(comp1, comp.astype(int))
			if str(comp.astype(int)) == str(comp1):
				angle_replace.append(np.radians(angles[idx]))
			if str(comp.astype(int)) == str(comp2):
				angle_replace.append(np.radians(angles[idx]))

	if replace is not None:
		rect_height = 0.06
		# angle_replace = sorted(angle_replace)
		# print(angle_replace)
		data_replace = data['-'.join(map(str, replace)) + '_replace']
		n_segments = len(data_replace)
		angles = np.linspace(angle_replace[1], angle_replace[0], n_segments + 1)
		# rect_height =
		bar_length = pm.distance_calculator(n, len(comp1)) * scaling_factor
		mol = x * bar_length
		radius = mol[-1] + y_bias + 0.06
		for idx in range(n_segments):
				# Midpoint angle for each rectangle
			angle = (angles[idx] + angles[idx + 1]) / 2
			theta_diff = angles[idx + 1] - angles[idx]

			color = cmap(norm(data_replace[idx]))

			rect = Rectangle(
				xy=(angle - theta_diff, radius - rect_height),  # (angle_start, radius_start)
				width=theta_diff,
				height=rect_height,
				facecolor=color,
				edgecolor='black',
				zorder=0,
				linewidth=1
			)
			ax.add_patch(rect)

		

	return total


def plot_text(composition, total, plot_params, ax, n, mol_gradation):
	angles = np.linspace(0, 2 * np.pi, len(total), endpoint=False)
	y_bias, scaling_factor, fontsize = plot_params
	x = np.linspace(0, 1, mol_gradation)
	for idx, comp in enumerate(total):
		angle_deg = np.degrees(angles[idx])
		rotation = angle_deg + 180 if 90 < angle_deg < 270 else angle_deg
		# angle_radians = np.radians(angles[idx])
		composition = np.array(composition)
		comp = comp.astype(int)
		name = '-'.join(composition[comp])
		radius = (x * pm.distance_calculator(n, 1) * scaling_factor)[-1]+ 1.1*y_bias

		if '-' in name:
			if 90 < angle_deg < 270:
				ax.text(
					angles[idx],
					radius,
					name,
					ha="right",
					va="center",
					color="black",
					rotation=rotation,
					rotation_mode="anchor",
					fontsize=fontsize,
					transform=ax.transData,
					fontname = "Helvetica"
				)
			else:

				ax.text(
					angles[idx],
					radius,
					name,
					ha="left",
					va="center",
					color="black",
					rotation=rotation,
					rotation_mode="anchor",
					fontsize=fontsize,
					transform=ax.transData,
					fontname = "Helvetica"
				)


		else:
			if 90 < angle_deg < 270:
				ax.text(
					angles[idx],
					radius,
					name,
					ha="right",
					va="center",
					color="#BB5566",
					fontweight="bold",
					rotation=rotation,
					rotation_mode="anchor",
					fontsize=fontsize,
					transform=ax.transData,
					fontname = "Helvetica"
				)
			else:

				ax.text(
					angles[idx],
					radius,
					name,
					ha="left",
					va="center",
					color="#BB5566",
					fontweight="bold",
					rotation=rotation,
					rotation_mode="anchor",
					fontsize=fontsize,
					transform=ax.transData,
					fontname = "Helvetica"
				)


def plot_circles_center(data, n, plot_grid, color_params, plot_params):
	fig, ax = plot_grid
	cmap, norm, line_colors = color_params
	y_bias, scaling_factor, fontsize = plot_params
	
	for N in range(1, n):
		draw_circle_in_polar(radius=pm.distance_calculator(n, N) * 0.5, ax=ax, y_bias=y_bias)
	
	point1 = data[list(data.keys())[0]][0]
	scatter_center(scatter=point1, ax=ax, cmap=cmap, norm=norm, y_bias=y_bias)

def main(composition, plot_grid, constraint_element_index, subset_idx, replacement, custom_data, is_custom, property_str, cbar_hide, cbar_ax, fontsize = 10):

	mol_gradation = 15
	n = len(composition)

	mol_dict = get_mol_grid(n, mol_gradation, constraint_element_index, replacement)

	if is_custom and custom_data is not None:
		data = custom_data
	elif is_custom and custom_data is None:
		print("Custom data is set to True but no data was provided. But take the mol_dict")
		return mol_dict
	else:
		data = get_data(mol_dict, property_str, composition)

	# plot_grid = set_plot_grid()
	color_param = color_params(*cbar_property(property_str, composition))
	plot_param = list(plot_params(y_bias=0.4))
	plot_param[2] = fontsize
	matplotlib.rcParams.update({'font.size': plot_param[2]})
	total = plot_bars(n, subset_idx, plot_grid, data, color_param, plot_param, replacement, mol_gradation)
	plot_text(composition, total, plot_param, plot_grid[1], n, mol_gradation)
	plot_circles_center(data, n, plot_grid, color_param, plot_param)
	if not cbar_hide:
		property_cbar(color_param[0], color_param[1], cbar_ax, plot_param[2], property_str)

	# composition_name = '-'.join(composition)
	# plt.savefig(f"plots/{composition_name}_{str(subset_idx)}_{str(replacement)}{property_str}.png", dpi=200,
	# 			bbox_inches='tight')


# if __name__ == "__main__":
#
# 	# mol_gradation = 15
# 	composition = ['A', 'B', 'C', 'D']
# 	constraint_element_index = None
# 	property_str = 'gibbs'
# 	subset_idx = [2]
# 	replacement = None
# 	custom_data = None
# 	is_custom = False
# 	plot_grid = set_plot_grid()
# 	main(composition, plot_grid, constraint_element_index, subset_idx, replacement, custom_data, is_custom, property_str)
	# n = len(composition)
	# constraint_element_index = None
	# property_str = 'gibbs'
	# subset_idx = [2]
	# replacement = None
	# custom_data = None
	# is_custom = False
	# mol_dict = get_mol_grid(n, mol_gradation, constraint_element_index, replacement)
	# # print(mol_dict)
	# # if adding custom data, use mol_dict for populating the data with same format as mol_dict
	# if is_custom and custom_data is not None:
	# 	data = custom_data
	# elif is_custom and custom_data is None:
	# 	raise ValueError("Custom data is set to True but no data was provided.")
	# else:
	# 	data = get_data(mol_dict, property_str, composition)
	# print(data['2-3'], data['0-1'])
	# print(data['0-2-3'], data['1'])
	# # print(data['2-3_replace'])
	# plot_grid = set_plot_grid()
	# color_param = color_params(*cbar_property(property_str, composition))
	# plot_param = plot_params(y_bias=0.4)
	# matplotlib.rcParams.update({'font.size': plot_param[2]})
	# total = plot_bars(n, subset_idx, plot_grid, data, color_param, plot_param, replacement)
	# plot_text(composition, total, plot_param, plot_grid[1])
	# plot_circles_center(data, n, plot_grid, color_param, plot_param)
	# property_cbar(color_param[0], color_param[1], plot_grid[1], plot_param[2], property_str)
	# composition_name = '-'.join(composition)
	# plt.savefig(f"plots/{composition_name}_{str(subset_idx)}_{str(replacement)}{property_str}.png", dpi = 200, bbox_inches='tight')
