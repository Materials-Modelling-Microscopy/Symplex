from itertools import combinations
import numpy as np
from typing import List, Iterable
import yaml


def load_properties_yaml():
	file_path = "/Users/mcube/Documents/Symplex/plotting_parameters.yaml"

	with open(file_path,'r') as input_file:
		return list(yaml.safe_load_all(input_file))


def create_high_sym_mol_grid(
		change_idx: List[int], x: Iterable, n: int, N: int
) -> np.ndarray:
	
	mol_list = []
	for mol in x:
		addition = np.zeros(n)
		addition += [(1 - mol) / n] * n
		print(mol, addition)
		for i in change_idx:
			addition[i] = (1 / N - 1 / n) * mol + (1 / n)
		mol_list.append(addition)
	
	return np.array(mol_list)

def create_central_point_sym_mol_grid(
		change_idx: List[int], x, n: int, N: int, central_point : Iterable
) -> np.ndarray:

	central_point = np.array(central_point)
	assert len(central_point) == n
	assert np.round(np.sum(central_point), 3) == 1.0

	end_point = np.zeros_like(central_point)
	end_point[change_idx] = np.round(1/N, 3)
	return np.round(np.array([(1 - t) * central_point + t * end_point for t in x]), 3)

def create_mol_grid_transmutation(
		transmutation_indice: List[int], n: int, x: Iterable
) -> np.ndarray:
	
	mols = []
	for i in x:
		subtract = np.zeros(n)
		subtract += 1 / (n - 1)
		subtract[transmutation_indice[0]] += -1 / (n - 1) + i / (n - 1)
		subtract[transmutation_indice[1]] -= i / (n - 1)
		mols.append(subtract)
	
	return np.array(mols)


def create_two_special_points(x, n: int, central_point: Iterable, special_point: Iterable
							  ) -> np.ndarray:

	central_point = np.array(central_point)
	special_point = np.array(special_point)
	assert len(central_point) == n
	assert np.round(np.sum(central_point), 3) == 1.0

	assert len(special_point) == n
	assert np.round(np.sum(special_point), 3) == 1.0

	return np.round(np.array([(1 - t) * central_point + t * special_point for t in x]), 3)


def find_indices(main_list: list, subset: list) -> list:
	
	indices = []
	for value in subset:
		try:
			index = main_list.index(
				value
			)  # Find the index of the value in the main list
			indices.append(index)
		except ValueError:
			indices.append(
				None
			)  # If the value is not found, append None or handle as needed
	return indices


def get_combinations(components, constraint_element_index):
	combinations_components = [list(combinations(components, i)) for i in range(1, len(components))]
	combinations_components = [item for sublist in combinations_components for item in sublist]
	constraint_element_index = 0 if constraint_element_index is None else constraint_element_index
	total = [list(i) for i in combinations_components if constraint_element_index in i]
	
	total_mirror = [sorted(tuple(set(components).difference(i))) for i in total]
	
	return total + total_mirror


def get_mol_grid(n, mol_gradation, constraint_element_index, replace):
	components = [i for i in range(n)]
	x = np.linspace(0, 1, mol_gradation)
	
	total = get_combinations(components, constraint_element_index)
	
	total_check = [''.join(map(str, i)) for i in total]
	assert len(total) == len(list(set(total_check)))
	
	mol_dict = {}
	for idx2, comp in enumerate(total):
		temp_i = comp
		comp_str = '-'.join(map(str, comp))
		mol_grid = create_high_sym_mol_grid(
			x=x, n=n, N=len(temp_i), change_idx=temp_i
		)
		mol_dict[comp_str] = mol_grid

	if replace:
		indices = replace
		comp_str = '-'.join(map(str, indices))
		mol_grid = create_mol_grid_transmutation(transmutation_indice=indices, n=n, x=x)
		mol_dict[comp_str + '_replace'] = mol_grid

	return mol_dict


def get_mol_grid_central(n, central_point, mol_gradation, constraint_element_index):
	components = [i for i in range(n)]
	x = np.linspace(0, 1, mol_gradation)

	total = get_combinations(components, constraint_element_index)

	total_check = [''.join(map(str, i)) for i in total]
	assert len(total) == len(list(set(total_check)))

	mol_dict = {}
	for idx2, comp in enumerate(total):
		temp_i = comp
		comp_str = '-'.join(map(str, comp))
		mol_grid = create_central_point_sym_mol_grid(
			central_point=central_point, x=x, n=n, N=len(temp_i), change_idx=temp_i
		)
		mol_dict[comp_str] = mol_grid


	return mol_dict

def get_mol_grid_special(n, central_point, special_points, mol_gradation):

	x = np.linspace(0, 1, mol_gradation)
	components = [i for i in range(n)]
	mol_dict = {}
	for idx2, comp in enumerate(special_points):
		temp_i = comp
		comp_str = '-'.join(f"{i}_{v}" for i, v in zip(components, temp_i))
		mol_grid = create_two_special_points(
			central_point=central_point, x=x, n=n, special_point=comp
		)
		mol_dict[comp_str] = mol_grid


	return mol_dict

print(get_mol_grid_special(n = 4, central_point = [0.25, 0.25, 0.25 ,0.25], special_points = [[0, 0.2, 0.8, 0], [0.2, 0.4, 0.1, 0.3]], mol_gradation=10))