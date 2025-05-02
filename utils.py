from itertools import combinations
import numpy as np
from typing import List, Iterable


#creation of various mol_grids

def create_high_sym_mol_grid(
		change_idx: List[int], x: Iterable, n: int, N: int
) -> np.ndarray:
	"""
	Create a high-symmetry mole fraction grid based on the provided indices.

	Args:
			change_idx (List[int]): Indices of the components to modify in the grid.
			x (Iterable): List of mole fractions to use for the grid.
			n (int): The total number of components in the system.
			N (int): The number of components to symmetrically modify.

	Returns:
			np.ndarray: A 2D array where each row represents a high-symmetry configuration
									of the components.
	"""
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
	"""
	Create a mole fraction grid with transmutation between two components.

	Args:
			transmutation_indice (List[int]): Indices of the components undergoing transmutation.
			n (int): The total number of components in the system.
			x ( Iterable): List of mole fractions for transmutation.

	Returns:
			np.ndarray: A 2D array where each row represents a configuration with
									transmuted mole fractions between two components.
	"""
	mols = []
	for i in x:
		subtract = np.zeros(n)
		subtract += 1 / (n - 1)
		subtract[transmutation_indice[0]] += -1 / (n - 1) + i / (n - 1)
		subtract[transmutation_indice[1]] -= i / (n - 1)
		mols.append(subtract)
	
	return np.array(mols)


def find_indices(main_list: list, subset: list) -> list:
	"""
	Finds the indices of the `subset` elements in the `main_list`. If a value in the subset is not found,
	returns `None` for that value's index.

	Args:
			main_list (list): The list to search within.
			subset (list): The list of values whose indices are to be found in the main list.

	Returns:
			list: A list of indices corresponding to the subset elements in the main list. If an element is not found, `None` is returned in its place.

	Example::

			main_list = ['a', 'b', 'c', 'd']
			subset = ['b', 'd', 'x']
			result = FancyListExtractions.find_indices(main_list, subset)
			# Output: [1, 3, None]
	"""
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
