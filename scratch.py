from typing import List, Iterable

from utils import get_mol_grid, create_high_sym_mol_grid, get_mol_grid_central
import numpy as np


mol_gradation = 15
n = 5
constraint_element_index = 0
replacement = False

# mol_dict = get_mol_grid(n, mol_gradation, constraint_element_index, replacement)
# print(mol_dict)
change_idx = [0, 1]
N =len(change_idx)

x = np.linspace(0, 1, mol_gradation)

# high_sym = create_high_sym_mol_grid(x = x, change_idx = change_idx, n = n, N= N)
# print(high_sym)


def create_central_point_sym_mol_grid(
		change_idx: List[int], x, n: int, N: int, central_point : Iterable
) -> np.ndarray:

	central_point = np.array(central_point)
	assert len(central_point) == n
	assert np.round(np.sum(central_point), 3) == 1.0

	end_point = np.zeros_like(central_point)
	end_point[change_idx] = np.round(1/N, 3)
	return np.round(np.array([(1 - t) * central_point + t * end_point for t in x]), 3)

# create_central_point_sym_mol_grid(change_idx, x, n, N, central_point=[0.2, 0.2, 0.2, 0.3, 0.1])
print(get_mol_grid_central(n = n, central_point=[0.2, 0.2, 0.2, 0.3, 0.1], mol_gradation=mol_gradation, constraint_element_index=constraint_element_index))
