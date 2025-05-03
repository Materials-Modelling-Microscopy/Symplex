import matplotlib.pyplot as plt
import numpy as np
from utils import load_properties_yaml

property_yaml = load_properties_yaml()[1]

def draw_circle_in_polar(radius: float, ax: plt.Axes, y_bias) -> None:

	theta = np.linspace(0, 2 * np.pi, 100)
	ax.plot(
		theta,
		[radius + y_bias] * len(theta),
		linewidth= property_yaml["circle_linewidth"],
		zorder=0,
		color=property_yaml["circle_color"],
		linestyle="--",
		alpha=property_yaml["circle_alpha"],
	)
	
def scatter_center(scatter: float, ax: plt.Axes, cmap, norm, y_bias) -> None:
	
	theta = np.linspace(0, 2 * np.pi, 100)
	ax.plot(
		theta,
		[y_bias - property_yaml["center_bias"]] * len(theta),
		linewidth=property_yaml["center_linewidth"],
		zorder=100,
		color=property_yaml["center_edgecolor"],
		linestyle="-",
		alpha=property_yaml["center_alpha"],
	)
	ax.fill(
		theta,
		[y_bias - property_yaml["center_bias"]] * len(theta),
		zorder=100,
		color=cmap(norm(scatter)),
		alpha=1,
	)