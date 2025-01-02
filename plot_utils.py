import matplotlib.pyplot as plt
import numpy as np


def draw_circle_in_polar(radius: float, ax: plt.Axes, y_bias) -> None:
	"""
	Draws a circle with the given radius on a polar plot.

	Args:
		radius (float): The radius of the circle.
		ax (plt.Axes): The polar plot axes on which the circle is drawn.
	"""
	theta = np.linspace(0, 2 * np.pi, 100)
	ax.plot(
		theta,
		[radius + y_bias] * len(theta),
		linewidth=0.8,
		zorder=0,
		color="black",
		linestyle="--",
		alpha=0.1,
	)
	
def scatter_center(scatter: float, ax: plt.Axes, cmap, norm, y_bias) -> None:
	"""
	Plots the central scatter point in a polar plot.

	Args:
		scatter (float): The value to normalize the color of the scatter point.
		ax (plt.Axes): The polar plot axes.
	"""
	
	theta = np.linspace(0, 2 * np.pi, 100)
	# self.draw_circle_in_polar(radius=self.y_bias, ax=ax)
	ax.plot(
		theta,
		[y_bias - 0.05] * len(theta),
		linewidth=1.5,
		zorder=100,
		color="black",
		linestyle="-",
		alpha=1,
	)
	ax.fill(
		theta,
		[y_bias - 0.05] * len(theta),
		zorder=100,
		color=cmap(norm(scatter)),
		alpha=1,
	)