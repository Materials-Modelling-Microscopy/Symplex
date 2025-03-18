import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import matplotlib
from main import main


matplotlib.rcParams['font.family'] = 'Helvetica'

fig = plt.figure(figsize=(7.48, 5))

gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 1.1])

ax1 = fig.add_subplot(gs[1], projection = 'polar')
img_path = "/Users/mcube/Library/CloudStorage/Box-Box/Pravan_John_shared/Projects/Vizualization/visualization_figures/Fig1a.jpg"  # Replace with your actual image file
img = mpimg.imread(img_path)

ax1.set_yticks([])
ax1.set_xticks([])
ax1.spines["polar"].set_visible(False)
ax1.grid(False)

composition = ['A', 'B', 'C', 'D']
constraint_element_index = None
property_str = 'gibbs'
subset_idx = None
replacement = None
custom_data = None
is_custom = False
plot_grid = fig, ax1
main(composition, plot_grid, constraint_element_index, subset_idx, replacement, custom_data, is_custom, property_str, cbar_hide=False, cbar_ax=ax1, fontsize = 12)

ax2 = fig.add_subplot(gs[0])
ax2.set_yticks([])
ax2.set_xticks([])
ax2.grid(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.imshow(img)

ax1.text(s = "(a)", fontsize = 12, y = 1.15, x = -0.05, transform=ax2.transAxes)
ax2.text(s = "(b)", fontsize = 12, y = 1.15, x = 1.08, transform=ax2.transAxes)

# plt.tight_layout()
plt.savefig("./fig1_revised.png", dpi = 200)
# plt.show()