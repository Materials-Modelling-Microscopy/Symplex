import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib
from main import main

matplotlib.rcParams['font.size'] = 11
# Create a figure with GridSpec
fig = plt.figure(figsize=(3.54, 9))

outer_gs = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[2, 2, 1.5])
inner_gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec = outer_gs[2], hspace = 0)


ax2 = fig.add_subplot(outer_gs[0], projection = 'polar')
ax1 = fig.add_subplot(outer_gs[1], projection = 'polar')
# ax2_0 = fig.add_subplot(gs[2])
polar_plots = [ax1, ax2]
# ax2_0.set_yticks([])
# ax2_0.set_xticks([])
# ax2_0.grid(False)
# ax2_0.spines['top'].set_visible(False)
# ax2_0.spines['bottom'].set_visible(False)
# ax2_0.spines['left'].set_visible(False)
# ax2_0.spines['right'].set_visible(False)

for polar_plot in polar_plots:
	polar_plot.set_yticks([])
	polar_plot.set_xticks([])
	polar_plot.spines["polar"].set_visible(False)
	polar_plot.grid(False)

plt.subplots_adjust(right = 0.99)
composition = ['A', 'B', 'C', 'D']
constraint_element_index = None
property_str = 'gibbs'

constraint_element_index = None
property_str = 'gibbs'
subset_idx = [2]
replacement = None
custom_data = None
is_custom = False
plot_grid = fig, ax1
main(composition, plot_grid, constraint_element_index, subset_idx, replacement, custom_data, is_custom, property_str, cbar_hide=False, cbar_ax = ax1, fontsize = 11)

constraint_element_index = None
property_str = 'gibbs'
subset_idx = [1,3]
replacement = [2,3]
custom_data = None
is_custom = False
plot_grid = fig, ax2
main(composition, plot_grid, constraint_element_index, subset_idx, replacement, custom_data, is_custom, property_str, cbar_hide=True, cbar_ax = ax2, fontsize = 11)



ax4 = fig.add_subplot(inner_gs[0])
ax5 = fig.add_subplot(inner_gs[1])
ax3 = fig.add_subplot(inner_gs[2])

# plt.subplots_adjust(hspace = 0.4)

# pos = ax2.get_position()  # Get current position
# ax2.set_position([pos.x0, pos.y0 - 0.05, pos.width, pos.height])
#
#

# ax5.set_xticks([])

AB = [9.417123906141093, 9.528975565517337, 9.866717349570619, 10.437047439907088, 11.251611051051574, 12.327808151691682, 13.690124250264986, 15.372269694646654, 17.420676917426757, 19.9004792824959, 22.90649724687649, 26.585737532669185, 31.19198129419523, 37.26543752088299, 47.36481195307054]
CD = [9.417123906141093, 9.528975565517309, 9.866717349570592, 10.437047439907102, 11.251611051051574, 12.327808151691656, 13.690124250264974, 15.372269694646626, 17.420676917426743, 19.900479282495887, 22.90649724687649, 26.585737532669185, 31.191981294195244, 37.265437520883005, 47.36481195307054]
ACD = [9.417123906141093, 7.818613866838311, 6.1531938026324475, 4.429627506475631, 2.658456460460326, 0.8525234828552358, -0.9722676633526495, -2.795807535955208, -4.591971846705858, -6.325582701966454, -7.946874935278195, -9.38021284379284, -10.49676793964109, -11.024668523755206, -9.682945111141649]
B = [9.417123906141093, 13.743073847468201, 17.260262955698185, 19.877822074503616, 21.541603734744633, 22.22473618053525, 21.923500354807196, 20.65678735439129, 18.468924564366105, 15.437018510925604, 11.68630530896765, 7.423093396628312, 3.016056474505084, -0.7347261151440243, 8.315000687570434e-09]
ABC_ABD = [13.650388222191667, 10.446944672004896, 9.664283322898648, 9.606509879748746, 9.925492324028145, 10.466821257623687, 11.151023357130002, 11.93865887598687, 12.817690023796679, 13.800154590957028, 14.925492324028136, 16.2731765464154, 17.997616656231965, 20.446944672004904, 25.31705488885834]

linewidth = 2
ms = 4
x0 = np.linspace(0, 1, 15)*0.5
x1 = x0 + 0.5
ax3.plot(x0, AB[::-1], marker = 'o', c = '#4477AA', linewidth = linewidth, ms = ms)
ax3.plot(x1, CD, marker = 'o', c = '#4477AA', linewidth = linewidth, ms = ms)
ax3.axvline(0.5, c = 'black', linestyle = '--', zorder = 0)

x4 = np.linspace(0, 1, 15)
ax4.plot(x4, ABC_ABD, marker = 'o', c = '#4477AA', linewidth = linewidth, ms = ms)

x2 = np.linspace(0, 1, 15)
x2 = x2/3
x3 = np.linspace(0.33, 1, 15)
ax5.plot(x2, ACD[::-1], marker = 'o', c = '#4477AA', linewidth = linewidth, ms = ms)
ax5.plot(x3, B, marker = 'o', c = '#4477AA', linewidth = linewidth, ms = ms)
ax5.axvline(0.33, c = 'black', linestyle = '--', zorder = 0)
ax5.set_ylabel('$G_{mix} (meV/atom)$')
ax4.set_yticks([15, 25])
ax3.set_yticks([25, 45])
ax5.set_xticks([])
ax3.set_xticks([])
ax4.set_xticks([])
ax3.text(s = 'AB', x = -0.04, y = 19)
ax3.text(s = 'CD', x = 0.96, y = 19)
ax3.text(s = 'ABCD', x = 0.425, y = 25, bbox=dict(facecolor='white', alpha=1, edgecolor= 'white'), zorder = 0)
ax4.text(s = 'ABD', x = -0.01, y = 18)
ax4.text(s = 'ABC', x = 0.9, y = 32)
ax5.text(s = 'ABCD', x = 0.34, y = -15, bbox=dict(facecolor='white', alpha=1, edgecolor= 'white'), zorder = 0)
ax5.text(s = 'B', x = 0.98, y = 5)
ax5.text(s = 'ACD', x = 0.0, y = 0)
ax3.set_ylim([-25, 60])
ax5.set_ylim([-25, 60])
ax4.set_ylim([-25, 60])
ax3.set_yticks([0, 40])
ax4.set_yticks([0, 40])
ax5.set_yticks([0, 40])

pos = ax3.get_position()  # Get current position
ax3.set_position([pos.x0 + 0.05, pos.y0, pos.width - 0.1, pos.height])
pos = ax5.get_position()  # Get current position
ax5.set_position([pos.x0 + 0.05, pos.y0, pos.width - 0.1, pos.height])
pos = ax4.get_position()  # Get current position
ax4.set_position([pos.x0 + 0.05, pos.y0, pos.width - 0.1, pos.height])

ax1.text(s = "(b)", fontsize = 11, y = 1.1, x = -0.58, transform=ax1.transAxes)
ax2.text(s = "(a)", fontsize = 11, y = 1.1, x = -0.35, transform=ax2.transAxes)
ax4.text(s = "(c)", fontsize = 11, y = 1.2, x = -0.2, transform=ax4.transAxes)

# plt.tight_layout()
plt.savefig("./fig2_revised.png", dpi = 200)
# plt.show()