import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.image as mpimg
from main import main
import pandas as pd
from matplotlib.colors import ListedColormap
import seaborn as sns

matplotlib.rcParams['font.size'] = 10
# Create a figure with GridSpec
fig = plt.figure(figsize=(7.48, 10))

outer_gs = gridspec.GridSpec(4, 2, figure=fig, width_ratios=[2, 1])

ax1 = fig.add_subplot(outer_gs[0:2, 0], projection = 'polar')
ax2 = fig.add_subplot(outer_gs[2:, 0], projection = 'polar')
ax3 = fig.add_subplot(outer_gs[0, 1])
ax4 = fig.add_subplot(outer_gs[1, 1])
ax5 = fig.add_subplot(outer_gs[2, 1])
ax6 = fig.add_subplot(outer_gs[3, 1])


ax1.set_yticks([])
ax1.set_xticks([])
ax1.spines["polar"].set_visible(False)
ax1.grid(False)

ax2.set_yticks([])
ax2.set_xticks([])
ax2.spines["polar"].set_visible(False)
ax2.grid(False)

composition = ['Cr',	'Mn',	'Fe',	'Co',	'Ni']
constraint_element_index = None
property_str = 'entropy'
subset_idx = None
replacement = None
custom_data = None
is_custom = False
plot_grid = fig, ax1
main(composition, plot_grid, constraint_element_index, subset_idx, replacement, custom_data, is_custom, property_str, cbar_hide=False, cbar_ax=ax1)

composition = ['Cr',	'Mn',	'Fe',	'Co',	'Ni']
constraint_element_index = None
property_str = 'density'
subset_idx = None
replacement = None
custom_data = None
is_custom = False
plot_grid = fig, ax2
main(composition, plot_grid, constraint_element_index, subset_idx, replacement, custom_data, is_custom, property_str, cbar_hide=False, cbar_ax=ax2)


df = pd.read_excel('/Users/mcube/Desktop/Projects/far_heaa/src/far_heaa/raymundo/data/1_Properties_in_Cantor_Space.xlsx')
# dfc = pd.read_excel('data/1_Cantor_benchmark.xlsx')
elements = ['Cr',	'Mn',	'Fe',	'Co',	'Ni']
props = ['Density Avg','FCC YS T C PRIOR','umap0','umap1']
prop = 'entropy'
composition = df[['Cr', 'Mn', 'Fe', 'Co', 'Ni']].to_numpy()
entropy = -np.sum(composition*np.log(composition + 1e-9), axis = 1)
df['entropy'] = entropy
df = df.sort_values(by=[prop], ascending=True)
dfp = df
cmap = plt.get_cmap('tab20c', 8)
cmap = ListedColormap(cmap(np.linspace(0.0, 0.8, 7))[:-1])
# ax3.scatter(df['umap0'],df['umap1'],c='grey')
ax3.scatter(dfp['umap0'],dfp['umap1'],c=dfp[prop],cmap=cmap, edgecolor = 'black', linewidth=0.5, s = 20)
ax3.grid(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.set_xticks([])
ax3.set_yticks([])

df = df.sort_values(by=[prop], ascending=False)
dfp = df
# ax4.scatter(df['umap0'],df['umap1'],c='grey')
ax4.scatter(dfp['umap0'],dfp['umap1'],c=dfp[prop],cmap=cmap, edgecolor='black', linewidth=0.5)
ax4.grid(False)
ax4.spines['top'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.set_xticks([])
ax4.set_yticks([])

prop = 'Density Avg'
df = df.sort_values(by=[prop],ascending=False)
ax5.scatter(df['umap0'],df['umap1'],c='grey')
dfp = df
ax5.scatter(dfp['umap0'],dfp['umap1'],c=dfp[prop],cmap='jet',alpha=.7)
ax5.grid(False)
ax5.spines['top'].set_visible(False)
ax5.spines['bottom'].set_visible(False)
ax5.spines['left'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.set_xticks([])
ax5.set_yticks([])



df = pd.read_excel('/Users/mcube/Desktop/Projects/far_heaa/src/far_heaa/raymundo/data/1_Properties_in_Cantor_Space.xlsx')
elements = ['Co','Cr','Fe','Mn','Ni']
df[elements] = (df[elements]*100).astype(int)
prop = 'Density Avg'
order=  list(map(int, range(0, 100, 20)))
# plt.figure(figsize=(2.86, 2.86))
boxprops = dict(facecolor='#ff7f0e', edgecolor='black', )
whiskerprops = dict(color='black')
medianprops = dict(color='black' )
capprops = dict(color='black')
sns.boxplot(x=df['Cr'], y=df[prop], showfliers=False, boxprops=boxprops, whiskerprops=whiskerprops,
				 medianprops=medianprops, capprops=capprops,order=order, ax = ax6)
# ax6.set_xticks([0, 20, 40, 60, 80, 100])


plt.subplots_adjust(hspace=0.5)
# plt.tight_layout()
#
pos = ax6.get_position()  # Get current position
ax6.set_position([pos.x0 , pos.y0 + 0.05, pos.width - 0.01, pos.height - 0.01])

pos = ax5.get_position()  # Get current position
ax5.set_position([pos.x0 - 0.025 , pos.y0 + 0.015, pos.width+0.015, pos.height + 0.015])

pos = ax4.get_position()  # Get current position
ax4.set_position([pos.x0 - 0.025 , pos.y0 + 0.015, pos.width+0.015, pos.height + 0.015])

pos = ax3.get_position()  # Get current position
ax3.set_position([pos.x0 - 0.025 , pos.y0 + 0.015, pos.width+0.015, pos.height + 0.015])

ax1.text(s = "(a)", fontsize = 11, y = 1.1, x = -0.35, transform=ax1.transAxes)
ax2.text(s = "(c)", fontsize = 11, y = 1.1, x = -0.35, transform=ax2.transAxes)
ax3.text(s = "(b)", fontsize = 11, y = 1.0, x = -0.2, transform=ax3.transAxes)
ax5.text(s = "(d)", fontsize = 11, y = 1.0, x = -0.2, transform=ax5.transAxes)
ax6.text(s = "(e)", fontsize = 11, y = 1.0, x = -0.3, transform=ax6.transAxes)

ax3.text(s = 'Co', x = 0.15, y = 0.99,  transform=ax3.transAxes)
ax3.text(s = 'Ni', x = 0.75, y = 1.0,  transform=ax3.transAxes)
ax3.text(s = 'Cr', x = 0.98, y = 0.35,  transform=ax3.transAxes)
ax3.text(s = 'Mn', x = 0.47, y = -0.08,  transform=ax3.transAxes)
ax3.text(s = 'Fe', x = -0.08, y = 0.35,  transform=ax3.transAxes)

ax4.text(s = 'Co', x = 0.15, y = 0.99,  transform=ax4.transAxes)
ax4.text(s = 'Ni', x = 0.75, y = 1.0,  transform=ax4.transAxes)
ax4.text(s = 'Cr', x = 0.98, y = 0.35,  transform=ax4.transAxes)
ax4.text(s = 'Mn', x = 0.47, y = -0.08,  transform=ax4.transAxes)
ax4.text(s = 'Fe', x = -0.08, y = 0.35,  transform=ax4.transAxes)

ax5.text(s = 'Co', x = 0.15, y = 0.99,  transform=ax5.transAxes)
ax5.text(s = 'Ni', x = 0.75, y = 1.0,  transform=ax5.transAxes)
ax5.text(s = 'Cr', x = 0.98, y = 0.35,  transform=ax5.transAxes)
ax5.text(s = 'Mn', x = 0.47, y = -0.08,  transform=ax5.transAxes)
ax5.text(s = 'Fe', x = -0.08, y = 0.35,  transform=ax5.transAxes)

plt.savefig("./fig4_revised.png", dpi = 200)