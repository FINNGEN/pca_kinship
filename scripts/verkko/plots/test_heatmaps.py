import heatmaps as hm
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

def test():
    n_points = 10000
    x_vals = np.exp(np.random.randn(n_points))
    y_vals = x_vals * 0 + np.exp(np.random.randn(n_points))

    fig, ax, cbar, im = hm.plot_counts(x_vals, y_vals)
    fname = os.path.join(os.path.dirname(__file__),
                         "gallery/heatmap_counts_example.pdf")
    plt.savefig(fname, format='pdf')

    fig, ax, cbar, im = hm.plot_prob_density(x_vals, y_vals)
    fname = os.path.join(os.path.dirname(__file__),
                         "gallery/heatmap_density_example.pdf")
    plt.savefig(fname, format='pdf')

    fig, ax, cbar, im = hm.plot_conditional_prob_density(x_vals, y_vals)
    fname = os.path.join(os.path.dirname(__file__),
                         "gallery/heatmap_cond_prob_example.pdf")
    plt.savefig(fname, format='pdf')


    #test subroutine plotting with gridspec
    f=plt.figure(figsize=(18/2.4,20/2.4))
    gs = gridspec.GridSpec(4,1)
    axes=[plt.subplot(gs[i]) for i in range(4)]

    fig, ax1, cbar, im = hm.plot_counts(x_vals,y_vals,ax=axes[0],use_gridspec=True)
    ax1.set_title('Counts')
    fig, ax2, cbar, im = hm.plot_prob_density(x_vals,y_vals,ax=axes[1],cmap='spring',use_gridspec=True)
    ax2.set_title('Prob density')
    fig, ax3, cbar, im = hm.plot_conditional_prob_density(x_vals,y_vals,ax=axes[2],cmap='jet',use_gridspec=True)
    ax3.set_title('Conditional prob density')

    #here the use_gridspec is false > the colorbar does not move when the grid layout is scaled
    fig, ax4, cbar, im = hm.plot_counts(x_vals,y_vals,ax=axes[3],use_gridspec=False)
    ax4.set_title('With use_gridspec=False')
    #update the layout, see how the last plot fails
    gs.update(left=0.09, right=0.98, bottom=0.05, top=0.95, hspace=0.75)

    fname = os.path.join(os.path.dirname(__file__),
                         "gallery/heatmap_subroutine_plot_example.pdf")
    plt.savefig(fname, format='pdf')


if __name__ == "__main__":
    test()
