"""
Plotting Utilities
==================

Helper utilities for consistent, publication-quality plotting styles
across all experiments.
"""

import matplotlib.pyplot as plt

def set_publication_style():
    """Sets a unified matplotlib style for all generated plots."""
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'lines.linewidth': 2,
        'figure.figsize': (10, 6),
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.format': 'png',
        'savefig.bbox': 'tight'
    })

def plot_tracking_error(t, errors, labels, title, filename):
    """Utility to plot tracking errors."""
    plt.figure()
    for err, label in zip(errors, labels):
        plt.plot(t, err, label=label)
    plt.yscale('log')
    plt.xlabel('Time (s)')
    plt.ylabel('||e(t)||')
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
