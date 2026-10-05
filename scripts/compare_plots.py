#!/usr/bin/env python3

######## Imports ########
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.animation as animation
import matplotlib.mlab as mlab

import numpy as np
import mcfacts.vis.LISA as li
import mcfacts.vis.PhenomA as pa
import pandas as pd
import os
import astropy.constants as const
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import poisson

# Grab those txt files
from importlib import resources as impresources
from mcfacts.vis import data
from mcfacts.vis import plotting
from mcfacts.vis import styles
from mcfacts.objects.snapshot import TxtSnapshotHandler, IniSnapshotHandler

mpl.rcParams['text.usetex'] = True

# Use the McFACTS plot style
#plt.style.use("mcfacts.vis.mplstyle.mcfacts_figures")
plt.style.use("mcfacts.vis.mplstyle.shawn_recoil_figs")

figsize = "apj_col"


######## Arg ########
def arg():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--plots-directory",
                        default=".",
                        type=str, help="directory to save plots")
    parser.add_argument("--sur-dir",
                        default="runs_sur",
                        type=str, help="output_dir of the surrogate run")
    parser.add_argument("--nosur-dir",
                        default="runs_nosur",
                        type=str, help="output_dir of the nosur run")
    parser.add_argument("--nosur-filter-dir",
                        default="runs_nosur_filter",
                        type=str, help="output_dir of the nosur filter run")
    parser.add_argument("--prec-dir",
                        default="runs_prec",
                        type=str, help="output_dir of the precession run")
    parser.add_argument("--fname-population",
                        default="population",
                        type=str, help="population snapshot file prefix in each run directory")
    opts = parser.parse_args()
    assert os.path.isdir(opts.sur_dir)
    assert os.path.isdir(opts.nosur_dir)
    assert os.path.isdir(opts.nosur_filter_dir)
    assert os.path.isdir(opts.prec_dir)
    return opts


# Old output_mergers_*.dat column layouts, so the column indices used below still apply
MERGER_COLS = ["galaxy_id", "bin_orb_a", "mass_final", "chi_eff", "spin_final",
               "spin_angle_final", "mass", "mass_2", "spin", "spin_2",
               "spin_angle", "spin_angle_2", "gen", "gen_2", "time_merged",
               "chi_p", "v_kick", "lum_shock", "lum_jet"]
LVK_COLS = ["galaxy_id", "time_merged", "bin_sep", lambda a: a.mass + a.mass_2, "bin_ecc",
            "gw_strain", "gw_freq", "gen", "gen_2"]
EMRI_COLS = ["galaxy_id", "time", "orb_a", "mass", "orb_ecc", "gw_strain", "gw_freq"]


def load_population(run_dir, file_name):
    """Load a run's population snapshot (as read by fiducial_plots.py) into
    (mergers, emris, lvk) arrays using the old output_mergers_*.dat column layouts.
    """
    agn_objects = TxtSnapshotHandler().load_cabinet(run_dir, file_name).agn_objects

    def columns(array_name, cols):
        if array_name not in agn_objects:  # load_cabinet skips empty files
            return np.empty((0, len(cols)))
        arr = agn_objects[array_name]
        return np.column_stack([np.asarray(col(arr) if callable(col) else getattr(arr, col), dtype=float)
                                for col in cols])

    return (columns("blackholes_merged", MERGER_COLS),
            columns("blackholes_emri", EMRI_COLS),
            columns("blackholes_lvk", LVK_COLS))


def linefunc(x, m):
    """Model for a line passing through (x,y) = (0,1).

    Function for a line used when fitting to the data.
    """
    return m * (x - 1)


def make_gen_masks(table, col1, col2):
    """Create masks for retrieving different sets of a merged or binary population based on generation.
    """
    # Column of generation data
    gen_obj1 = table[:, col1]
    gen_obj2 = table[:, col2]

    # Masks for hierarchical generations
    # g1 : all 1g-1g objects
    # g2 : 2g-1g and 2g-2g objects
    # g3 : >=3g-Ng (first object at least 3rd gen; second object any gen)
    # Pipe operator (|) = logical OR. (&)= logical AND.
    g1_mask = (gen_obj1 == 1) & (gen_obj2 == 1)
    g2_mask = ((gen_obj1 == 2) | (gen_obj2 == 2)) & ((gen_obj1 <= 2) & (gen_obj2 <= 2))
    gX_mask = (gen_obj1 >= 3) | (gen_obj2 >= 3)

    return g1_mask, g2_mask, gX_mask


######## Main ########
def main():
    # plt.style.use('seaborn-v0_8-poster')

    # Load data from output files
    opts = arg()

    sur_mergers, sur_emris, sur_lvk = load_population(opts.sur_dir, opts.fname_population)
    nosur_mergers, nosur_emris, nosur_lvk = load_population(opts.nosur_dir, opts.fname_population)
    nosur_filter_mergers, _, _ = load_population(opts.nosur_filter_dir, opts.fname_population)
    prec_mergers, prec_emris, prec_lvk = load_population(opts.prec_dir, opts.fname_population)
    
    # Exclude all rows with NaNs or zeros in the final mass column
    sur_merger_nan_mask = (np.isfinite(sur_mergers[:, 2])) & (sur_mergers[:, 2] != 0)
    sur_mergers = sur_mergers[sur_merger_nan_mask]
    
    # Exclude all rows with NaNs or zeros in the final mass column
    nosur_merger_nan_mask = (np.isfinite(nosur_mergers[:, 2])) & (nosur_mergers[:, 2] != 0)
    nosur_mergers = nosur_mergers[nosur_merger_nan_mask]

    # Exclude all rows with NaNs or zeros in the final mass column
    nosur_filter_merger_nan_mask = (np.isfinite(nosur_filter_mergers[:, 2])) & (nosur_filter_mergers[:, 2] != 0)
    nosur_filter_mergers = nosur_filter_mergers[nosur_filter_merger_nan_mask]

    # Exclude all rows with NaNs or zeros in the final mass column
    prec_merger_nan_mask = (np.isfinite(prec_mergers[:, 2])) & (prec_mergers[:, 2] != 0)
    prec_mergers = prec_mergers[prec_merger_nan_mask]

    sur_merger_g1_mask, sur_merger_g2_mask, sur_merger_gX_mask = make_gen_masks(sur_mergers, 12, 13)
    nosur_merger_g1_mask, nosur_merger_g2_mask, nosur_merger_gX_mask = make_gen_masks(nosur_mergers, 12, 13)
    nosur_filter_merger_g1_mask, nosur_filter_merger_g2_mask, nosur_filter_merger_gX_mask = make_gen_masks(nosur_filter_mergers, 12, 13)
    prec_merger_g1_mask, prec_merger_g2_mask, prec_merger_gX_mask = make_gen_masks(prec_mergers, 12, 13)


    # Ensure no union between sets
    assert all(sur_merger_g1_mask & sur_merger_g2_mask) == 0
    assert all(sur_merger_g1_mask & sur_merger_gX_mask) == 0
    assert all(sur_merger_g2_mask & sur_merger_gX_mask) == 0

    assert all(nosur_merger_g1_mask & nosur_merger_g2_mask) == 0
    assert all(nosur_merger_g1_mask & nosur_merger_gX_mask) == 0
    assert all(nosur_merger_g2_mask & nosur_merger_gX_mask) == 0

    assert all(nosur_filter_merger_g1_mask & nosur_filter_merger_g2_mask) == 0
    assert all(nosur_filter_merger_g1_mask & nosur_filter_merger_gX_mask) == 0
    assert all(nosur_filter_merger_g2_mask & nosur_filter_merger_gX_mask) == 0

    assert all(prec_merger_g1_mask & prec_merger_g2_mask) == 0
    assert all(prec_merger_g1_mask & prec_merger_gX_mask) == 0
    assert all(prec_merger_g2_mask & prec_merger_gX_mask) == 0

    # Ensure no elements are missed
    assert all(sur_merger_g1_mask | sur_merger_g2_mask | sur_merger_gX_mask) == 1
    assert all(nosur_merger_g1_mask | nosur_merger_g2_mask | nosur_merger_gX_mask) == 1
    assert all(nosur_filter_merger_g1_mask | nosur_filter_merger_g2_mask | nosur_filter_merger_gX_mask) == 1
    assert all(prec_merger_g1_mask | prec_merger_g2_mask | prec_merger_gX_mask) == 1

    print("success")

    # ========================================
    # NOSUR - Number of Mergers vs Mass
    # ========================================

    # Plot intial and final mass distributions
    # figsize=plotting.set_size(figsize)
    fig, ax = plt.subplots(2, 1, figsize=plotting.set_size(figsize), sharex=True)

    nosur = ax[0]
    sur = ax[1]

    # Plot intial and final mass distributions
    counts, bins = np.histogram(nosur_mergers[:, 2])
    # plt.hist(bins[:-1], bins, weights=counts)
    bins = np.arange(int(nosur_mergers[:, 2].min()), int(nosur_mergers[:, 2].max()) + 2, 1)

    # masking generation values for histogram data
    nosur_hist_data = [nosur_mergers[:, 2][nosur_merger_g1_mask], nosur_mergers[:, 2][nosur_merger_g2_mask],
                       nosur_mergers[:, 2][nosur_merger_gX_mask]]
    hist_label = ['1g-1g', '2g-1g or 2g-2g', r'$\geq$3g-Ng']
    hist_color = [styles.color_gen1, styles.color_gen2, styles.color_genX]

    nosur.hist(nosur_hist_data, bins=bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label,
               stacked=True)
    nosur_hist_data_int = list(map(int, nosur_hist_data[0]))
    mode, count = stats.mode(nosur_hist_data_int, axis=None, keepdims=False)

    nosur.set_ylabel(r'Number of Mergers' '\n' r'McFACTSsur', fontsize=5, wrap=True)
    nosur.set(
        yticks=(np.linspace(int(count * 0.30), int(count * 0.90), 2)),
        # ylim=(-5,max(counts)),
        xlabel=(r'Remnant Mass [$M_\odot$]'),
        xscale=('log'),
        xticks=(np.geomspace(int(nosur_mergers[:, 2].min()), int(nosur_mergers[:, 2].max()), 5).astype(int))
    )

    nosur.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:.0f}'))
    nosur.xaxis.set_minor_formatter(mticker.NullFormatter())
    nosur.tick_params(axis='x', direction='in', which='both')
    # plt.grid(True, color='gray', ls='dashed')
    nosur.yaxis.grid(True, color='gray', ls='dashed')

    # svf_ax = plt.gca()
    # svf_ax.set_axisbelow(True)
    nosur.tick_params(axis='x', direction='in', which='both')
    # nosur.grid(True, color='gray', ls='dashed', alpha=0.4)
    nosur.yaxis.grid(True, color='gray', ls='dashed')

    if figsize == 'apj_col':
        nosur.legend(fontsize=5)
    elif figsize == 'apj_page':
        nosur.legend()

    # ========================================
    # SUR - Number of Mergers vs Mass
    # ========================================

    counts, bins = np.histogram(sur_mergers[:, 2])
    # plt.hist(bins[:-1], bins, weights=counts)
    bins = np.arange(int(sur_mergers[:, 2].min()), int(sur_mergers[:, 2].max()) + 2, 1)

    sur_hist_data = [sur_mergers[:, 2][sur_merger_g1_mask], sur_mergers[:, 2][sur_merger_g2_mask],
                     sur_mergers[:, 2][sur_merger_gX_mask]]
    hist_label = ['1g-1g', '2g-1g or 2g-2g', r'$\geq$3g-Ng']
    hist_color = [styles.color_gen1, styles.color_gen2, styles.color_genX]

    sur.hist(sur_hist_data, bins=bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label,
             stacked=True)
    sur_hist_data_int = list(map(int, sur_hist_data[0]))
    mode, count = stats.mode(sur_hist_data_int, axis=None, keepdims=False)

    sur.set_ylabel(r'Number of Mergers' '\n' r'McFACTSsur', fontsize=5, wrap=True)
    sur.set(
        yticks=(np.linspace(int(count * 0.30), int(count * 0.90), 2)),
        # ylim=(-5,max(counts)),
        xlabel=(r'Remnant Mass [$M_\odot$]'),
        xscale=('log'),
        xticks=(np.geomspace(int(sur_mergers[:, 2].min()), int(sur_mergers[:, 2].max()), 5).astype(int))
    )

    sur.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:.0f}'))
    sur.xaxis.set_minor_formatter(mticker.NullFormatter())
    sur.tick_params(axis='x', direction='in', which='both')
    # plt.grid(True, color='gray', ls='dashed')
    sur.yaxis.grid(True, color='gray', ls='dashed')

    # svf_ax = sur.gca()
    # svf_ax.set_axisbelow(True)
    # plt.xticks(np.geomspace(20, 200, 5).astype(int))

    # if figsize == 'apj_col':
    #    sur.legend(fontsize=5)
    # elif figsize == 'apj_page':
    #    sur.legend()

    # sur.savefig(opts.plots_directory + r"/merger_remnant_mass.png", format='png')
    plt.savefig(opts.plots_directory + r"/merger_remnant_mass.png", format='png')
    # plt.show()

    # ========================================
    # SUR - Merger Mass vs Radius
    # ========================================

    # Retrieve the migration trap radius used in run
    trap_radius = IniSnapshotHandler().load_settings(opts.sur_dir, "settings").disk_radius_trap

    # plt.title('Migration Trap influence')
    for i in range(len(sur_mergers[:, 1])):
        if sur_mergers[i, 1] < 10.0:
            sur_mergers[i, 1] = 10.0

    # Separate generational subpopulations
    sur_gen1_orb_a = sur_mergers[:, 1][sur_merger_g1_mask]
    sur_gen2_orb_a = sur_mergers[:, 1][sur_merger_g2_mask]
    sur_genX_orb_a = sur_mergers[:, 1][sur_merger_gX_mask]
    sur_gen1_mass = sur_mergers[:, 2][sur_merger_g1_mask]
    sur_gen2_mass = sur_mergers[:, 2][sur_merger_g2_mask]
    sur_genX_mass = sur_mergers[:, 2][sur_merger_gX_mask]

    fig, ax = plt.subplots(1, 2, figsize=(5, 2), sharey=True, layout='constrained',
                           gridspec_kw={'wspace': 0, 'hspace': 0})
    sur = ax[1]
    nosur = ax[0]

    sur.scatter(sur_gen1_orb_a, sur_gen1_mass,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolors="none",
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    sur.scatter(sur_gen2_orb_a, sur_gen2_mass,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolors="none",
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    sur.scatter(sur_genX_orb_a, sur_genX_mass,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolors="none",
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.axvline(trap_radius, color='k', linestyle='--', zorder=0,
                label=f'Trap Radius = {trap_radius:.0f} ' + r'$R_g$')

    # plt.text(650, 602, 'Migration Trap', rotation='vertical', size=18, fontweight='bold')
    # sur.set_ylabel(r'Remnant Mass [$M_\odot$]')
    sur.set_xlabel(r'Radius$_{sur}$ [$R_g$]')
    sur.set_xscale('log')
    sur.set_yscale('log')

    if figsize == 'apj_col':
        sur.legend(fontsize=5)
    elif figsize == 'apj_page':
        sur.legend()

    sur.set_ylim(18, 1000)

    # svf_ax = sur.gca()
    # svf_ax.set_axisbelow(True)
    sur.grid(True, color='gray', ls='dashed')
    # sur.savefig(opts.plots_directory + "/merger_mass_v_radius.png", format='png')
    # sur.close()

    # ========================================
    # NOSUR - Merger Mass vs Radius
    # ========================================

    # Read the log file
    # log_data = ReadLog(opts.fname_log)

    # Retrieve the migration trap radius used in run
    # trap_radius = log_data["disk_radius_trap"]

    # plt.title('Migration Trap influence')
    for i in range(len(nosur_mergers[:, 1])):
        if nosur_mergers[i, 1] < 10.0:
            nosur_mergers[i, 1] = 10.0

    # Separate generational subpopulations
    nosur_gen1_orb_a = nosur_mergers[:, 1][nosur_merger_g1_mask]
    nosur_gen2_orb_a = nosur_mergers[:, 1][nosur_merger_g2_mask]
    nosur_genX_orb_a = nosur_mergers[:, 1][nosur_merger_gX_mask]
    nosur_gen1_mass = nosur_mergers[:, 2][nosur_merger_g1_mask]
    nosur_gen2_mass = nosur_mergers[:, 2][nosur_merger_g2_mask]
    nosur_genX_mass = nosur_mergers[:, 2][nosur_merger_gX_mask]

    nosur.scatter(nosur_gen1_orb_a, nosur_gen1_mass,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolors="none",
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    nosur.scatter(nosur_gen2_orb_a, nosur_gen2_mass,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolors="none",
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    nosur.scatter(nosur_genX_orb_a, nosur_genX_mass,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolors="none",
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    nosur.axvline(trap_radius, color='k', linestyle='--', zorder=0,
                  label=f'Trap Radius = {trap_radius:.0f} ' + r'$R_g$')

    # plt.text(650, 602, 'Migration Trap', rotation='vertical', size=18, fontweight='bold')
    nosur.set_ylabel(r'Remnant Mass [$M_\odot$]')
    nosur.set_xlabel(r'Radius$_{nosur}$ [$R_g$]')
    nosur.set_xscale('log')
    nosur.set_yscale('log')

    # if figsize == 'apj_col':
    #    nosur.legend(fontsize=5)
    # elif figsize == 'apj_page':
    #    nosur.legend()

    nosur.set_ylim(18, 1000)

    # svf_ax = nosur.gca()
    # svf_ax.set_axisbelow(True)
    nosur.grid(True, color='gray', ls='dashed')
    plt.savefig(opts.plots_directory + "/merger_mass_v_radius.png", format='png')
    # plt.show()

    # ========================================
    # SUR - q vs Chi Effective
    # ========================================

    # retrieve component masses and mass ratio
    m1 = np.zeros(sur_mergers.shape[0])
    m2 = np.zeros(sur_mergers.shape[0])
    mass_ratio = np.zeros(sur_mergers.shape[0])
    for i in range(sur_mergers.shape[0]):
        if sur_mergers[i, 6] < sur_mergers[i, 7]:
            m1[i] = sur_mergers[i, 7]
            m2[i] = sur_mergers[i, 6]
            mass_ratio[i] = sur_mergers[i, 6] / sur_mergers[i, 7]
        else:
            mass_ratio[i] = sur_mergers[i, 7] / sur_mergers[i, 6]
            m1[i] = sur_mergers[i, 6]
            m2[i] = sur_mergers[i, 7]
        # mass_ratio[i] = 1 / mass_ratio[i]

    # (q,X_eff) Figure details here:
    # Want to highlight higher generation mergers on this plot
    sur_chi_eff = sur_mergers[:, 3]

    # Get 1g-1g population
    sur_gen1_chi_eff = sur_chi_eff[sur_merger_g1_mask]
    sur_gen1_mass_ratio = mass_ratio[sur_merger_g1_mask]
    # 2g-1g and 2g-2g population
    sur_gen2_chi_eff = sur_chi_eff[sur_merger_g2_mask]
    sur_gen2_mass_ratio = mass_ratio[sur_merger_g2_mask]
    # >=3g-Ng population (i.e., N=1,2,3,4,...)
    sur_genX_chi_eff = sur_chi_eff[sur_merger_gX_mask]
    sur_genX_mass_ratio = mass_ratio[sur_merger_gX_mask]
    # all 2+g mergers; H = hierarchical
    sur_genH_chi_eff = sur_chi_eff[(sur_merger_g2_mask + sur_merger_gX_mask)]
    sur_genH_mass_ratio = mass_ratio[(sur_merger_g2_mask + sur_merger_gX_mask)]

    # points for plotting line fit
    x = np.linspace(-1, 1, num=2)

    # fit the hierarchical mergers (any binaries with 2+g) to a line passing through 0,1
    # popt contains the model parameters, pcov the covariances
    # poptHigh, pcovHigh = curve_fit(linefunc, high_gen_mass_ratio, high_gen_chi_eff)

    # plot the 1g-1g population
    fig, ax = plt.subplots(1, 2, figsize=(5, 3), sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0})
    # ax2 = fig.add_subplot(111)
    sur = ax[1]
    nosur = ax[0]

    # GW231123 error bar additions - mainly used for paper and thesis
    # comment section out if not wanting error bars
    chi = .3  # +0.2 -0.4

    m1 = 137.0  # +23 -18
    m2 = 101.0  # +22 -51

    q = m2 / m1

    # upper
    e_q_up = 1.0 / m1 * (22.0 ** 2.0 + (q * 23.0) ** 2.0) ** (0.5)

    # lower
    e_q_low = 1.0 / m1 * ((-51.0) ** 2.0 + (q * (-18.0)) ** 2.0) ** (0.5)

    e_q = [[e_q_low], [e_q_up]]
    e_chi = [[0.4], [0.2]]

    # GW231123 addition to plot
    sur.errorbar(chi, q, xerr=e_chi, yerr=e_q, ecolor='c', label="GW231123")

    # end comment for 231123 section

    # 1g-1g mergers
    sur.scatter(sur_gen1_chi_eff, sur_gen1_mass_ratio,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(sur_gen2_chi_eff, sur_gen2_mass_ratio,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(sur_genX_chi_eff, sur_genX_mass_ratio,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolor='none',
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    if len(sur_genH_chi_eff) > 0:
        poptHier, pcovHier = curve_fit(linefunc, sur_genH_mass_ratio, sur_genH_chi_eff)
        errHier = np.sqrt(np.diag(pcovHier))[0]
        # plot the line fitting the hierarchical mergers
        sur.plot(linefunc(x, *poptHier), x,
                 ls='dashed',
                 lw=1,
                 color='gray',
                 zorder=3,
                 label=r'$d\chi/dq(\geq$2g)=' +
                       f'{poptHier[0]:.2f}' +
                       r'$\pm$' + f'{errHier:.2f}'
                 )
        #         #  alpha=linealpha,

    if len(sur_chi_eff) > 0:
        poptAll, pcovAll = curve_fit(linefunc, mass_ratio, sur_chi_eff)
        errAll = np.sqrt(np.diag(pcovAll))[0]
        sur.plot(linefunc(x, *poptAll), x,
                 ls='solid',
                 lw=1,
                 color='black',
                 zorder=3,
                 label=r'$d\chi/dq$(all)=' +
                       f'{poptAll[0]:.2f}' +
                       r'$\pm$' + f'{errAll:.2f}'
                 )
        #  alpha=linealpha,

    # sur.set_ylabel(r'$q = M_2 / M_1$')  # ($M_1 > M_2$)')
    sur.set_xlabel(r'$\chi_{\rm eff}^{sur}$')
    sur.set_ylim(0, 1)
    sur.set_xlim(-1, 1)
    sur.set_axisbelow = True

    # if figsize == 'apj_col':
    #    sur.legend(loc='lower left', fontsize=5)
    # elif figsize == 'apj_page':
    #    sur.legend(loc='lower left')

    sur.grid('on', color='gray', ls='dotted')
    # plt.savefig(opts.plots_directory + "/q_chi_eff_231123.png", format='png')  # ,dpi=600)

    # ========================================
    # NOSUR - q vs Chi Effective
    # ========================================

    # retrieve component masses and mass ratio
    m1 = np.zeros(nosur_mergers.shape[0])
    m2 = np.zeros(nosur_mergers.shape[0])
    mass_ratio = np.zeros(nosur_mergers.shape[0])
    for i in range(nosur_mergers.shape[0]):
        if nosur_mergers[i, 6] < nosur_mergers[i, 7]:
            m1[i] = nosur_mergers[i, 7]
            m2[i] = nosur_mergers[i, 6]
            mass_ratio[i] = nosur_mergers[i, 6] / nosur_mergers[i, 7]
        else:
            mass_ratio[i] = nosur_mergers[i, 7] / nosur_mergers[i, 6]
            m1[i] = nosur_mergers[i, 6]
            m2[i] = nosur_mergers[i, 7]

    # (q,X_eff) Figure details here:
    # Want to highlight higher generation mergers on this plot
    nosur_chi_eff = nosur_mergers[:, 3]

    # Get 1g-1g population
    nosur_gen1_chi_eff = nosur_chi_eff[nosur_merger_g1_mask]
    nosur_gen1_mass_ratio = mass_ratio[nosur_merger_g1_mask]
    # 2g-1g and 2g-2g population
    nosur_gen2_chi_eff = nosur_chi_eff[nosur_merger_g2_mask]
    nosur_gen2_mass_ratio = mass_ratio[nosur_merger_g2_mask]
    # >=3g-Ng population (i.e., N=1,2,3,4,...)
    nosur_genX_chi_eff = nosur_chi_eff[nosur_merger_gX_mask]
    nosur_genX_mass_ratio = mass_ratio[nosur_merger_gX_mask]
    # all 2+g mergers; H = hierarchical
    nosur_genH_chi_eff = nosur_chi_eff[(nosur_merger_g2_mask + nosur_merger_gX_mask)]
    nosur_genH_mass_ratio = mass_ratio[(nosur_merger_g2_mask + nosur_merger_gX_mask)]

    # points for plotting line fit
    x = np.linspace(-1, 1, num=2)

    # fit the hierarchical mergers (any binaries with 2+g) to a line passing through 0,1
    # popt contains the model parameters, pcov the covariances
    # poptHigh, pcovHigh = curve_fit(linefunc, high_gen_mass_ratio, high_gen_chi_eff)

    # plot the 1g-1g population
    # fig = plt.figure(figsize=(plotting.set_size(figsize)[0], 2.8))
    # ax2 = fig.add_subplot(111)
    # 1g-1g mergers

    # adding in GW231123 error bars for NOsur plot
    nosur.errorbar(chi, q, xerr=e_chi, yerr=e_q, ecolor='c', label="GW231123")

    nosur.scatter(nosur_gen1_chi_eff, nosur_gen1_mass_ratio,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolor='none',
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    # plot the 2g+ mergers
    nosur.scatter(nosur_gen2_chi_eff, nosur_gen2_mass_ratio,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolor='none',
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    # plot the 3g+ mergers
    nosur.scatter(nosur_genX_chi_eff, nosur_genX_mass_ratio,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolor='none',
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    if len(nosur_genH_chi_eff) > 0:
        poptHier, pcovHier = curve_fit(linefunc, nosur_genH_mass_ratio, nosur_genH_chi_eff)
        errHier = np.sqrt(np.diag(pcovHier))[0]
        # plot the line fitting the hierarchical mergers
        nosur.plot(linefunc(x, *poptHier), x,
                   ls='dashed',
                   lw=1,
                   color='gray',
                   zorder=3,
                   label=r'$d\chi/dq(\geq$2g)=' +
                         f'{poptHier[0]:.2f}' +
                         r'$\pm$' + f'{errHier:.2f}'
                   )
        #         #  alpha=linealpha,

    if len(nosur_chi_eff) > 0:
        poptAll, pcovAll = curve_fit(linefunc, mass_ratio, nosur_chi_eff)
        errAll = np.sqrt(np.diag(pcovAll))[0]
        nosur.plot(linefunc(x, *poptAll), x,
                   ls='solid',
                   lw=1,
                   color='black',
                   zorder=3,
                   label=r'$d\chi/dq$(all)=' +
                         f'{poptAll[0]:.2f}' +
                         r'$\pm$' + f'{errAll:.2f}'
                   )
        #  alpha=linealpha,

    nosur.set(
        ylabel=(r'$q = M_2 / M_1$'),  # ($M_1 > M_2$)'),
        xlabel=(r'$\chi_{\rm eff}^{nosur}$'),
        ylim=(0, 1),
        xlim=(-1, 1),
        axisbelow=True
    )

    if figsize == 'apj_col':
        nosur.legend(loc='best', fontsize=4)
    elif figsize == 'apj_page':
        nosur.legend(loc='best')

    nosur.grid('on', color='gray', ls='dotted')
    plt.savefig(opts.plots_directory + "/q_chi_eff.png", format='png')  # ,dpi=600)
    # plt.show()

    # ========================================
    # NOSUR - q vs Chi Effective                   | STAND ALONE PLOT |
    # ========================================

    # retrieve component masses and mass ratio
    m1 = np.zeros(nosur_mergers.shape[0])
    m2 = np.zeros(nosur_mergers.shape[0])
    mass_ratio = np.zeros(nosur_mergers.shape[0])
    for i in range(nosur_mergers.shape[0]):
        if nosur_mergers[i, 6] < nosur_mergers[i, 7]:
            m1[i] = nosur_mergers[i, 7]
            m2[i] = nosur_mergers[i, 6]
            mass_ratio[i] = nosur_mergers[i, 6] / nosur_mergers[i, 7]
        else:
            mass_ratio[i] = nosur_mergers[i, 7] / nosur_mergers[i, 6]
            m1[i] = nosur_mergers[i, 6]
            m2[i] = nosur_mergers[i, 7]

    # (q,X_eff) Figure details here:
    # Want to highlight higher generation mergers on this plot
    nosur_chi_eff = nosur_mergers[:, 3]

    # Get 1g-1g population
    nosur_gen1_chi_eff = nosur_chi_eff[nosur_merger_g1_mask]
    nosur_gen1_mass_ratio = mass_ratio[nosur_merger_g1_mask]
    # 2g-1g and 2g-2g population
    nosur_gen2_chi_eff = nosur_chi_eff[nosur_merger_g2_mask]
    nosur_gen_mass_ratio = mass_ratio[nosur_merger_g2_mask]
    # >=3g-Ng population (i.e., N=1,2,3,4,...)
    nosur_genX_chi_eff = nosur_chi_eff[nosur_merger_gX_mask]
    nosur_genX_mass_ratio = mass_ratio[nosur_merger_gX_mask]
    # all 2+g mergers; H = hierarchical
    nosur_genH_chi_eff = nosur_chi_eff[(nosur_merger_g2_mask + nosur_merger_gX_mask)]
    nosur_genH_mass_ratio = mass_ratio[(nosur_merger_g2_mask + nosur_merger_gX_mask)]

    # points for plotting line fit
    x = np.linspace(-1, 1, num=2)

    # fit the hierarchical mergers (any binaries with 2+g) to a line passing through 0,1
    # popt contains the model parameters, pcov the covariances
    # poptHigh, pcovHigh = curve_fit(linefunc, high_gen_mass_ratio, high_gen_chi_eff)

    fig, ax = plt.subplots(1, 1, figsize=(4, 3), sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0})
    # ax2 = fig.add_subplot(111)
    # sur = ax[1]
    nosur = ax

    # 1g-1g mergers
    nosur.scatter(nosur_gen1_chi_eff, nosur_gen1_mass_ratio,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolor='none',
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    # plot the 2g+ mergers
    nosur.scatter(nosur_gen2_chi_eff, nosur_gen2_mass_ratio,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolor='none',
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    # plot the 3g+ mergers
    nosur.scatter(nosur_genX_chi_eff, nosur_genX_mass_ratio,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolor='none',
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    if len(nosur_genH_chi_eff) > 0:
        poptHier, pcovHier = curve_fit(linefunc, nosur_genH_mass_ratio, nosur_genH_chi_eff)
        errHier = np.sqrt(np.diag(pcovHier))[0]
        # plot the line fitting the hierarchical mergers
        nosur.plot(linefunc(x, *poptHier), x,
                   ls='dashed',
                   lw=1,
                   color='gray',
                   zorder=3,
                   label=r'$d\chi/dq(\geq$2g)=' +
                         f'{poptHier[0]:.2f}' +
                         r'$\pm$' + f'{errHier:.2f}'
                   )
        #         #  alpha=linealpha,

    if len(nosur_chi_eff) > 0:
        poptAll, pcovAll = curve_fit(linefunc, mass_ratio, nosur_chi_eff)
        errAll = np.sqrt(np.diag(pcovAll))[0]
        nosur.plot(linefunc(x, *poptAll), x,
                   ls='solid',
                   lw=1,
                   color='black',
                   zorder=3,
                   label=r'$d\chi/dq$(all)=' +
                         f'{poptAll[0]:.2f}' +
                         r'$\pm$' + f'{errAll:.2f}'
                   )
        #  alpha=linealpha,

    nosur.set_ylabel(r'$q = M_2 / M_1$')  # ($M_1 > M_2$)')
    nosur.set_xlabel(r'$\chi_{\rm eff}^{nosur}$')
    nosur.set_ylim(0, 1)
    nosur.set_xlim(-1, 1)
    nosur.set_axisbelow = True

    if figsize == 'apj_col':
        nosur.legend(loc='best', fontsize=4)
    elif figsize == 'apj_page':
        nosur.legend(loc='best')

    nosur.grid('on', color='gray', ls='dotted')
    # plt.savefig(opts.plots_directory + "/q_chi_eff_nosur.png", format='png')  # ,dpi=600)
    # plt.show()

    # ========================================
    # SUR - Chi Effective vs Disk Radius
    # ========================================

    # plot the 1g-1g population
    fig, ax = plt.subplots(1, 2, figsize=(5, 3), sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0})
    # ax2 = fig.add_subplot(111)
    sur = ax[1]
    nosur = ax[0]

    # 1g-1g mergers
    sur.scatter(sur_gen1_orb_a, sur_gen1_chi_eff,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(sur_gen2_orb_a, sur_gen2_chi_eff,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(sur_genX_orb_a, sur_genX_chi_eff,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolor='none',
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set_xlabel(r'$Radius_{sur} [R_g]$')
    # sur.set_ylabel(r'$\chi_{\rm eff}$')
    sur.set_xscale('log')
    sur.set_ylim(-0.4, 1)
    sur.set_axisbelow = True

    # if figsize == 'apj_col':
    #    sur.legend(loc='lower left', fontsize=4)
    # elif figsize == 'apj_page':
    #    sur.legend(loc='lower left')

    sur.grid('on', color='gray', ls='dotted')
    # plt.savefig(opts.plots_directory + "/q_chi_eff.png", format='png')  # ,dpi=600)

    # ========================================
    # NOSUR - Chi Effective vs Disk Radius
    # ========================================

    # 1g-1g mergers
    nosur.scatter(nosur_gen1_orb_a, nosur_gen1_chi_eff,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolor='none',
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    # plot the 2g+ mergers
    nosur.scatter(nosur_gen2_orb_a, nosur_gen2_chi_eff,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolor='none',
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    # plot the 3g+ mergers
    nosur.scatter(nosur_genX_orb_a, nosur_genX_chi_eff,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolor='none',
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    nosur.set_xlabel(r'$Radius_{nosur} [R_g]$')
    nosur.set_ylabel(r'$\chi_{\rm eff}$')
    nosur.set_xscale('log')
    nosur.set_ylim(-0.4, 1)
    nosur.set_axisbelow = True

    if figsize == 'apj_col':
        nosur.legend(loc='lower left', fontsize=5)
    elif figsize == 'apj_page':
        nosur.legend(loc='lower left')

    nosur.grid('on', color='gray', ls='dotted')
    plt.savefig(opts.plots_directory + "/chi_eff_radius.png", format='png')  # ,dpi=600)

    # ========================================
    # SUR - Chi_p vs Disk Radius
    # ========================================

    # Can break out higher mass Chi_p events as test/illustration.
    # Set up default arrays for high mass BBH (>40Msun say) to overplot vs chi_p.
    sur_chi_p = sur_mergers[:, 15]
    sur_gen1_chi_p = sur_chi_p[sur_merger_g1_mask]
    sur_gen2_chi_p = sur_chi_p[sur_merger_g2_mask]
    sur_genX_chi_p = sur_chi_p[sur_merger_gX_mask]

    fig, ax = plt.subplots(1, 2, figsize=(5, 3), sharey=True)
    # ax1 = fig.add_subplot(111)

    sur = ax[1]
    nosur = ax[0]

    sur.scatter(sur_gen1_orb_a, sur_gen1_chi_p,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g')

    # plot the 2g+ mergers
    sur.scatter(sur_gen2_orb_a, sur_gen2_chi_p,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g')

    # plot the 3g+ mergers
    sur.scatter(sur_genX_orb_a, sur_genX_chi_p,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolor='none',
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng')

    sur.axvline(trap_radius, color='k', linestyle='--', zorder=0,
                label=f'Trap Radius = {trap_radius:.0f} ' + r'$R_g$')

    # plt.title("In-plane effective Spin vs. Merger radius")
    sur.set(
        # ylabel=r'$\chi_{\rm p}$',
        xlabel=r'Radius$_{sur}$ [$R_g$]',
        xscale='log',
        ylim=(0, 1),
        axisbelow=True)

    sur.grid(True, color='gray', ls='dashed')

    if figsize == 'apj_col':
        sur.legend(fontsize=5, loc='upper right')
    elif figsize == 'apj_page':
        sur.legend()

    # plt.savefig(opts.plots_directory + "/r_chi_p.png", format='png')
    # plt.close()

    # plt.figure()
    # index = 2
    # mode = 10
    # pareto = (np.random.pareto(index, 1000) + 1) * mode

    # x = np.linspace(1,100)
    # p = index*mode**index / x**(index+1)

    # # count, bins, _ = plt.hist(pareto, 100)
    # plt.plot(x, p)
    # plt.xlim(0,100)
    # plt.show()

    # ========================================
    # NOSUR - Chi_p vs Disk Radius
    # ========================================

    # Can break out higher mass Chi_p events as test/illustration.
    # Set up default arrays for high mass BBH (>40Msun say) to overplot vs chi_p.
    nosur_chi_p = nosur_mergers[:, 15]
    nosur_gen1_chi_p = nosur_chi_p[nosur_merger_g1_mask]
    nosur_gen2_chi_p = nosur_chi_p[nosur_merger_g2_mask]
    nosur_genX_chi_p = nosur_chi_p[nosur_merger_gX_mask]

    # fig = plt.figure(figsize=plotting.set_size(figsize))
    # ax1 = fig.add_subplot(111)

    nosur.scatter(nosur_gen1_orb_a, nosur_gen1_chi_p,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolor='none',
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g')

    # plot the 2g+ mergers
    nosur.scatter(nosur_gen2_orb_a, nosur_gen2_chi_p,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolor='none',
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g')

    # plot the 3g+ mergers
    nosur.scatter(nosur_genX_orb_a, nosur_genX_chi_p,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolor='none',
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng')

    nosur.axvline(trap_radius, color='k', linestyle='--', zorder=0,
                  label=f'Trap Radius = {trap_radius:.0f} ' + r'$R_g$')

    # plt.title("In-plane effective Spin vs. Merger radius")
    nosur.set(
        ylabel=r'$\chi_{\rm p}$',
        xlabel=r'Radius$_{nosur}$ [$R_g$]',
        xscale='log',
        ylim=(0, 1),
        axisbelow=True)

    nosur.grid(True, color='gray', ls='dashed')

    # if figsize == 'apj_col':
    #    nosur.legend(fontsize=5, loc='upper left')
    # elif figsize == 'apj_page':
    #    nosur.legend()

    plt.savefig(opts.plots_directory + "/chi_p_radius.png", format='png')
    # plt.show()

    # plt.figure()
    # index = 2
    # mode = 10
    # pareto = (np.random.pareto(index, 1000) + 1) * mode

    # x = np.linspace(1,100)
    # p = index*mode**index / x**(index+1)

    # # count, bins, _ = plt.hist(pareto, 100)
    # plt.plot(x, p)
    # plt.xlim(0,100)
    # plt.show()

    # ========================================
    # SUR - Time of Merger vs Remnant Mass
    # ========================================

    sur_all_time = sur_mergers[:, 14]
    sur_gen1_time = sur_all_time[sur_merger_g1_mask]
    sur_gen2_time = sur_all_time[sur_merger_g2_mask]
    sur_genX_time = sur_all_time[sur_merger_gX_mask]

    fig, ax = plt.subplots(1, 2, figsize=plotting.set_size(figsize), sharey=True)
    # ax3 = fig.add_subplot(111)

    sur = ax[1]
    nosur = ax[0]

    # plt.title("Time of Merger after AGN Onset")
    # ax3.scatter(mergers[:,14]/1e6, mergers[:,2], s=pointsize_merge_time, color='darkolivegreen')
    sur.scatter(sur_gen1_time / 1e6, sur_gen1_mass,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(sur_gen2_time / 1e6, sur_gen2_mass,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(sur_genX_time / 1e6, sur_genX_mass,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolor='none',
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set(
        xlabel=r'Time$_{sur}$ [Myr]',
        # ylabel=r'Remnant Mass [$M_\odot$]',
        yscale="log",
        axisbelow=True
    )

    plt.grid(True, color='gray', ls='dashed')

    # if figsize == 'apj_col':
    #    sur.legend(fontsize=5)
    # elif figsize == 'apj_page':
    #    sur.legend()

    # plt.savefig(opts.plots_directory + '/time_of_merger.png', format='png')
    # plt.close()

    # ========================================
    # NOSUR - Time of Merger vs Remnant Mass
    # ========================================

    nosur_all_time = nosur_mergers[:, 14]
    nosur_gen1_time = nosur_all_time[nosur_merger_g1_mask]
    nosur_gen2_time = nosur_all_time[nosur_merger_g2_mask]
    nosur_genX_time = nosur_all_time[nosur_merger_gX_mask]

    # fig = plt.figure(figsize=plotting.set_size(figsize))
    # ax3 = fig.add_subplot(111)

    # plt.title("Time of Merger after AGN Onset")
    # ax3.scatter(mergers[:,14]/1e6, mergers[:,2], s=pointsize_merge_time, color='darkolivegreen')
    nosur.scatter(nosur_gen1_time / 1e6, nosur_gen1_mass,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolor='none',
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    # plot the 2g+ mergers
    nosur.scatter(nosur_gen2_time / 1e6, nosur_gen2_mass,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolor='none',
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    # plot the 3g+ mergers
    nosur.scatter(nosur_genX_time / 1e6, nosur_genX_mass,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolor='none',
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    nosur.set(
        xlabel=r'Time$_{nosur}$ [Myr]',
        ylabel=r'Remnant Mass [$M_\odot$]',
        yscale="log",
        axisbelow=True
    )

    plt.grid(True, color='gray', ls='dashed')

    if figsize == 'apj_col':
        nosur.legend(fontsize=5, loc='upper left')
    elif figsize == 'apj_page':
        nosur.legend()

    plt.savefig(opts.plots_directory + '/time_of_merger.png', format='png')
    # plt.show()

    # ========================================
    # SUR - Mass 1 vs Mass 2
    # ========================================

    # Sort Objects into Mass 1 and Mass 2 by generation
    sur_mass_mask_g1 = sur_mergers[sur_merger_g1_mask, 6] > sur_mergers[sur_merger_g1_mask, 7]
    sur_gen1_mass_1 = np.zeros(np.sum(sur_merger_g1_mask))
    sur_gen1_mass_1[sur_mass_mask_g1] = sur_mergers[sur_merger_g1_mask, 6][sur_mass_mask_g1]
    sur_gen1_mass_1[~sur_mass_mask_g1] = sur_mergers[sur_merger_g1_mask, 7][~sur_mass_mask_g1]
    sur_gen1_mass_2 = np.zeros(np.sum(sur_merger_g1_mask))
    sur_gen1_mass_2[~sur_mass_mask_g1] = sur_mergers[sur_merger_g1_mask, 6][~sur_mass_mask_g1]
    sur_gen1_mass_2[sur_mass_mask_g1] = sur_mergers[sur_merger_g1_mask, 7][sur_mass_mask_g1]

    sur_mass_mask_g2 = sur_mergers[sur_merger_g2_mask, 6] > sur_mergers[sur_merger_g2_mask, 7]
    sur_gen2_mass_1 = np.zeros(np.sum(sur_merger_g2_mask))
    sur_gen2_mass_1[sur_mass_mask_g2] = sur_mergers[sur_merger_g2_mask, 6][sur_mass_mask_g2]
    sur_gen2_mass_1[~sur_mass_mask_g2] = sur_mergers[sur_merger_g2_mask, 7][~sur_mass_mask_g2]
    sur_gen2_mass_2 = np.zeros(np.sum(sur_merger_g2_mask))
    sur_gen2_mass_2[~sur_mass_mask_g2] = sur_mergers[sur_merger_g2_mask, 6][~sur_mass_mask_g2]
    sur_gen2_mass_2[sur_mass_mask_g2] = sur_mergers[sur_merger_g2_mask, 7][sur_mass_mask_g2]

    sur_mass_mask_gX = sur_mergers[sur_merger_gX_mask, 6] > sur_mergers[sur_merger_gX_mask, 7]
    sur_genX_mass_1 = np.zeros(np.sum(sur_merger_gX_mask))
    sur_genX_mass_1[sur_mass_mask_gX] = sur_mergers[sur_merger_gX_mask, 6][sur_mass_mask_gX]
    sur_genX_mass_1[~sur_mass_mask_gX] = sur_mergers[sur_merger_gX_mask, 7][~sur_mass_mask_gX]
    sur_genX_mass_2 = np.zeros(np.sum(sur_merger_gX_mask))
    sur_genX_mass_2[~sur_mass_mask_gX] = sur_mergers[sur_merger_gX_mask, 6][~sur_mass_mask_gX]
    sur_genX_mass_2[sur_mass_mask_gX] = sur_mergers[sur_merger_gX_mask, 7][sur_mass_mask_gX]

    # Check that there aren't any zeros remaining.
    assert (sur_gen1_mass_1 > 0).all()
    assert (sur_gen1_mass_2 > 0).all()
    assert (sur_gen2_mass_1 > 0).all()
    assert (sur_gen2_mass_2 > 0).all()
    assert (sur_genX_mass_1 > 0).all()
    assert (sur_genX_mass_2 > 0).all()

    pointsize_m1m2 = 5
    fig, ax = plt.subplots(1, 2, figsize=(5, 3), sharey=True)
    # ax4 = fig.add_subplot(111)

    sur = ax[1]
    nosur = ax[0]

    # plt.scatter(m1, m2, s=pointsize_m1m2, color='k')
    sur.scatter(sur_gen1_mass_1, sur_gen1_mass_2,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(sur_gen2_mass_1, sur_gen2_mass_2,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(sur_genX_mass_1, sur_genX_mass_2,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolor='none',
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set(
        xlabel=r'$M_1^{sur}$ [$M_\odot$]',
        xscale='log',
        yscale='log',
        axisbelow=(True),
        # aspect=('equal')
    )

    # sur.legend(fontsize=5)

    # plt.grid(True, color='gray', ls='dotted')
    # plt.savefig(opts.plots_directory + '/m1m2.png', format='png')
    # plt.show()

    # ========================================
    # NOSUR - Mass 1 vs Mass 2
    # ========================================

    # Sort Objects into Mass 1 and Mass 2 by generation
    nosur_mass_mask_g1 = nosur_mergers[nosur_merger_g1_mask, 6] > nosur_mergers[nosur_merger_g1_mask, 7]
    nosur_gen1_mass_1 = np.zeros(np.sum(nosur_merger_g1_mask))
    nosur_gen1_mass_1[nosur_mass_mask_g1] = nosur_mergers[nosur_merger_g1_mask, 6][nosur_mass_mask_g1]
    nosur_gen1_mass_1[~nosur_mass_mask_g1] = nosur_mergers[nosur_merger_g1_mask, 7][~nosur_mass_mask_g1]
    nosur_gen1_mass_2 = np.zeros(np.sum(nosur_merger_g1_mask))
    nosur_gen1_mass_2[~nosur_mass_mask_g1] = nosur_mergers[nosur_merger_g1_mask, 6][~nosur_mass_mask_g1]
    nosur_gen1_mass_2[nosur_mass_mask_g1] = nosur_mergers[nosur_merger_g1_mask, 7][nosur_mass_mask_g1]

    nosur_mass_mask_g2 = nosur_mergers[nosur_merger_g2_mask, 6] > nosur_mergers[nosur_merger_g2_mask, 7]
    nosur_gen2_mass_1 = np.zeros(np.sum(nosur_merger_g2_mask))
    nosur_gen2_mass_1[nosur_mass_mask_g2] = nosur_mergers[nosur_merger_g2_mask, 6][nosur_mass_mask_g2]
    nosur_gen2_mass_1[~nosur_mass_mask_g2] = nosur_mergers[nosur_merger_g2_mask, 7][~nosur_mass_mask_g2]
    nosur_gen2_mass_2 = np.zeros(np.sum(nosur_merger_g2_mask))
    nosur_gen2_mass_2[~nosur_mass_mask_g2] = nosur_mergers[nosur_merger_g2_mask, 6][~nosur_mass_mask_g2]
    nosur_gen2_mass_2[nosur_mass_mask_g2] = nosur_mergers[nosur_merger_g2_mask, 7][nosur_mass_mask_g2]

    nosur_mass_mask_gX = nosur_mergers[nosur_merger_gX_mask, 6] > nosur_mergers[nosur_merger_gX_mask, 7]
    nosur_genX_mass_1 = np.zeros(np.sum(nosur_merger_gX_mask))
    nosur_genX_mass_1[nosur_mass_mask_gX] = nosur_mergers[nosur_merger_gX_mask, 6][nosur_mass_mask_gX]
    nosur_genX_mass_1[~nosur_mass_mask_gX] = nosur_mergers[nosur_merger_gX_mask, 7][~nosur_mass_mask_gX]
    nosur_genX_mass_2 = np.zeros(np.sum(nosur_merger_gX_mask))
    nosur_genX_mass_2[~nosur_mass_mask_gX] = nosur_mergers[nosur_merger_gX_mask, 6][~nosur_mass_mask_gX]
    nosur_genX_mass_2[nosur_mass_mask_gX] = nosur_mergers[nosur_merger_gX_mask, 7][nosur_mass_mask_gX]

    # Check that there aren't any zeros remaining.
    assert (nosur_gen1_mass_1 > 0).all()
    assert (nosur_gen1_mass_2 > 0).all()
    assert (nosur_gen2_mass_1 > 0).all()
    assert (nosur_gen2_mass_2 > 0).all()
    assert (nosur_genX_mass_1 > 0).all()
    assert (nosur_genX_mass_2 > 0).all()

    pointsize_m1m2 = 5
    # fig = plt.figure(figsize=plotting.set_size(figsize))
    # ax4 = fig.add_subplot(111)

    # plt.scatter(m1, m2, s=pointsize_m1m2, color='k')
    nosur.scatter(nosur_gen1_mass_1, nosur_gen1_mass_2,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolor='none',
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    # plot the 2g+ mergers
    nosur.scatter(nosur_gen2_mass_1, nosur_gen2_mass_2,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolor='none',
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    # plot the 3g+ mergers
    nosur.scatter(nosur_genX_mass_1, nosur_genX_mass_2,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolor='none',
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    nosur.set(
        xlabel=r'$M_1^{nosur}$ [$M_\odot$]',
        ylabel=r'$M_2$ [$M_\odot$]',
        xscale='log',
        yscale='log',
        axisbelow=(True),
        # aspect=('equal')
    )
    # ylabel=r'$M_2$ [$M_\odot$]',

    nosur.legend(fontsize=5)

    # plt.grid(True, color='gray', ls='dotted')
    plt.savefig(opts.plots_directory + '/m1m2.png', format='png')
    # plt.show()

    # ========================================
    # NOSUR - Mass 1 vs Mass 2                       | STAND ALONE PLOTS |
    # ========================================

    # Sort Objects into Mass 1 and Mass 2 by generation
    nosur_mass_mask_g1 = nosur_mergers[nosur_merger_g1_mask, 6] > nosur_mergers[nosur_merger_g1_mask, 7]
    nosur_gen1_mass_1 = np.zeros(np.sum(nosur_merger_g1_mask))
    nosur_gen1_mass_1[nosur_mass_mask_g1] = nosur_mergers[nosur_merger_g1_mask, 6][nosur_mass_mask_g1]
    nosur_gen1_mass_1[~nosur_mass_mask_g1] = nosur_mergers[nosur_merger_g1_mask, 7][~nosur_mass_mask_g1]
    nosur_gen1_mass_2 = np.zeros(np.sum(nosur_merger_g1_mask))
    nosur_gen1_mass_2[~nosur_mass_mask_g1] = nosur_mergers[nosur_merger_g1_mask, 6][~nosur_mass_mask_g1]
    nosur_gen1_mass_2[nosur_mass_mask_g1] = nosur_mergers[nosur_merger_g1_mask, 7][nosur_mass_mask_g1]

    nosur_mass_mask_g2 = nosur_mergers[nosur_merger_g2_mask, 6] > nosur_mergers[nosur_merger_g2_mask, 7]
    nosur_gen2_mass_1 = np.zeros(np.sum(nosur_merger_g2_mask))
    nosur_gen2_mass_1[nosur_mass_mask_g2] = nosur_mergers[nosur_merger_g2_mask, 6][nosur_mass_mask_g2]
    nosur_gen2_mass_1[~nosur_mass_mask_g2] = nosur_mergers[nosur_merger_g2_mask, 7][~nosur_mass_mask_g2]
    nosur_gen2_mass_2 = np.zeros(np.sum(nosur_merger_g2_mask))
    nosur_gen2_mass_2[~nosur_mass_mask_g2] = nosur_mergers[nosur_merger_g2_mask, 6][~nosur_mass_mask_g2]
    nosur_gen2_mass_2[nosur_mass_mask_g2] = nosur_mergers[nosur_merger_g2_mask, 7][nosur_mass_mask_g2]

    nosur_mass_mask_gX = nosur_mergers[nosur_merger_gX_mask, 6] > nosur_mergers[nosur_merger_gX_mask, 7]
    nosur_genX_mass_1 = np.zeros(np.sum(nosur_merger_gX_mask))
    nosur_genX_mass_1[nosur_mass_mask_gX] = nosur_mergers[nosur_merger_gX_mask, 6][nosur_mass_mask_gX]
    nosur_genX_mass_1[~nosur_mass_mask_gX] = nosur_mergers[nosur_merger_gX_mask, 7][~nosur_mass_mask_gX]
    nosur_genX_mass_2 = np.zeros(np.sum(nosur_merger_gX_mask))
    nosur_genX_mass_2[~nosur_mass_mask_gX] = nosur_mergers[nosur_merger_gX_mask, 6][~nosur_mass_mask_gX]
    nosur_genX_mass_2[nosur_mass_mask_gX] = nosur_mergers[nosur_merger_gX_mask, 7][nosur_mass_mask_gX]

    # Check that there aren't any zeros remaining.
    assert (nosur_gen1_mass_1 > 0).all()
    assert (nosur_gen1_mass_2 > 0).all()
    assert (nosur_gen2_mass_1 > 0).all()
    assert (nosur_gen2_mass_2 > 0).all()
    assert (nosur_genX_mass_1 > 0).all()
    assert (nosur_genX_mass_2 > 0).all()

    pointsize_m1m2 = 5
    # fig = plt.figure(figsize=plotting.set_size(figsize))
    fig, ax = plt.subplots(1, 1, figsize=plotting.set_size(figsize), sharey=True)
    # ax4 = fig.add_subplot(111)

    # plt.scatter(m1, m2, s=pointsize_m1m2, color='k')
    ax.scatter(nosur_gen1_mass_1, nosur_gen1_mass_2,
               s=styles.markersize_gen1,
               marker=styles.marker_gen1,
               edgecolor=styles.color_gen1,
               facecolor='none',
               alpha=styles.markeralpha_gen1,
               label='1g-1g'
               )

    # plot the 2g+ mergers
    ax.scatter(nosur_gen2_mass_1, nosur_gen2_mass_2,
               s=styles.markersize_gen2,
               marker=styles.marker_gen2,
               edgecolor=styles.color_gen2,
               facecolor='none',
               alpha=styles.markeralpha_gen2,
               label='2g-1g or 2g-2g'
               )

    # plot the 3g+ mergers
    ax.scatter(nosur_genX_mass_1, nosur_genX_mass_2,
               s=styles.markersize_genX,
               marker=styles.marker_genX,
               edgecolor=styles.color_genX,
               facecolor='none',
               alpha=styles.markeralpha_genX,
               label=r'$\geq$3g-Ng'
               )

    ax.set(
        xlabel=r'$M_1$ [$M_\odot$] - (nosur)',
        ylabel=r'$M_2$ [$M_\odot$]',
        xscale='log',
        yscale='log',
        axisbelow=(True),
        # aspect=('equal')
    )
    # ylabel=r'$M_2$ [$M_\odot$]',

    ax.legend(fontsize=5)

    # plt.grid(True, color='gray', ls='dotted')
    # plt.savefig(opts.plots_directory + '/m1m2_nosur.png', format='png')
    # plt.show()

    # ===============================
    ### SUR - kick velocity histogram ###
    # ===============================
    fig, ax = plt.subplots(1, 2, figsize=plotting.set_size(figsize), constrained_layout=True, sharey=True)

    sur = ax[1]
    nosur = ax[0]

    # make your bins...
    kick_bins = np.logspace(np.log10(sur_mergers[:, 16].min()), np.log10(sur_mergers[:, 16].max() + 10), 50)

    sur_hist_data = [sur_mergers[:, 16][sur_merger_g1_mask], sur_mergers[:, 16][sur_merger_g2_mask],
                     sur_mergers[:, 16][sur_merger_gX_mask]]
    hist_label = ['1g-1g', '2g-1g or 2g-2g', r'$\geq$3g-Ng']
    hist_color = [styles.color_gen1, styles.color_gen2, styles.color_genX]

    # plot the distribution of mergers as a function of generation
    sur.hist(sur_hist_data, bins=kick_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label,
             stacked=True)
    # sur.set_ylabel(r'n')
    sur.set_xlabel(r'v$_{kick}^{sur}$ [km/s]')
    sur.set_xscale('log')

    # if figsize == 'apj_col':
    #    sur.legend(fontsize=5)
    # elif figsize == 'apj_page':
    #    sur.legend()

    # plt.title(r"Distribution of v$_{kick}$")
    # sur.grid(True, color='gray', ls='dashed')
    # plt.savefig(opts.plots_directory + "/v_kick_distribution.png", format='png')
    # plt.close()

    # ===============================
    ### NOSUR - kick velocity histogram ###
    # ===============================
    # fig = plt.figure(figsize=plotting.set_size(figsize))

    # make your bins...
    kick_bins = np.logspace(np.log10(nosur_mergers[:, 16].min()), np.log10(nosur_mergers[:, 16].max() + 10), 50)

    nosur_hist_data = [nosur_mergers[:, 16][nosur_merger_g1_mask], nosur_mergers[:, 16][nosur_merger_g2_mask],
                       nosur_mergers[:, 16][nosur_merger_gX_mask]]
    hist_label = ['1g-1g', '2g-1g or 2g-2g', r'$\geq$3g-Ng']
    hist_color = [styles.color_gen1, styles.color_gen2, styles.color_genX]

    # plot the distribution of mergers as a function of generation
    nosur.hist(nosur_hist_data, bins=kick_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label,
               stacked=True)
    nosur.set_ylabel(r'n')
    nosur.set_xlabel(r'v$_{kick}^{nosur}$ [km/s]')
    nosur.set_xscale('log')

    if figsize == 'apj_col':
        nosur.legend(fontsize=4)
    elif figsize == 'apj_page':
        nosur.legend()

    # plt.title(r"Distribution of v$_{kick}$")
    # plt.grid(True, color='gray', ls='dashed')
    plt.savefig(opts.plots_directory + "/vkick_distribution.png", format='png')
    # plt.show()

    # ===============================
    ### SUR - Kick Velocity vs Radius ###
    # ===============================

    sur_all_kick = sur_mergers[:, 16]
    sur_gen1_vkick = sur_all_kick[sur_merger_g1_mask]
    sur_gen2_vkick = sur_all_kick[sur_merger_g2_mask]
    sur_genX_vkick = sur_all_kick[sur_merger_gX_mask]

    # figsize is hardcoded here. don't change, shrink everything illegibly
    fig, axs = plt.subplots(nrows=2, ncols=2, sharey=True, figsize=(4, 3),
                            gridspec_kw={'width_ratios': [3, 1], 'wspace': 0, 'hspace': 0}, constrained_layout=True)

    sur_scatter = axs[1][0]
    sur_hist = axs[1][1]
    nosur_scatter = axs[0][0]
    nosur_hist = axs[0][1]

    # plot 1g-1g mergers
    sur_scatter.scatter(sur_gen1_orb_a, sur_gen1_vkick,
                        s=styles.markersize_gen1,
                        marker=styles.marker_gen1,
                        edgecolor=styles.color_gen1,
                        facecolors="none",
                        alpha=styles.markeralpha_gen1,
                        label='1g-1g'
                        )

    # plot 2g-mg mergers
    sur_scatter.scatter(sur_gen2_orb_a, sur_gen2_vkick,
                        s=styles.markersize_gen2,
                        marker=styles.marker_gen2,
                        edgecolor=styles.color_gen2,
                        facecolors="none",
                        alpha=styles.markeralpha_gen2,
                        label='2g-1g or 2g-2g'
                        )

    # plot 3g-ng mergers
    sur_scatter.scatter(sur_genX_orb_a, sur_genX_vkick,
                        s=styles.markersize_genX,
                        marker=styles.marker_genX,
                        edgecolor=styles.color_genX,
                        facecolors="none",
                        alpha=styles.markeralpha_genX,
                        label=r'$\geq$3g-Ng'
                        )

    # plot trap radius
    trap_radius = 700
    sur_scatter.axvline(trap_radius, color='k', linestyle='--', zorder=0,
                        label=f'Trap Radius = {trap_radius} ' + r'$R_g$')

    # plotting escape velocity [km/s]
    max_radius = np.max(sur_gen1_orb_a) * (const.G.value * 1.e8 * const.M_sun.value / const.c.value ** 2)
    r_kms = np.array(np.linspace(3e3, max_radius, 100))
    r_rg = np.array(np.linspace(5e2, np.max(sur_gen1_orb_a), 100))

    # SMBH = 1e8
    v_esc8 = np.sqrt((2 * const.G.value * 1.e8 * const.M_sun.value) / r_kms) / 1e3
    sur_scatter.plot(r_rg[-99:], v_esc8[-99:], label='$v_{esc}$, SMBH = 1e8', color='teal')
    sur_scatter.set_xlim(3e2, 7e4)
    # SMBH = 1e7
    # v_esc7 = np.sqrt((2 * const.G.value * 1.e7 * const.M_sun.value)/ r_kms) / 1e3
    # sur1.plot(r_rg[-99:], v_esc7[-99:], label='$v_{esc}$, SMBH = 1e7', color='green')
    # sur1.set_xlim(3e2, 7e4)
    # SMBH = 1e6
    # v_esc6 = np.sqrt((2 * const.G.value * 1.e6 * const.M_sun.value)/ r_kms) / 1e3
    # sur1.plot(r_rg[-99:], v_esc6[-99:], label='$v_{esc}$, SMBH = 1e6', color='orange')
    # sur1.set_xlim(3e2, 7e4)
    # SMBH = 1e5
    # v_esc5 = np.sqrt((2 * const.G.value * 1.e5 * const.M_sun.value)/ r_kms) / 1e3
    # sur1.plot(r_rg[-99:], v_esc5[-99:], label='$v_{esc}$, SMBH = 1e5', color='red')
    # sur1.set_xlim(3e2, 7e4)

    # plotting keplarian velocity [km/s]
    # v_kep = np.sqrt((const.G.value * 1.e8 * const.M_sun.value)/ r_kms) / 1e3
    # ur1.plot(r_rg[-99:], v_kep[-99:], label='SMBH Keplarian Velocity', color='red')

    # configure scatter plot
    sur_scatter.set(
        ylabel=(r'$v_{kick}^{sur}$ [km/s]'),
        xlabel=(r'Radius [$R_g$]'),
        xscale=('log'),
        yscale=('log'),
        # xlim=(3e2, 7e4)
    )

    sur_scatter.grid(True, color='gray', ls='dashed')
    if figsize == 'apj_col':
        sur_scatter.legend(fontsize=3, loc='lower right')
    elif figsize == 'apj_page':
        sur_scatter.legend()

    # calculate mean kick velocity for all mergers
    sur_mean_kick = np.mean(sur_mergers[:, 16])

    kick_bins = np.logspace(np.log10(sur_mergers[:, 16].min()), np.log10(sur_mergers[:, 16].max() + 10), 50)

    # configure histogram
    sur_hist.grid(True, color='gray', ls='dashed')
    # sur_hist_data = [sur_mergers[:, 4][sur_merger_g1_mask], sur_mergers[:, 4][sur_merger_g2_mask], sur_mergers[:, 4][sur_merger_gX_mask]]
    sur_hist.hist(sur_hist_data, bins=kick_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                  label=hist_label, stacked=True, orientation='horizontal')
    sur_hist.axhline(sur_mean_kick, color='black', linewidth=1, linestyle='dashdot',
                     label=r'$\langle v_{kick}\rangle $ =' + f"{sur_mean_kick:.2f}")
    sur_hist.yaxis.tick_right()

    sur_hist_data_int = list(
        map(int, sur_hist_data[0] / 25))  # dividing by 25 to offset the data points for the kick counts
    mode, spin_count_sur = stats.mode(sur_hist_data_int, axis=None, keepdims=False)
    sur_hist.set(
        xlabel=r'n',
        xlim=[0, int(spin_count_sur * 1.30)],
        # xticks=[300, 1000],
        xticks=np.linspace(int(spin_count_sur * 0.30), int(spin_count_sur * 0.90), 2)
    )

    plt.setp(sur_scatter.get_yticklabels(), visible=True)
    plt.setp(sur_hist.get_yticklabels(), visible=False)

    if figsize == 'apj_col':
        sur_hist.legend(fontsize=3, loc='lower right')
    elif figsize == 'apj_page':
        sur_hist.legend()

    # plt.title(r"v$_{kick} vs. semi-major axis with distribution of v$_{kick}$")
    # plt.tight_layout()
    # plt.savefig(opts.plots_directory + '/v_kick_vs_radius.png', format='png')
    # plt.close()

    # ===============================
    ### NOSUR - Kick Velocity vs Radius ###
    # ===============================

    nosur_all_kick = nosur_mergers[:, 16]
    nosur_gen1_vkick = nosur_all_kick[nosur_merger_g1_mask]
    nosur_gen2_vkick = nosur_all_kick[nosur_merger_g2_mask]
    nosur_genX_vkick = nosur_all_kick[nosur_merger_gX_mask]

    # figsize is hardcoded here. don't change, shrink everything illegibly
    # fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(5.5,3), gridspec_kw={'width_ratios': [3, 1], 'wspace':0, 'hspace':0})

    # plot 1g-1g mergers
    nosur_scatter.scatter(nosur_gen1_orb_a, nosur_gen1_vkick,
                          s=styles.markersize_gen1,
                          marker=styles.marker_gen1,
                          edgecolor=styles.color_gen1,
                          facecolors="none",
                          alpha=styles.markeralpha_gen1,
                          label='1g-1g'
                          )

    # plot 2g-mg mergers
    nosur_scatter.scatter(nosur_gen2_orb_a, nosur_gen2_vkick,
                          s=styles.markersize_gen2,
                          marker=styles.marker_gen2,
                          edgecolor=styles.color_gen2,
                          facecolors="none",
                          alpha=styles.markeralpha_gen2,
                          label='2g-1g or 2g-2g'
                          )

    # plot 3g-ng mergers
    nosur_scatter.scatter(nosur_genX_orb_a, nosur_genX_vkick,
                          s=styles.markersize_genX,
                          marker=styles.marker_genX,
                          edgecolor=styles.color_genX,
                          facecolors="none",
                          alpha=styles.markeralpha_genX,
                          label=r'$\geq$3g-Ng'
                          )

    # plot trap radius
    trap_radius = 700
    nosur_scatter.axvline(trap_radius, color='k', linestyle='--', zorder=0,
                          label=f'Trap Radius = {trap_radius} ' + r'$R_g$')

    # configure scatter plot
    nosur_scatter.set(
        ylabel=(r'$v_{kick}^{nosur}$ [km/s]'),
        xlabel=(r'Radius [$R_g$]'),
        xscale=('log'),
        yscale=('log'),
        xlim=(3e2, 7e4),

    )
    nosur_scatter.grid(True, color='gray', ls='dashed')

    # SMBH = 1e8
    nosur_scatter.plot(r_rg[-99:], v_esc8[-99:], label='$v_{esc}$, SMBH = 1e8', color='teal')
    # SMBH = 1e7
    # nosur1.plot(r_rg[-99:], v_esc7[-99:], label='$v_{esc}$, SMBH = 1e7', color='green')
    # SMBH = 1e6
    # nosur1.plot(r_rg[-99:], v_esc6[-99:], label='$v_{esc}$, SMBH = 1e6', color='orange')
    # MBH = 1e5
    # nosur1.plot(r_rg[-99:], v_esc5[-99:], label='$v_{esc}$, SMBH = 1e5', color='red')

    # if figsize == 'apj_col':
    #    nosur1.legend(fontsize=4, loc='upper left')
    # elif figsize == 'apj_page':
    #    nosur1.legend()

    # nosur1.set_xlim(5e2, np.max(sur_gen1_orb_a))

    # if figsize == 'apj_col':
    #    nosur1.legend(fontsize=5, loc='upper left')
    # elif figsize == 'apj_page':
    #    nosur1.legend()

    # calculate mean kick velocity for all mergers
    nosur_mean_kick = np.mean(nosur_mergers[:, 16])

    kick_bins = np.logspace(np.log10(nosur_mergers[:, 16].min()), np.log10(nosur_mergers[:, 16].max() + 10), 50)

    # configure histogram
    nosur_hist.grid(True, color='gray', ls='dashed')
    # sur_hist_data = [sur_mergers[:, 4][sur_merger_g1_mask], sur_mergers[:, 4][sur_merger_g2_mask], sur_mergers[:, 4][sur_merger_gX_mask]]
    nosur_hist.hist(nosur_hist_data, bins=kick_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                    label=hist_label, stacked=True, orientation='horizontal')
    nosur_hist.axhline(nosur_mean_kick, color='black', linewidth=1, linestyle='dashdot',
                       label=r'$\langle v_{kick}\rangle $ =' + f"{nosur_mean_kick:.2f}")
    nosur_hist.yaxis.tick_right()

    # setting graph tick marks to line up with paired graph
    nosur_hist_data_int = list(
        map(int, sur_hist_data[0] / 25))  # dividing by 25 to offset the data points for the kick counts
    mode, spin_count_sur = stats.mode(sur_hist_data_int, axis=None, keepdims=False)
    nosur_hist.set(
        xlabel=r'n',
        xlim=[0, int(spin_count_sur * 1.30)],
        # xticks=[300, 1000],
        xticks=np.linspace(int(spin_count_sur * 0.30), int(spin_count_sur * 0.90), 2)
    )

    plt.setp(nosur_scatter.get_yticklabels(), visible=True)
    plt.setp(nosur_hist.get_yticklabels(), visible=False)

    if figsize == 'apj_col':
        nosur_hist.legend(fontsize=3, loc='lower right')
    elif figsize == 'apj_page':
        nosur_hist.legend()

    # plt.title(r"v$_{kick} vs. semi-major axis with distribution of v$_{kick}$")
    plt.tight_layout()
    plt.savefig(opts.plots_directory + '/vkick_vs_radius.png', format='png')
    # plt.savefig(opts.plots_directory + '/v_kick_vs_radius_nosur.png', format='png')
    # plt.show()

    # ===============================
    ### NOSUR - Kick Velocity vs Radius ### STAND ALONE PLOT
    # ===============================

    nosur_all_kick = nosur_mergers[:, 16]
    nosur_gen1_vkick = nosur_all_kick[nosur_merger_g1_mask]
    nosur_gen2_vkick = nosur_all_kick[nosur_merger_g2_mask]
    nosur_genX_vkick = nosur_all_kick[nosur_merger_gX_mask]

    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(4, 3),
                            gridspec_kw={'width_ratios': [3, 1], 'wspace': 0, 'hspace': 0})
    kvel1 = axs[0]
    kvel2 = axs[1]
    # figsize is hardcoded here. don't change, shrink everything illegibly
    # fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(5.5,3), gridspec_kw={'width_ratios': [3, 1], 'wspace':0, 'hspace':0})

    # plot 1g-1g mergers
    kvel1.scatter(nosur_gen1_orb_a, nosur_gen1_vkick,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolors="none",
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    # plot 2g-mg mergers
    kvel1.scatter(nosur_gen2_orb_a, nosur_gen2_vkick,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolors="none",
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    # plot 3g-ng mergers
    kvel1.scatter(nosur_genX_orb_a, nosur_genX_vkick,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolors="none",
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    # plot trap radius
    trap_radius = 700
    kvel1.axvline(trap_radius, color='k', linestyle='--', zorder=0,
                  label=f'Trap Radius = {trap_radius} ' + r'$R_g$')

    # configure scatter plot
    kvel1.set_ylabel(r'$v_{kick}^{nosur}$ [km/s]')
    kvel1.set_xlabel(r'Radius [$R_g$]')
    kvel1.set_xscale('log')
    kvel1.set_yscale('log')
    kvel1.set_xlim(3e2, 7e4)
    kvel1.grid(True, color='gray', ls='dashed')
    if figsize == 'apj_col':
        kvel1.legend(fontsize=5, loc='lower left')
    elif figsize == 'apj_page':
        kvel1.legend()

    # calculate mean kick velocity for all mergers
    nosur_mean_kick = np.mean(nosur_mergers[:, 16])

    kick_bins = np.logspace(np.log10(nosur_mergers[:, 16].min()), np.log10(nosur_mergers[:, 16].max() + 10), 50)

    # configure histogram
    kvel2.hist(nosur_hist_data, bins=kick_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label,
               stacked=True, orientation='horizontal')
    kvel2.axhline(nosur_mean_kick, color='black', linewidth=1, linestyle='dashdot',
                  label=r'$\langle v_{kick}\rangle $ =' + f"{nosur_mean_kick:.2f}")
    kvel2.grid(True, color='gray', ls='dashed')
    kvel2.set_yscale('log')
    kvel2.yaxis.tick_right()
    kvel2.set_xlim(0, 500)
    kvel2.set_xlabel(r'n')
    kvel2.set_xticks([100, 400])

    if figsize == 'apj_col':
        kvel2.legend(fontsize=4, loc='lower right')
    elif figsize == 'apj_page':
        kvel2.legend()

    # plt.title(r"v$_{kick} vs. semi-major axis with distribution of v$_{kick}$")
    plt.tight_layout()
    # plt.savefig(opts.plots_directory + '/v_kick_vs_radius_nosur.png', format='png')
    # plt.show()

    # ===============================
    ### SUR - kick velocity histogram across disk radius ###
    # ===============================
    fig, ax = plt.subplots(1, 2, figsize=plotting.set_size(figsize), sharey=True)

    sur = ax[1]
    nosur = ax[0]

    # make your bins...
    sur_radius_bins = np.logspace(np.log10(sur_mergers[:, 1].min()), np.log10(sur_mergers[:, 1].max() + 10), 50)

    sur_hist_data = [sur_mergers[:, 1][sur_merger_g1_mask], sur_mergers[:, 1][sur_merger_g2_mask],
                     sur_mergers[:, 1][sur_merger_gX_mask]]
    hist_label = ['1g-1g', '2g-1g or 2g-2g', r'$\geq$3g-Ng']
    hist_color = [styles.color_gen1, styles.color_gen2, styles.color_genX]

    # plot the distribution of mergers as a function of generation
    sur.hist(sur_hist_data, bins=sur_radius_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
             label=hist_label, stacked=True)
    sur.set_xlabel(r'Radius$_{sur}$ [$R_g$]')
    # sur.set_ylabel(r'No. of mergers')
    sur.set_xscale('log')

    if figsize == 'apj_col':
        sur.legend(fontsize=6)
    elif figsize == 'apj_page':
        sur.legend()

    # plt.title(r"Distribution of v$_{kick}$")
    # sur.grid(True, color='gray', ls='dashed')
    # plt.savefig(opts.plots_directory + "/v_kick_distribution.png", format='png')
    # plt.close()

    # ===============================
    ### NOSUR - kick velocity histogram across disk radius ###
    # ===============================
    # fig = plt.figure(figsize=plotting.set_size(figsize))

    # make your bins...
    nosur_radius_bins = np.logspace(np.log10(nosur_mergers[:, 1].min()), np.log10(nosur_mergers[:, 1].max() + 10), 50)

    nosur_hist_data = [nosur_mergers[:, 1][nosur_merger_g1_mask], nosur_mergers[:, 1][nosur_merger_g2_mask],
                       nosur_mergers[:, 1][nosur_merger_gX_mask]]
    hist_label = ['1g-1g', '2g-1g or 2g-2g', r'$\geq$3g-Ng']
    hist_color = [styles.color_gen1, styles.color_gen2, styles.color_genX]

    # plot the distribution of mergers as a function of generation
    nosur.hist(nosur_hist_data, bins=nosur_radius_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
               label=hist_label, stacked=True)
    nosur.set_xlabel(r'Radius$_{nosur}$ [$R_g$]')
    nosur.set_ylabel(r'No. of mergers')
    nosur.set_xscale('log')

    if figsize == 'apj_col':
        nosur.legend(fontsize=6)
    elif figsize == 'apj_page':
        nosur.legend()

    # plt.title(r"Distribution of v$_{kick}$")
    # plt.grid(True, color='gray', ls='dashed')
    plt.savefig(opts.plots_directory + "/vkick_dist_radius.png", format='png')
    # plt.show()

    # ========================================
    # Spin vs Kick Velocity
    # ========================================

    plot = plt.figure(figsize=(12, 4), constrained_layout=False)

    # ======= NO SURROGATE =========
    nosur = plot.add_gridspec(nrows=1, ncols=4, left=0.05, right=0.2625, wspace=0)
    nosur_scatter = plot.add_subplot(nosur[0, :-1])
    nosur_hist = plot.add_subplot(nosur[0, 3], sharey=nosur_scatter)

    nosur_spin = nosur_mergers[:, 4]
    nosur_gen1_spin = nosur_spin[nosur_merger_g1_mask]
    nosur_gen2_spin = nosur_spin[nosur_merger_g2_mask]
    nosur_genX_spin = nosur_spin[nosur_merger_gX_mask]

    nosur_all_kick = nosur_mergers[:, 16]
    nosur_gen1_vkick = nosur_all_kick[nosur_merger_g1_mask]
    nosur_gen2_vkick = nosur_all_kick[nosur_merger_g2_mask]
    nosur_genX_vkick = nosur_all_kick[nosur_merger_gX_mask]

    # some random data
    x = nosur_all_kick
    y = nosur_spin

    def scatter_hist(x, y, ax, ax_histx, ax_histy):
        # no labels
        # ax_histx.tick_params(axis="x", labelbottom=False)
        # ax_histy.tick_params(axis="y", labelleft=False)

        # the scatter plot:
        ax.scatter(x, y)

        # now determine nice limits by hand:
        spin_bins = np.logspace(np.log10(nosur_mergers[:, 4].min()), np.log10(nosur_mergers[:, 4].max()), 50)
        kick_bins = np.logspace(np.log10(nosur_mergers[:, 16].min()), np.log10(nosur_mergers[:, 16].max()), 50)

        # bins = np.arange(-lim, lim + binwidth, binwidth)
        nosur_spin_hist_data = [nosur_mergers[:, 4][nosur_merger_g1_mask], nosur_mergers[:, 4][nosur_merger_g2_mask],
                                nosur_mergers[:, 4][nosur_merger_gX_mask]]
        nosur_kick_hist_data = [nosur_mergers[:, 16][nosur_merger_g1_mask], nosur_mergers[:, 16][nosur_merger_g2_mask],
                                nosur_mergers[:, 16][nosur_merger_gX_mask]]

        # setting grid lines to the plots
        ax.grid(True, color='gray', ls='dashed', alpha=0.4)
        ax_histx.grid(True, color='gray', ls='dashed', alpha=0.4)
        ax_histy.grid(True, color='gray', ls='dashed', alpha=0.4)

        ax_histx.hist(nosur_kick_hist_data, bins=kick_bins, color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label,
                      stacked=False)
        ax_histy.hist(nosur_spin_hist_data, bins=spin_bins, color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label,
                      stacked=True, orientation='horizontal')

        ax.set(
            xlabel=r'$v_{kick}$ [km/s]',
            ylabel=r'$a_{final}$',
            xscale="log",
            axisbelow=True,
            xlim=([1.1e0, 2e3]),
            ylim=(0.21, 1.01)
            # , title=r'$a_{final}^{nosur}$',
        )
        ax_histx.set(
            xscale='log',
            ylabel=r'n'
        )

    fig, axs = plt.subplot_mosaic([['histx', '.'],
                                   ['scatter', 'histy']],
                                  figsize=(6, 6),
                                  width_ratios=(4, 1), height_ratios=(1, 4),
                                  layout='constrained')
    scatter_hist(x, y, axs['scatter'], axs['histx'], axs['histy'])
    plt.savefig(opts.plots_directory + "/vkick_test_fig.png", format='png')

    nosur_scatter.patch.set_linewidth(3)
    nosur_hist.patch.set_linewidth(3)

    # color for presentations
    # nosur_scatter.patch.set_edgecolor('#D81B60')
    # nosur_hist.patch.set_edgecolor('#D81B60')

    # plot the 1g-1g mergers
    nosur_scatter.scatter(nosur_gen1_vkick, nosur_gen1_spin,
                          s=styles.markersize_gen1,
                          marker=styles.marker_gen1,
                          edgecolor=styles.color_gen1,
                          facecolor='none',
                          alpha=styles.markeralpha_gen1,
                          label='1g-1g'
                          )

    # plot the 2g+ mergers
    nosur_scatter.scatter(nosur_gen2_vkick, nosur_gen2_spin,
                          s=styles.markersize_gen2,
                          marker=styles.marker_gen2,
                          edgecolor=styles.color_gen2,
                          facecolor='none',
                          alpha=styles.markeralpha_gen2,
                          label='2g-1g or 2g-2g'
                          )

    # plot the 3g+ mergers
    nosur_scatter.scatter(nosur_genX_vkick, nosur_genX_spin,
                          s=styles.markersize_genX,
                          marker=styles.marker_genX,
                          edgecolor=styles.color_genX,
                          facecolor='none',
                          alpha=styles.markeralpha_genX,
                          label=r'$\geq$3g-Ng'
                          )

    nosur_scatter.grid(True, color='gray', ls='dashed', alpha=0.4)
    nosur_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        ylabel=r'$a_{final}$',
        title=r'$a_{final}^{nosur}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    spin_bins = np.logspace(np.log10(sur_mergers[:, 4].min()), np.log10(sur_mergers[:, 4].max()), 50)

    nosur_hist.grid(True, color='gray', ls='dashed', alpha=0.4)
    nosur_hist_data = [nosur_mergers[:, 4][nosur_merger_g1_mask], nosur_mergers[:, 4][nosur_merger_g2_mask],
                       nosur_mergers[:, 4][nosur_merger_gX_mask]]
    nosur_hist_entries, nosur_hist_bin_edges, _ = nosur_hist.hist(nosur_hist_data, bins=spin_bins, align='left',
                                                               color=hist_color, alpha=0.9, rwidth=0.8,
                                                               label=hist_label, stacked=True, orientation='horizontal')
    nosur_hist_data_int = list(map(int, nosur_hist_data[0] * 100))
    mode, spin_count_nosur = stats.mode(nosur_hist_data_int, axis=None, keepdims=False)

    # get poisson deviated random numbers
    # data = np.random.poisson(2, 1000)
    # merger generations are seperated to fit each poisson
    nosur_hist_data_gen1 = nosur_mergers[:, 4][nosur_merger_g1_mask]
    nosur_hist_data_gen2 = nosur_mergers[:, 4][nosur_merger_g2_mask]
    nosur_hist_data_genX = nosur_mergers[:, 4][nosur_merger_gX_mask]

    # the bins should be of integer width, because poisson is an integer distribution
    # bins = np.arange(11) - 0.5
    # entries, bin_edges, patches = plt.hist(data, bins=bins, density=True, label='Data')
    # nosur_hist_bins_gen1 = np.logspace(np.log10(nosur_mergers[:, 4][nosur_merger_g1_mask].min()),
    #                                    np.log10(nosur_mergers[:, 4][nosur_merger_g1_mask].max()), 50)
    # nosur_hist_bins_gen2 = np.logspace(np.log10(nosur_mergers[:, 4][nosur_merger_g2_mask].min()),
    #                                    np.log10(nosur_mergers[:, 4][nosur_merger_g2_mask].max()), 50)
    # nosur_hist_bins_genX = np.logspace(np.log10(nosur_mergers[:, 4][nosur_merger_gX_mask].min()),
    #                                    np.log10(nosur_mergers[:, 4][nosur_merger_gX_mask].max()), 50)
    #
    # nosur_hist.grid(True, color='gray', ls='dashed', alpha=0.4)
    # nosur_hist_data = [nosur_mergers[:, 4][nosur_merger_g1_mask], nosur_mergers[:, 4][nosur_merger_g2_mask],
    #                    nosur_mergers[:, 4][nosur_merger_gX_mask]]
    # nosur_hist_entries, nosur_hist_bin_edges, _ = nosur_hist.hist(nosur_hist_data, bins=spin_bins, align='left',
    #                                                            color=hist_color, alpha=0.9, rwidth=0.8,
    #                                                            label=hist_label, stacked=True, orientation='horizontal')
    # nosur_hist_data_int = list(map(int, nosur_hist_data[0] * 100))
    # mode, spin_count_nosur = stats.mode(nosur_hist_data_int, axis=None, keepdims=False)
    #
    # # calculate bin centers
    # bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    #
    # def fit_function(k, lamb):
    #     '''poisson function, parameter lamb is the fit parameter'''
    #     return poisson.pmf(k, lamb)
    #
    # # fit with curve_fit
    # parameters, cov_matrix = curve_fit(fit_function, bin_centers, entries)
    #
    # # plot poisson-deviation with fitted parameter
    # x_plot = np.arange(0, 15)
    #
    # plt.plot(
    #     x_plot,
    #     fit_function(x_plot, *parameters),
    #     marker='o', linestyle='',
    #     label='Fit result',
    # )
    # plt.legend()
    # plt.show()

    nosur_hist.yaxis.tick_right()
    nosur_hist.set(
        xlabel=r'n',
        xlim=[0, int(spin_count_nosur * 1.80)],
        # xticks=[300, 1000],
        xticks=np.linspace(int(spin_count_nosur * 0.50), int(spin_count_nosur * 1.30), 2)
    )
    plt.setp(nosur_hist.get_yticklabels(), visible=False)

    if figsize == 'apj_col':
        nosur_scatter.legend(fontsize=5, loc='lower left')
    elif figsize == 'apj_page':
        nosur_scatter.legend()

    # ======= NO SURROGATE W/FILTER =========
    nosur_filter = plot.add_gridspec(nrows=1, ncols=4, left=0.275, right=0.4875, wspace=0)
    nosur_filter_scatter = plot.add_subplot(nosur_filter[0, :-1])
    nosur_filter_hist = plot.add_subplot(nosur_filter[0, 3], sharey=nosur_filter_scatter)

    nosur_filter_spin = nosur_filter_mergers[:, 4]
    nosur_filter_gen1_spin = nosur_filter_spin[nosur_filter_merger_g1_mask]
    nosur_filter_gen2_spin = nosur_filter_spin[nosur_filter_merger_g2_mask]
    nosur_filter_genX_spin = nosur_filter_spin[nosur_filter_merger_gX_mask]

    nosur_filter_all_kick = nosur_filter_mergers[:, 16]
    nosur_filter_gen1_vkick = nosur_filter_all_kick[nosur_filter_merger_g1_mask]
    nosur_filter_gen2_vkick = nosur_filter_all_kick[nosur_filter_merger_g2_mask]
    nosur_filter_genX_vkick = nosur_filter_all_kick[nosur_filter_merger_gX_mask]

    nosur_filter_scatter.patch.set_linewidth(3)
    nosur_filter_hist.patch.set_linewidth(3)

    # color for presentations
    # nosur_filter_scatter.patch.set_edgecolor('#1E88E5')
    # nosur_filter_hist.patch.set_edgecolor('#1E88E5')

    # plot the 1g-1g mergers
    nosur_filter_scatter.scatter(nosur_filter_gen1_vkick, nosur_filter_gen1_spin,
                                 s=styles.markersize_gen1,
                                 marker=styles.marker_gen1,
                                 edgecolor=styles.color_gen1,
                                 facecolor='none',
                                 alpha=styles.markeralpha_gen1,
                                 label='1g-1g'
                                 )

    # plot the 2g+ mergers
    nosur_filter_scatter.scatter(nosur_filter_gen2_vkick, nosur_filter_gen2_spin,
                                 s=styles.markersize_gen2,
                                 marker=styles.marker_gen2,
                                 edgecolor=styles.color_gen2,
                                 facecolor='none',
                                 alpha=styles.markeralpha_gen2,
                                 label='2g-1g or 2g-2g'
                                 )

    # plot the 3g+ mergers
    nosur_filter_scatter.scatter(nosur_filter_genX_vkick, nosur_filter_genX_spin,
                                 s=styles.markersize_genX,
                                 marker=styles.marker_genX,
                                 edgecolor=styles.color_genX,
                                 facecolor='none',
                                 alpha=styles.markeralpha_genX,
                                 label=r'$\geq$3g-Ng'
                                 )

    nosur_filter_scatter.grid(True, color='gray', ls='dashed', alpha=0.4)
    nosur_filter_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        title=r'$a_{final}^{nosur\_filter}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )
    plt.setp(nosur_filter_scatter.get_yticklabels(), visible=False)

    spin_bins = np.logspace(np.log10(nosur_mergers[:, 4].min()), np.log10(nosur_mergers[:, 4].max()), 50)

    nosur_filter_hist.grid(True, color='gray', ls='dashed', alpha=0.4)
    nosur_filter_hist_data = [nosur_filter_mergers[:, 4][nosur_filter_merger_g1_mask],
                              nosur_filter_mergers[:, 4][nosur_filter_merger_g2_mask],
                              nosur_filter_mergers[:, 4][nosur_filter_merger_gX_mask]]
    nosur_filter_hist.hist(nosur_filter_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9,
                           rwidth=0.8, label=hist_label, stacked=True, orientation='horizontal')
    nosur_filter_hist_data_int = list(map(int, nosur_filter_hist_data[0] * 100))
    mode, spin_count_nosur_filter = stats.mode(nosur_filter_hist_data_int, axis=None, keepdims=False)

    nosur_filter_hist.set(
        xlabel=r'n',
        xlim=[0, int(spin_count_nosur_filter * 1.80)],
        # xticks=[300, 1000],
        xticks=np.linspace(int(spin_count_nosur * 0.50), int(spin_count_nosur * 1.30), 2)
    )
    plt.setp(nosur_filter_hist.get_yticklabels(), visible=False)

    # ======= PRECESSION =========
    prec = plot.add_gridspec(nrows=1, ncols=4, left=0.5, right=0.7125, wspace=0)
    prec_scatter = plot.add_subplot(prec[0, :-1])
    prec_hist = plot.add_subplot(prec[0, 3], sharey=prec_scatter)

    prec_spin = prec_mergers[:, 4]
    prec_gen1_spin = prec_spin[prec_merger_g1_mask]
    prec_gen2_spin = prec_spin[prec_merger_g2_mask]
    prec_genX_spin = prec_spin[prec_merger_gX_mask]

    prec_all_kick = prec_mergers[:, 16]
    prec_gen1_vkick = prec_all_kick[prec_merger_g1_mask]
    prec_gen2_vkick = prec_all_kick[prec_merger_g2_mask]
    prec_genX_vkick = prec_all_kick[prec_merger_gX_mask]

    prec_scatter.patch.set_linewidth(3)
    prec_hist.patch.set_linewidth(3)

    # color for presentations
    # prec_scatter.patch.set_edgecolor('#E4AC04')
    # prec_hist.patch.set_edgecolor("#E4AC04")

    prec_scatter.scatter(prec_gen1_vkick, prec_gen1_spin,
                         s=styles.markersize_gen1,
                         marker=styles.marker_gen1,
                         edgecolor=styles.color_gen1,
                         facecolor='none',
                         alpha=styles.markeralpha_gen1,
                         label='1g-1g'
                         )

    # plot the 2g+ mergers
    prec_scatter.scatter(prec_gen2_vkick, prec_gen2_spin,
                         s=styles.markersize_gen2,
                         marker=styles.marker_gen2,
                         edgecolor=styles.color_gen2,
                         facecolor='none',
                         alpha=styles.markeralpha_gen2,
                         label='2g-1g or 2g-2g'
                         )

    # plot the 3g+ mergers
    prec_scatter.scatter(prec_genX_vkick, prec_genX_spin,
                         s=styles.markersize_genX,
                         marker=styles.marker_genX,
                         edgecolor=styles.color_genX,
                         facecolor='none',
                         alpha=styles.markeralpha_genX,
                         label=r'$\geq$3g-Ng'
                         )

    prec_scatter.grid(True, color='gray', ls='dashed', alpha=0.4)
    prec_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        title=r'$a_{final}^{prec}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    spin_bins = np.logspace(np.log10(prec_mergers[:, 4].min()), np.log10(prec_mergers[:, 4].max()), 50)

    prec_hist.grid(True, color='gray', ls='dashed', alpha=0.4)
    prec_hist_data = [prec_mergers[:, 4][prec_merger_g1_mask], prec_mergers[:, 4][prec_merger_g2_mask],
                      prec_mergers[:, 4][prec_merger_gX_mask]]
    prec_hist.hist(prec_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                   label=hist_label, stacked=True, orientation='horizontal')
    prec_hist.yaxis.tick_right()
    prec_hist_data_int = list(map(int, prec_hist_data[0] * 100))
    mode, spin_count_prec = stats.mode(prec_hist_data_int, axis=None, keepdims=False)

    prec_hist.set(
        xlabel=r'n',
        xlim=[0, int(spin_count_prec * 1.80)],
        # xticks=[300, 1000],
        xticks=np.linspace(int(spin_count_nosur * 0.50), int(spin_count_nosur * 1.30), 2)
    )
    plt.setp(prec_scatter.get_yticklabels(), visible=False)
    plt.setp(prec_hist.get_yticklabels(), visible=False)
    # ======= SURROGATE =========
    sur = plot.add_gridspec(nrows=1, ncols=4, left=0.725, right=0.95, wspace=0)
    sur_scatter = plot.add_subplot(sur[0, :-1])
    sur_hist = plot.add_subplot(sur[0, 3], sharey=sur_scatter)

    sur_scatter.patch.set_linewidth(3)
    sur_hist.patch.set_linewidth(3)

    # color for presentations
    # sur_scatter.patch.set_edgecolor("#006D5B")
    # sur_hist.patch.set_edgecolor('#006D5B')

    sur_spin = sur_mergers[:, 4]
    sur_gen1_spin = sur_spin[sur_merger_g1_mask]
    sur_gen2_spin = sur_spin[sur_merger_g2_mask]
    sur_genX_spin = sur_spin[sur_merger_gX_mask]

    sur_all_kick = sur_mergers[:, 16]
    sur_gen1_vkick = sur_all_kick[sur_merger_g1_mask]
    sur_gen2_vkick = sur_all_kick[sur_merger_g2_mask]
    sur_genX_vkick = sur_all_kick[sur_merger_gX_mask]

    # plot the 1g-1g mergers
    sur_scatter.scatter(sur_gen1_vkick, sur_gen1_spin,
                        s=styles.markersize_gen1,
                        marker=styles.marker_gen1,
                        edgecolor=styles.color_gen1,
                        facecolor='none',
                        alpha=styles.markeralpha_gen1,
                        label='1g-1g'
                        )

    # plot the 2g+ mergers
    sur_scatter.scatter(sur_gen2_vkick, sur_gen2_spin,
                        s=styles.markersize_gen2,
                        marker=styles.marker_gen2,
                        edgecolor=styles.color_gen2,
                        facecolor='none',
                        alpha=styles.markeralpha_gen2,
                        label='2g-1g or 2g-2g'
                        )

    # plot the 3g+ mergers
    sur_scatter.scatter(sur_genX_vkick, sur_genX_spin,
                        s=styles.markersize_genX,
                        marker=styles.marker_genX,
                        edgecolor=styles.color_genX,
                        facecolor='none',
                        alpha=styles.markeralpha_genX,
                        label=r'$\geq$3g-Ng'
                        )

    sur_scatter.grid(True, color='gray', ls='dashed', alpha=0.4)
    sur_scatter.set(
        xlabel=r'$v_{\textrm{kick}}$ [km/s]',
        title=r'$a_{final}^{sur}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    spin_bins = np.logspace(np.log10(sur_mergers[:, 4].min()), np.log10(sur_mergers[:, 4].max()), 50)

    sur_hist.grid(True, color='gray', ls='dashed', alpha=0.4)
    sur_hist_data = [sur_mergers[:, 4][sur_merger_g1_mask], sur_mergers[:, 4][sur_merger_g2_mask],
                     sur_mergers[:, 4][sur_merger_gX_mask]]
    sur_hist.hist(sur_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                  label=hist_label, stacked=True, orientation='horizontal')
    sur_hist.yaxis.tick_right()
    sur_hist_data_int = list(map(int, sur_hist_data[0] * 100))
    mode, spin_count_sur = stats.mode(sur_hist_data_int, axis=None, keepdims=False)

    sur_hist.set(
        xlabel=r'n',
        xlim=[0, int(spin_count_sur * 1.80)],
        # xticks=[300, 1000],
        xticks=np.linspace(int(spin_count_nosur * 0.50), int(spin_count_nosur * 1.30), 2)
    )

    plt.setp(sur_scatter.get_yticklabels(), visible=False)
    plt.setp(sur_hist.get_yticklabels(), visible=False)

    plt.tight_layout()
    # plt.show()
    plt.savefig(opts.plots_directory + '/spin_vkick.png', format='png')
    # plt.close()

    # ========================================
    # Spin vs Kick Velocity (q < 1/8)
    # ========================================

    plot = plt.figure(figsize=(12, 4), constrained_layout=False)

    # ======= NO SURROGATE =========
    nosur = plot.add_gridspec(nrows=1, ncols=4, left=0.05, right=0.275, wspace=0)
    nosur_scatter = plot.add_subplot(nosur[0, :-1])
    nosur_hist = plot.add_subplot(nosur[0, 3], sharey=nosur_scatter)

    new_nosur_gen1_vkick, new_nosur_gen2_vkick, new_nosur_genX_vkick = [], [], []
    new_nosur_gen1_spin, new_nosur_gen2_spin, new_nosur_genX_spin = [], [], []

    for i in range(len(nosur_gen1_vkick)):
        if nosur_gen1_mass_ratio[i] < (1.0 / 8.0):
            new_nosur_gen1_vkick.append(nosur_gen1_vkick[i])
            new_nosur_gen1_spin.append(nosur_gen1_spin[i])
    for i in range(len(nosur_gen2_vkick)):
        if nosur_gen2_mass_ratio[i] < (1.0 / 8.0):
            new_nosur_gen2_vkick.append(nosur_gen2_vkick[i])
            new_nosur_gen2_spin.append(nosur_gen2_spin[i])
    for i in range(len(nosur_genX_vkick)):
        if nosur_genX_mass_ratio[i] < (1.0 / 8.0):
            new_nosur_genX_vkick.append(nosur_genX_vkick[i])
            new_nosur_genX_spin.append(nosur_genX_spin[i])

    # plot the 1g-1g mergers
    nosur_scatter.scatter(new_nosur_gen1_vkick, new_nosur_gen1_spin,
                          s=styles.markersize_gen1,
                          marker=styles.marker_gen1,
                          edgecolor=styles.color_gen1,
                          facecolor='none',
                          alpha=styles.markeralpha_gen1,
                          label='1g-1g'
                          )

    # plot the 2g+ mergers
    nosur_scatter.scatter(new_nosur_gen2_vkick, new_nosur_gen2_spin,
                          s=styles.markersize_gen2,
                          marker=styles.marker_gen2,
                          edgecolor=styles.color_gen2,
                          facecolor='none',
                          alpha=styles.markeralpha_gen2,
                          label='2g-1g or 2g-2g'
                          )

    # plot the 3g+ mergers
    nosur_scatter.scatter(new_nosur_genX_vkick, new_nosur_genX_spin,
                          s=styles.markersize_genX,
                          marker=styles.marker_genX,
                          edgecolor=styles.color_genX,
                          facecolor='none',
                          alpha=styles.markeralpha_genX,
                          label=r'$\geq$3g-Ng'
                          )

    nosur_hist.grid(True, color='gray', ls='dashed')
    # nosur_hist_data = [new_nosur_gen1_spin, new_nosur_gen2_spin, new_nosur_genX_spin]
    nosur_hist_data = [nosur_mergers[:, 4][nosur_merger_g1_mask], nosur_mergers[:, 4][nosur_merger_g2_mask],
                       nosur_mergers[:, 4][nosur_merger_gX_mask]]
    nosur_hist.hist(nosur_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                    label=hist_label, stacked=True, orientation='horizontal')
    nosur_hist.yaxis.tick_right()
    nosur_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(nosur_hist.get_yticklabels(), visible=False)
    nosur_scatter.grid(True, color='gray', ls='dashed')
    nosur_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        ylabel=r'$a_{final}$',
        title=r'$a_{final}^{nosur}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    spin_bins = np.logspace(np.log10(sur_mergers[:, 4].min()), np.log10(sur_mergers[:, 4].max()), 50)

    if figsize == 'apj_col':
        nosur_scatter.legend(fontsize=5, loc='lower left')
    elif figsize == 'apj_page':
        nosur_scatter.legend()

    # ======= NO SURROGATE W/FILTER =========
    nosur_filter = plot.add_gridspec(nrows=1, ncols=4, left=0.725, right=0.95, wspace=0)
    nosur_filter_scatter = plot.add_subplot(nosur_filter[0, :-1])
    nosur_filter_hist = plot.add_subplot(nosur_filter[0, 3], sharey=nosur_filter_scatter)

    nosur_filter_spin = nosur_filter_mergers[:, 4]
    nosur_filter_gen1_spin = nosur_filter_spin[nosur_filter_merger_g1_mask]
    nosur_filter_gen2_spin = nosur_filter_spin[nosur_filter_merger_g2_mask]
    nosur_filter_genX_spin = nosur_filter_spin[nosur_filter_merger_gX_mask]

    nosur_filter_m1 = np.zeros(nosur_filter_mergers.shape[0])
    nosur_filter_m2 = np.zeros(nosur_filter_mergers.shape[0])
    nosur_filter_mass_ratio = np.zeros(nosur_filter_mergers.shape[0])

    for i in range(nosur_filter_mergers.shape[0]):
        if nosur_filter_mergers[i, 6] < nosur_filter_mergers[i, 7]:
            nosur_filter_m1[i] = nosur_filter_mergers[i, 7]
            nosur_filter_m2[i] = nosur_filter_mergers[i, 6]
            nosur_filter_mass_ratio[i] = nosur_filter_mergers[i, 6] / nosur_filter_mergers[i, 7]
        else:
            nosur_filter_mass_ratio[i] = nosur_filter_mergers[i, 7] / nosur_filter_mergers[i, 6]
            nosur_filter_m1[i] = nosur_filter_mergers[i, 6]
            nosur_filter_m2[i] = nosur_filter_mergers[i, 7]

    # (q,X_eff) Figure details here:
    # Want to highlight higher generation mergers on this plot
    # nosur_filter_chi_eff = nosur_filter_mergers[:, 3]

    # Get 1g-1g population
    # nosur_gen1_chi_eff = sur_chi_eff[sur_merger_g1_mask]
    nosur_filter_gen1_mass_ratio = nosur_filter_mass_ratio[nosur_filter_merger_g1_mask]
    # 2g-1g and 2g-2g population
    # nosur_gen2_chi_eff = sur_chi_eff[sur_merger_g2_mask]
    nosur_filter_gen2_mass_ratio = nosur_filter_mass_ratio[nosur_filter_merger_g2_mask]
    # >=3g-Ng population (i.e., N=1,2,3,4,...)
    # nosur_genX_chi_eff = sur_chi_eff[sur_merger_gX_mask]
    nosur_filter_genX_mass_ratio = nosur_filter_mass_ratio[nosur_filter_merger_gX_mask]

    new_nosur_filter_gen1_vkick, new_nosur_filter_gen2_vkick, new_nosur_filter_genX_vkick = [], [], []
    new_nosur_filter_gen1_spin, new_nosur_filter_gen2_spin, new_nosur_filter_genX_spin = [], [], []

    for i in range(len(nosur_filter_gen1_mass_ratio)):
        if nosur_filter_gen1_mass_ratio[i] < (1.0 / 8.0):
            new_nosur_filter_gen1_vkick.append(nosur_filter_gen1_vkick[i])
            new_nosur_filter_gen1_spin.append(nosur_filter_gen1_spin[i])
    for i in range(len(nosur_filter_gen2_mass_ratio)):
        if nosur_filter_gen2_mass_ratio[i] < (1.0 / 8.0):
            new_nosur_filter_gen2_vkick.append(nosur_filter_gen2_vkick[i])
            new_nosur_filter_gen2_spin.append(nosur_filter_gen2_spin[i])
    for i in range(len(nosur_filter_genX_mass_ratio)):
        if nosur_filter_genX_mass_ratio[i] < (1.0 / 8.0):
            new_nosur_filter_genX_vkick.append(nosur_filter_genX_vkick[i])
            new_nosur_filter_genX_spin.append(nosur_filter_genX_spin[i])

    # plot the 1g-1g mergers
    nosur_filter_scatter.scatter(new_nosur_filter_gen1_vkick, new_nosur_filter_gen1_spin,
                                 s=styles.markersize_gen1,
                                 marker=styles.marker_gen1,
                                 edgecolor=styles.color_gen1,
                                 facecolor='none',
                                 alpha=styles.markeralpha_gen1,
                                 label='1g-1g'
                                 )

    # plot the 2g+ mergers
    nosur_filter_scatter.scatter(new_nosur_filter_gen2_vkick, new_nosur_filter_gen2_spin,
                                 s=styles.markersize_gen2,
                                 marker=styles.marker_gen2,
                                 edgecolor=styles.color_gen2,
                                 facecolor='none',
                                 alpha=styles.markeralpha_gen2,
                                 label='2g-1g or 2g-2g'
                                 )

    # plot the 3g+ mergers
    nosur_filter_scatter.scatter(new_nosur_filter_genX_vkick, new_nosur_filter_genX_spin,
                                 s=styles.markersize_genX,
                                 marker=styles.marker_genX,
                                 edgecolor=styles.color_genX,
                                 facecolor='none',
                                 alpha=styles.markeralpha_genX,
                                 label=r'$\geq$3g-Ng'
                                 )

    nosur_filter_scatter.grid(True, color='gray', ls='dashed')
    nosur_filter_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        title=r'$a_{final}^{nosur\_filter}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    spin_bins = np.logspace(np.log10(sur_mergers[:, 4].min()), np.log10(sur_mergers[:, 4].max()), 50)

    nosur_filter_hist.grid(True, color='gray', ls='dashed')
    nosur_filter_hist_data = [nosur_filter_mergers[:, 4][nosur_filter_merger_g1_mask],
                              nosur_filter_mergers[:, 4][nosur_filter_merger_g2_mask],
                              nosur_filter_mergers[:, 4][nosur_filter_merger_gX_mask]]
    nosur_filter_hist.hist(nosur_filter_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9,
                           rwidth=0.8, label=hist_label, stacked=True, orientation='horizontal')
    nosur_filter_hist.yaxis.tick_right()
    nosur_filter_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(nosur_filter_hist.get_yticklabels(), visible=False)

    # ======= SURROGATE =========
    sur = plot.add_gridspec(nrows=1, ncols=4, left=0.275, right=0.5, wspace=0)
    sur_scatter = plot.add_subplot(sur[0, :-1])
    sur_hist = plot.add_subplot(sur[0, 3], sharey=sur_scatter)

    sur_spin = sur_mergers[:, 4]
    sur_gen1_spin = sur_spin[sur_merger_g1_mask]
    sur_gen2_spin = sur_spin[sur_merger_g2_mask]
    sur_genX_spin = sur_spin[sur_merger_gX_mask]

    sur_all_kick = sur_mergers[:, 16]
    sur_gen1_vkick = sur_all_kick[sur_merger_g1_mask]
    sur_gen2_vkick = sur_all_kick[sur_merger_g2_mask]
    sur_genX_vkick = sur_all_kick[sur_merger_gX_mask]

    sur_m1 = np.zeros(sur_mergers.shape[0])
    sur_m2 = np.zeros(sur_mergers.shape[0])
    sur_mass_ratio = np.zeros(sur_mergers.shape[0])

    for i in range(sur_mergers.shape[0]):
        if sur_mergers[i, 6] < sur_mergers[i, 7]:
            sur_m1[i] = sur_mergers[i, 7]
            sur_m2[i] = sur_mergers[i, 6]
            sur_mass_ratio[i] = sur_mergers[i, 6] / sur_mergers[i, 7]
        else:
            sur_m1[i] = sur_mergers[i, 6]
            sur_m2[i] = sur_mergers[i, 6]
            sur_mass_ratio[i] = sur_mergers[i, 7] / sur_mergers[i, 6]

    new_sur_gen1_vkick, new_sur_gen2_vkick, new_sur_genX_vkick = [], [], []
    new_sur_gen1_spin, new_sur_gen2_spin, new_sur_genX_spin = [], [], []

    for i in range(len(sur_gen1_mass_ratio)):
        if sur_gen1_mass_ratio[i] < (1.0 / 8.0):
            new_sur_gen1_vkick.append(sur_gen1_vkick[i])
            new_sur_gen1_spin.append(sur_gen1_spin[i])
    for i in range(len(sur_gen2_mass_ratio)):
        if sur_gen2_mass_ratio[i] < (1.0 / 8.0):
            new_sur_gen2_vkick.append(sur_gen2_vkick[i])
            new_sur_gen2_spin.append(sur_gen2_spin[i])
    for i in range(len(sur_genX_mass_ratio)):
        if sur_genX_mass_ratio[i] < (1.0 / 8.0):
            new_sur_genX_vkick.append(sur_genX_vkick[i])
            new_sur_genX_spin.append(sur_genX_spin[i])

    # plot the 1g-1g mergers
    sur_scatter.scatter(new_sur_gen1_vkick, new_sur_gen1_spin,
                        s=styles.markersize_gen1,
                        marker=styles.marker_gen1,
                        edgecolor=styles.color_gen1,
                        facecolor='none',
                        alpha=styles.markeralpha_gen1,
                        label='1g-1g'
                        )

    # plot the 2g+ mergers
    sur_scatter.scatter(new_sur_gen2_vkick, new_sur_gen2_spin,
                        s=styles.markersize_gen2,
                        marker=styles.marker_gen2,
                        edgecolor=styles.color_gen2,
                        facecolor='none',
                        alpha=styles.markeralpha_gen2,
                        label='2g-1g or 2g-2g'
                        )

    # plot the 3g+ mergers
    sur_scatter.scatter(new_sur_genX_vkick, new_sur_genX_spin,
                        s=styles.markersize_genX,
                        marker=styles.marker_genX,
                        edgecolor=styles.color_genX,
                        facecolor='none',
                        alpha=styles.markeralpha_genX,
                        label=r'$\geq$3g-Ng'
                        )

    sur_scatter.grid(True, color='gray', ls='dashed')
    sur_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        title=r'$a_{final}^{sur}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    sur_hist.grid(True, color='gray', ls='dashed')
    sur_hist_data = [sur_mergers[:, 4][sur_merger_g1_mask], sur_mergers[:, 4][sur_merger_g2_mask],
                     sur_mergers[:, 4][sur_merger_gX_mask]]
    sur_hist.hist(sur_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                  label=hist_label, stacked=True, orientation='horizontal')
    sur_hist.yaxis.tick_right()
    sur_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(sur_scatter.get_yticklabels(), visible=False)
    plt.setp(sur_hist.get_yticklabels(), visible=False)

    # ======= PRECESSION =========
    prec = plot.add_gridspec(nrows=1, ncols=4, left=0.5, right=0.725, wspace=0)
    prec_scatter = plot.add_subplot(prec[0, :-1])
    prec_hist = plot.add_subplot(prec[0, 3], sharey=prec_scatter)

    prec_spin = prec_mergers[:, 4]
    prec_gen1_spin = prec_spin[prec_merger_g1_mask]
    prec_gen2_spin = prec_spin[prec_merger_g2_mask]
    prec_genX_spin = prec_spin[prec_merger_gX_mask]

    prec_all_kick = prec_mergers[:, 16]
    prec_gen1_vkick = prec_all_kick[prec_merger_g1_mask]
    prec_gen2_vkick = prec_all_kick[prec_merger_g2_mask]
    prec_genX_vkick = prec_all_kick[prec_merger_gX_mask]

    prec_m1 = np.zeros(prec_mergers.shape[0])
    prec_m2 = np.zeros(prec_mergers.shape[0])
    prec_mass_ratio = np.zeros(prec_mergers.shape[0])

    for i in range(prec_mergers.shape[0]):
        if prec_mergers[i, 6] < prec_mergers[i, 7]:
            prec_m1[i] = prec_mergers[i, 7]
            prec_m2[i] = prec_mergers[i, 6]
            prec_mass_ratio[i] = prec_mergers[i, 6] / prec_mergers[i, 7]
        else:
            prec_m1[i] = prec_mergers[i, 6]
            prec_m2[i] = prec_mergers[i, 7]
            prec_mass_ratio[i] = prec_mergers[i, 7] / prec_mergers[i, 6]

    # Get 1g-1g population
    # nosur_gen1_chi_eff = sur_chi_eff[sur_merger_g1_mask]
    prec_gen1_mass_ratio = prec_mass_ratio[prec_merger_g1_mask]
    # 2g-1g and 2g-2g population
    # nosur_gen2_chi_eff = sur_chi_eff[sur_merger_g2_mask]
    prec_gen2_mass_ratio = prec_mass_ratio[prec_merger_g2_mask]
    # >=3g-Ng population (i.e., N=1,2,3,4,...)
    # nosur_genX_chi_eff = sur_chi_eff[sur_merger_g
    prec_genX_mass_ratio = prec_mass_ratio[prec_merger_gX_mask]

    new_prec_gen1_vkick, new_prec_gen2_vkick, new_prec_genX_vkick = [], [], []
    new_prec_gen1_spin, new_prec_gen2_spin, new_prec_genX_spin = [], [], []

    for i in range(len(prec_gen1_mass_ratio)):
        if prec_gen1_mass_ratio[i] < (1.0 / 8.0):
            new_prec_gen1_vkick.append(prec_gen1_vkick[i])
            new_prec_gen1_spin.append(prec_gen1_spin[i])
    for i in range(len(prec_gen2_mass_ratio)):
        if prec_gen2_mass_ratio[i] < (1.0 / 8.0):
            new_prec_gen2_vkick.append(prec_gen2_vkick[i])
            new_prec_gen2_spin.append(prec_gen2_spin[i])
    for i in range(len(prec_genX_mass_ratio)):
        if prec_genX_mass_ratio[i] < (1.0 / 8.0):
            new_prec_genX_vkick.append(prec_genX_vkick[i])
            new_prec_genX_spin.append(prec_genX_spin[i])

    prec_scatter.scatter(new_prec_gen1_vkick, new_prec_gen1_spin,
                         s=styles.markersize_gen1,
                         marker=styles.marker_gen1,
                         edgecolor=styles.color_gen1,
                         facecolor='none',
                         alpha=styles.markeralpha_gen1,
                         label='1g-1g'
                         )

    # plot the 2g+ mergers
    prec_scatter.scatter(new_prec_gen2_vkick, new_prec_gen2_spin,
                         s=styles.markersize_gen2,
                         marker=styles.marker_gen2,
                         edgecolor=styles.color_gen2,
                         facecolor='none',
                         alpha=styles.markeralpha_gen2,
                         label='2g-1g or 2g-2g'
                         )

    # plot the 3g+ mergers
    prec_scatter.scatter(new_prec_genX_vkick, new_prec_genX_spin,
                         s=styles.markersize_genX,
                         marker=styles.marker_genX,
                         edgecolor=styles.color_genX,
                         facecolor='none',
                         alpha=styles.markeralpha_genX,
                         label=r'$\geq$3g-Ng'
                         )

    prec_scatter.grid(True, color='gray', ls='dashed')
    prec_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        title=r'$a_{final}^{prec}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    prec_hist.grid(True, color='gray', ls='dashed')
    prec_hist_data = [prec_mergers[:, 4][prec_merger_g1_mask], prec_mergers[:, 4][prec_merger_g2_mask],
                      prec_mergers[:, 4][prec_merger_gX_mask]]
    prec_hist.hist(prec_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                   label=hist_label, stacked=True, orientation='horizontal')
    prec_hist.yaxis.tick_right()
    prec_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(prec_scatter.get_yticklabels(), visible=False)
    plt.setp(prec_hist.get_yticklabels(), visible=False)

    plt.tight_layout()
    plt.savefig(opts.plots_directory + '/spin_vkick_8below.png', format='png')

    # ========================================
    # Spin vs Kick Velocity (2g2g vs 2g1g)
    # ========================================

    plot = plt.figure(figsize=(12, 4), constrained_layout=False)

    # ======= NO SURROGATE =========
    nosur = plot.add_gridspec(nrows=1, ncols=4, left=0.05, right=0.2625, wspace=0)
    nosur_scatter = plot.add_subplot(nosur[0, :-1])
    nosur_hist = plot.add_subplot(nosur[0, 3], sharey=nosur_scatter)

    nosur_spin = nosur_mergers[:, 4]
    nosur_gen1_spin = nosur_spin[nosur_merger_g1_mask]
    nosur_gen2_spin = nosur_spin[nosur_merger_g2_mask]
    nosur_genX_spin = nosur_spin[nosur_merger_gX_mask]

    nosur_all_kick = nosur_mergers[:, 16]
    nosur_gen1_vkick = nosur_all_kick[nosur_merger_g1_mask]
    nosur_gen2_vkick = nosur_all_kick[nosur_merger_g2_mask]
    nosur_genX_vkick = nosur_all_kick[nosur_merger_gX_mask]

    # plot the 1g-1g mergers
    nosur_scatter.scatter(nosur_gen1_vkick, nosur_gen1_spin,
                          s=styles.markersize_gen1,
                          marker=styles.marker_gen1,
                          edgecolor=styles.color_gen1,
                          facecolor='none',
                          alpha=styles.markeralpha_gen1,
                          label='1g-1g'
                          )

    # plot the 2g+ mergers
    nosur_scatter.scatter(nosur_gen2_vkick, nosur_gen2_spin,
                          s=styles.markersize_gen2,
                          marker=styles.marker_gen2,
                          edgecolor=styles.color_gen2,
                          facecolor='none',
                          alpha=styles.markeralpha_gen2,
                          label='2g-1g or 2g-2g'
                          )

    # plot the 3g+ mergers
    nosur_scatter.scatter(nosur_genX_vkick, nosur_genX_spin,
                          s=styles.markersize_genX,
                          marker=styles.marker_genX,
                          edgecolor=styles.color_genX,
                          facecolor='none',
                          alpha=styles.markeralpha_genX,
                          label=r'$\geq$3g-Ng'
                          )

    nosur_scatter.grid(True, color='gray', ls='dashed')
    nosur_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        ylabel=r'$a_{final}$',
        title=r'$a_{final}^{nosur}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    spin_bins = np.logspace(np.log10(sur_mergers[:, 4].min()), np.log10(sur_mergers[:, 4].max()), 50)

    nosur_hist.grid(True, color='gray', ls='dashed')
    nosur_hist_data = [nosur_mergers[:, 4][nosur_merger_g1_mask], nosur_mergers[:, 4][nosur_merger_g2_mask],
                       nosur_mergers[:, 4][nosur_merger_gX_mask]]
    nosur_hist.hist(nosur_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                    label=hist_label, stacked=True, orientation='horizontal')
    nosur_hist.yaxis.tick_right()
    nosur_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(nosur_hist.get_yticklabels(), visible=False)

    if figsize == 'apj_col':
        nosur_scatter.legend(fontsize=5, loc='lower left')
    elif figsize == 'apj_page':
        nosur_scatter.legend()

    # ======= NO SURROGATE W/FILTER =========
    nosur_filter = plot.add_gridspec(nrows=1, ncols=4, left=0.275, right=0.4875, wspace=0)
    nosur_filter_scatter = plot.add_subplot(nosur_filter[0, :-1])
    nosur_filter_hist = plot.add_subplot(nosur_filter[0, 3], sharey=nosur_filter_scatter)

    nosur_filter_spin = nosur_filter_mergers[:, 4]
    nosur_filter_gen1_spin = nosur_filter_spin[nosur_filter_merger_g1_mask]
    nosur_filter_gen2_spin = nosur_filter_spin[nosur_filter_merger_g2_mask]
    nosur_filter_genX_spin = nosur_filter_spin[nosur_filter_merger_gX_mask]

    nosur_filter_all_kick = nosur_filter_mergers[:, 16]
    nosur_filter_gen1_vkick = nosur_filter_all_kick[nosur_filter_merger_g1_mask]
    nosur_filter_gen2_vkick = nosur_filter_all_kick[nosur_filter_merger_g2_mask]
    nosur_filter_genX_vkick = nosur_filter_all_kick[nosur_filter_merger_gX_mask]

    # plot the 1g-1g mergers
    nosur_filter_scatter.scatter(nosur_filter_gen1_vkick, nosur_filter_gen1_spin,
                                 s=styles.markersize_gen1,
                                 marker=styles.marker_gen1,
                                 edgecolor=styles.color_gen1,
                                 facecolor='none',
                                 alpha=styles.markeralpha_gen1,
                                 label='1g-1g'
                                 )

    # plot the 2g+ mergers
    nosur_filter_scatter.scatter(nosur_filter_gen2_vkick, nosur_filter_gen2_spin,
                                 s=styles.markersize_gen2,
                                 marker=styles.marker_gen2,
                                 edgecolor=styles.color_gen2,
                                 facecolor='none',
                                 alpha=styles.markeralpha_gen2,
                                 label='2g-1g or 2g-2g'
                                 )

    # plot the 3g+ mergers
    nosur_filter_scatter.scatter(nosur_filter_genX_vkick, nosur_filter_genX_spin,
                                 s=styles.markersize_genX,
                                 marker=styles.marker_genX,
                                 edgecolor=styles.color_genX,
                                 facecolor='none',
                                 alpha=styles.markeralpha_genX,
                                 label=r'$\geq$3g-Ng'
                                 )

    nosur_filter_scatter.grid(True, color='gray', ls='dashed')
    nosur_filter_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        title=r'$a_{final}^{nosur\_filter}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    spin_bins = np.logspace(np.log10(sur_mergers[:, 4].min()), np.log10(sur_mergers[:, 4].max()), 50)

    nosur_filter_hist.grid(True, color='gray', ls='dashed')
    nosur_filter_hist_data = [nosur_filter_mergers[:, 4][nosur_filter_merger_g1_mask],
                              nosur_filter_mergers[:, 4][nosur_filter_merger_g2_mask],
                              nosur_filter_mergers[:, 4][nosur_filter_merger_gX_mask]]
    nosur_filter_hist.hist(nosur_filter_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9,
                           rwidth=0.8, label=hist_label, stacked=True, orientation='horizontal')
    nosur_filter_hist.yaxis.tick_right()
    nosur_filter_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(nosur_filter_hist.get_yticklabels(), visible=False)

    # ======= SURROGATE =========
    sur = plot.add_gridspec(nrows=1, ncols=4, left=0.725, right=0.95, wspace=0)
    sur_scatter = plot.add_subplot(sur[0, :-1])
    sur_hist = plot.add_subplot(sur[0, 3], sharey=sur_scatter)

    sur_spin = sur_mergers[:, 4]
    sur_gen1_spin = sur_spin[sur_merger_g1_mask]
    sur_gen2_spin = sur_spin[sur_merger_g2_mask]
    sur_genX_spin = sur_spin[sur_merger_gX_mask]

    sur_all_kick = sur_mergers[:, 16]
    sur_gen1_vkick = sur_all_kick[sur_merger_g1_mask]
    sur_gen2_vkick = sur_all_kick[sur_merger_g2_mask]
    sur_genX_vkick = sur_all_kick[sur_merger_gX_mask]
    '''
    sur_gen2_vkick_2g1g, sur_gen2_spin_2g1g = [], []
    sur_gen2_spin_2g2g, sur_gen2_vkick_2g2g = [], []

    for i in range(len(sur_mergers)):
        if (sur_mergers[i, 12] == 1.0) and (sur_mergers[i, 13] == 2.0) and (sur_mergers[i, 4] in sur_gen2_spin) and (sur_mergers[i, 16] in sur_gen2_vkick):
            sur_gen2_vkick_2g1g.append(sur_mergers[i, 16])
            sur_gen2_spin_2g1g.append(sur_mergers[i, 4])
        elif (sur_mergers[i, 12] == 2.0) and (sur_mergers[i, 13] == 1.0) and (sur_mergers[i, 4] in sur_gen2_spin) and (sur_mergers[i, 16] in sur_gen2_vkick):
            sur_gen2_vkick_2g1g.append(sur_mergers[i, 16])
            sur_gen2_spin_2g1g.append(sur_mergers[i, 4])
        elif (sur_mergers[i, 12] == 2.0) and (sur_mergers[i, 13] == 2.0) and (sur_mergers[i, 4] in sur_gen2_spin) and (sur_mergers[i, 16] in sur_gen2_vkick):
            sur_gen2_vkick_2g2g.append(sur_mergers[i, 16])
            sur_gen2_spin_2g2g.append(sur_mergers[i, 4])'''

    new_sur_gen1_vkick_above, new_sur_gen2_vkick_above, new_sur_genX_vkick_above = [], [], []
    new_sur_gen1_spin_above, new_sur_gen2_spin_above, new_sur_genX_spin_above = [], [], []

    # Categorizing the mass ratio cutoff values - ABOVE
    for i in range(len(sur_gen1_mass_ratio)):
        if sur_gen1_mass_ratio[i] > 0.25:
            new_sur_gen1_vkick_above.append(sur_gen1_vkick[i])
            new_sur_gen1_spin_above.append(sur_gen1_spin[i])
    for i in range(len(sur_gen2_mass_ratio)):
        if sur_gen2_mass_ratio[i] > 0.25:
            new_sur_gen2_vkick_above.append(sur_gen2_vkick[i])
            new_sur_gen2_spin_above.append(sur_gen2_spin[i])
    for i in range(len(sur_genX_mass_ratio)):
        if sur_genX_mass_ratio[i] > 0.25:
            new_sur_genX_vkick_above.append(sur_genX_vkick[i])
            new_sur_genX_spin_above.append(sur_genX_spin[i])

    # Categorizing the mass ratio cutoff values - BELOW
    new_sur_gen1_vkick_below, new_sur_gen2_vkick_below, new_sur_genX_vkick_below = [], [], []
    new_sur_gen1_spin_below, new_sur_gen2_spin_below, new_sur_genX_spin_below = [], [], []
    for i in range(len(sur_gen1_mass_ratio)):
        if sur_gen1_mass_ratio[i] < 0.25:
            new_sur_gen1_vkick_below.append(sur_gen1_vkick[i])
            new_sur_gen1_spin_below.append(sur_gen1_spin[i])
    for i in range(len(sur_gen2_mass_ratio)):
        if sur_gen2_mass_ratio[i] < 0.25:
            new_sur_gen2_vkick_below.append(sur_gen2_vkick[i])
            new_sur_gen2_spin_below.append(sur_gen2_spin[i])
    for i in range(len(sur_genX_mass_ratio)):
        if sur_genX_mass_ratio[i] < 0.25:
            new_sur_genX_vkick_below.append(sur_genX_vkick[i])
            new_sur_genX_spin_below.append(sur_genX_spin[i])

    # plot the 1g-1g mergers
    sur_scatter.scatter(sur_gen1_vkick, sur_gen1_spin,
                        s=styles.markersize_gen1,
                        marker=styles.marker_gen1,
                        edgecolor=styles.color_gen1,
                        facecolor='none',
                        alpha=styles.markeralpha_gen1,
                        label='1g-1g'
                        )

    # plot the 2g+ mergers
    '''sur_scatter.scatter(sur_gen2_vkick_2g2g, sur_gen2_spin_2g2g,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor='blue',
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    sur_scatter.scatter(sur_gen2_vkick_2g1g, sur_gen2_spin_2g1g,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor='black',
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )'''
    # sur_scatter.scatter(sur_gen2_vkick, sur_gen2_spin,
    #            s=styles.markersize_gen2,
    #            marker=styles.marker_gen2,
    #            edgecolor=styles.color_gen2,
    #            facecolor='none',
    #            alpha=styles.markeralpha_gen2,
    #            label='2g-mg'
    #            )

    sur_scatter.scatter(new_sur_gen2_vkick_above, new_sur_gen2_spin_above,
                        s=styles.markersize_gen2,
                        marker=styles.marker_gen2,
                        edgecolor='black',
                        # edgecolor='#FF7E00',
                        facecolor='none',
                        alpha=styles.markeralpha_gen2,
                        label='2g-1g or 2g-2g (q > 0.25)'
                        )

    sur_scatter.scatter(new_sur_gen2_vkick_below, new_sur_gen2_spin_below,
                        s=styles.markersize_gen2,
                        marker=styles.marker_gen2,
                        edgecolor='blue',
                        # edgecolor='#FF00FE',
                        facecolor='none',
                        alpha=styles.markeralpha_gen2,
                        label='2g-1g or 2g-2g (q < 0.25)'
                        )

    # plot the 3g+ mergers
    # sur_scatter.scatter(sur_genX_vkick, sur_genX_spin,
    #            s=styles.markersize_genX,
    #            marker=styles.marker_genX,
    #            edgecolor=styles.color_genX,
    #            facecolor='none',
    #            alpha=styles.markeralpha_genX,
    #            label='3g-ng'
    #            )

    sur_scatter.scatter(new_sur_genX_vkick_above, new_sur_genX_spin_above,
                        s=styles.markersize_genX,
                        marker=styles.marker_genX,
                        edgecolor='#00800A',
                        facecolor='none',
                        alpha=styles.markeralpha_genX,
                        label=r'$\geq$3g-Ng (q > 0.25)'
                        )
    sur_scatter.scatter(new_sur_genX_vkick_below, new_sur_genX_spin_below,
                        s=styles.markersize_genX,
                        marker=styles.marker_genX,
                        edgecolor='#FF00FE',
                        facecolor='none',
                        alpha=styles.markeralpha_genX,
                        label=r'$\geq$3g-Ng (q < 0.25)'
                        )

    sur_scatter.grid(True, color='gray', ls='dashed')
    sur_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        title=r'$a_{final}^{sur}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    sur_hist.grid(True, color='gray', ls='dashed')
    sur_hist_data = [sur_mergers[:, 4][sur_merger_g1_mask], sur_mergers[:, 4][sur_merger_g2_mask],
                     sur_mergers[:, 4][sur_merger_gX_mask]]
    sur_hist.hist(sur_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                  label=hist_label, stacked=True, orientation='horizontal')
    sur_hist.yaxis.tick_right()
    sur_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(sur_scatter.get_yticklabels(), visible=False)
    plt.setp(sur_hist.get_yticklabels(), visible=False)

    if figsize == 'apj_col':
        sur_scatter.legend(fontsize=4, loc='lower left')
    elif figsize == 'apj_page':
        sur_scatter.legend()

    # ======= PRECESSION =========
    prec = plot.add_gridspec(nrows=1, ncols=4, left=0.5, right=0.7125, wspace=0)
    prec_scatter = plot.add_subplot(prec[0, :-1])
    prec_hist = plot.add_subplot(prec[0, 3], sharey=prec_scatter)

    prec_spin = prec_mergers[:, 4]
    prec_gen1_spin = prec_spin[prec_merger_g1_mask]
    prec_gen2_spin = prec_spin[prec_merger_g2_mask]
    prec_genX_spin = prec_spin[prec_merger_gX_mask]

    prec_all_kick = prec_mergers[:, 16]
    prec_gen1_vkick = prec_all_kick[prec_merger_g1_mask]
    prec_gen2_vkick = prec_all_kick[prec_merger_g2_mask]
    prec_genX_vkick = prec_all_kick[prec_merger_gX_mask]

    prec_scatter.scatter(prec_gen1_vkick, prec_gen1_spin,
                         s=styles.markersize_gen1,
                         marker=styles.marker_gen1,
                         edgecolor=styles.color_gen1,
                         facecolor='none',
                         alpha=styles.markeralpha_gen1,
                         label='1g-1g'
                         )

    # plot the 2g+ mergers
    prec_scatter.scatter(prec_gen2_vkick, prec_gen2_spin,
                         s=styles.markersize_gen2,
                         marker=styles.marker_gen2,
                         edgecolor=styles.color_gen2,
                         facecolor='none',
                         alpha=styles.markeralpha_gen2,
                         label='2g-1g or 2g-2g'
                         )

    # plot the 3g+ mergers
    prec_scatter.scatter(prec_genX_vkick, prec_genX_spin,
                         s=styles.markersize_genX,
                         marker=styles.marker_genX,
                         edgecolor=styles.color_genX,
                         facecolor='none',
                         alpha=styles.markeralpha_genX,
                         label=r'$\geq$3g-Ng'
                         )

    prec_scatter.grid(True, color='gray', ls='dashed')
    prec_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        title=r'$a_{final}^{prec}$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    prec_hist.grid(True, color='gray', ls='dashed')
    prec_hist_data = [prec_mergers[:, 4][prec_merger_g1_mask], prec_mergers[:, 4][prec_merger_g2_mask],
                      prec_mergers[:, 4][prec_merger_gX_mask]]
    prec_hist.hist(prec_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                   label=hist_label, stacked=True, orientation='horizontal')
    prec_hist.yaxis.tick_right()
    prec_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(prec_scatter.get_yticklabels(), visible=False)
    plt.setp(prec_hist.get_yticklabels(), visible=False)

    plt.tight_layout()
    # plt.show()
    plt.savefig(opts.plots_directory + '/spin_vkick_q_cutoff.png', format='png')
    # plt.close()

    # ========================================
    # Spin vs Kick Velocity (2g2g vs 2g1g) Surrogate Standalone Plot
    # ========================================

    plot = plt.figure(figsize=plotting.set_size(figsize), constrained_layout=False)

    sur = plot.add_gridspec(nrows=1, ncols=4, left=0.15, right=0.85, wspace=0)
    sur_scatter = plot.add_subplot(sur[0, :-1])
    sur_hist = plot.add_subplot(sur[0, -1], sharey=sur_scatter)

    sur_spin = sur_mergers[:, 4]
    sur_gen1_spin = sur_spin[sur_merger_g1_mask]
    sur_gen2_spin = sur_spin[sur_merger_g2_mask]
    sur_genX_spin = sur_spin[sur_merger_gX_mask]

    sur_all_kick = sur_mergers[:, 16]
    sur_gen1_vkick = sur_all_kick[sur_merger_g1_mask]
    sur_gen2_vkick = sur_all_kick[sur_merger_g2_mask]
    sur_genX_vkick = sur_all_kick[sur_merger_gX_mask]

    sur_gen2_vkick_2g1g, sur_gen2_spin_2g1g = [], []
    sur_gen2_spin_2g2g, sur_gen2_vkick_2g2g = [], []

    # sorting the mergers into the generation values (2g-2g, 2g-1g)
    for i in range(len(sur_mergers)):
        if (sur_mergers[i, 12] == 1.0) and (sur_mergers[i, 13] == 2.0) and (sur_mergers[i, 4] in sur_gen2_spin) and (
                sur_mergers[i, 16] in sur_gen2_vkick):
            sur_gen2_vkick_2g1g.append(sur_mergers[i, 16])
            sur_gen2_spin_2g1g.append(sur_mergers[i, 4])
        elif (sur_mergers[i, 12] == 2.0) and (sur_mergers[i, 13] == 1.0) and (sur_mergers[i, 4] in sur_gen2_spin) and (
                sur_mergers[i, 16] in sur_gen2_vkick):
            sur_gen2_vkick_2g1g.append(sur_mergers[i, 16])
            sur_gen2_spin_2g1g.append(sur_mergers[i, 4])
        elif (sur_mergers[i, 12] == 2.0) and (sur_mergers[i, 13] == 2.0) and (sur_mergers[i, 4] in sur_gen2_spin) and (
                sur_mergers[i, 16] in sur_gen2_vkick):
            sur_gen2_vkick_2g2g.append(sur_mergers[i, 16])
            sur_gen2_spin_2g2g.append(sur_mergers[i, 4])

    # setting the mass ratio limit above and below a specified value (q=0.6)
    new_sur_gen1_vkick_above, new_sur_gen2_vkick_above, new_sur_genX_vkick_above = [], [], []
    new_sur_gen1_spin_above, new_sur_gen2_spin_above, new_sur_genX_spin_above = [], [], []

    for i in range(len(sur_gen1_mass_ratio)):
        if sur_gen1_mass_ratio[i] > 0.6:
            new_sur_gen1_vkick_above.append(sur_gen1_vkick[i])
            new_sur_gen1_spin_above.append(sur_gen1_spin[i])
    for i in range(len(sur_gen2_mass_ratio)):
        if sur_gen2_mass_ratio[i] > 0.6:
            new_sur_gen2_vkick_above.append(sur_gen2_vkick[i])
            new_sur_gen2_spin_above.append(sur_gen2_spin[i])
    for i in range(len(sur_genX_mass_ratio)):
        if sur_genX_mass_ratio[i] > 0.6:
            new_sur_genX_vkick_above.append(sur_genX_vkick[i])
            new_sur_genX_spin_above.append(sur_genX_spin[i])

    new_sur_gen1_vkick_below, new_sur_gen2_vkick_below, new_sur_genX_vkick_below = [], [], []
    new_sur_gen1_spin_below, new_sur_gen2_spin_below, new_sur_genX_spin_below = [], [], []
    for i in range(len(sur_gen1_mass_ratio)):
        if sur_gen1_mass_ratio[i] < 0.6:
            new_sur_gen1_vkick_below.append(sur_gen1_vkick[i])
            new_sur_gen1_spin_below.append(sur_gen1_spin[i])
    for i in range(len(sur_gen2_mass_ratio)):
        if sur_gen2_mass_ratio[i] < 0.6:
            new_sur_gen2_vkick_below.append(sur_gen2_vkick[i])
            new_sur_gen2_spin_below.append(sur_gen2_spin[i])
    for i in range(len(sur_genX_mass_ratio)):
        if sur_genX_mass_ratio[i] < 0.6:
            new_sur_genX_vkick_below.append(sur_genX_vkick[i])
            new_sur_genX_spin_below.append(sur_genX_spin[i])

    # plot the 1g-1g mergers
    sur_scatter.scatter(sur_gen1_vkick, sur_gen1_spin,
                        s=styles.markersize_gen1,
                        marker=styles.marker_gen1,
                        edgecolor=styles.color_gen1,
                        facecolor='none',
                        alpha=styles.markeralpha_gen1,
                        label='1g-1g'
                        )

    # plot the 2g+ mergers
    '''sur_scatter.scatter(sur_gen2_vkick_2g2g, sur_gen2_spin_2g2g,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor='blue',
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    sur_scatter.scatter(sur_gen2_vkick_2g1g, sur_gen2_spin_2g1g,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor='black',
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )'''

    sur_scatter.scatter(new_sur_gen2_vkick_above, new_sur_gen2_spin_above,
                        s=styles.markersize_gen2,
                        marker=styles.marker_gen2,
                        edgecolor='black',
                        facecolor='none',
                        alpha=styles.markeralpha_gen2,
                        label='2g-1g or 2g-2g'
                        )

    sur_scatter.scatter(new_sur_gen2_vkick_below, new_sur_gen2_spin_below,
                        s=styles.markersize_gen2,
                        marker=styles.marker_gen2,
                        edgecolor='blue',
                        facecolor='none',
                        alpha=styles.markeralpha_gen2,
                        label='2g-1g or 2g-2g'
                        )

    # plot the 3g+ mergers
    sur_scatter.scatter(sur_genX_vkick, sur_genX_spin,
                        s=styles.markersize_genX,
                        marker=styles.marker_genX,
                        edgecolor=styles.color_genX,
                        facecolor='none',
                        alpha=styles.markeralpha_genX,
                        label=r'$\geq$3g-Ng'
                        )

    sur_scatter.grid(True, color='gray', ls='dashed')
    sur_scatter.set(
        xlabel=r'$v_{kick}$ [km/s]',
        ylabel=r'$a_{final}$',
        title=r'$a_{final}^{sur} \ [q_{cutoff}=0.6]$',
        xscale="log",
        axisbelow=True,
        xlim=([1.1e0, 2e3]),
        ylim=(0.21, 1.01)
    )

    sur_hist.grid(True, color='gray', ls='dashed')
    sur_hist_data = [sur_mergers[:, 4][sur_merger_g1_mask], sur_mergers[:, 4][sur_merger_g2_mask],
                     sur_mergers[:, 4][sur_merger_gX_mask]]
    sur_hist.hist(sur_hist_data, bins=spin_bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8,
                  label=hist_label, stacked=True, orientation='horizontal')
    sur_hist.yaxis.tick_right()
    sur_hist_data_int = list(map(int, sur_hist_data[0] * 100))
    mode, spin_count_sur = stats.mode(sur_hist_data_int, axis=None, keepdims=False)

    sur_hist.set(
        xlabel=r'n',
        xlim=[0, 1200],
        xticks=[300, 1000]
    )
    plt.setp(sur_scatter.get_yticklabels(), visible=True)
    plt.setp(sur_hist.get_yticklabels(), visible=False)

    sur_hist.set(
        xlabel=r'n',
        xlim=[0, int(spin_count_sur * 1.30)],
        # xticks=[300, 1000],
        xticks=np.linspace(int(spin_count_sur * 0.30), int(spin_count_sur * 0.90), 2)
    )
    if figsize == 'apj_col':
        sur_scatter.legend(fontsize=5, loc='lower left')
    elif figsize == 'apj_page':
        sur_scatter.legend()

    plt.tight_layout()
    # plt.show()
    plt.savefig(opts.plots_directory + '/spin_vkick_2g2g_sur.png', format='png')
    # plt.close()

    # ========================================
    # SUR - Q vs Spin
    # ========================================

    sur_spin = sur_mergers[:, 4]
    sur_gen1_spin = sur_spin[sur_merger_g1_mask]
    sur_gen2_spin = sur_spin[sur_merger_g2_mask]
    sur_genX_spin = sur_spin[sur_merger_gX_mask]

    sur_all_kick = sur_mergers[:, 16]
    sur_gen1_vkick = sur_all_kick[sur_merger_g1_mask]
    sur_gen2_vkick = sur_all_kick[sur_merger_g2_mask]
    sur_genX_vkick = sur_all_kick[sur_merger_gX_mask]

    fig, ax = plt.subplots(1, 2, figsize=(4.5, 2.5), sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0})
    # ax3 = fig.add_subplot(111)

    sur = ax[1]
    nosur = ax[0]

    # plot the 1g-1g mergers
    sur.scatter(sur_gen1_spin, sur_gen1_mass_ratio,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(sur_gen2_spin, sur_gen2_mass_ratio,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(sur_genX_spin, sur_genX_mass_ratio,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolor='none',
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set(
        xlabel='Spin (sur)',
        ylabel='Mass Ratio [q]',
        # xscale="log",
        axisbelow=True,
        # xlim=([2e0,4e3])
    )

    sur.grid(True, color='gray', ls='dashed')

    if figsize == 'apj_col':
        sur.legend(fontsize=6, loc='upper left')
    elif figsize == 'apj_page':
        sur.legend()

    # plt.savefig(opts.plots_directory + '/time_of_merger.png', format='png')
    # plt.close()

    # ========================================
    # NOSUR - Q vs Spin
    # ========================================

    nosur_spin = nosur_mergers[:, 4]
    nosur_gen1_spin = nosur_spin[nosur_merger_g1_mask]
    nosur_gen2_spin = nosur_spin[nosur_merger_g2_mask]
    nosur_genX_spin = nosur_spin[nosur_merger_gX_mask]

    nosur_all_kick = nosur_mergers[:, 16]
    nosur_gen1_vkick = nosur_all_kick[nosur_merger_g1_mask]
    nosur_gen2_vkick = nosur_all_kick[nosur_merger_g2_mask]
    nosur_genX_vkick = nosur_all_kick[nosur_merger_gX_mask]

    # plot the 1g-1g mergers
    nosur.scatter(nosur_gen1_spin, nosur_gen1_mass_ratio,
                  s=styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolor='none',
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    # plot the 2g+ mergers
    nosur.scatter(nosur_gen2_spin, nosur_gen2_mass_ratio,
                  s=styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolor='none',
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    # plot the 3g+ mergers
    nosur.scatter(nosur_genX_spin, nosur_genX_mass_ratio,
                  s=styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolor='none',
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    nosur.set(
        xlabel='Spin (nosur)',
        ylabel='Mass Ratio [q]',
        # xscale="log",
        axisbelow=True,
        # xlim=([2e0,4e3])
    )

    nosur.grid(True, color='gray', ls='dashed')

    if figsize == 'apj_col':
        nosur.legend(fontsize=6, loc='upper left')
    elif figsize == 'apj_page':
        nosur.legend()

    plt.savefig(opts.plots_directory + '/q_spin.png', format='png')
    # plt.close()

    # ========================================
    # SUR - Spin 2 vs Final Spin
    # ========================================

    sur_spin = sur_mergers[:, 4]
    sur_gen1_spin = sur_spin[sur_merger_g1_mask]
    sur_gen2_spin = sur_spin[sur_merger_g2_mask]
    sur_genX_spin = sur_spin[sur_merger_gX_mask]

    sur_spin2 = sur_mergers[:, 9]
    sur_gen1_spin2 = sur_spin2[sur_merger_g1_mask]
    sur_gen2_spin2 = sur_spin2[sur_merger_g2_mask]
    sur_genX_spin2 = sur_spin2[sur_merger_gX_mask]

    fig, ax = plt.subplots(1, 2, figsize=(5, 3), sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0})
    # ax3 = fig.add_subplot(111)

    sur = ax[1]
    nosur = ax[0]

    # plot the 1g-1g mergers
    sur.scatter(sur_gen1_spin, sur_gen1_spin2,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(sur_gen2_spin, sur_gen2_spin2,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(sur_genX_spin, sur_genX_spin2,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolor='none',
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set(
        xlabel=r'$a_{final}^{sur}$',
        # xscale="log",
        xlim=(0.2, 1.04),
        ylim=(-0.05, 1.04),
        axisbelow=True,
        # xlim=([2e0,4e3])
    )

    sur.grid(True, color='gray', ls='dashed')

    # if figsize == 'apj_col':
    #   sur.legend(fontsize=5, loc='upper left')
    # elif figsize == 'apj_page':
    #    sur.legend()

    # plt.savefig(opts.plots_directory + '/time_of_merger.png', format='png')
    # plt.close()

    # # ========================================
    # # NOSUR - Spin 2 vs Final Spin
    # # ========================================

    # nosur_spin = nosur_mergers[:, 4]
    # nosur_gen1_spin = nosur_spin[nosur_merger_g1_mask]
    # nosur_gen2_spin = nosur_spin[nosur_merger_g2_mask]
    # nosur_genX_spin = nosur_spin[nosur_merger_gX_mask]

    # nosur_spin2 = nosur_mergers[:, 9]
    # nosur_gen1_spin2 = nosur_spin2[nosur_merger_g1_mask]
    # nosur_gen2_spin2 = nosur_spin2[nosur_merger_g2_mask]
    # nosur_genX_spin2 = nosur_spin2[nosur_merger_gX_mask]

    # # plot the 1g-1g mergers
    # nosur.scatter(nosur_gen1_spin, nosur_gen1_spin2,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # nosur.scatter(nosur_gen2_spin, nosur_gen2_spin2,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # nosur.scatter(nosur_genX_spin, nosur_genX_spin2,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # nosur.set(
    #     xlabel=r'$a_{final}^{nosur}$',
    #     ylabel=r'$a_2$',
    #     xlim=(0.2, 1.04),
    #     ylim=(-0.05, 1.04),
    #     #xscale="log",
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # nosur.grid(True, color='gray', ls='dashed')

    # if figsize == 'apj_col':
    #     nosur.legend(fontsize=5, loc='upper left')
    # elif figsize == 'apj_page':
    #     nosur.legend()

    # plt.savefig(opts.plots_directory + '/spin2_spin_final.png', format='png')
    # #plt.close()

    # # ========================================
    # # SUR - Spin 2 vs Spin 1
    # # ========================================

    # sur_spin1 = sur_mergers[:, 8]
    # sur_gen1_spin1 = sur_spin1[sur_merger_g1_mask]
    # sur_gen2_spin1 = sur_spin1[sur_merger_g2_mask]
    # sur_genX_spin1 = sur_spin1[sur_merger_gX_mask]

    # sur_spin2 = sur_mergers[:, 9]
    # sur_gen1_spin2 = sur_spin2[sur_merger_g1_mask]
    # sur_gen2_spin2 = sur_spin2[sur_merger_g2_mask]
    # sur_genX_spin2 = sur_spin2[sur_merger_gX_mask]

    # fig, ax = plt.subplots(1, 2, figsize=(5, 3), sharey=True, gridspec_kw={'wspace':0, 'hspace':0})
    # #ax3 = fig.add_subplot(111)

    # sur = ax[1]
    # nosur = ax[0]

    # # plot the 1g-1g mergers
    # sur.scatter(sur_gen1_spin1, sur_gen1_spin2,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # sur.scatter(sur_gen2_spin1, sur_gen2_spin2,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # sur.scatter(sur_genX_spin1, sur_genX_spin2,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # sur.set(
    #     xlabel=r'$a_{final}^{sur}$',
    #     #xscale="log",
    #     xlim=(0, 1.02),
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # sur.grid(True, color='gray', ls='dashed')

    # if figsize == 'apj_col':
    #     sur.legend(fontsize=5, loc='lower left')
    # elif figsize == 'apj_page':
    #     sur.legend()

    # #plt.savefig(opts.plots_directory + '/time_of_merger.png', format='png')
    # #plt.close()

    # # ========================================
    # # NOSUR - Spin 2 vs Spin 1
    # # ========================================

    # nosur_spin1 = nosur_mergers[:, 8]
    # nosur_gen1_spin1 = nosur_spin1[nosur_merger_g1_mask]
    # nosur_gen2_spin1 = nosur_spin1[nosur_merger_g2_mask]
    # nosur_genX_spin1 = nosur_spin1[nosur_merger_gX_mask]

    # nosur_spin2 = nosur_mergers[:, 9]
    # nosur_gen1_spin2 = nosur_spin2[nosur_merger_g1_mask]
    # nosur_gen2_spin2 = nosur_spin2[nosur_merger_g2_mask]
    # nosur_genX_spin2 = nosur_spin2[nosur_merger_gX_mask]

    # # plot the 1g-1g mergers
    # nosur.scatter(nosur_gen1_spin1, nosur_gen1_spin2,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # nosur.scatter(nosur_gen2_spin1, nosur_gen2_spin2,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # nosur.scatter(nosur_genX_spin1, nosur_genX_spin2,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # nosur.set(
    #     xlabel=r'$a_1^{nosur}$',
    #     ylabel=r'$a_2^{nosur}$',
    #     xlim=(0, 1.02),
    #     #xscale="log",
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # nosur.grid(True, color='gray', ls='dashed')

    # if figsize == 'apj_col':
    #     nosur.legend(fontsize=5, loc='lower left')
    # elif figsize == 'apj_page':
    #     nosur.legend()

    # plt.savefig(opts.plots_directory + '/spin2_spin1_final.png', format='png')
    # #plt.close()

    # # ===============================
    # ### Spin_Final vs. Spin Angle Final
    # # ===============================
    # all_spin_angle = sur_mergers[:, 5]
    # gen1_spin_angle = all_spin_angle[sur_merger_g1_mask]
    # gen2_spin_angle = all_spin_angle[sur_merger_g2_mask]
    # genX_spin_angle = all_spin_angle[sur_merger_gX_mask]
    # fig = plt.figure(figsize=plotting.set_size(figsize))
    # ax3 = fig.add_subplot(111)
    # # plt.title("Time of Merger after AGN Onset")
    # # ax3.scatter(mergers[:,14]/1e6, mergers[:,2], s=pointsize_merge_time, color='darkolivegreen')
    # ax3.scatter(gen1_spin_angle,sur_gen1_spin,
    #                 s=styles.markersize_gen1,
    #                 marker=styles.marker_gen1,
    #                 edgecolor=styles.color_gen1,
    #                 facecolor='none',
    #                 alpha=styles.markeralpha_gen1,
    #                 label='1g-1g'
    #                 )
    #     # plot the 2g+ mergers
    # ax3.scatter(gen2_spin_angle, sur_gen2_spin,
    #                 s=styles.markersize_gen2,
    #                 marker=styles.marker_gen2,
    #                 edgecolor=styles.color_gen2,
    #                 facecolor='none',
    #                 alpha=styles.markeralpha_gen2,
    #                 label='2g-1g or 2g-2g'
    #                 )
    #     # plot the 3g+ mergers
    # ax3.scatter(genX_spin_angle, sur_genX_spin,
    #                 s=styles.markersize_genX,
    #                 marker=styles.marker_genX,
    #                 edgecolor=styles.color_genX,
    #                 facecolor='none',
    #                 alpha=styles.markeralpha_genX,
    #                 label=r'$\geq$3g-Ng'
    #                 )
    # ax3.set(
    #         xlabel=r'Spin angle',
    #         ylabel=r'a$_{\mathrm{remnant}}$',
    #         #xscale="log",
    #         #yscale="log",
    #         axisbelow=True
    #     )
    # plt.grid(True, color='gray', ls='dashed')
    # plt.savefig(opts.plots_directory + '/spin_vs_angle.png', format='png')

    # # ========================================
    # # NOSUR - Spin Angle Distributions
    # # ========================================

    # # Plot final spin distributions
    # fig, ax = plt.subplots(2, 1, figsize=plotting.set_size(figsize), sharex=True)

    # sur = ax[1]
    # nosur = ax[0]

    # counts, bins = np.histogram(nosur_mergers[:, 5])
    # bins = np.arange(int(nosur_mergers[:, 5].min()), int(nosur_mergers[:, 5].max()), 0.1)

    # nosur_hist_data = [nosur_mergers[:, 5][nosur_merger_g1_mask], nosur_mergers[:, 5][nosur_merger_g2_mask], nosur_mergers[:, 5][nosur_merger_gX_mask]]
    # hist_label = ['1g-1g', '2g-1g or 2g-2g', r'$\geq$3g-Ng']
    # hist_color = [styles.color_gen1, styles.color_gen2, styles.color_genX]

    # nosur.hist(nosur_hist_data, bins=bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label, stacked=True)

    # nosur.set_ylabel('Spin Angle Distribution - (nosur)', fontsize=5, wrap=True)
    # nosur.set_xlabel(r'Final Spin Angle')
    # #nosur.set_xscale('log')

    # if figsize == 'apj_col':
    #     nosur.legend(fontsize=5)
    # elif figsize == 'apj_page':
    #     nosur.legend()

    # # ========================================
    # # SUR - Spin Angle Distributions
    # # ========================================

    # # Plot final spin distribution
    # counts, bins = np.histogram(sur_mergers[:, 5])
    # bins = np.arange(int(sur_mergers[:, 5].min()), int(sur_mergers[:, 5].max()), 0.1)

    # sur_hist_data = [sur_mergers[:, 5][sur_merger_g1_mask], sur_mergers[:, 5][sur_merger_g2_mask], sur_mergers[:, 5][sur_merger_gX_mask]]
    # hist_label = ['1g-1g', '2g-1g or 2g-2g', r'$\geq$3g-Ng']
    # hist_color = [styles.color_gen1, styles.color_gen2, styles.color_genX]

    # sur.hist(nosur_hist_data, bins=bins, align='left', color=hist_color, alpha=0.9, rwidth=0.8, label=hist_label, stacked=True)

    # sur.set_ylabel('Spin Angle Distribution - (sur)', fontsize=5, wrap=True)
    # sur.set_xlabel(r'Final Spin Angle')
    # #nosur.set_xscale('log')

    # if figsize == 'apj_col':
    #     sur.legend(fontsize=5)
    # elif figsize == 'apj_page':
    #     sur.legend()

    # #nosur.savefig(opts.plots_directory + r"/merger_remnant_mass.png", format='png')
    # #plt.savefig(opts.plots_directory + r"/spin_final_dist.png", format='png')

    # # ========================================
    # # SUR - Q vs Kick Velocity
    # # ========================================

    # fig, ax = plt.subplots(1, 2, figsize=(4.5,2.5), sharey=True, gridspec_kw={'wspace':0, 'hspace':0})
    # #ax3 = fig.add_subplot(111)
    # sur = ax[0]
    # nosur = ax[1]
    # # plot the 1g-1g mergers
    # sur.scatter(sur_gen1_vkick, sur_gen1_mass_ratio,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # sur.scatter(sur_gen2_vkick, sur_gen2_mass_ratio,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # sur.scatter(sur_genX_vkick, sur_genX_mass_ratio,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # sur.set(
    #     xlabel=r'$v_{kick}^{sur}$',
    #     ylabel='Mass Ratio [q]',
    #     xscale="log",
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # #plt.savefig(opts.plots_directory + '/q_vkick.png', format='png')
    # #plt.close()

    # # ========================================
    # # NOSUR - Q vs Kick Velocity
    # # ========================================

    # # plot the 1g-1g mergers
    # nosur.scatter(nosur_gen1_vkick, nosur_gen1_mass_ratio,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # nosur.scatter(nosur_gen2_vkick, nosur_gen2_mass_ratio,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # nosur.scatter(nosur_genX_vkick, nosur_genX_mass_ratio,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # nosur.set(
    #     xlabel=r'$v_{kick}^{nosur}$',
    #     ylabel='Mass Ratio [q]',
    #     xscale="log",
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # nosur.grid(True, color='gray', ls='dashed')
    # if figsize == 'apj_col':
    #     nosur.legend(fontsize=5, loc='lower left')
    # elif figsize == 'apj_page':
    #     nosur.legend()

    # plt.savefig(opts.plots_directory + '/q_vkick.png', format='png')
    # #plt.close()

    # # ========================================
    # # NOSUR - Final Spin vs Mass
    # # ========================================

    # nosur_spin = nosur_mergers[:, 4]
    # nosur_gen1_spin = nosur_spin[nosur_merger_g1_mask]
    # nosur_gen2_spin = nosur_spin[nosur_merger_g2_mask]
    # nosur_genX_spin = nosur_spin[nosur_merger_gX_mask]

    # nosur_gen1_mass = nosur_mergers[:, 2][nosur_merger_g1_mask]
    # nosur_gen2_mass = nosur_mergers[:, 2][nosur_merger_g2_mask]
    # nosur_genX_mass = nosur_mergers[:, 2][nosur_merger_gX_mask]

    # # plot the 1g-1g mergers
    # nosur.scatter(nosur_gen1_mass, nosur_gen1_spin,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # nosur.scatter(nosur_gen2_mass, nosur_gen2_spin,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # nosur.scatter(nosur_genX_mass, nosur_genX_spin,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # nosur.set(
    #     xlabel=r'$mass$',
    #     ylabel=r'$a_{final}^{nosur}$',
    #     #xlim=(0.38, 1.02),
    #     #xscale="log",
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # nosur.grid(True, color='gray', ls='dashed')

    # if figsize == 'apj_col':
    #     nosur.legend(fontsize=5, loc='upper left')
    # elif figsize == 'apj_page':
    #     nosur.legend()

    # #plt.savefig(opts.plots_directory + '/spin_mass.png', format='png')
    # #plt.close()

    # # ========================================
    # # SUR - Final Spin vs Mass
    # # ========================================

    # sur_spin = sur_mergers[:, 4]
    # sur_gen1_spin = sur_spin[sur_merger_g1_mask]
    # sur_gen2_spin = sur_spin[sur_merger_g2_mask]
    # sur_genX_spin = sur_spin[sur_merger_gX_mask]

    # sur_gen1_mass = sur_mergers[:, 2][sur_merger_g1_mask]
    # sur_gen2_mass = sur_mergers[:, 2][sur_merger_g2_mask]
    # sur_genX_mass = sur_mergers[:, 2][sur_merger_gX_mask]

    # # plot the 1g-1g mergers
    # sur.scatter(sur_gen1_mass, sur_gen1_spin,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # sur.scatter(sur_gen2_mass, sur_gen2_spin,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # sur.scatter(sur_genX_mass, sur_genX_spin,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # sur.set(
    #     xlabel=r'$mass$',
    #     ylabel=r'$a_{final}^{sur}$',
    #     #xlim=(0.38, 1.02),
    #     #xscale="log",
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # sur.grid(True, color='gray', ls='dashed')

    # if figsize == 'apj_col':
    #     sur.legend(fontsize=5, loc='upper left')
    # elif figsize == 'apj_page':
    #     sur.legend()

    # plt.savefig(opts.plots_directory + '/spin_mass.png', format='png')
    # #plt.close()

    # # ========================================
    # # SUR - Q vs Final Mass
    # # ========================================

    # fig, ax = plt.subplots(1, 2, figsize=(4.5,2.5), sharey=True, gridspec_kw={'wspace':0, 'hspace':0})
    # #ax3 = fig.add_subplot(111)
    # sur = ax[0]
    # nosur = ax[1]
    # # plot the 1g-1g mergers
    # sur.scatter(sur_gen1_mass, sur_gen1_mass_ratio,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # sur.scatter(sur_gen2_mass, sur_gen2_mass_ratio,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # sur.scatter(sur_genX_mass, sur_genX_mass_ratio,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # sur.set(
    #     xlabel='Final Mass (sur)',
    #     ylabel='Mass Ratio [q]',
    #     xscale="log",
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # #plt.savefig(opts.plots_directory + '/q_vkick.png', format='png')
    # #plt.close()

    # # ========================================
    # # NOSUR - Q vs Kick Velocity
    # # ========================================

    # # plot the 1g-1g mergers
    # nosur.scatter(nosur_gen1_mass, nosur_gen1_mass_ratio,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label='1g-1g'
    #             )

    # # plot the 2g+ mergers
    # nosur.scatter(nosur_gen2_mass, nosur_gen2_mass_ratio,
    #             s=styles.markersize_gen2,
    #             marker=styles.marker_gen2,
    #             edgecolor=styles.color_gen2,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen2,
    #             label='2g-1g or 2g-2g'
    #             )

    # # plot the 3g+ mergers
    # nosur.scatter(nosur_genX_mass, nosur_genX_mass_ratio,
    #             s=styles.markersize_genX,
    #             marker=styles.marker_genX,
    #             edgecolor=styles.color_genX,
    #             facecolor='none',
    #             alpha=styles.markeralpha_genX,
    #             label=r'$\geq$3g-Ng'
    #             )

    # nosur.set(
    #     xlabel='Final mass (nosur)',
    #     ylabel='Mass Ratio [q]',
    #     xscale="log",
    #     axisbelow=True,
    #     #xlim=([2e0,4e3])
    # )

    # nosur.grid(True, color='gray', ls='dashed')
    # if figsize == 'apj_col':
    #     nosur.legend(fontsize=5, loc='lower left')
    # elif figsize == 'apj_page':
    #     nosur.legend()

    # plt.savefig(opts.plots_directory + '/q_mass.png', format='png')
    # #plt.close()

    # ========================================
    # Spin v Kick - Highlighting the 10+10 mergers
    # ========================================

    sur_spin1 = sur_mergers[:, 8]
    sur_gen1_spin1 = sur_spin1[sur_merger_g1_mask]
    sur_gen2_spin1 = sur_spin1[sur_merger_g2_mask]
    sur_genX_spin1 = sur_spin1[sur_merger_gX_mask]

    sur_spin2 = sur_mergers[:, 9]
    sur_gen1_spin2 = sur_spin2[sur_merger_g1_mask]
    sur_gen2_spin2 = sur_spin2[sur_merger_g2_mask]
    sur_genX_spin2 = sur_spin2[sur_merger_gX_mask]

    sur_spin = sur_mergers[:, 4]
    sur_gen1_spin = sur_spin[sur_merger_g1_mask]
    sur_gen2_spin = sur_spin[sur_merger_g2_mask]
    sur_genX_spin = sur_spin[sur_merger_gX_mask]

    sur_all_kick = sur_mergers[:, 16]
    sur_gen1_vkick = sur_all_kick[sur_merger_g1_mask]
    sur_gen2_vkick = sur_all_kick[sur_merger_g2_mask]
    sur_genX_vkick = sur_all_kick[sur_merger_gX_mask]

    sur_spin_angle1 = sur_mergers[:, 10]
    sur_gen1_spin_angle1 = sur_spin_angle1[sur_merger_g1_mask]
    sur_gen2_spin_angle1 = sur_spin_angle1[sur_merger_g2_mask]
    sur_genX_spin_angle1 = sur_spin_angle1[sur_merger_gX_mask]

    sur_spin_angle2 = sur_mergers[:, 11]
    sur_gen1_spin_angle2 = sur_spin_angle2[sur_merger_g1_mask]
    sur_gen2_spin_angle2 = sur_spin_angle2[sur_merger_g2_mask]
    sur_genX_spin_angle2 = sur_spin_angle2[sur_merger_gX_mask]

    new_sur_gen1_mass_final, new_sur_gen2_mass_final, new_sur_genX_mass_final = [], [], []
    new_sur_gen1_spin_final, new_sur_gen2_spin_final, new_sur_genX_spin_final = [], [], []

    new_sur_gen1_vkick_70, new_sur_gen2_vkick_70, new_sur_genX_vkick_70 = [], [], []
    new_sur_gen1_vkick_70200, new_sur_gen2_vkick_70200, new_sur_genX_vkick_70200 = [], [], []
    new_sur_gen1_vkick_200, new_sur_gen2_vkick_200, new_sur_genX_vkick_200 = [], [], []

    new_sur_gen1_spin_angle_1_70, new_sur_gen2_spin_angle_1_70, new_sur_genX_spin_angle_1_70 = [], [], []
    new_sur_gen1_spin_angle_2_70, new_sur_gen2_spin_angle_2_70, new_sur_genX_spin_angle_2_70 = [], [], []

    new_sur_gen1_spin_angle_1_70200, new_sur_gen2_spin_angle_1_70200, new_sur_genX_spin_angle_1_70200 = [], [], []
    new_sur_gen1_spin_angle_2_70200, new_sur_gen2_spin_angle_2_70200, new_sur_genX_spin_angle_2_70200 = [], [], []

    new_sur_gen1_spin_angle_1_200, new_sur_gen2_spin_angle_1_200, new_sur_genX_spin_angle_1_200 = [], [], []
    new_sur_gen1_spin_angle_2_200, new_sur_gen2_spin_angle_2_200, new_sur_genX_spin_angle_2_200 = [], [], []

    new_sur_gen1_spin_angle_1_all, new_sur_gen2_spin_angle_1_all, new_sur_genX_spin_angle_1_all = [], [], []
    new_sur_gen1_spin_angle_2_all, new_sur_gen2_spin_angle_2_all, new_sur_genX_spin_angle_2_all = [], [], []

    new_sur_gen1_vkick_all = []

    # Categorizing the masses that are less than 15 Msun and the spins that are chi > 0.1

    # Gen 1 mergers
    for i in range(len(sur_gen1_mass_1)):
        if sur_gen1_mass_1[i] < 15.00 and sur_gen1_mass_2[i] < 15.00 and sur_gen1_spin1[i] > 0.1 and sur_gen1_spin2[
            i] > 0.1:
            new_sur_gen1_mass_final.append(sur_gen1_mass[i])
            new_sur_gen1_spin_final.append(sur_gen1_spin[i])

            new_sur_gen1_spin_angle_1_all.append(sur_gen1_spin_angle1[i])
            new_sur_gen1_spin_angle_2_all.append(sur_gen1_spin_angle2[i])

            new_sur_gen1_vkick_all.append(sur_gen1_vkick[i])

            if sur_gen1_vkick[i] < 70.0:
                new_sur_gen1_spin_angle_1_70.append(sur_gen1_spin_angle1[i])
                new_sur_gen1_spin_angle_2_70.append(sur_gen1_spin_angle2[i])
                new_sur_gen1_vkick_70.append(sur_gen1_vkick[i])
            if sur_gen1_vkick[i] > 70.0 and sur_gen1_vkick[i] < 200.0:
                new_sur_gen1_spin_angle_1_70200.append(sur_gen1_spin_angle1[i])
                new_sur_gen1_spin_angle_2_70200.append(sur_gen1_spin_angle2[i])
                new_sur_gen1_vkick_70200.append(sur_gen1_vkick[i])
            if sur_gen1_vkick[i] > 200.0:
                new_sur_gen1_spin_angle_1_200.append(sur_gen1_spin_angle1[i])
                new_sur_gen1_spin_angle_2_200.append(sur_gen1_spin_angle2[i])
                new_sur_gen1_vkick_200.append(sur_gen1_vkick[i])

    # Gen 2 mergers
    for i in range(len(sur_gen2_mass_1)):
        if sur_gen2_mass_1[i] < 15.00 and sur_gen2_mass_2[i] < 15.00 and sur_gen2_spin1[i] > 0.1 and sur_gen2_spin2[
            i] > 0.1:
            new_sur_gen2_mass_final.append(sur_gen2_mass[i])
            new_sur_gen2_spin_final.append(sur_gen2_spin[i])

            new_sur_gen2_spin_angle_1_all.append(sur_gen2_spin_angle1[i])
            new_sur_gen2_spin_angle_2_all.append(sur_gen2_spin_angle2[i])

            if sur_gen2_vkick[i] < 70.0:
                new_sur_gen2_spin_angle_1_70.append(sur_gen2_spin_angle1[i])
                new_sur_gen2_spin_angle_2_70.append(sur_gen2_spin_angle2[i])
                new_sur_gen2_vkick_70.append(sur_gen2_vkick[i])
            if sur_gen1_vkick[i] > 70.0 and sur_gen2_vkick[i] < 200.0:
                new_sur_gen2_spin_angle_1_70200.append(sur_gen2_spin_angle1[i])
                new_sur_gen2_spin_angle_2_70200.append(sur_gen2_spin_angle2[i])
                new_sur_gen2_vkick_70200.append(sur_gen2_vkick[i])
            if sur_gen1_vkick[i] > 200.0:
                new_sur_gen2_spin_angle_1_200.append(sur_gen2_spin_angle1[i])
                new_sur_gen2_spin_angle_2_200.append(sur_gen2_spin_angle2[i])
                new_sur_gen2_vkick_200.append(sur_gen2_vkick[i])

    # Gen X mergers
    for i in range(len(sur_genX_mass_1)):
        if sur_genX_mass_1[i] < 15.00 and sur_genX_mass_2[i] < 15.00 and sur_genX_spin1[i] > 0.1 and sur_genX_spin2[
            i] > 0.1:
            new_sur_genX_mass_final.append(sur_genX_mass[i])
            new_sur_genX_spin_final.append(sur_genX_spin[i])

            new_sur_genX_spin_angle_1_all.append(sur_genX_spin_angle1[i])
            new_sur_genX_spin_angle_2_all.append(sur_genX_spin_angle2[i])

            if sur_genX_vkick[i] < 70.0:
                new_sur_genX_spin_angle_1_70.append(sur_genX_spin_angle1[i])
                new_sur_genX_spin_angle_2_70.append(sur_genX_spin_angle2[i])
                new_sur_genX_vkick_70.append(sur_genX_vkick[i])
            if sur_gen1_vkick[i] > 70.0 and sur_genX_vkick[i] < 200.0:
                new_sur_genX_spin_angle_1_70200.append(sur_genX_spin_angle1[i])
                new_sur_genX_spin_angle_2_70200.append(sur_genX_spin_angle2[i])
                new_sur_genX_vkick_70200.append(sur_genX_vkick[i])
            if sur_gen1_vkick[i] > 200.0:
                new_sur_genX_spin_angle_1_200.append(sur_genX_spin_angle1[i])
                new_sur_genX_spin_angle_2_200.append(sur_genX_spin_angle2[i])
                new_sur_genX_vkick_200.append(sur_genX_vkick[i])

    # ============= 70 km/s =============

    fig = plt.figure(figsize=(4, 4))
    sur = fig.add_subplot()

    # plot the 1g-1g mergers
    sur.scatter(new_sur_gen1_spin_angle_1_70, new_sur_gen1_spin_angle_2_70,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(new_sur_gen2_spin_angle_1_70, new_sur_gen2_spin_angle_2_70,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(new_sur_genX_spin_angle_1_70, new_sur_genX_spin_angle_2_70,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolors="none",
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set_xlabel('Spin 1')
    sur.set_ylabel('Spin 2')
    sur.set_title(r'Kick Velocity $> 70km/s$')

    if figsize == 'apj_col':
        sur.legend(fontsize=4, loc='upper left')
    elif figsize == 'apj_page':
        sur.legend()

    plt.grid(True, color='gray', ls='dashed')
    plt.savefig(opts.plots_directory + '/spin_angle1_spin_angle2_70.png', format='png')

    # ============= 70 - 200 km/s =============

    fig = plt.figure(figsize=(4, 4))
    sur = fig.add_subplot()

    # plot the 1g-1g mergers
    sur.scatter(new_sur_gen1_spin_angle_1_70200, new_sur_gen1_spin_angle_2_70200,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(new_sur_gen2_spin_angle_1_70200, new_sur_gen2_spin_angle_2_70200,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(new_sur_genX_spin_angle_1_70200, new_sur_genX_spin_angle_2_70200,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolors="none",
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set_xlabel('Spin 1')
    sur.set_ylabel('Spin 2')
    sur.set_title(r'Kick Velocity $> 70km/s$ and $< 200km/s$')

    if figsize == 'apj_col':
        sur.legend(fontsize=4, loc='upper left')
    elif figsize == 'apj_page':
        sur.legend()

    plt.grid(True, color='gray', ls='dashed')
    plt.savefig(opts.plots_directory + '/spin_angle1_spin_angle2_70200.png', format='png')

    # ============= > 200km/s =============

    fig = plt.figure(figsize=(4, 4))
    sur = fig.add_subplot()

    # plot the 1g-1g mergers
    sur.scatter(new_sur_gen1_spin_angle_1_200, new_sur_gen1_spin_angle_2_200,
                s=styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    # plot the 2g+ mergers
    sur.scatter(new_sur_gen2_spin_angle_1_200, new_sur_gen2_spin_angle_2_200,
                s=styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    # plot the 3g+ mergers
    sur.scatter(new_sur_genX_spin_angle_1_200, new_sur_genX_spin_angle_2_200,
                s=styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolors="none",
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set_xlabel('Spin 1')
    sur.set_ylabel('Spin 2')
    sur.set_title(r'Kick Velocity $> 200km/s$')

    if figsize == 'apj_col':
        sur.legend(fontsize=4, loc='upper left')
    elif figsize == 'apj_page':
        sur.legend()

    plt.grid(True, color='gray', ls='dashed')
    plt.savefig(opts.plots_directory + '/spin_angle1_spin_angle2_200.png', format='png')

    # plt.xlabel(r'$v_{kick}$ [km/s]')
    # plt.ylabel(r'$a_{final}$')
    # plt.title(r'$a_{final}^{sur}$')
    # plt.xscale("log")
    # #plt.axisbelow=True
    # plt.xlim(([1.1e0,2e3]))
    # plt.ylim((0.21, 1.01))

    # mean = np.mean(new_sur_gen1_vkick_final)
    # median = np.median(new_sur_gen1_vkick_final)
    # std = np.std(new_sur_gen1_vkick_final)
    # size = np.size(new_sur_gen1_vkick_final)

    # print("Mean: ", mean)
    # print("Median: ", median)
    # print("Std: ", std)
    # print("Size: ", size)

    # plt.vlines(mean, ymin=0, ymax=1, label='Mean', color='black')
    # plt.vlines(median, ymin=0, ymax=1, label='Median', color='green')
    # plt.vlines(std, ymin=0, ymax=1, label='Std', color='pink')

    # if figsize == 'apj_col':
    #         plt.legend(fontsize=4)
    # elif figsize == 'apj_page':
    #         plt.legend()

    # plt.grid(True, color='gray', ls='dashed', alpha=0.4)
    # plt.savefig(opts.plots_directory + '/spin_vkick_1010.png', format='png')

    # ================================================================================
    # SUR - Animating Spin Angle 1 vs Spin Angle 2 across Kick Velocities
    # ================================================================================

    from matplotlib.widgets import Button, Slider

    fig = plt.figure(figsize=(4, 4))
    sur = fig.add_subplot()
    fig.subplots_adjust(left=0.25, bottom=0.25)

    def sp1_to_sp2(spang1):
        return new_sur_gen1_spin_angle_2_all

    # plot mergers
    scatter = sur.scatter(new_sur_gen1_spin_angle_1_all, new_sur_gen1_spin_angle_2_all,
                          s=styles.markersize_gen1,
                          marker=styles.marker_gen1,
                          edgecolor=styles.color_gen1,
                          facecolor='none',
                          alpha=styles.markeralpha_gen1
                          )

    axfreq = fig.add_axes((0.25, 0.1, 0.65, 0.03))
    freq_slider = Slider(
        ax=axfreq,
        label='Kick Velocity',
        valmin=np.min(new_sur_gen1_vkick_all),
        valmax=np.max(new_sur_gen1_vkick_all),
        valinit=100,
    )

    def update(val):
        mask = new_sur_gen1_vkick_all <= freq_slider.val
        filtered = np.column_stack([new_sur_gen1_spin_angle_1_all[mask], new_sur_gen1_spin_angle_2_all[mask]])
        scatter.set_offsets(filtered)

        scatter.set_array(new_sur_gen1_vkick_all[mask])
        fig.canvas.draw_idle()

        print(f"Showing {mask.sum()} of {len(new_sur_gen1_vkick_all)} points")

    freq_slider.on_changed(update)

    resetax = fig.add_axes((0.8, 0.025, 0.1, 0.04))
    button = Button(resetax, 'Reset', hovercolor='0.975')

    def reset(event):
        freq_slider.reset()

    button.on_clicked(reset)

    plt.show()

    # from matplotlib.animation import PillowWriter

    # sp1, sp2, vkick = zip(*sorted(zip(new_sur_gen1_spin_angle_1_all, new_sur_gen1_spin_angle_2_all, new_sur_gen1_vkick_all)))

    # # plot mergers
    # scatter = sur.scatter(sp1, sp2,
    #             s=styles.markersize_gen1,
    #             marker=styles.marker_gen1,
    #             edgecolor=styles.color_gen1,
    #             facecolor='none',
    #             alpha=styles.markeralpha_gen1,
    #             label=vkick
    #             )

    # def update(frame):
    #     # for each frame, update the data stored on each artist.
    #     x = new_sur_gen1_spin_angle_1_all[:frame]
    #     y = new_sur_gen1_spin_angle_2_all[:frame]
    #     # update the scatter plot:
    #     data = np.stack([x, y]).T
    #     scatter.set_offsets(data)
    #     return (scatter)

    # filename = opts.plots_directory + "/spin1_spin2_vkick.gif"
    # ani = animation.FuncAnimation(fig=fig, func=update, frames=60, interval=10)
    # ani.save(filename=filename, writer="pillow")

    # filename = opts.plots_directory + "/spin1_spin2_vkick.gif"
    # writer = PillowWriter(fps=30)

    # l = plt.plot([],[], 'k-')

    # def sp1_to_sp2(spang1):
    #     return new_sur_gen1_spin_angle_2_all

    # xlist, ylist = [], []

    # with writer.saving(sur, filename, 100):
    #     for spin_ang1 in new_sur_gen1_spin_angle_1_all:
    #         xlist.append(spin_ang1)
    #         ylist.append(sp1_to_sp2(spin_ang1))

    #         l.set_data(xlist, ylist)
    #         writer.grab_frame()

    # fig, ax = plt.subplots()
    # rng = np.random.default_rng(19680801)
    # data = np.array([20, 20, 20, 20])
    # x = np.array([1, 2, 3, 4])

    # artists = []
    # colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:purple']
    # for i in range(len(new_sur_gen1_vki)):
    #     data += rng.integers(low=0, high=10, size=data.shape)
    #     container = ax.barh(x, data, color=colors)
    #     artists.append(container)
    #
    # ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=400)
    # plt.show()

    # ========================================
    # SUR - LVK and LISA Strain vs Freq
    # ========================================

    # Read LIGO O3 sensitivity data (https://git.ligo.org/sensitivity-curves/o3-sensitivity-curves)
    H1 = impresources.files(data) / 'O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_strain_asd.txt'
    L1 = impresources.files(data) / 'O3-L1-C01_CLEAN_SUB60HZ-1262141640.0_sensitivity_strain_asd.txt'

    # Adjust sep according to your delimiter (e.g., '\t' for tab-delimited files)
    dfh1 = pd.read_csv(H1, sep='\t', header=None)  # Use header=None if the file doesn't contain header row
    dfl1 = pd.read_csv(L1, sep='\t', header=None)

    # Access columns as df[0], df[1], ...
    f_H1 = dfh1[0]
    h_H1 = dfh1[1]

    # H - hanford
    # L - Ligvston

    # Using https://github.com/eXtremeGravityInstitute/LISA_Sensitivity/blob/master/LISA.py
    # Create LISA object
    sur_lisa = li.LISA()

    #   lisa_freq is the frequency (x-axis) being created
    #   lisa_sn is the sensitivity curve of LISA
    sur_lisa_freq = np.logspace(np.log10(1.0e-5), np.log10(1.0e0), 1000)
    sur_lisa_sn = sur_lisa.Sn(sur_lisa_freq)

    # Create figure and ax
    fig, svf_ax = plt.subplots(1, 2, figsize=(plotting.set_size(figsize)[0], 2.9), sharey=True)

    sur = svf_ax[1]
    nosur = svf_ax[0]

    sur.set_xlabel(r'f [Hz]')  # , fontsize=20, labelpad=10)
    sur.set_ylabel(r'${\rm h}_{\rm char}$')  # , fontsize=20, labelpad=10)
    # ax.tick_params(axis='both', which='major', labelsize=20)

    sur.set_xlim(0.5e-7, 1.0e+4)
    sur.set_ylim(1.0e-28, 1.0e-15)

    # ----------Finding the rows in which EMRIs signals are either identical or zeroes and removing them----------
    sur_identical_rows_emris = np.where(sur_emris[:, 5] == sur_emris[:, 6])
    sur_zero_rows_emris = np.where(sur_emris[:, 6] == 0)
    sur_emris = np.delete(sur_emris, sur_identical_rows_emris, 0)
    # emris = np.delete(emris,zero_rows_emris,0)
    sur_emris[~np.isfinite(sur_emris)] = 1.e-40

    # ----------Finding the rows in which LVKs signals are either identical or zeroes and removing them----------
    sur_identical_rows_lvk = np.where(sur_lvk[:, 5] == sur_lvk[:, 6])
    sur_zero_rows_lvk = np.where(sur_lvk[:, 6] == 0)
    sur_lvk = np.delete(sur_lvk, sur_identical_rows_lvk, 0)
    # lvk = np.delete(lvk,zero_rows_lvk,0)
    sur_lvk[~np.isfinite(sur_lvk)] = 1.e-40

    sur_lvk_g1_mask, sur_lvk_g2_mask, sur_lvk_gX_mask = make_gen_masks(sur_lvk, 7, 8)

    sur_lvk_g1 = sur_lvk[sur_lvk_g1_mask]
    sur_lvk_g2 = sur_lvk[sur_lvk_g2_mask]
    sur_lvk_gX = sur_lvk[sur_lvk_gX_mask]

    # ----------Setting the values for the EMRIs and LVKs signals and inverting them----------
    sur_inv_freq_emris = 1 / sur_emris[:, 6]
    # inv_freq_lvk = 1/lvk[:,6]
    # ma_freq_emris = np.ma.where(freq_emris == 0)
    # ma_freq_lvk = np.ma.where(freq_lvk == 0)
    # indices_where_zeros_emris = np.where(freq_emris = 0.)
    # freq_emris = freq_emris[freq_emris !=0]
    # freq_lvk = freq_lvk[freq_lvk !=0]

    # inv_freq_emris = 1.0/ma_freq_emris
    # inv_freq_lvk = 1.0/ma_freq_lvk
    # timestep =1.e4yr
    timestep = 1.e4
    sur_strain_per_freq_emris = sur_emris[:, 5] * sur_inv_freq_emris / timestep

    sur_strain_per_freq_lvk_g1 = sur_lvk_g1[:, 5] * (1 / sur_lvk_g1[:, 6]) / timestep
    sur_strain_per_freq_lvk_g2 = sur_lvk_g2[:, 5] * (1 / sur_lvk_g2[:, 6]) / timestep
    sur_strain_per_freq_lvk_gX = sur_lvk_gX[:, 5] * (1 / sur_lvk_gX[:, 6]) / timestep

    # plot the characteristic detector strains
    sur.loglog(sur_lisa_freq, np.sqrt(sur_lisa_freq * sur_lisa_sn),
               label='LISA Sensitivity',
               #   color='darkred',
               zorder=0)

    sur.loglog(f_H1, h_H1,
               label='LIGO O3, H1 Sensitivity',
               #   color='darkblue',
               zorder=0)

    sur.scatter(sur_emris[:, 6], sur_strain_per_freq_emris,
                s=0.4 * styles.markersize_gen1,
                alpha=styles.markeralpha_gen1
                )

    sur.scatter(sur_lvk_g1[:, 6], sur_strain_per_freq_lvk_g1,
                s=0.4 * styles.markersize_gen1,
                marker=styles.marker_gen1,
                edgecolor=styles.color_gen1,
                facecolor='none',
                alpha=styles.markeralpha_gen1,
                label='1g-1g'
                )

    sur.scatter(sur_lvk_g2[:, 6], sur_strain_per_freq_lvk_g2,
                s=0.4 * styles.markersize_gen2,
                marker=styles.marker_gen2,
                edgecolor=styles.color_gen2,
                facecolor='none',
                alpha=styles.markeralpha_gen2,
                label='2g-1g or 2g-2g'
                )

    sur.scatter(sur_lvk_gX[:, 6], sur_strain_per_freq_lvk_gX,
                s=0.4 * styles.markersize_genX,
                marker=styles.marker_genX,
                edgecolor=styles.color_genX,
                facecolor='none',
                alpha=styles.markeralpha_genX,
                label=r'$\geq$3g-Ng'
                )

    sur.set_yscale('log')
    sur.set_xscale('log')

    # ax.loglog(f_L1, h_L1,label = 'LIGO O3, L1 Sensitivity') # plot the characteristic strain
    # ax.loglog(f_gw,h,color ='black', label='GW150914')

    if figsize == 'apj_col':
        plt.legend(fontsize=7, loc="best")
    elif figsize == 'apj_page':
        plt.legend(loc="upper right")

    sur.set_xlabel(r'$\nu_{\rm GW}$ [Hz]')  # , fontsize=20, labelpad=10)
    sur.set_ylabel(r'$h_{\rm char}/\nu_{\rm GW}$')  # , fontsize=20, labelpad=10)

    # plt.savefig(opts.plots_directory + './gw_strain.png', format='png')
    # plt.show()

    # ========================================
    # NOSUR - LVK and LISA Strain vs Freq
    # ========================================

    # Read LIGO O3 sensitivity data (https://git.ligo.org/sensitivity-curves/o3-sensitivity-curves)
    H1 = impresources.files(data) / 'O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_strain_asd.txt'
    L1 = impresources.files(data) / 'O3-L1-C01_CLEAN_SUB60HZ-1262141640.0_sensitivity_strain_asd.txt'

    # Adjust sep according to your delimiter (e.g., '\t' for tab-delimited files)
    dfh1 = pd.read_csv(H1, sep='\t', header=None)  # Use header=None if the file doesn't contain header row
    dfl1 = pd.read_csv(L1, sep='\t', header=None)

    # Access columns as df[0], df[1], ...
    f_H1 = dfh1[0]
    h_H1 = dfh1[1]

    # H - hanford
    # L - Ligvston

    # Using https://github.com/eXtremeGravityInstitute/LISA_Sensitivity/blob/master/LISA.py
    # Create LISA object
    nosur_lisa = li.LISA()

    #   lisa_freq is the frequency (x-axis) being created
    #   lisa_sn is the sensitivity curve of LISA
    nosur_lisa_freq = np.logspace(np.log10(1.0e-5), np.log10(1.0e0), 1000)
    nosur_lisa_sn = nosur_lisa.Sn(nosur_lisa_freq)

    # Create figure and ax
    # fig, svf_ax = plt.subplots(1, 2, figsize=(plotting.set_size(figsize)[0], 2.9))

    nosur.set_xlabel(r'f [Hz]')  # , fontsize=20, labelpad=10)
    # nosur.set_ylabel(r'${\rm h}_{\rm char}$')  # , fontsize=20, labelpad=10)
    # ax.tick_params(axis='both', which='major', labelsize=20)

    nosur.set_xlim(0.5e-7, 1.0e+4)
    nosur.set_ylim(1.0e-28, 1.0e-15)

    # ----------Finding the rows in which EMRIs signals are either identical or zeroes and removing them----------
    nosur_identical_rows_emris = np.where(nosur_emris[:, 5] == nosur_emris[:, 6])
    nosur_zero_rows_emris = np.where(nosur_emris[:, 6] == 0)
    nosur_emris = np.delete(nosur_emris, nosur_identical_rows_emris, 0)
    # emris = np.delete(emris,zero_rows_emris,0)
    nosur_emris[~np.isfinite(nosur_emris)] = 1.e-40

    # ----------Finding the rows in which LVKs signals are either identical or zeroes and removing them----------
    nosur_identical_rows_lvk = np.where(nosur_lvk[:, 5] == nosur_lvk[:, 6])
    nosur_zero_rows_lvk = np.where(nosur_lvk[:, 6] == 0)
    nosur_lvk = np.delete(nosur_lvk, nosur_identical_rows_lvk, 0)
    # lvk = np.delete(lvk,zero_rows_lvk,0)
    nosur_lvk[~np.isfinite(nosur_lvk)] = 1.e-40

    nosur_lvk_g1_mask, nosur_lvk_g2_mask, nosur_lvk_gX_mask = make_gen_masks(nosur_lvk, 7, 8)

    nosur_lvk_g1 = nosur_lvk[nosur_lvk_g1_mask]
    nosur_lvk_g2 = nosur_lvk[nosur_lvk_g2_mask]
    nosur_lvk_gX = nosur_lvk[nosur_lvk_gX_mask]

    # ----------Setting the values for the EMRIs and LVKs signals and inverting them----------
    nosur_inv_freq_emris = 1 / nosur_emris[:, 6]
    # inv_freq_lvk = 1/lvk[:,6]
    # ma_freq_emris = np.ma.where(freq_emris == 0)
    # ma_freq_lvk = np.ma.where(freq_lvk == 0)
    # indices_where_zeros_emris = np.where(freq_emris = 0.)
    # freq_emris = freq_emris[freq_emris !=0]
    # freq_lvk = freq_lvk[freq_lvk !=0]

    # inv_freq_emris = 1.0/ma_freq_emris
    # inv_freq_lvk = 1.0/ma_freq_lvk
    # timestep =1.e4yr
    timestep = 1.e4
    nosur_strain_per_freq_emris = nosur_emris[:, 5] * nosur_inv_freq_emris / timestep

    nosur_strain_per_freq_lvk_g1 = nosur_lvk_g1[:, 5] * (1 / nosur_lvk_g1[:, 6]) / timestep
    nosur_strain_per_freq_lvk_g2 = nosur_lvk_g2[:, 5] * (1 / nosur_lvk_g2[:, 6]) / timestep
    nosur_strain_per_freq_lvk_gX = nosur_lvk_gX[:, 5] * (1 / nosur_lvk_gX[:, 6]) / timestep

    # plot the characteristic detector strains
    nosur.loglog(nosur_lisa_freq, np.sqrt(nosur_lisa_freq * nosur_lisa_sn),
                 label='LISA Sensitivity',
                 #   color='darkred',
                 zorder=0)

    nosur.loglog(f_H1, h_H1,
                 label='LIGO O3, H1 Sensitivity',
                 #   color='darkblue',
                 zorder=0)

    nosur.scatter(nosur_emris[:, 6], nosur_strain_per_freq_emris,
                  s=0.4 * styles.markersize_gen1,
                  alpha=styles.markeralpha_gen1
                  )

    nosur.scatter(nosur_lvk_g1[:, 6], nosur_strain_per_freq_lvk_g1,
                  s=0.4 * styles.markersize_gen1,
                  marker=styles.marker_gen1,
                  edgecolor=styles.color_gen1,
                  facecolor='none',
                  alpha=styles.markeralpha_gen1,
                  label='1g-1g'
                  )

    nosur.scatter(nosur_lvk_g2[:, 6], nosur_strain_per_freq_lvk_g2,
                  s=0.4 * styles.markersize_gen2,
                  marker=styles.marker_gen2,
                  edgecolor=styles.color_gen2,
                  facecolor='none',
                  alpha=styles.markeralpha_gen2,
                  label='2g-1g or 2g-2g'
                  )

    nosur.scatter(nosur_lvk_gX[:, 6], nosur_strain_per_freq_lvk_gX,
                  s=0.4 * styles.markersize_genX,
                  marker=styles.marker_genX,
                  edgecolor=styles.color_genX,
                  facecolor='none',
                  alpha=styles.markeralpha_genX,
                  label=r'$\geq$3g-Ng'
                  )

    nosur.set_yscale('log')
    nosur.set_xscale('log')

    # ax.loglog(f_L1, h_L1,label = 'LIGO O3, L1 Sensitivity') # plot the characteristic strain
    # ax.loglog(f_gw,h,color ='black', label='GW150914')

    if figsize == 'apj_col':
        plt.legend(fontsize=6, loc="best")
    elif figsize == 'apj_page':
        plt.legend(loc="upper right")

    nosur.set_xlabel(r'$\nu_{\rm GW}$ [Hz]')  # , fontsize=20, labelpad=10)
    # nosur.set_ylabel(r'$h_{\rm char}/\nu_{\rm GW}$')  # , fontsize=20, labelpad=10)

    plt.savefig(opts.plots_directory + '/gw_strain.png', format='png')
    # plt.show()


######## Execution ########
if __name__ == "__main__":
    main()
