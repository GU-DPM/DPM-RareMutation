from DPM_lib import os, plt, mpl, math, itertools, PdfPages  # ticker
from DPM_miscellaneous import DPM_miscellaneous_treatment_change_time
from DPM_generate import DPM_generate_pho_cumulative
from DPM_constant import *

""" This script plots the DPM (Dynamic Precision Medicine) model simulation results."""
if not PLT_INTERACTIVE:
    mpl.use('pdf')  # 'pdf'
else:
    # mpl.use('TkAgg')
    mpl.use('MacOSX')
print('Plot backend is ', mpl.get_backend())


def DPM_plot_all(Num_drug, strategy, Simduration, Limit_molecular_detection, Limit_radiologicdetection, Limit_mortality, pathsave, savename):
    X_all = dict()
    for i, (i_key, i_strategy) in enumerate(strategy.items()):
        if i_strategy is not None:
            title = i_key
            DPM_plot_1strategy(Num_drug, i_strategy, Simduration, Limit_molecular_detection, Limit_radiologicdetection, Limit_mortality, title,
                               pathsave)
            t, X, _ = i_strategy
            i_strategy = (t, X.sum(axis=0))
            X_all[i_key] = i_strategy
        else:
            X_all[i_key] = None
    DPM_plot_allstrategy(Num_drug, X_all, Simduration, Limit_radiologicdetection, Limit_mortality, pathsave)
    return


def DPM_plot_1strategy(Num_drug, strategy, Simduration, Limit_moleculardetection, Limit_radiologicdetection, Limit_mortality, title, pathsave):
    # Plot simulation result of one strategy.
    assert Num_drug == 2
    color_X = ('#45936D', '#edc165', '#009ee9', '#b95559')
    label_X = LABEL_X_2DRUG
    color_drug = ('#ED1C24', '#0000FF')
    label_drug = LABEL_2DRUG
    legend_order = LEGEND_ORDER_2DRUG
    legend_colnum = LEGEND_COLNUM_2DRUG

    t, X, d = strategy
    Xtotal = np.sum(X, axis=0)
    # X = np.vstack((Xtotal, X))
    diffpts = DPM_miscellaneous_treatment_change_time(d)

    # plt.rcParams['figure.dpi'] = FIG_DPI
    plt.rcParams['font.size'] = 14
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(9, 5)
    ax = plt.axes()
    for i_type in range(X.shape[0]):
        X1 = X[i_type, :]
        # X1[X1 < MINX] = MINX
        plt.plot(t/T_NORM, X1, color=color_X[i_type], label=label_X[i_type])

    ymaxval = 16  # math.ceil(math.log(np.max(X), LOGBASE)) + YINCREASE
    # xmaxval = t[-1]/T_NORM

    X_Limit_moleculardetection = Xtotal * Limit_moleculardetection
    plt.plot(t/T_NORM, X_Limit_moleculardetection, color='k', linestyle='dotted')
    if ymaxval > math.log(Limit_radiologicdetection, LOGBASE):
        plt.axhline(y=Limit_radiologicdetection, color='k', linestyle='dashdot')
        plt.text(np.max(t)/T_NORM/2, 2 * Limit_radiologicdetection, 'Limit of radiologic detection', fontsize=TEXTFONT)
    if ymaxval > math.log(Limit_mortality, LOGBASE):
        plt.axhline(y=Limit_mortality, color='k', linestyle='dashed')
        plt.text(np.max(t)/T_NORM/2, 2 * Limit_mortality, 'Limit of mortality', fontsize=TEXTFONT)

    diffpts = np.append(diffpts, (t[-1]))
    diffpts = diffpts.astype(int)
    for i_treat in range(len(diffpts)-1):
        d_i = d[:, diffpts[i_treat]:diffpts[i_treat+1]]
        d_i = np.unique(d_i, axis=1)
        d_i_sum = np.sum(d_i)
        y_basal = ymaxval - YINCREASE
        for i_drug in range(Num_drug):
            if d_i[i_drug] != 0:
                plt.fill_between(t[diffpts[i_treat]:diffpts[i_treat + 1]]/T_NORM, LOGBASE ** y_basal,
                                 LOGBASE ** (y_basal + YINCREASE * d_i[i_drug]/d_i_sum), color=color_drug[i_drug])
                plt.hlines(y=LOGBASE ** (y_basal + YINCREASE * d_i[i_drug] / d_i_sum), xmin=t[diffpts[i_treat]] / T_NORM,
                           xmax=t[diffpts[i_treat + 1] - 1]/T_NORM, color='black')
                # if i_treat != 0:
                plt.vlines(x=t[diffpts[i_treat]]/T_NORM, ymin=LOGBASE ** (ymaxval-YINCREASE), ymax=LOGBASE ** ymaxval, color='black')
                y_basal = y_basal + YINCREASE * d_i[i_drug]/d_i_sum

    # Plot artificial label.
    y_basal = ymaxval - YINCREASE
    d_lab = np.full(Num_drug, 1/Num_drug, dtype=float)
    t_lab = np.array(range(-int(t[-1]), -int(t[-1]/2), 1))
    for i in range(Num_drug):
        plt.fill_between(t_lab.tolist(), LOGBASE ** y_basal, LOGBASE ** (y_basal + YINCREASE * d_lab[i]), color=color_drug[i],
                         label=label_drug[i])
        y_basal = y_basal + YINCREASE * d_lab[i]
    # plt.plot(t_lab, t_lab * Limit_moleculardetection, color='k', linestyle='dotted', label='Limit of molecular detection')

    plt.xlabel('Time[days]')
    plt.ylabel('Number of cells')
    plt.title(title)
    handles, labels = plt.gca().get_legend_handles_labels()
    leg = plt.legend([handles[idx] for idx in legend_order], [labels[idx] for idx in legend_order],
                     bbox_to_anchor=LEGEND_BOXANCHOR, loc='upper center', ncol=legend_colnum)
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(0)
    ax.set_yscale('log', base=LOGBASE)
    # plt.xlim([XMINVAL, xmaxval + XINCREASE])
    plt.xlim([XMINVAL, Simduration + XINCREASE])
    plt.ylim([YMINVAL, LOGBASE ** ymaxval])
    ylocmaj = mpl.ticker.LogLocator(base=LOGBASE, numticks=len(range(-2, int(ymaxval), 2)))
    ax.yaxis.set_major_locator(ylocmaj)
    ylocmin = mpl.ticker.LogLocator(base=LOGBASE, subs=np.arange(0, LOGBASE) * 1/LOGBASE, numticks=100)
    ax.yaxis.set_minor_locator(ylocmin)
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_yticks([10**-2, 10**0, 10**3, 10**6, 10**9, 10 ** 13])
    ax.set_xticks(np.arange(0, math.ceil(Simduration/STEPSIZE_DEFAULT_VAL)+1)*STEPSIZE_DEFAULT_VAL)
    ax.grid()
    plt.tight_layout()
    plt.show()
    return


def DPM_plot_allstrategy(Num_drug, X_all, Simduration, Limit_radiologicdetection, Limit_mortality, pathsave):
    assert Num_drug == 2
    # plt.rcParams['figure.dpi'] = FIG_DPI
    color = ['r', 'b', 'b']
    linestyle = ['solid', 'solid', 'dashed']
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(9, 5)
    ax = plt.axes()
    Xmax = 0
    tmax = 0
    for i, (i_key, i_strategy) in enumerate(X_all.items()):
        if i_strategy is not None:
            t, i_X_all = i_strategy
            i_X_all[i_X_all < MINX] = MINX
            Xmax = max(Xmax, max(i_X_all))
            tmax = max(tmax, max(t))
            plt.plot(t/T_NORM, i_X_all, color=color[i], linestyle=linestyle[i], label=i_key)

    ymaxval = 16  # math.ceil(math.log(Xmax, LOGBASE)) + YINCREASE/2
    # xmaxval = tmax/T_NORM
    if ymaxval > math.log(Limit_radiologicdetection, LOGBASE):
        plt.axhline(y=Limit_radiologicdetection, color='k', linestyle='dashdot')
        plt.text(tmax/T_NORM/2, 2 * Limit_radiologicdetection, 'Limit of radiologic detection', fontsize=TEXTFONT+1)
    if ymaxval > math.log(Limit_mortality, LOGBASE):
        plt.axhline(y=Limit_mortality, color='k', linestyle='dashed')
        plt.text(tmax/T_NORM/2, 2 * Limit_mortality, 'Limit of mortality', fontsize=TEXTFONT+1)

    plt.xlabel('Time[Months]')
    plt.ylabel('Number of cells')
    plt.title('Total cell number')
    handles, labels = plt.gca().get_legend_handles_labels()
    leg = plt.legend(handles, labels, bbox_to_anchor=LEGEND_BOXANCHOR, loc='upper center', ncol=LEGEND_COLNUM_XALL-1)
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(0)
    ax.set_yscale('log', base=LOGBASE)
    plt.xlim([XMINVAL, Simduration + XINCREASE])
    plt.ylim([LOGBASE ** YMINVAL_XALL, LOGBASE ** ymaxval])
    ylocmaj = mpl.ticker.LogLocator(base=LOGBASE, numticks=len(range(-2, int(ymaxval), 2)))
    ax.yaxis.set_major_locator(ylocmaj)
    ylocmin = mpl.ticker.LogLocator(base=LOGBASE, subs=np.arange(0, LOGBASE) * 1 / LOGBASE, numticks=100)
    ax.yaxis.set_minor_locator(ylocmin)
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xticks(np.arange(0, math.ceil(Simduration/STEPSIZE_DEFAULT_VAL)+1)*STEPSIZE_DEFAULT_VAL)
    ax.set_yticks([10**-2, 10**0, 10**3, 10**6, 10**9, 10 ** 13])
    plt.show()
    plt.tight_layout()
    ax.grid()
    return


def DPM_plot_KM(data0, data, par, titlestr):
    if 'p' in par.keys():
        p_str = DPM_plot_pstr(par['p'])
    else:
        p_str = ''

    plt.rcParams['font.size'] = 20
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(12.25, 9)
    ax = plt.axes()
    l1, = plt.plot(data0['t'], data0['val'], color=par['color'][0])
    ax.fill_between(data0['t'], data0['interval_lower'], data0['interval_upper'],
                    color=par['color'][0], alpha=.3)
    l2, = plt.plot(data['t'], data['val'], color=par['color'][1])
    ax.fill_between(data['t'], data['interval_lower'], data['interval_upper'],
                    color=par['color'][1], alpha=.3)
    l3, = plt.plot([], [], ' ')
    l4, = plt.plot([], [], ' ')
    l5, = plt.plot([], [], ' ')
    if p_str:
        l6, = plt.plot([], [], ' ')
    else:
        l6 = None
    plt.title(f"{int(par['totalnum'])} of patients ({titlestr})")
    plt.xticks(list(range(0, par['duration']+par['xtick step'], par['xtick step'])))
    plt.yticks(np.arange(0, 1.2, .2).tolist())
    plt.xlabel('Survival time[days]')
    plt.ylabel('Fraction of surviving patients')
    if p_str:
        plt.legend([l1, l2, l3, l4, l5, l6],
                   [par['0'], par['1'],
                    f"median survival time of {par['0']} (days): {int(data0['median_survival'])}",
                    f"median survival time of {par['1']} (days): {int(data['median_survival'])}",
                    f"hazard ratio {par['1']} over {par['0']}: {round(par['hzr'], 3)}",
                    p_str],
                   loc='upper right', frameon=False)
    else:
        plt.legend([l1, l2, l3, l4, l5],
                   [par['0'], par['1'],
                    f"median survival time of {par['0']} (days): {int(data0['median_survival'])}",
                    f"median survival time of {par['1']} (days): {int(data['median_survival'])}",
                    f"hazard ratio {par['1']} over {par['0']}: {round(par['hzr'], 3)}"],
                   loc='upper right', frameon=False)
    return


def DPM_plot_KM_multi(km, par, titlestr):
    plt.rcParams['font.size'] = 14
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(12.25, 9)
    ax = plt.axes()
    keys = km.keys()
    legh = []
    legstr = []
    color = COLOR_KM_MULTI
    for i, i_key in enumerate(keys):
        i_km = km[i_key]
        l_1i, = plt.plot(i_km['t'], i_km['val'], color=color[i])
        ax.fill_between(i_km['t'], i_km['interval_lower'], i_km['interval_upper'], color=color[i], alpha=.3)
        l_2i, = plt.plot([], [], ' ')
        legh.extend([l_1i, l_2i])
        legstr.extend([i_key, f"median survival time of {i_key} (days): {int(i_km['median_survival'])}"])
        # l_i, = plt.plot([], [], ' ')
        # legh.append(l_i)
        # legstr.extend([f"# of patients {par['totalnum']:,}, ({titlestr})"])

    plt.title(f"{int(par['totalnum'])} of patients ({titlestr})")
    plt.xticks(list(range(0, par['duration']+par['xtick step'], par['xtick step'])))
    plt.yticks(np.arange(0, 1.2, .2).tolist())
    plt.xlabel('Survival time[days]')
    plt.ylabel('Fraction of surviving patients')
    plt.legend(legh, legstr, loc='upper right', frameon=False)
    return


def DPM_plot_KM_multi2(km, par, titlestr):
    # mpl.use('TkAgg')
    plt.rcParams['font.size'] = 18
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(16, 10)
    ax = plt.axes()
    keys = list(km.keys())
    color_mis = ['k', 'g', 'r', 'c', 'm', 'y', 'blueviolet', 'maroon', 'sienna', 'darkgoldenrod', 'darkolivegreen', 'forestgreen', 'darkcyan',
                 'dodgerblue', 'darkviolet', 'deeppink']
    color_nomis = 'b'
    lytype = ['-', '--', ':', '-.']
    num_mis = len(MISSPECIFICATION_LOD)

    legh, legstr, median_survival, totalnum, flag_samenum = [], [], [], [], False
    key_order = []
    for i in range(int(len(keys[:-1])/num_mis)):
        key_order.extend(keys[:-1][i::int(len(keys[:-1])/num_mis)])
    key_order.extend(keys[-1:])
    keys = key_order
    for i, i_key in enumerate(keys):
        i_km = km[i_key]
        if i_key is keys[-1]:
            linewidth, color, linestyle = 2, color_nomis, '-'
        else:
            linewidth, color, linestyle = 1, color_mis[int(np.floor(i / num_mis))], lytype[i % num_mis]

        l_1i, = plt.plot(i_km['t'], i_km['val'], color=color, linestyle=linestyle, linewidth=linewidth)
        ax.fill_between(i_km['t'], i_km['interval_lower'], i_km['interval_upper'], color=color, alpha=.3)
        median_survival.append(i_km['median_survival'])
        totalnum.append(i_km['num'])
        legh.extend([l_1i])
        legstr.append(i_key + f", ms (d): {int(i_km['median_survival'])}" + f", num: {int(i_km['num'])}")

    if all(i == totalnum[0] for i in totalnum):
        flag_samenum = True
        legstr = [i_legstr[:i_legstr.rindex(', num')] for i_legstr in legstr]

    sort_index = np.argsort(median_survival)[::-1]
    plt.title(f"{int(totalnum[0])} of patients ({titlestr})") if flag_samenum else plt.title(f"({titlestr})")
    plt.xticks(list(range(0, par['duration']+par['xtick step'], par['xtick step'])))
    plt.yticks(np.arange(0, 1.2, .2).tolist())
    plt.ylim([0, 1.02])
    plt.xlim([-50, 2000])
    plt.xlabel('Survival time[days]')
    plt.ylabel('Fraction of surviving patients')
    plt.legend([legh[i] for i in sort_index], [legstr[i] for i in sort_index], loc='upper right', frameon=False, ncol=2)

    if PLT_INTERACTIVE:
        plt.show()
    return


def DPM_plot_KM_multi_sel(km, par, LOD, titlestr):
    # mpl.use('TkAgg')
    plt.rcParams['font.size'] = 18
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(16, 10)
    ax = plt.axes()
    color = ['k', 'b']
    lytype = ['-', '--', ':', '-.']
    linewidth = 1

    legh, legstr, median_survival, totalnum, flag_samenum = [], [], [], [], False
    for i, i_strategy in enumerate(km.keys()):
        i_km = km[i_strategy]
        for j, j_mis in enumerate(MISSPECIFICATION_LOD):
            j_key = LOD + ' ' + str(j_mis)
            l_1i, = plt.plot(i_km[j_key]['t'], i_km[j_key]['val'], color=color[i], linestyle=lytype[j], linewidth=linewidth)
            # ax.fill_between(i_km[j_key]['t'], i_km[j_key]['interval_lower'], i_km[j_key]['interval_upper'], color=color, alpha=.3)
            median_survival.append(i_km[j_key]['median_survival'])
            totalnum.append(i_km[j_key]['num'])
            legh.extend([l_1i])
            legstr.append(i_strategy + ' ' + str(j_mis) + f", ms (d): {int(i_km[j_key]['median_survival'])}" +
                          f", num: {int(i_km[j_key]['num'])}")

    if all(i == totalnum[0] for i in totalnum):
        flag_samenum = True
        legstr = [i_legstr[:i_legstr.rindex(', num')] for i_legstr in legstr]

    plt.title(f"{int(totalnum[0])} of patients ({titlestr})") if flag_samenum else plt.title(f"({titlestr})")
    plt.xticks(list(range(0, par['duration']+par['xtick step'], par['xtick step'])))
    plt.yticks(np.arange(0, 1.2, .2).tolist())
    plt.ylim([0, 1.02])
    plt.xlim([-50, 2000])
    plt.xlabel('Survival time[days]')
    plt.ylabel('Fraction of surviving patients')
    plt.legend(legh, legstr, loc='upper right', frameon=False, ncol=1)
    plt.tight_layout()

    if PLT_INTERACTIVE:
        plt.show()
    return


def DPM_plot_LOD_multi(median_survial, LOD, titlestr, ylablestr='Median survival time (Days)', ylim=None):
    Strategy_name = list(median_survial[0].keys())
    LOD[-1] = 1e-9
    LOD = [float(i_LOD) for i_LOD in LOD]
    color = ['r', 'b', 'blueviolet', 'g', 'c', 'm', 'y', 'k']
    lsty = ['-', '--', ':', '-.']

    plt.rcParams['font.size'] = 21
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(23, 13)
    legh, legstr = [], []
    for i, i_strategy in enumerate(Strategy_name):
        for j, j_str in enumerate(MISSPECIFICATION_LOD_STR):
            j_idx = [i for i in range(len(median_survial[j][i_strategy]))]
            j_sort_index = np.argsort(np.array(LOD)[j_idx])
            j_l, = plt.plot(np.array(LOD)[j_idx][j_sort_index],
                            np.array(median_survial[j][i_strategy])[j_idx][j_sort_index], color=color[i],
                            linestyle=lsty[j], linewidth=2, marker='d', markersize=7)
            legh.append(j_l)
            legstr.append(f"set by {j_str}: {i_strategy}")

    plt.title(titlestr)
    plt.xlabel('LOD')
    plt.ylabel(ylablestr)
    plt.xscale("log")
    plt.xticks(LOD)
    if ylim is not None:
        plt.ylim(ylim)
    plt.legend(legh, legstr, bbox_to_anchor=(0.5, 1.17), loc='upper center', frameon=False, ncol=4)
    return


def DPM_plot_hz_ratio(hz_ratio, LOD, titlestr):
    LOD[-1] = 1e-9
    LOD = [float(i_LOD) for i_LOD in LOD]
    color = ['r', 'b', 'blueviolet', 'g', 'c', 'm', 'y', 'k']

    plt.rcParams['font.size'] = 21
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(23, 13)
    legh, legstr = [], []

    for i, i_str in enumerate(MISSPECIFICATION_LOD_STR):
        i_idx = [i for i in range(len(hz_ratio[i]))]
        i_sort_index = np.argsort(np.array(LOD)[i_idx])
        i_val = [list(hz_ratio[i][j].values()) for j in range(len(hz_ratio[i]))]
        i_val = list(itertools.chain(*i_val))

        i_l, = plt.plot(np.array(LOD)[i_idx][i_sort_index], np.array(i_val)[i_sort_index],
                        color=color[i], linestyle='-', linewidth=2, marker='d', markersize=7)

        legh.append(i_l)
        legstr.append(f"set by {i_str}: {list(hz_ratio[i][0].keys())[0]}")

    plt.title(titlestr)
    plt.xlabel('LOD')
    plt.ylabel('Hazard ratio')
    plt.xscale("log")
    plt.xticks(LOD)
    plt.legend(legh, legstr, loc='upper left', frameon=False)
    return


def DPM_plot_hz(hazard, LOD, titlestr):
    Strategy_name = list(hazard[0].keys())
    LOD[-1] = 1e-9
    LOD = [float(i_LOD) for i_LOD in LOD]
    color = ['r', 'b', 'blueviolet', 'g', 'c', 'm', 'y', 'k']
    lsty = ['-', '--', ':', '-.']

    plt.rcParams['font.size'] = 21
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(23, 13)
    legh, legstr = [], []

    def DPM_plot_hz_1(hazard_i, strategy_name_i):
        idx = [ii for ii in range(len(hazard_i[strategy_name_i]))]

        val = [hazard_i[i_strategy][ii]['mean'].values[0] for ii in range(len(hazard_i[i_strategy]))]
        val_lowerci = [hazard_i[i_strategy][ii]['ci'].values[0] for ii in range(len(hazard_i[i_strategy]))]
        val_upperci = [hazard_i[i_strategy][ii]['ci'].values[1] for ii in range(len(hazard_i[i_strategy]))]

        val_lowerci = np.subtract(np.array(val), np.array(val_lowerci))
        val_upperci = np.subtract(np.array(val_upperci), np.array(val))

        return idx, val, val_lowerci, val_upperci

    for i, i_strategy in enumerate(Strategy_name):
        for j, j_str in enumerate(MISSPECIFICATION_LOD_STR):
            j_idx, j_val, j_val_lowerci, j_val_upperci = DPM_plot_hz_1(hazard[j], i_strategy)
            j_sort_index = np.argsort(np.array(LOD)[j_idx])
            j_l = plt.errorbar(np.array(LOD)[j_idx][j_sort_index], np.array(j_val)[j_sort_index],
                               yerr=np.array([j_val_lowerci, j_val_upperci]).tolist(), color=color[i],
                               linestyle=lsty[j], linewidth=2, marker='d', markersize=7)

            legh.append(j_l)
            legstr.append(f"set by {j_str}: {i_strategy}")

    plt.title(titlestr)
    plt.xlabel('LOD')
    plt.ylabel('Hazard')
    plt.xscale("log")
    plt.xticks(LOD)
    plt.legend(legh, legstr, bbox_to_anchor=(0.5, 1.17), loc='upper center', frameon=False, ncol=4)

    return


def DPM_plot_pho_cumulative(pathsave):
    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    color = ['k', 'r', 'b', 'blueviolet', 'g', 'c', 'y']
    # color_mis = ['k', 'g', 'r', 'c', 'm', 'y', 'blueviolet', 'maroon', 'sienna', 'darkgoldenrod', 'darkolivegreen', 'forestgreen', 'darkcyan',
    #              'dodgerblue', 'darkviolet', 'deeppink']
    pho_lower = [7.1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    xticks = [7e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]
    yticks = np.power(10, np.arange(-7, 1, 1), dtype=float)
    k = 7.1e-7
    pho_min = k
    pho_max = 0.5
    for i, i_pho_lower in enumerate(pho_lower):
        num = 50
        step = np.power(10, np.linspace(np.log10(i_pho_lower), np.log10(pho_max), num))
        val = []
        for i_pho_upper in step:
            i_val = DPM_generate_pho_cumulative(i_pho_lower, i_pho_upper, k, pho_min, pho_max)
            val.append(i_val)

        plt.plot(step[1:], val[1:], label=i_pho_lower, color=color[i])
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel('rho_upper')
        plt.xticks(xticks)
        plt.yticks(list(yticks[1:]))
    plt.legend(loc='lower left', frameon=False, ncol=1)
    plt.ylabel('Values')
    plt.title('Fraction of neutral genome sites at different ranges of VAF')
    pathsave = os.path.join(pathsave, f"VAF.pdf")
    plt.savefig(pathsave, format='pdf', bbox_inches='tight')
    plt.close()
    return


def DPM_plot_contour(data_ref, data, par, titlestr):
    duration, binsize, name = par['duration'], par['binsize'], par['name']

    data_ref = [i_val if i_val <= duration else duration for i_val in data_ref]
    data = [i_val if i_val <= duration else duration for i_val in data]
    h, x, y = np.histogram2d(data_ref, data, bins=[np.arange(0, duration + binsize, binsize), np.arange(0, duration + binsize, binsize)],
                             range=[[0, duration + binsize], [0, duration + binsize]])
    h = h.T
    # h = h/np.sum(h)
    cmap = plt.colormaps['bwr']
    levels = mpl.ticker.MaxNLocator(nbins=math.ceil(np.log10(h).max())).tick_values(0, np.log10(h).max())
    norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # levels = np.array([2e-5, 1e-4, 1e-3, 1e-2, 1e-1])
    xc = (x[:-1] + x[1:]) / 2.0  # x-axis bin centers
    yc = (y[:-1] + y[1:]) / 2.0
    X, Y = np.meshgrid(xc, yc)
    plt.rcParams['font.size'] = 14
    plt.rcParams['figure.figsize'] = (13, 9)
    plt.figure()
    # fig = plt.gcf()
    # ax = plt.axes()
    # cp = ax.contour(X, Y, h, levels=levels, cmap='bwr')  # , colors='black'
    # ax.set_facecolor('white')
    # if 'id_sigbetter' in par.keys():
    #     idx = par['id_sigbetter']
    #     id_data_ref, id_data = [data_ref[i] for i in idx], [data[i] for i in idx]
    #     plt.scatter(id_data_ref, id_data, s=5, c='r', label=f'{name[1]} significantly better')
    # if 'id_sigworse' in par.keys():
    #     idx = par['id_sigworse']
    #     id_data_ref, id_data = [data_ref[i] for i in idx], [data[i] for i in idx]
    #     plt.scatter(id_data_ref, id_data, s=5, c='k', label=f'{name[0]} significantly better')
    im = plt.pcolormesh(X, Y, np.log10(h), shading='nearest', edgecolors='black', linewidth=0.5, norm=norm, cmap=cmap)
    cbar = plt.colorbar(im)
    plt.xlim([0, duration])
    plt.ylim([0, duration])
    plt.xticks(list(range(0, duration+binsize, par['xtick step'])))
    plt.yticks(list(range(0, duration+binsize, par['xtick step'])))
    plt.title(f'{int(np.sum(h))} of patients ({titlestr})')
    plt.xlabel(name[0] + ' [days]')
    plt.ylabel(name[1] + ' [days]')
    cbar.set_label('log10(# of patientss) ', rotation=270, labelpad=30)
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.08), frameon=False, ncol=2)
    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=1.0)
    return


def DPM_plot_density2d_surface(data_ref, data, par):
    duration, binsize, name = par['duration'], par['binsize'], par['name']
    data_ref = [i_val if i_val <= duration else duration for i_val in data_ref]
    data = [i_val if i_val <= duration else duration for i_val in data]
    h, x, y = np.histogram2d(data_ref, data, bins=[np.arange(0, duration + binsize, binsize), np.arange(0, duration + binsize, binsize)],
                             range=[[0, duration], [0, duration]])  # density=True
    h = h.T
    xc = (x[:-1] + x[1:]) / 2.0  # x-axis bin centers
    yc = (y[:-1] + y[1:]) / 2.0  # y-axis bin centers
    X, Y = np.meshgrid(xc, yc)

    cmap = plt.colormaps['bwr']
    levels = mpl.ticker.MaxNLocator(nbins=math.ceil(np.log10(h).max())).tick_values(0, np.log10(h).max())
    norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    plt.rcParams['font.size'] = 14
    plt.rcParams['figure.figsize'] = (13, 9)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    im = ax.plot_surface(X, Y, np.log10(h), rcount=np.inf, ccount=np.inf, vmin=np.log10(h).min(), vmax=np.log10(h).max(), norm=norm, cmap=cmap)
    plt.colorbar(im, location='left')
    plt.xticks(list(range(0, duration+binsize, par['xtick step'])))
    plt.yticks(list(range(0, duration+binsize, par['xtick step'])))
    ax.tick_params(axis='both', which='major', pad=5)
    ax.zaxis.set_ticks(np.arange(0, np.log10(h).max(), 1))
    ax.tick_params(axis='z', which='major', pad=15)
    plt.title(f'Total {int(np.sum(h))} of patients')
    plt.xlabel(name[0] + ' [days]', labelpad=20)
    plt.ylabel(name[1] + ' [days]', labelpad=20)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel('log10(# of patients)', labelpad=20, rotation=270)
    plt.tight_layout()
    # xlm = (-26.25, 628.25)
    # ylm = (-26.25, 628.25)
    # zlm = (0.0, 0.012463285612752565)
    azim = -50
    elev = 45
    ax.view_init(elev=elev, azim=azim)
    return


def DPM_plot_misspec_pop(ax, ind, hist, hist_pdf, hist_loguni, bins, i_cell,  legend_loc='upper left', legend_fontsize=9):
    ax.stairs(hist, bins, color='r', lw=1, label='sim')
    ax.stairs(hist_pdf, bins, color='k', lw=1, alpha=0.6, label='pdf')
    ax.stairs(hist_loguni, bins, color='b', lw=1, alpha=0.6, label='loguni')
    ax.set_xscale('log', base=LOGBASE)
    ax.set_ylabel('Probability')
    ax.set_xlim((1e-9, bins[-1]))
    # ax.set_ylim((0, 1))
    ax.set_title(f"{i_cell} cell, Step at {ind}", fontsize=8)

    ax.legend(loc=legend_loc, prop={'size': legend_fontsize}, frameon=False, ncol=1)
    return


def DPM_plot_pstr(p):
    if p > 0.05:
        p_str = 'n.s.'
    elif 0.01 < p <= 0.05:
        p_str = '*p ≤ 0.05'
    elif 0.001 < p <= 0.01:
        p_str = '**p ≤ 0.01'
    elif 0.0001 < p <= 0.01:
        p_str = '***p ≤ 0.001'
    else:
        p_str = '****p ≤ 0.0001'
    return p_str


def DPM_plot_df2pdf(df, filename, color_, numpages=(1, 1), pagesize=(11, 8.5)):
    def DPM_plot_df2pdf_1(df_, pagesize_):
        # draw as table
        alternating_colors = [color_[0] * len(df_.columns)]
        alternating_colors = alternating_colors * len(df_)
        # plt.rc('font', size=30)
        fig_, ax = plt.subplots(figsize=pagesize)
        ax.axis('tight')
        ax.axis('off')
        the_table = ax.table(cellText=df_.values,
                             rowLabels=None,
                             colLabels=df_.columns,
                             rowColours=color_[1]*len(df_),
                             colColours=color_[1]*len(df_.columns),
                             cellColours=alternating_colors,
                             loc='center')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)
        the_table.scale(1, 4)
        return fig_

    with PdfPages(filename) as pdf:
        nh, nv = numpages
        rows_per_page = len(df)//nh
        cols_per_page = len(df.columns)//nv
        for i in range(0, nh):
            for j in range(0, nv):
                page = df.iloc[(i*rows_per_page):min((i+1)*rows_per_page, len(df)),
                               (j*cols_per_page):min((j+1)*cols_per_page, len(df.columns))]
                fig = DPM_plot_df2pdf_1(page, pagesize)
                if nh > 1 or nv > 1:
                    # Add part/page number at bottom-center of page
                    fig.text(0.5, 0.5/pagesize[0],
                             'Part-{}x{}: Page-{}'.format(i+1, j+1, i*nv + j + 1), ha='center', fontsize=9)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close()
    return
