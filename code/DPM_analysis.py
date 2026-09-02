from DPM_lib import deepcopy, plt, pd, KaplanMeierFitter, CoxPHFitter, statistics, ExponentialFitter, tabulate, itertools, \
    os, np, bz2, re, functools, pickle, tqdm, random, sns, supervenn, rv_continuous, loguniform, pearsonr, alt
from DPM_plot import DPM_plot_KM, DPM_plot_KM_multi, DPM_plot_KM_multi2, DPM_plot_KM_multi_sel, DPM_plot_contour, DPM_plot_LOD_multi, \
    DPM_plot_hz_ratio, DPM_plot_hz, DPM_plot_misspec_pop
from DPM_constant import *
from DPM_generate import DPM_generate_pdf
from DPM_miscellaneous import DPM_miscellaneous_fillful


def DPM_analysis_stat(stoptime, Strategy_name, Simduration, bool_select=None, titlestr='All'):
    num_betterall, num_sigbetterall, num_sigbetter, num_sigworse, paramID_sigbetter, paramID_sigworse = \
        dict(), dict(), dict(), dict(), dict(), dict()
    if bool_select is not None:
        stoptime = {item: list(itertools.compress(value, bool_select)) for (item, value) in stoptime.items()}
    df_stoptime = pd.DataFrame.from_dict(stoptime, dtype=int)
    for i_strategy in Strategy_name:
        df_stoptime[i_strategy] = df_stoptime[i_strategy].apply(lambda x: x if x < Simduration else Simduration)
    survival_median = df_stoptime.median()[Strategy_name].astype(int)
    survival_mean = df_stoptime.mean()[Strategy_name].apply(np.ceil).astype(int)
    survival_5y = df_stoptime[Strategy_name].apply(lambda x: x >= Simduration).apply(sum) / df_stoptime.shape[0] * 100
    survival_5y = survival_5y.round(2)
    for i, i_strategy in enumerate(df_stoptime[Strategy_name]):
        sub_Strategy_name = [ii for ii in Strategy_name if ii is not i_strategy]
        # No. of cases strategy numerically better than all others
        num_betterall[i_strategy] = \
            np.array(df_stoptime[i_strategy] > df_stoptime[sub_Strategy_name].max(axis=1)).sum()
        # No. of cases strategy significantly better than all others
        num_sigbetterall[i_strategy] = \
            np.array((df_stoptime[i_strategy] > df_stoptime[sub_Strategy_name].max(axis=1) + 8 * 7) &
                     (df_stoptime[i_strategy] > 1.25 * df_stoptime[sub_Strategy_name].max(axis=1))).sum()
        i_sigbetter, i_sigworse, i_paramID_sigbetter, i_paramID_sigworse = dict(), dict(), dict(), dict()
        i_stoptime = df_stoptime[i_strategy]
        i_sigbetter[i_strategy], i_sigworse[i_strategy], i_paramID_sigbetter[i_strategy], i_paramID_sigworse[i_strategy] = \
            'N.A.', 'N.A.', 'N.A.', 'N.A.'
        for i_substrategy in sub_Strategy_name:
            bool_sigbetter = (df_stoptime[i_substrategy] > i_stoptime + 8 * 7) & (df_stoptime[i_substrategy] > 1.25 * i_stoptime)
            i_sigbetter[i_substrategy] = np.array(bool_sigbetter).sum()
            i_paramID_sigbetter[i_substrategy] = bool_sigbetter.index[bool_sigbetter]

            bool_sigworse = (i_stoptime > df_stoptime[i_substrategy] + 8 * 7) & (i_stoptime > 1.25 * df_stoptime[i_substrategy])
            i_sigworse[i_substrategy] = np.array(bool_sigworse).sum()
            i_paramID_sigworse[i_substrategy] = bool_sigworse.index[bool_sigworse]

        i_sigbetter = {k: i_sigbetter[k] for k in Strategy_name}
        i_sigworse = {k: i_sigworse[k] for k in Strategy_name}
        i_paramID_sigbetter = {k: i_paramID_sigbetter[k] for k in Strategy_name}
        i_paramID_sigworse = {k: i_paramID_sigworse[k] for k in Strategy_name}

        num_sigbetter[i_strategy] = i_sigbetter
        num_sigworse[i_strategy] = i_sigworse
        paramID_sigbetter[i_strategy] = i_paramID_sigbetter
        paramID_sigworse[i_strategy] = i_paramID_sigworse

    data = list()
    data.append(['Median survival, day'] + list(survival_median.values))
    data.append(['Mean survival, day'] + list(survival_mean.values))
    data.append(['Survival at 5y, %'] + list(survival_5y.values))
    data.append(['No. of cases strategy numerically better than all others'] + list(num_betterall.values()))
    data.append(['No. of cases strategy significantly better than all others'] + list(num_sigbetterall.values()))
    # for i_strategy in Strategy_name:
    #     i_num = num_sigbetter[i_strategy]
    #     data.append([f'No. of cases significantly better than {i_strategy}'] + list(i_num.values()))
    # data.append(['Median survival, day'] + list(survival_median.values))
    # data.append(['Mean survival, day'] + list(survival_mean.values))
    # data.append(['Survival at 5y, %'] + list(survival_5y.values))
    # data.append(['No. of cases strategy numerically better than all others'] + list(num_betterall.values()))
    # data.append(['No. of cases strategy significantly better than all others'] + list(num_sigbetterall.values()))
    for i_strategy in Strategy_name:
        i_num_sigbetter, i_num_sigworse = num_sigbetter[i_strategy], num_sigworse[i_strategy]
        data.append([f'No. of cases significantly better than {i_strategy}'] + list(i_num_sigbetter.values()))
        data.append([f'No. of cases significantly worse than {i_strategy}'] + list(i_num_sigworse.values()))

    col_names = [f'Patients ({titlestr})']
    col_names.extend(Strategy_name)
    print(tabulate(data, headers=col_names, tablefmt='rst', numalign='left'))

    # df_data = pd.DataFrame(data, columns=[f'Patients ({titlestr})']+Strategy_name)

    # Significantly better means at least 8 week of absolute improvement
    par = {'duration': Simduration, 'binsize': 14, 'xtick step': 300}
    for i in list(itertools.combinations(Strategy_name, 2)):
        if 'strategy0' in i:
            ref = stoptime['strategy0']
            treat = stoptime[list(filter(lambda x: x != 'strategy0', i))[0]]
            name = ('strategy0', list(filter(lambda x: x != 'strategy0', i))[0])
        else:
            ref, treat = stoptime[i[0]], stoptime[i[1]]
            name = i
        # i_ID_sigbetter = paramID_sigbetter[name[0]][name[1]]
        # i_ID_sigworese = paramID_sigworse[name[0]][name[1]]
        par['name'] = name
        # par['id_sigbetter'] = i_ID_sigbetter
        # par['id_sigworse'] = i_ID_sigworese
        DPM_plot_contour(ref, treat, par, titlestr)
        # DPM_plot_density2d_surface(ref, treat, par)

    km = dict()
    for i_strategy in Strategy_name:
        i_stoptime = stoptime[i_strategy]
        i_stoptime = [i_val if i_val <= Simduration else Simduration + 1 for i_val in i_stoptime]
        km[i_strategy] = DPM_analysis_KM(i_stoptime, Simduration)

    # plot all KM
    par = {'duration': Simduration, 'xtick step': 300, 'totalnum': df_stoptime.shape[0]}
    DPM_plot_KM_multi(km, par, titlestr)

    p = DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration)
    hz = dict()
    for i in list(itertools.combinations(Strategy_name, 2)):
        if 'strategy0' in i:
            ref = stoptime['strategy0']
            treat = stoptime[list(filter(lambda x: x != 'strategy0', i))[0]]
            name = ('strategy0', list(filter(lambda x: x != 'strategy0', i))[0])
        else:
            ref, treat = stoptime[i[0]], stoptime[i[1]]
            name = i
        hz[name] = DPM_analysis_HZ(ref, treat, Simduration)

    for i in hz.keys():
        idx_p = [idx for idx, i_val in enumerate(p.keys()) if set(i) == set(i_val)][0]
        km_ref = km[i[0]]
        km_treat = km[i[1]]
        par = {'color': ['k', 'b'], '0': i[0], '1': i[1], 'hzr': hz[i], 'p': p[list(p.keys())[idx_p]],
               'duration': Simduration, 'xtick step': 300, 'totalnum': df_stoptime.shape[0]}
        DPM_plot_KM(km_ref, km_treat, par, titlestr)
        plt.close('all')
    return


def DPM_analysis_misspec_pop(para, pop, mispopcorr_mis, info, LOD, pathsave, pathload):
    class MutationFraction(rv_continuous):
        def _pdf(self, x, k, LOD, const):
            return (1.0/const)*DPM_generate_pdf(x, k, LOD)

    def DPM_analysis_pop_1(pop_, mutationFraction_distribution_, LOD_, pathsave_):
        bins = np.insert(np.logspace(-9, np.log10(i_LOD), num=int(1e3)), 0, 0)
        for i, i_cell in enumerate(ALL_POSSIBLE_CELLTYPE_2DRUG):
            plt.rcParams['font.size'] = 9
            fig, ax = plt.subplots(6, 7)
            fig.set_size_inches(15, 10)
            plt.tight_layout(pad=0.9, h_pad=0.9, w_pad=0.9)
            pearson_pdf = list()
            pearson_loguni = list()
            with tqdm(total=len(pop_), ncols=150) as pbar:
                for j, j_pop_ in enumerate(pop_):
                    j_ax = ax[j // 7, j % 7]
                    j_pop_ = np.asarray(j_pop_)
                    j_total = np.sum(j_pop_, axis=1)
                    j_pop_ = j_pop_[:, i]/j_total
                    ind = np.flatnonzero(j_pop_ < LOD_)
                    j_pop_ = j_pop_[ind]

                    pdf_samples = mutationFraction_distribution_.rvs(k=MUTATION_RATE, LOD=i_LOD, const=1, size=len(j_pop_))
                    loguni_samples = loguniform.rvs(MUTATION_RATE, i_LOD, size=len(j_pop_))

                    hist, _ = np.histogram(j_pop_, bins=bins, weights=[1/len(j_pop_)]*len(j_pop_))
                    hist_pdf, _ = np.histogram(pdf_samples, bins, weights=[1/len(j_pop_)]*len(j_pop_))
                    hist_loguni, _ = np.histogram(loguni_samples, bins, weights=[1/len(j_pop_)]*len(j_pop_))

                    # plt.stairs(hist, bins, color='r', lw=1, label='sim')
                    # plt.stairs(hist_pdf, bins, color='r', lw=1, label='sim')
                    # plt.stairs(hist_loguni, bins, color='r', lw=1, label='sim')

                    pearson_pdf.append(pearsonr(hist, hist_pdf)[0])
                    pearson_loguni.append(pearsonr(hist, hist_loguni)[0])

                    DPM_plot_misspec_pop(j_ax, j, hist, hist_pdf, hist_loguni, bins, i_cell, legend_fontsize=8)
                    pbar.update(1)

                ax[-1, -1].plot(np.arange(len(pearson_pdf)), pearson_pdf, color='k', lw=1, label='pdf')
                ax[-1, -1].plot(np.arange(len(pearson_loguni)), pearson_loguni, color='b', lw=1, label='loguni')
                ax[-1, -1].set_ylabel('Pearson corr')
                ax[-1, -1].set_xlabel('Step')
                ax[-1, -1].legend(loc='best', prop={'size': 8}, frameon=False, ncol=1)
            i_pathsave_ = os.path.join(pathsave_, f"{i_cell} cell.pdf")
            plt.savefig(i_pathsave_, format='pdf', bbox_inches='tight')
            plt.close('all')
        return

    paramID_sigbetter = info['total']['paramID sigbetter']
    pathsave = os.path.join(pathsave, 'pop')
    misspecification_LOD = list(mispopcorr_mis[0].keys())
    strategy_name = list(pop.keys())
    celltype = ['R1pop', 'R2pop']
    miscorr = {i_strategy: {i_mis: {i_celltype: [] for i_celltype in celltype} for i_mis in misspecification_LOD} for i_strategy in strategy_name}
    for i, i_LOD in enumerate(LOD):
        i_mispopcorr = mispopcorr_mis[i]
        i_LOD = float(i_LOD)
        mutationFraction_distribution = MutationFraction(name='mutationFraction_distribution', a=MUTATION_RATE, b=i_LOD)
        # i_pdf = mutationFraction_distribution.pdf(x=bins, k=MUTATION_RATE, LOD=i_LOD, const=1)
        # i_loguni = loguniform.pdf(bins, MUTATION_RATE, i_LOD)
        for i_strategy in strategy_name:
            i_pop = pop[i_strategy]
            # DPM_analysis_pop_1(i_pop, mutationFraction_distribution, i_LOD, i_pathsave)
            for i_mis in misspecification_LOD:
                for i_celltype in celltype:
                    miscorr[i_strategy][i_mis][i_celltype].append(i_mispopcorr[i_mis][i_strategy][i_celltype])

    color = ['r', 'b']
    lsty = ['-', '--']
    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    for i, i_strategy in enumerate(strategy_name):
        plt.rcParams['font.size'] = 16
        plt.figure()
        fig = plt.gcf()
        fig.set_size_inches(13, 8)
        legh, legstr = [], []
        for j, j_mis in enumerate(misspecification_LOD):
            for k, k_celltype in enumerate(celltype):
                j_l, = plt.plot(LOD, miscorr[i_strategy][j_mis][k_celltype], color=color[k], linestyle=lsty[j], marker='s',
                                linewidth=2, markersize=7)
                legh.append(j_l)
                legstr.append(f"set by {j_mis}: {k_celltype}")
        plt.title(i_strategy)
        plt.xlabel('LOD')
        plt.ylabel('pearson correlation coefficient')
        plt.legend(legh, legstr, loc='upper right', frameon=False)
        plt.tight_layout()
        i_pathsave = os.path.join(pathsave, i_strategy+'.pdf')
        plt.savefig(i_pathsave, format='pdf')
        plt.close()

    # para_df = pd.DataFrame.from_dict(para)
    # para_df['X'] = para_df.loc[:, ['Spop', 'R1pop', 'R2pop', 'R12pop']].sum(axis=1)
    # print('{:.5e}'.format(para_df['X'].min()))
    # print('{:.5e}'.format(para_df['X'].max()))
    # para_df['Spop per'] = para_df['Spop'] / para_df['X']
    # para_df['R1pop per'] = para_df['R1pop'] / para_df['X']
    # para_df['R2pop per'] = para_df['R2pop'] / para_df['X']
    # para_df['R12pop per'] = para_df['R12pop'] / para_df['X']
    # assert set(para_df['R12pop per']) == {0}
    # para_df['R1pop per'].min()
    # para_df['R1pop per'].max()
    # para_df['R2pop per'].min()
    # para_df['R2pop per'].max()
    return


def DPM_analysis_misspec_dosage(stoptime, dosage, LOD, filename_pattern, pathloadmis, pathsave):
    def DPM_analysis_misspec_dosage_1(val, dose_):
        if dose_ == '(1.0,0.0)':
            if np.isnan(val):
                val = 1
            else:
                val += 1
        elif dose_ == '(0.0,1.0)':
            if np.isnan(val):
                val = -1
            else:
                val -= 1
        elif dose_ == '(0.5,0.5)':
            if np.isnan(val):
                val = 0
            else:
                val += 0
        return val

    def DPM_analysis_misspec_dosage_2(stoptime_, stoptime_mis_, dosage_, dosage_mis_, num_, dosage_info_):
        first_ = NUM_STEP
        totaldiff_ = 0
        dosage_info_['num'][i_LOD][i_mis] += 1   # record number of misworse, miseven and misbetter instances
        dosage_info_['survial'][i_LOD][i_mis].append(abs(stoptime_ - stoptime_mis_))

        for m in range(len(dosage_)):
            dosage_info_['num_true_ateachstep'][i_LOD][i_mis][m] += 1
            dosage_info_['dose_true'][i_LOD][i_mis][m] = \
                DPM_analysis_misspec_dosage_1(dosage_info_['dose_true'][i_LOD][i_mis][m], dosage_[m])
            assert_ = abs(np.nanmax(dosage_info_['dose_true'][i_LOD][i_mis]/dosage_info_['num_true_ateachstep'][i_LOD][i_mis])) <= 1
            assert assert_

        for m in range(len(dosage_mis_)):
            # if treated by drug
            if dosage_mis_[m] != '(0.0,0.0)':
                dosage_info_['num_mis_ateachstep_withdrug'][i_LOD][i_mis][m] += 1
                dosage_info_['dose_mis'][i_LOD][i_mis][m] = \
                    DPM_analysis_misspec_dosage_1(dosage_info_['dose_mis'][i_LOD][i_mis][m], dosage_mis_[m])
                assert_ = abs(np.nanmax(dosage_info_['dose_mis'][i_LOD][i_mis] /
                                        dosage_info_['num_mis_ateachstep_withdrug'][i_LOD][i_mis])) <= 1
                assert assert_

        for m in range(0, num_):
            assert dosage_mis_[m] in set_dose and dosage_[m] in set_dose
            dosage_info_['num_mis_ateachstep'][i_LOD][i_mis][m] += 1
            if dosage_mis_[m] != dosage_[m]:
                totaldiff_ += 1
                if m < first_:
                    first_ = m
            else:
                ''' if same, dosage_mis_[m] != '(0.0,0.0)' '''
                # dosage_info_['num_samedose'][i_LOD][i_mis][m] += 1
                pass
        for m in range(0, num_):
            assert np.all(np.diff(dosage_info_['num_samedose'][i_LOD][i_mis]) <= 0)
            if dosage_mis_[m] == dosage_[m]:
                dosage_info_['num_samedose'][i_LOD][i_mis][m] += 1
            else:
                break
        percent_samedose_ = (dosage_info_['num_samedose'][i_LOD][i_mis]/dosage_info_['num_mis_ateachstep'][i_LOD][i_mis] * 100)

        dosage_info_['firstdiff'][i_LOD][i_mis].append(first_)
        dosage_info_['numdiff'][i_LOD][i_mis].append(totaldiff_)

        return dosage_info_

    def DPM_analysis_misspec_dosage_3(df_, name_):
        df_ = df_.stack().reset_index()
        df_.columns = ['c1', 'c2', 'values']
        df_['DF'] = name_
        return df_

    pathsave = os.path.join(pathsave, 'dosage')
    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    filename = os.path.join(pathsave, 'dose_info.pckl')
    strategy_name = list(stoptime.keys())
    strategy_name.remove('paramID')
    misspecification_LOD = MISSPECIFICATION_LOD_STR
    if not os.path.exists(filename):
        set_dose = {'(1.0,0.0)', '(0.0,1.0)', '(0.5,0.5)', '(0.0,0.0)'}
        result_list = {i_LOD: {i_mis: [] for i_mis in misspecification_LOD} for i_LOD in LOD}
        result_zero = {i_LOD: {i_mis: 0 for i_mis in misspecification_LOD} for i_LOD in LOD}
        result_array = {i_LOD: {i_mis: np.zeros((NUM_STEP,)) for i_mis in misspecification_LOD} for i_LOD in LOD}
        result_array_nan = np.empty((NUM_STEP,))
        result_array_nan[:] = np.nan
        result_array_nan = {i_LOD: {i_mis: deepcopy(result_array_nan) for i_mis in misspecification_LOD} for i_LOD in LOD}

        result = {'firstdiff': deepcopy(result_list),                             # 1
                  'numdiff': deepcopy(result_list),                               # 2
                  'survial': deepcopy(result_list),                               # 3
                  'num': deepcopy(result_zero),                                   # 4
                  'num_samedose': deepcopy(result_array),                         # 5
                  'dose_mis': deepcopy(result_array_nan),                         # 6
                  'num_mis_ateachstep_withdrug': deepcopy(result_array),          # 7
                  'num_mis_ateachstep': deepcopy(result_array),                   # 8
                  'dose_true': deepcopy(result_array_nan),                        # 9
                  'num_true_ateachstep': deepcopy(result_array)}                  # 10
        dosage_info = dict(zip(strategy_name, [None] * len(strategy_name)))
        for i_strategy in strategy_name:
            i_dosage_info = {'total': deepcopy(result),
                             'misworse': deepcopy(result),
                             'misbetter': deepcopy(result),
                             'miseven': deepcopy(result),
                             'misworsesig': deepcopy(result),
                             'misbettersig': deepcopy(result)}
            for i_LOD in LOD:
                i_pathloadmis = os.path.join(pathloadmis, i_LOD, filename_pattern + '_' + i_LOD)
                i_filename_stoptime_mis = i_pathloadmis + '_stoptime.pckl'
                i_filename_dosage_mis = i_pathloadmis + '_dosage.pckl'
                with bz2.BZ2File(i_filename_stoptime_mis, 'rb') as f_:
                    i_stoptime_mis = pickle.load(f_)
                with bz2.BZ2File(i_filename_dosage_mis, 'rb') as f_:
                    i_dosage_mis = pickle.load(f_)
                for i_mis in misspecification_LOD:
                    j_stoptime_mis = i_stoptime_mis[i_mis]
                    j_dosage_mis = i_dosage_mis[i_mis]
                    with tqdm(total=len(stoptime['paramID']), ncols=150, desc=f'{i_LOD},{i_mis}') as pbar:
                        for k, i_paramID in enumerate(stoptime['paramID']):
                            assert j_stoptime_mis['paramID'][k] == i_paramID == j_dosage_mis['paramID'][k]
                            k_stoptime_mis = j_stoptime_mis[i_strategy][k]
                            k_stoptme = stoptime[i_strategy][k]

                            k_dosage = dosage[i_strategy][k].split(';')
                            k_dosage = list(filter(lambda x: x != '-1', k_dosage))

                            k_dosage_mis = j_dosage_mis[i_strategy][k].split(';')
                            k_dosage_mis = list(filter(lambda x: x != '-1', k_dosage_mis))
                            k_num = min([len(k_dosage), len(k_dosage_mis)])

                            k_key = 'total'
                            i_dosage_info[k_key] = \
                                DPM_analysis_misspec_dosage_2(k_stoptme, k_stoptime_mis, k_dosage, k_dosage_mis, k_num,
                                                              i_dosage_info[k_key])

                            k_flag_misworsesig = k_stoptme > k_stoptime_mis + 30 * 2 & k_stoptme > 1.25 * k_stoptime_mis
                            k_flag_misbetterig = k_stoptime_mis > k_stoptme + 30 * 2 & k_stoptime_mis > 1.25 * k_stoptme

                            if k_stoptme > k_stoptime_mis:
                                i_dosage_info['misworse'] = \
                                    DPM_analysis_misspec_dosage_2(k_stoptme, k_stoptime_mis, k_dosage, k_dosage_mis, k_num,
                                                                  i_dosage_info['misworse'])
                                if k_flag_misworsesig:
                                    i_dosage_info['misworsesig'] = \
                                        DPM_analysis_misspec_dosage_2(k_stoptme, k_stoptime_mis, k_dosage, k_dosage_mis, k_num,
                                                                      i_dosage_info['misworsesig'])
                            elif k_stoptme < k_stoptime_mis:
                                i_dosage_info['misbetter'] = \
                                    DPM_analysis_misspec_dosage_2(k_stoptme, k_stoptime_mis, k_dosage, k_dosage_mis, k_num,
                                                                  i_dosage_info['misbetter'])
                                if k_flag_misbetterig:
                                    i_dosage_info['misbettersig'] = \
                                        DPM_analysis_misspec_dosage_2(k_stoptme, k_stoptime_mis, k_dosage, k_dosage_mis, k_num,
                                                                      i_dosage_info['misbettersig'])
                            elif k_stoptme == k_stoptime_mis:
                                i_dosage_info['miseven'] = \
                                    DPM_analysis_misspec_dosage_2(k_stoptme, k_stoptime_mis, k_dosage, k_dosage_mis, k_num,
                                                                  i_dosage_info['miseven'])

                            pbar.update(1)
                        assert i_dosage_info['misworse']['num'][i_LOD][i_mis] + \
                               i_dosage_info['misbetter']['num'][i_LOD][i_mis] + \
                               i_dosage_info['miseven']['num'][i_LOD][i_mis] == len(stoptime['paramID'])
                        dosage_info[i_strategy] = i_dosage_info

        with bz2.BZ2File(filename, 'wb') as f:
            pickle.dump(dosage_info, f)
    else:
        with bz2.BZ2File(filename, 'rb') as f:
            dosage_info = pickle.load(f)

        pathsave_survivalhist = os.path.join(pathsave, 'survivalhist')
        if not os.path.exists(pathsave_survivalhist):
            os.makedirs(pathsave_survivalhist)
        pathsave_firstdiff = os.path.join(pathsave, 'firstdiff')
        if not os.path.exists(pathsave_firstdiff):
            os.makedirs(pathsave_firstdiff)
        pathsave_numdiff = os.path.join(pathsave, 'numdiff')
        if not os.path.exists(pathsave_numdiff):
            os.makedirs(pathsave_numdiff)
        pathsave_percent_samedose = os.path.join(pathsave, 'percentsamedose')
        if not os.path.exists(pathsave_percent_samedose):
            os.makedirs(pathsave_percent_samedose)
        pathsave_dose_mis = os.path.join(pathsave, 'dosemis')
        if not os.path.exists(pathsave_dose_mis):
            os.makedirs(pathsave_dose_mis)
        pathsave_dose_true = os.path.join(pathsave, 'dosetrue')
        if not os.path.exists(pathsave_dose_true):
            os.makedirs(pathsave_dose_true)

        num_total = len(stoptime['paramID'])
        binsize_survialhist = 25
        xtick_survivalhist = 600
        bins_survivalhist = np.arange(-SIMDURATION_DEFAULT_VAL, SIMDURATION_DEFAULT_VAL+binsize_survialhist, binsize_survialhist)

        binsize_firstdiff = 4
        bins_firstdiffhist = np.arange(0, NUM_STEP + binsize_firstdiff, binsize_firstdiff)

        maxstep_percent_samedose = 10
        vmin_percent_samedose = 60

        misresult_type = ['misworse', 'misbetter', 'miseven', 'misworsesig', 'misbettersig', 'total']  #
        rowname_percent_samedose = [','.join(pair) for pair in itertools.product(LOD, misspecification_LOD)]

        num_LOD = len(LOD)
        num_misresult_type = len(misresult_type)
        num_mis = len(misspecification_LOD)
        for i_strategy in strategy_name:
            i_dosage_info = dosage_info[i_strategy]
            i_df_misworse = pd.DataFrame(np.zeros((num_LOD, num_mis)), index=LOD, columns=misspecification_LOD)
            i_df_misworsesig = pd.DataFrame(np.zeros((num_LOD, num_mis)), index=LOD, columns=misspecification_LOD)
            i_df_misbetter = pd.DataFrame(np.zeros((num_LOD, num_mis)), index=LOD, columns=misspecification_LOD)
            i_df_misbettersig = pd.DataFrame(np.zeros((num_LOD, num_mis)), index=LOD, columns=misspecification_LOD)
            i_df_miseven = pd.DataFrame(np.zeros((num_LOD, num_mis)), index=LOD, columns=misspecification_LOD)

            dose_mis = dict(zip(misresult_type, [np.zeros((num_mis*num_LOD, NUM_STEP)) for _ in range(num_misresult_type)]))
            dose_true = dict(zip(misresult_type, [np.zeros((num_mis*num_LOD, NUM_STEP)) for _ in range(num_misresult_type)]))
            percent_samedose = dict(zip(misresult_type, [np.zeros((num_mis*num_LOD, NUM_STEP)) for _ in range(num_misresult_type)]))
            firstdiff = dict(zip(misresult_type, [np.zeros((num_mis*num_LOD, NUM_STEP)) for _ in range(num_misresult_type)]))
            numdiff = dict(zip(misresult_type, [np.zeros((num_mis*num_LOD, NUM_STEP)) for _ in range(num_misresult_type)]))

            for i, i_LOD in enumerate(LOD):
                i_firstdiff = dict(zip(misresult_type, [pd.DataFrame(columns=['val', 'type']) for _ in range(num_misresult_type)]))
                i_numdiff = dict(zip(misresult_type, [pd.DataFrame(columns=['val', 'type']) for _ in range(num_misresult_type)]))

                i_percent_samedose = dict(zip(misresult_type, [np.zeros((num_mis, NUM_STEP)) for _ in range(num_misresult_type)]))
                i_dose_mis = dict(zip(misresult_type, [np.zeros((num_mis, NUM_STEP)) for _ in range(num_misresult_type)]))
                i_dose_true = dict(zip(misresult_type, [np.zeros((num_mis, NUM_STEP)) for _ in range(num_misresult_type)]))

                for j, j_mis in enumerate(misspecification_LOD):
                    i_df_misworsesig.loc[i_LOD][j_mis] = i_dosage_info['misworsesig']['num'][i_LOD][j_mis]/num_total * 100
                    i_df_misworse.loc[i_LOD][j_mis] = i_dosage_info['misworse']['num'][i_LOD][j_mis]/num_total * 100
                    # i_df_misworse.loc[i_LOD][j_mis] = i_df_misworse.loc[i_LOD][j_mis] - i_df_misworsesig.loc[i_LOD][j_mis]

                    i_df_misbettersig.loc[i_LOD][j_mis] = i_dosage_info['misbettersig']['num'][i_LOD][j_mis]/num_total * 100
                    i_df_misbetter.loc[i_LOD][j_mis] = i_dosage_info['misbetter']['num'][i_LOD][j_mis]/num_total * 100
                    # i_df_misbetter.loc[i_LOD][j_mis] = i_df_misbetter.loc[i_LOD][j_mis] - i_df_misbettersig.loc[i_LOD][j_mis]

                    i_df_miseven.loc[i_LOD][j_mis] = i_dosage_info['miseven']['num'][i_LOD][j_mis]/num_total * 100

                    assert all(v == 0 for v in i_dosage_info['miseven']['survial'][i_LOD][j_mis])
                    i_survival_misworse = -np.array(i_dosage_info['misworse']['survial'][i_LOD][j_mis])
                    i_survival_misbetter = np.array(i_dosage_info['misbetter']['survial'][i_LOD][j_mis])

                    j_survival_misworse = pd.DataFrame(list(zip(i_survival_misworse, ['misworse'] * len(i_survival_misworse))),
                                                       columns=['val', 'type'])
                    j_survival_misbetter = pd.DataFrame(list(zip(i_survival_misbetter, ['misbetter'] * len(i_survival_misbetter))),
                                                        columns=['val', 'type'])
                    j_survival = pd.concat([j_survival_misworse, j_survival_misbetter])

                    for i_type in misresult_type:
                        j_firstdiff = i_dosage_info[i_type]['firstdiff'][i_LOD][j_mis]
                        j_firstdiff_hist, _ = np.histogram(j_firstdiff, bins=np.arange(NUM_STEP+1))
                        firstdiff[i_type][i * len(misspecification_LOD) + j, :] = j_firstdiff_hist/sum(j_firstdiff_hist) * 100
                        j_firstdiff = pd.DataFrame(list(zip(j_firstdiff, [j_mis] * len(j_firstdiff))), columns=['val', 'type'])
                        i_firstdiff[i_type] = pd.concat([i_firstdiff[i_type], j_firstdiff])

                        j_numdiff = i_dosage_info[i_type]['numdiff'][i_LOD][j_mis]
                        j_numdiff_hist, _ = np.histogram(j_numdiff, bins=np.arange(NUM_STEP+1))
                        numdiff[i_type][i * len(misspecification_LOD) + j, :] = j_numdiff_hist/sum(j_numdiff_hist) * 100
                        j_numdiff = pd.DataFrame(list(zip(j_numdiff, [j_mis] * len(j_numdiff))), columns=['val', 'type'])
                        i_numdiff[i_type] = pd.concat([i_numdiff[i_type], j_numdiff])

                        j_percent_samedose = (i_dosage_info[i_type]['num_samedose'][i_LOD][j_mis] /
                                              i_dosage_info[i_type]['num_mis_ateachstep'][i_LOD][j_mis] * 100)

                        assert np.all(np.diff(i_dosage_info[i_type]['num_samedose'][i_LOD][j_mis]) <= 0)
                        assert np.all(np.diff(i_dosage_info[i_type]['num_mis_ateachstep'][i_LOD][j_mis]) <= 0)

                        i_percent_samedose[i_type][j, :] = j_percent_samedose
                        percent_samedose[i_type][i*len(misspecification_LOD)+j, :] = j_percent_samedose

                        j_dose_mis = (i_dosage_info[i_type]['dose_mis'][i_LOD][j_mis] /
                                      i_dosage_info[i_type]['num_mis_ateachstep_withdrug'][i_LOD][j_mis])
                        i_dose_mis[i_type][j, :] = j_dose_mis
                        dose_mis[i_type][i*len(misspecification_LOD)+j, :] = j_dose_mis

                        j_dose_true = (i_dosage_info[i_type]['dose_true'][i_LOD][j_mis] /
                                       i_dosage_info[i_type]['num_true_ateachstep'][i_LOD][j_mis])
                        i_dose_true[i_type][j, :] = j_dose_true
                        dose_true[i_type][i*len(misspecification_LOD)+j, :] = j_dose_true

                    ''' hist plot of survial difference '''
                    ##
                    sns.histplot(data=j_survival, bins=bins_survivalhist, x='val', hue='type', palette=['#ff6f69', '#96ceb4'],
                                 stat='percent', linewidth=0, common_norm=False, legend=False)
                    plt.xlabel('survial difference')
                    plt.yscale('log')
                    plt.ylim([10**-2, 10**2])
                    plt.xlim([-SIMDURATION_DEFAULT_VAL-binsize_survialhist, SIMDURATION_DEFAULT_VAL+binsize_survialhist])
                    # plt.axvline(x=0, color='k', linestyle='-', linewidth=0.4)
                    plt.title(f'{i_strategy},{i_LOD},{j_mis}')
                    plt.xticks(np.arange(-SIMDURATION_DEFAULT_VAL, SIMDURATION_DEFAULT_VAL+xtick_survivalhist, xtick_survivalhist))
                    plt.yticks(10**np.arange(-2, 3, 1).astype(float))
                    plt.gcf().set_size_inches((5.5, 2.5))
                    plt.tight_layout()
                    plt.savefig(os.path.join(pathsave_survivalhist, f'{i_strategy},{i_LOD},{j_mis}' + '.pdf'), dpi=FIG_DPI)
                    plt.close()
                    ##
                for i_type in misresult_type:
                    '''hist plot of step index of first difference'''
                    i_pathsave_firstdiff = os.path.join(pathsave_firstdiff, i_type)
                    if not os.path.exists(i_pathsave_firstdiff):
                        os.makedirs(i_pathsave_firstdiff)
                    sns.histplot(data=i_firstdiff[i_type], bins=bins_firstdiffhist, x='val', hue='type', multiple='dodge',
                                 stat='percent', linewidth=0, legend=True, common_norm=False, shrink=0.8)
                    plt.xticks(bins_firstdiffhist)
                    plt.yticks(np.arange(0, 110, 10))
                    plt.xlabel('Step index')
                    plt.title(f'{i_strategy},{i_LOD},{i_type}')
                    plt.gcf().set_size_inches((10, 4))
                    plt.tight_layout()
                    plt.savefig(os.path.join(i_pathsave_firstdiff, f'{i_strategy},{i_LOD}' + '.pdf'), dpi=FIG_DPI)
                    plt.close()

                    '''hist plot of total number of steps under different drug treatments'''
                    i_pathsave_numdiff = os.path.join(pathsave_numdiff, i_type)
                    if not os.path.exists(i_pathsave_numdiff):
                        os.makedirs(i_pathsave_numdiff)
                    sns.histplot(data=i_numdiff[i_type], bins=bins_firstdiffhist, x='val', hue='type', multiple='dodge',
                                 stat='percent', linewidth=0, legend=True, common_norm=False, shrink=0.8)
                    plt.xticks(bins_firstdiffhist)
                    plt.yticks(np.arange(0, 110, 10))
                    plt.xlabel('Number of steps')
                    plt.title(f'{i_strategy},{i_LOD},{i_type}')
                    plt.gcf().set_size_inches((10, 4))
                    plt.tight_layout()
                    plt.savefig(os.path.join(i_pathsave_numdiff, f'{i_strategy},{i_LOD}' + '.pdf'), dpi=FIG_DPI)
                    plt.close()

                    '''Percentage of identical doses at each step'''
                    i_pathsave_percent_samedose = os.path.join(pathsave_percent_samedose, i_type)
                    if not os.path.exists(i_pathsave_percent_samedose):
                        os.makedirs(i_pathsave_percent_samedose)
                    sns.heatmap(i_percent_samedose[i_type][:, :maxstep_percent_samedose], annot=True, fmt='.1f', linewidth=.5,
                                vmin=vmin_percent_samedose, vmax=100, xticklabels=np.arange(1, maxstep_percent_samedose+1),
                                yticklabels=misspecification_LOD, cmap='rocket')
                    plt.xlabel(f'Number of steps,{i_strategy},{i_LOD},{i_type}')
                    plt.gcf().set_size_inches((8, 3))
                    plt.tight_layout()
                    plt.show()
                    plt.savefig(os.path.join(i_pathsave_percent_samedose, f'{i_strategy},{i_LOD}' + '.pdf'), dpi=FIG_DPI)
                    plt.close()

                    '''Misspecified drug doses at each step'''
                    i_pathsave_dose_mis = os.path.join(pathsave_dose_mis, i_type)
                    if not os.path.exists(i_pathsave_dose_mis):
                        os.makedirs(i_pathsave_dose_mis)
                    sns.heatmap(i_dose_mis[i_type], annot=True, fmt='.1f', linewidth=.5,
                                vmin=-1, vmax=1, xticklabels=np.arange(1, NUM_STEP+1),
                                yticklabels=misspecification_LOD, cmap='Spectral')
                    plt.xlabel(f'Number of steps,{i_strategy},{i_LOD},{i_type}')
                    plt.gcf().set_size_inches((18, 3))
                    plt.tight_layout()
                    plt.savefig(os.path.join(i_pathsave_dose_mis, f'{i_strategy},{i_LOD}' + '.pdf'), dpi=FIG_DPI)
                    plt.close()

                    '''True drug doses at each step'''
                    i_pathsave_dose_true = os.path.join(pathsave_dose_true, i_type)
                    if not os.path.exists(i_pathsave_dose_true):
                        os.makedirs(i_pathsave_dose_true)
                    sns.heatmap(i_dose_true[i_type], annot=True, fmt='.1f', linewidth=.5,
                                vmin=-1, vmax=1, xticklabels=np.arange(1, NUM_STEP+1),
                                yticklabels=misspecification_LOD, cmap='Spectral')
                    plt.xlabel(f'Number of steps,{i_strategy},{i_LOD},{i_type}')
                    plt.gcf().set_size_inches((18, 3))
                    plt.tight_layout()
                    plt.savefig(os.path.join(i_pathsave_dose_true, f'{i_strategy},{i_LOD}' + '.pdf'), dpi=FIG_DPI)
                    plt.close()

            rowname = [','.join([j_s, i_s]) for i_s in misspecification_LOD for j_s in list(reversed(LOD))]
            order = [i for item2 in rowname for i, item1 in enumerate(rowname_percent_samedose) if item1 in item2]
            for i_type in misresult_type:
                '''first difference'''
                i_pathsave_firstdiff = os.path.join(pathsave_firstdiff, i_type)
                if not os.path.exists(i_pathsave_firstdiff):
                    os.makedirs(i_pathsave_firstdiff)
                sns.heatmap(firstdiff[i_type][order, :], annot=False, fmt='.1f', linewidth=.5,
                            vmin=np.nanmin((firstdiff[i_type])),
                            vmax=np.nanmax((firstdiff[i_type])),
                            xticklabels=np.arange(1, NUM_STEP+1),
                            yticklabels=rowname, cmap='Spectral')
                plt.xlabel(f'Number of steps,{i_strategy}, {i_type}')
                plt.gcf().set_size_inches((15, 7))
                plt.tight_layout()
                plt.savefig(os.path.join(i_pathsave_firstdiff, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
                plt.close()

                '''total number of steps under different drug treatments'''
                i_pathsave_numdiff = os.path.join(pathsave_numdiff, i_type)
                if not os.path.exists(i_pathsave_numdiff):
                    os.makedirs(i_pathsave_numdiff)
                sns.heatmap(numdiff[i_type][order, :], annot=False, fmt='.1f', linewidth=.5,
                            vmin=np.nanmin((numdiff[i_type])),
                            vmax=np.nanmax((numdiff[i_type])),
                            xticklabels=np.arange(1, NUM_STEP+1),
                            yticklabels=rowname, cmap='Spectral')
                plt.xlabel(f'Number of steps,{i_strategy}, {i_type}')
                plt.gcf().set_size_inches((15, 7))
                plt.tight_layout()
                plt.savefig(os.path.join(i_pathsave_numdiff, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
                plt.close()

                '''Percentage of identical doses at each step'''
                i_pathsave_percent_samedose = os.path.join(pathsave_percent_samedose, i_type)
                if not os.path.exists(i_pathsave_percent_samedose):
                    os.makedirs(i_pathsave_percent_samedose)
                i_percent_samedose = percent_samedose[i_type][order, :maxstep_percent_samedose]
                sns.heatmap(i_percent_samedose, annot=True, fmt='.1f', linewidth=.5,
                            vmin=np.nanmin(i_percent_samedose),
                            vmax=np.nanmax(i_percent_samedose),
                            xticklabels=np.arange(1, maxstep_percent_samedose+1),
                            yticklabels=rowname, cmap='Spectral')
                plt.xlabel(f'Number of steps,{i_strategy}, {i_type}')
                plt.gcf().set_size_inches((12, 8))
                plt.tight_layout()
                plt.savefig(os.path.join(i_pathsave_percent_samedose, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
                plt.close()

                vmin = min([np.round(np.nanmin((dose_mis[i_type])), 1), np.round(np.nanmin((dose_true[i_type])), 1)])
                vmax = max([np.round(np.nanmax((dose_mis[i_type])), 1), np.round(np.nanmax((dose_true[i_type])), 1)])
                ''''Misspecified drug doses at each step'''
                i_pathsave_dose_mis = os.path.join(pathsave_dose_mis, i_type)
                if not os.path.exists(i_pathsave_dose_mis):
                    os.makedirs(i_pathsave_dose_mis)
                sns.heatmap(dose_mis[i_type][order, :], annot=False, fmt='.1f', linewidth=.5, vmin=vmin, vmax=vmax,
                            xticklabels=np.arange(1, NUM_STEP+1),
                            yticklabels=rowname, cmap='Spectral')
                plt.xlabel(f'Number of steps,{i_strategy}, {i_type}')
                plt.gcf().set_size_inches((13, 8))
                plt.tight_layout()
                plt.savefig(os.path.join(i_pathsave_dose_mis, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
                plt.close()

                '''True drug doses at each step'''
                i_pathsave_dose_true = os.path.join(pathsave_dose_true, i_type)
                if not os.path.exists(i_pathsave_dose_true):
                    os.makedirs(i_pathsave_dose_true)
                sns.heatmap(dose_true[i_type][order, :], annot=False, fmt='.1f', linewidth=.5, vmin=vmin, vmax=vmax,
                            xticklabels=np.arange(1, NUM_STEP+1),
                            yticklabels=rowname, cmap='Spectral')
                plt.xlabel(f'Number of steps,{i_strategy},{i_type}')
                plt.gcf().set_size_inches((13, 8))
                plt.tight_layout()
                plt.savefig(os.path.join(i_pathsave_dose_true, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
                plt.close()
            ''' total '''
            '''first difference'''
            sns.heatmap(firstdiff['total'][order, :], annot=False, fmt='.1f', linewidth=.5,
                        vmin=np.nanmin((firstdiff['total'])),
                        vmax=np.nanmax((firstdiff['total'])),
                        xticklabels=np.arange(1, NUM_STEP + 1),
                        yticklabels=rowname, cmap='Spectral')
            plt.xlabel(f'Number of steps,{i_strategy}')
            plt.gcf().set_size_inches((15, 7))
            plt.tight_layout()
            plt.savefig(os.path.join(pathsave_firstdiff, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
            plt.close()

            '''total number of steps under different drug treatments'''
            sns.heatmap(numdiff['total'][order, :], annot=False, fmt='.1f', linewidth=.5,
                        vmin=np.nanmin((numdiff['total'])),
                        vmax=np.nanmax((numdiff['total'])),
                        xticklabels=np.arange(1, NUM_STEP + 1),
                        yticklabels=rowname, cmap='Spectral')
            plt.xlabel(f'Number of steps,{i_strategy}')
            plt.gcf().set_size_inches((15, 7))
            plt.tight_layout()
            plt.savefig(os.path.join(pathsave_numdiff, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
            plt.close()

            '''Percentage of identical doses at each step'''
            i_percent_samedose = percent_samedose['total'][order, :maxstep_percent_samedose]
            sns.heatmap(i_percent_samedose, annot=True, fmt='.1f', linewidth=.5,
                        vmin=np.nanmin(i_percent_samedose),
                        vmax=np.nanmax(i_percent_samedose),
                        xticklabels=np.arange(1, maxstep_percent_samedose + 1),
                        yticklabels=rowname, cmap='Spectral')
            plt.xlabel(f'Number of steps,{i_strategy}')
            plt.gcf().set_size_inches((12, 8))
            plt.tight_layout()
            plt.savefig(os.path.join(pathsave_percent_samedose, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
            plt.close()

            '''True drug doses at each step'''
            assert np.all(dose_true['total'] == dose_true['total'][0], axis=1).all()

            ''''True and Misspecified drug doses at each step'''
            dose = np.insert(dose_mis['total'][order, :], 0, dose_true['total'][0, :], axis=0)
            sns.heatmap(dose, annot=False, fmt='.1f', linewidth=.5,
                        vmin=np.min(dose),
                        vmax=np.max(dose),
                        xticklabels=np.arange(1, NUM_STEP + 1),
                        yticklabels=['True']+rowname, cmap='Spectral')
            plt.xlabel(f'Number of steps,{i_strategy}')
            plt.gcf().set_size_inches((17, 5))
            plt.tight_layout()
            plt.savefig(os.path.join(pathsave_dose_mis, f'{i_strategy}' + '.pdf'), dpi=FIG_DPI)
            plt.close()

            # stacked bar plot
            alt.renderers.enable('browser')
            i_df_misworse_ = DPM_analysis_misspec_dosage_3(i_df_misworse, 'worse')
            # i_df_misworsesig_ = DPM_analysis_misspec_dosage_3(i_df_misworsesig, 'worsesig')
            i_df_misbetter_ = DPM_analysis_misspec_dosage_3(i_df_misbetter, 'better')
            # i_df_misbettersig_ = DPM_analysis_misspec_dosage_3(i_df_misbettersig, 'bettersig')
            i_df_miseven_ = DPM_analysis_misspec_dosage_3(i_df_miseven, 'even')
            i_df = pd.concat([i_df_misworse_, i_df_misbetter_, i_df_miseven_])
            Chart = (alt.Chart(i_df).mark_bar().encode(x=alt.X('c2:N', title=None, sort=['0', 'max', 'loguni', 'pdf']),
                                                       y=alt.Y('sum(values):Q', axis=alt.Axis(grid=False, title=None)),
                                                       column=alt.Column('c1:N', title=None, sort='descending'),
                                                       color=alt.Color('DF:N', scale=alt.Scale(
                                                           range=['#96ceb4', '#ffcc5c', '#ff8c69']))).
                     configure_view(strokeOpacity=0))
            Chart.show()
    return


def DPM_analysis_misspec_example(para, stoptime, dosage, pop, filename_pattern, pathloadmis, pathsave, LOD='1e-03', mis='0'):
    i_pathloadmis = os.path.join(pathloadmis, LOD, filename_pattern + '_' + LOD)
    i_filename_stoptime_mis = i_pathloadmis + '_stoptime.pckl'
    i_filename_dosage_mis = i_pathloadmis + '_dosage.pckl'
    with bz2.BZ2File(i_filename_stoptime_mis, 'rb') as f_:
        stoptime_mis = pickle.load(f_)
        stoptime_mis = stoptime_mis[mis]
    with bz2.BZ2File(i_filename_dosage_mis, 'rb') as f_:
        dosage_mis = pickle.load(f_)[mis]

    Strategy_name = list(stoptime.keys())
    Strategy_name.remove('paramID')
    # ind_sigbetter = [i for i, val in enumerate(stoptime['paramID']) if val in paramID_sigbetter]

    arr1 = np.array(stoptime[Strategy_name[0]]) > STEPSIZE_DEFAULT_VAL * 3
    arr2 = np.array(stoptime[Strategy_name[0]]) < STEPSIZE_DEFAULT_VAL * 21
    arr3 = np.array(stoptime[Strategy_name[1]]) >= SIMDURATION_DEFAULT_VAL
    arr4 = np.array(stoptime[Strategy_name[1]]) == SIMDURATION_DEFAULT_VAL + STEPSIZE_DEFAULT_VAL
    arr5 = np.array(stoptime_mis[Strategy_name[1]]) < STEPSIZE_DEFAULT_VAL * 31
    arr6 = np.array(stoptime_mis[Strategy_name[1]]) > np.array(stoptime[Strategy_name[0]]) + STEPSIZE_DEFAULT_VAL * 2

    ind = functools.reduce(np.logical_and, [arr1, arr2, arr3, arr4, arr5, arr6])
    paramID_sel = list(itertools.compress(para['paramID'], ind))

    ind_sel = []
    popt0_sel = []
    stoptime_sel = []
    dosage_strategy0_sel = []
    dosage_strategy2_sel = []
    dosage_strategy2_mis_sel = []
    diff_strategy2_mis_strategy0_sel = []
    with tqdm(total=len(paramID_sel), ncols=100, desc=' ') as pbar:
        for i_paramID in paramID_sel:
            i = stoptime['paramID'].index(i_paramID)
            assert pop['paramID'][i] == dosage['paramID'][i] == stoptime_mis['paramID'][i] == dosage_mis['paramID'][i] == i_paramID

            i_dosage_strategy0 = dosage[Strategy_name[0]][i].split(';')[:2]
            flag1 = set(i_dosage_strategy0) == {'(1.0,0.0)'}
            i_dosage_strategy2 = dosage[Strategy_name[1]][i].split(';')[:1]
            flag2 = set(i_dosage_strategy2) == {'(0.0,1.0)'}
            i_dosage_strategy2_mis = dosage_mis[Strategy_name[1]][i].split(';')[:1]
            flag3 = set(i_dosage_strategy2_mis) == {'(1.0,0.0)'}

            i_pop_t0 = np.array([para['Spop'][i], para['R1pop'][i], para['R2pop'][i], para['R12pop'][i]]).reshape(-1, 1)
            i_pop = pop[Strategy_name[0]][i].split(';')
            i_pop = list(filter(lambda x: x != '-1', i_pop))
            i_pop_np = np.zeros((len(ALL_POSSIBLE_CELLTYPE_2DRUG), len(i_pop)))
            for j in range(len(i_pop)):
                i_pop_np[:, j] = np.array(re.sub('[()]', '', i_pop[j]).split(','), dtype=float)
            i_pop_np = np.concatenate((i_pop_t0, i_pop_np), axis=1)
            flag4 = i_pop_np[:, -1].argmax() == i_pop_np.shape[0] - 1
            flag5 = max(i_pop_np[1:, 0]) / sum(i_pop_t0) < 0.01

            if flag1 and flag2 and flag3 and flag4 and flag5:
                ind_sel.append(i)
                popt0_sel.append(tuple(i_pop_np[:, 0]/sum(i_pop_t0)))
                stoptime_sel.append((stoptime[Strategy_name[0]][i], stoptime[Strategy_name[1]][i],
                                     stoptime_mis[Strategy_name[0]][i], stoptime_mis[Strategy_name[1]][i]))
                dosage_strategy0_sel.append(i_dosage_strategy0)
                dosage_strategy2_sel.append(i_dosage_strategy2)
                dosage_strategy2_mis_sel.append(i_dosage_strategy2_mis)
                diff_strategy2_mis_strategy0_sel.append(stoptime_mis[Strategy_name[1]][i] - stoptime[Strategy_name[0]][i])
            pbar.update(1)
    ind_max = np.array(diff_strategy2_mis_strategy0_sel).argmax()
    ind_sel = ind_sel[ind_max]
    popt0_sel = popt0_sel[ind_max]
    para_sel = {key: value[ind_sel] for key, value in para.items()}

    para_sel = DPM_miscellaneous_fillful(para_sel)

    pathsave = os.path.join(pathsave, 'example')
    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    with bz2.BZ2File(os.path.join(pathsave, 'para_sel.pckl'), 'wb') as f:
        pickle.dump(para_sel, f)
    return


def DPM_analysis_misspec(km, info, km_mis, info_mis, LOD, setname, Strategy_name, Simduration, pathsave, plot):
    km_mis = dict(zip(MISSPECIFICATION_LOD_STR, km_mis))
    info_mis = dict(zip(MISSPECIFICATION_LOD_STR, info_mis))

    LOD_plot = LOD_LIST
    LOD_all = LOD + ['nomis']

    keys = []
    for i in MISSPECIFICATION_LOD_STR:
        keys.extend([i_LOD + ' ' + i for i_LOD in LOD if i_LOD in LOD_plot])
    keys.extend([LOD_all[-1]])

    for i, i_LOD in enumerate(LOD):
        for i_setname in setname:
            if plot:
                path_folder = os.path.join(pathsave, i_setname)
                if not os.path.exists(path_folder):
                    os.makedirs(path_folder)
                # True
                filename = os.path.join(path_folder, 'nomis.pdf')
                if not os.path.exists(filename):
                    titlestr = i_setname + ', ' + LOD_all[-1]
                    DPM_analysis_misspec_plot(Strategy_name, Simduration, km[i_setname], titlestr)
                    # plt.savefig(filename, dpi=FIG_DPI)
                    plt.close('all')
                for i_mis in MISSPECIFICATION_LOD_STR:
                    if i_LOD in LOD_plot:
                        filename = os.path.join(path_folder, 'LOD ' + i_LOD + ' ' + i_mis + '.pdf')
                        titlestr = i_setname + ', set by ' + i_mis + ', LOD: ' + i_LOD
                        DPM_analysis_misspec_plot(Strategy_name, Simduration, km_mis[i_mis][i][i_setname], titlestr)
                        # plt.savefig(filename, dpi=FIG_DPI)
                        plt.close('all')

    for i_setname in setname:
        result = {i_strategy: [] for i_strategy in Strategy_name}
        median_survial, hazard, firstd2, movenum_d2, numdchange = [], [], [], [], []
        for _ in range(len(MISSPECIFICATION_LOD_STR)):
            median_survial.append(deepcopy(result))
            hazard.append(deepcopy(result))
            firstd2.append(deepcopy(result))
            movenum_d2.append(deepcopy(result))
            numdchange.append(deepcopy(result))

        km_set_i = km[i_setname]
        info_i_set = info[i_setname]

        km_mis_set_i, hz_ratio_set_i, km_mis_set_i_num, info_mis_set_i = [], [], [], []
        for i, (key_i, km_mis_i) in enumerate(km_mis.items()):
            km_mis_i_set_i = [i_val[i_setname] for i_val in km_mis_i]
            km_mis_set_i.append(km_mis_i_set_i)

            info_mis_i = info_mis[key_i]
            info_mis_i_set_i = [i_val[i_setname] for i_val in info_mis_i]
            info_mis_set_i.append(info_mis_i_set_i)

            hz_ratio_mis_i_set_i = [i_km_mis_i_set_i['hz_ratio'] for i_km_mis_i_set_i in km_mis_i_set_i]
            hz_ratio_mis_i_set_i.append(km_set_i['hz_ratio'])
            hz_ratio_set_i.append(hz_ratio_mis_i_set_i)

            km_mis_i_set_i_num = [i_km_mis_i_set_i['num'] for i_km_mis_i_set_i in km_mis_i_set_i]
            km_mis_set_i_num.append(km_mis_i_set_i_num)

        path_folder = os.path.join(pathsave, i_setname)
        if not os.path.exists(path_folder):
            os.makedirs(path_folder)

        sel_LOD = '1e-02'
        sel_km = {i_key: None for i_key in Strategy_name}
        for i_strategy in Strategy_name:
            km_i_set_strategy, km_i_set_strategy_plot, info_i_set_strategy = [], [], []
            for i in range(len(MISSPECIFICATION_LOD_STR)):
                i_km_i_set_strategy = [i_val[i_strategy] for i_val in km_mis_set_i[i]]
                km_i_set_strategy.append(i_km_i_set_strategy)
                i_km_i_set_strategy_plot = [i_km_i_set_strategy[ii] for ii in range(len(i_km_i_set_strategy)) if LOD[i] in LOD_plot]
                km_i_set_strategy_plot.append(i_km_i_set_strategy_plot)
                info_i_set_strategy.append([i_info_0_i_set[i_strategy] for i_info_0_i_set in info_mis_set_i[i]])

            ms_i_set_strategy, hazard_i_set_strategy, firstd2_i_set_strategy, movenum_d2_i_set_strategy, \
                ave_numdrugchange_i_set_strategy = [], [], [], [], []

            for i in range(len(MISSPECIFICATION_LOD_STR)):
                i_ms_i_set_strategy = [i_km_i_set_strategy['median_survival'] for i_km_i_set_strategy in km_i_set_strategy[i]]
                i_ms_i_set_strategy.append(km_set_i[i_strategy]['median_survival'])
                ms_i_set_strategy.append(i_ms_i_set_strategy)
                median_survial[i][i_strategy] = i_ms_i_set_strategy

                i_hazard_i_set_strategy = [i_km_i_set_strategy['hazard'] for i_km_i_set_strategy in km_i_set_strategy[i]]
                i_hazard_i_set_strategy.append(km_set_i[i_strategy]['hazard'])
                hazard_i_set_strategy.append(i_hazard_i_set_strategy)
                hazard[i][i_strategy] = i_hazard_i_set_strategy

                i_firstd2_i_set_strategy = [i_info_i_set_strategy['first drug 2'] for i_info_i_set_strategy in info_i_set_strategy[i]]
                i_firstd2_i_set_strategy.append(info_i_set[i_strategy]['first drug 2'])
                firstd2_i_set_strategy.append(i_firstd2_i_set_strategy)
                firstd2[i][i_strategy] = i_firstd2_i_set_strategy

                i_movenum_d2_i_set_strategy = [i_info_i_set_strategy['average move number'] for i_info_i_set_strategy in info_i_set_strategy[i]]
                i_movenum_d2_i_set_strategy.append(info_i_set[i_strategy]['average move number'])
                movenum_d2_i_set_strategy.append(i_movenum_d2_i_set_strategy)
                movenum_d2[i][i_strategy] = i_movenum_d2_i_set_strategy

                i_ave_numdrugchange_i_set_strategy = [i_info_i_set_strategy['drug changes'] for i_info_i_set_strategy in info_i_set_strategy[i]]
                i_ave_numdrugchange_i_set_strategy.append(info_i_set[i_strategy]['drug changes'])
                ave_numdrugchange_i_set_strategy.append(i_ave_numdrugchange_i_set_strategy)
                numdchange[i][i_strategy] = i_ave_numdrugchange_i_set_strategy

            i_km = {i_key: None for i_key in keys}
            for i, i_LOD in enumerate(LOD):
                if i_LOD in LOD_plot:
                    for ii, i_mis in enumerate(MISSPECIFICATION_LOD_STR):
                        i_key = i_LOD + ' ' + i_mis
                        i_km[i_key] = {**km_i_set_strategy_plot[ii][i], **{'num': km_mis_set_i_num[ii][i]}}

            i_km[keys[-1]] = {**km_set_i[i_strategy], 'num': km_set_i['num']}

            titlestr = i_setname + ', ' + i_strategy
            par = {'duration': Simduration, 'xtick step': 300}
            DPM_plot_KM_multi2(i_km, par, titlestr)
            # plt.savefig(os.path.join(path_folder, i_strategy + '.pdf'), dpi=FIG_DPI)
            plt.close('all')

            key_select = list(i_km.keys())
            key_select = [i_key for i_key in key_select if sel_LOD in i_key]
            sel_km[i_strategy] = {key: i_km[key] for key in key_select}

        par = {'duration': Simduration, 'xtick step': 300}
        DPM_plot_KM_multi_sel(sel_km, par, sel_LOD, i_setname+f',{sel_LOD}')
        # plt.savefig(os.path.join(path_folder, f'km_sel.pdf'), dpi=FIG_DPI)
        plt.close('all')

        titlestr = i_setname
        DPM_plot_LOD_multi(median_survial, LOD_all, titlestr)
        # plt.savefig(os.path.join(path_folder, 'median survival.pdf'), dpi=FIG_DPI)
        plt.close('all')

        DPM_plot_hz_ratio(hz_ratio_set_i, LOD_all, titlestr)
        # plt.savefig(os.path.join(path_folder, 'hz ratio.pdf'), dpi=FIG_DPI)
        plt.close('all')

        DPM_plot_hz(hazard, LOD_all, titlestr)
        # plt.savefig(os.path.join(path_folder, 'hz.pdf'), dpi=FIG_DPI)
        plt.close('all')

        name = 'Average drug 2 introduction timestep'
        DPM_plot_LOD_multi(firstd2, LOD_all, titlestr, ylablestr=name, ylim=(2.6, 5.2))
        plt.savefig(os.path.join(path_folder, name+'.pdf'), dpi=FIG_DPI)
        plt.close('all')

        name = 'Average drug 2 intensity'
        DPM_plot_LOD_multi(numdchange, LOD_all, titlestr, ylablestr=name, ylim=(0.14, 0.26))
        plt.savefig(os.path.join(path_folder, name+'.pdf'), dpi=FIG_DPI)
        plt.close('all')

        name = 'Drug 2 intensity weighted average'
        DPM_plot_LOD_multi(movenum_d2, LOD_all, titlestr, ylablestr=name, ylim=(5.2, 7.8))
        plt.savefig(os.path.join(path_folder, name+'.pdf'), dpi=FIG_DPI)
        plt.close('all')
    return


def DPM_analysis_misspec_par(para, info, info_mis, LOD, pathsave, pathload):
    def DPM_analysis_misspec_par_1(i_value, i_paraval):
        num = [i_value.count(ii) for ii in i_paraval]
        assert sum(num) == len(i_value)
        return [i_num / len(i_value) * 100 for i_num in num]

    def DPM_analysis_misspec_par_2(name):
        bins = [0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.9, 1]
        ind = para_df[name] == 0
        ind = np.flatnonzero(ind)
        para_df.loc[ind, name] = 0
        for ii in range(len(bins)-1):
            ind = (bins[ii] < para_df[name]) & (para_df[name] <= bins[ii+1])
            ind = np.flatnonzero(ind)
            para_df.loc[ind, name] = bins[ii+1]

    paramID = para['paramID']
    # paramID_first2diff = [i for i, _ in enumerate(para['paramID']) if setind['first 2 diff'][i]]

    path_folder = os.path.join(pathsave, 'set')
    if not os.path.exists(path_folder):
        os.makedirs(path_folder)
    if not os.path.exists(os.path.join(pathsave, 'par')):
        os.makedirs(os.path.join(pathsave, 'par'))
    try:
        with bz2.BZ2File(os.path.join(pathload, 'para_index.pckl'), 'rb') as f:
            paramID_sigbetter, paramID_nosigbetter, paramID_mis_losesig, paramID_mis_gainsig = pickle.load(f)
    except FileNotFoundError:
        paramID_sigbetter, paramID_nosigbetter = [], []
        paramID_mis_losesig = {i_mis: {i_LOD: [] for i_LOD in LOD} for i_mis in MISSPECIFICATION_LOD_STR}
        paramID_mis_gainsig = deepcopy(paramID_mis_losesig)
        with tqdm(total=len(paramID), ncols=150) as pbar:
            for i, i_paramID in enumerate(paramID):
                if i_paramID in info['total']['paramID sigbetter']:
                    paramID_sigbetter.append(i)
                    for j, j_mis in enumerate(MISSPECIFICATION_LOD_STR):
                        for k, k_LOD in enumerate(LOD):
                            jk_info_mis = info_mis[j][k]['total']['paramID sigbetter']
                            if i_paramID not in jk_info_mis:
                                paramID_mis_losesig[j_mis][k_LOD].append(i)
                else:
                    paramID_nosigbetter.append(i)
                    for j, j_mis in enumerate(MISSPECIFICATION_LOD_STR):
                        for k, k_LOD in enumerate(LOD):
                            jk_info_mis = info_mis[j][k]['total']['paramID sigbetter']
                            if i_paramID in jk_info_mis:
                                paramID_mis_gainsig[j_mis][k_LOD].append(i)
                pbar.update(1)
        assert set(paramID_sigbetter + paramID_nosigbetter) == set(range(len(paramID)))
        with bz2.BZ2File(os.path.join(pathload, 'para_index.pckl'), 'wb') as f:
            pickle.dump((paramID_sigbetter, paramID_nosigbetter, paramID_mis_losesig, paramID_mis_gainsig), f)

    # barplot
    df = pd.DataFrame(columns=['num sigbetter', 'MISSPECIFICATION', 'LOD'])
    df_nomis = pd.DataFrame({'num sigbetter': [len(info['total']['paramID sigbetter'])]*len(MISSPECIFICATION_LOD_STR),
                             'MISSPECIFICATION': MISSPECIFICATION_LOD_STR, 'LOD': ['nomis']*len(MISSPECIFICATION_LOD_STR)})
    df = pd.concat([df, df_nomis], ignore_index=True)
    for i, i_mis in enumerate(MISSPECIFICATION_LOD_STR):
        for j, j_LOD in enumerate(LOD):
            i_df = pd.DataFrame({'num sigbetter': [len(info_mis[i][j]['total']['paramID sigbetter'])],
                                 'MISSPECIFICATION': [i_mis], 'LOD': [j_LOD]})
            df = pd.concat([df, i_df], ignore_index=True)

    color = ['r', 'g', 'm', 'b']
    my_pal = dict(zip(MISSPECIFICATION_LOD_STR, color))
    fig, ax = plt.subplots()
    sns.barplot(df, x='LOD', y='num sigbetter', hue='MISSPECIFICATION',  gap=0.15, ax=ax, palette=my_pal)
    fig.set_size_inches(12.5, 5)
    ax.set_ylim([0, 1200000])
    ax.ticklabel_format(style='plain', axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(path_folder, 'numsig.pdf'), dpi=FIG_DPI)
    plt.close()

    # para['X'] = [sum([para['Spop'][i], para['R1pop'][i], para['R2pop'][i], para['R12pop'][i]]) for i in range(len(para['Spop']))]
    # para_df = pd.DataFrame.from_dict(para)
    # para_df['Spop'] = para_df['Spop']/para_df['X']
    # para_df['R1pop'] = para_df['R1pop']/para_df['X']
    # para_df['R2pop'] = para_df['R2pop']/para_df['X']
    # para_df['R12pop'] = para_df['R12pop']/para_df['X']

    '''supervenn plot'''
    # set_total = set(range(len(para['paramID'])))
    # set_sigbetter = set(paramID_sigbetter)
    # for i_mis in MISSPECIFICATION_LOD_STR:
    #     sets = [set_total, set_sigbetter]
    #     labels = ['total', 'sigbetter']
    #     for i_LOD in LOD:
    #         sets.append(set(paramID_mis_losesig[i_mis][i_LOD]))
    #         labels.append('lose sigbetter ' + i_LOD)
    #     plt.figure(figsize=(18, 8))
    #     supervenn(sets, labels, widths_minmax_ratio=0.1, col_annotations_area_height=1.2, side_plots='right', chunks_ordering='minimize gaps')
    #     plt.title('mis ' + i_mis, fontsize=18)
    #     plt.ylabel('Sets', fontsize=18)
    #     plt.xlabel('Items', fontsize=18)
    #     plt.tight_layout()
    #     plt.savefig(os.path.join(path_folder, 'set ' + i_mis + '.pdf'), dpi=FIG_DPI)
    #     plt.close('all')
    #
    # ind_par = 1
    # DPM_analysis_misspec_par_2('Spop'), DPM_analysis_misspec_par_2('R1pop'), DPM_analysis_misspec_par_2('R2pop')
    #
    # result = {'total': None, 'sigbetter': None, **dict(zip(LOD, [None]*len(LOD)))}
    # val = {i_mis: {i_para: deepcopy(result) for i_para in list(para.keys())[ind_par:]} for i_mis in MISSPECIFICATION_LOD_STR}
    # for i_mis in MISSPECIFICATION_LOD_STR:
    #     for i, i_para in enumerate(list(para.keys())[ind_par:]):
    #         if i_para in ['X', 'Spop', 'R1pop', 'R2pop', 'R12pop', 'g0_S', 'T.R1..S.', 'T.R2..S.']:
    #
    #             i_paraval_sorted = sorted(set(para_df[i_para]))
    #             val[i_mis][i_para]['total'] = DPM_analysis_misspec_par_1(para_df[i_para].values.tolist(), i_paraval_sorted)
    #
    #             i_paraval_sigbetter = para_df.iloc[paramID_sigbetter, para_df.columns.get_loc(i_para)].values.tolist()
    #             val[i_mis][i_para]['sigbetter'] = DPM_analysis_misspec_par_1(i_paraval_sigbetter, i_paraval_sorted)
    #
    #             for i_LOD in LOD:
    #                 i_parval_mis = para_df.iloc[paramID_mis_losesig[i_mis][i_LOD], para_df.columns.get_loc(i_para)].values.tolist()
    #                 val[i_mis][i_para][i_LOD] = DPM_analysis_misspec_par_1(i_parval_mis, i_paraval_sorted)
    #         else:
    #             val[i_mis][i_para]['total'] = para[i_para]
    #             val[i_mis][i_para]['sigbetter'] = para_df.iloc[paramID_sigbetter, para_df.columns.get_loc(i_para)].values.tolist()
    #             for i_LOD in LOD:
    #                 val[i_mis][i_para][i_LOD] = para_df.iloc[paramID_mis_losesig[i_mis][i_LOD], para_df.columns.get_loc(i_para)].values.tolist()
    #
    # color = ['r', 'b', 'g', 'c', 'y', 'k', 'm', 'blueviolet']
    # title = ['Spop', 'R1pop', 'R2pop', 'R12pop', 'g0', 'Sa.S.D1', 'Sa.S.D2', 'Sa.R1.D1', 'Sa.R2.D2', 'T.StoR1', 'T.StoR2']
    # x = list(result.keys())
    # for i_mis in MISSPECIFICATION_LOD_STR:
    #     for i, i_para in enumerate(list(para.keys())[ind_par:]):
    #         print(i_para)
    #         if i_para == 'X':
    #             continue
    #         plt.rcParams['font.size'] = 21
    #         plt.figure()
    #         fig = plt.gcf()
    #         fig.set_size_inches(24, 11)
    #         i_val = val[i_mis][i_para]
    #         if i_para in ['Spop', 'R1pop', 'R2pop', 'R12pop', 'g0_S', 'T.R1..S.', 'T.R2..S.']:
    #             i_paraval_sorted = sorted(set(para_df[i_para]))
    #             if len(i_paraval_sorted) > len(color):
    #                 color = ["#" + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)]) for i in range(len(i_paraval_sorted))]
    #             i_val = pd.DataFrame.from_dict(i_val)
    #             baseline = np.zeros(i_val.shape[1])
    #
    #             for j in range(i_val.shape[0]):
    #                 i_row = i_val.iloc[j, :]
    #                 plt.plot(x, i_row, color=color[j], linewidth=2, marker='d', markersize=7)
    #                 # plt.bar(x, i_row.values, bottom=baseline, color=color[j])
    #                 # baseline += i_row.values
    #
    #             plt.ylim([-10, 80])
    #             plt.yticks(list(range(0, 80, 10)))
    #             plt.ylabel('percentage')
    #             plt.legend(['{:0.2e}'.format(i) for i in i_paraval_sorted], loc='upper center', frameon=False, ncol=int(len(i_paraval_sorted)))
    #         else:
    #             df = pd.DataFrame(columns=['group', 'value'])
    #             for key, value in i_val.items():
    #                 i_df = pd.DataFrame({'group': np.repeat(key, len(value)), 'value': value})
    #                 df = pd.concat([df, i_df])
    #
    #             df = pd.DataFrame(df.to_dict('records'))
    #             sns.violinplot(x='group', y='value', data=df, order=i_val.keys())
    #             # plt.xticks(labels=i_val.keys())
    #             # plt.yscale('log', base=10)
    #             # plt.ylim([1e-15, 1e3])
    #
    #         plt.title(title[i] + '  ' + i_mis)
    #         plt.xlabel('Group')
    #         plt.show()
    #         plt.savefig(os.path.join(pathsave, 'par', title[i] + '  ' + i_mis + '.pdf'), dpi=FIG_DPI)
    #         plt.close('all')
    return


def DPM_analysis_misspec_plot(Strategy_name, Simduration, data, titlestr):
    plt.ioff()
    km = {key: data[key] for key in Strategy_name}
    hz, p = data['hz_ratio'], data['p']
    for i in hz.keys():
        idx_p = [idx for idx, i_val in enumerate(p.keys()) if set(i) == set(i_val)][0]
        km_ref = km[i[0]]
        km_treat = km[i[1]]
        par = {'color': ['k', 'b'], '0': i[0], '1': i[1], 'hzr': hz[i], 'p': p[list(p.keys())[idx_p]], 'duration': Simduration,
               'xtick step': 300, 'totalnum': data['num']}
        DPM_plot_KM(km_ref, km_treat, par, titlestr)
        if PLT_INTERACTIVE:
            plt.show()
        # plt.close('all')
    return


def DPM_analysis_hazard_ratio(stoptime, Strategy_name, Simduration):
    p = DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration)
    hz = dict()
    for i in list(itertools.combinations(Strategy_name, 2)):
        if 'strategy0' in i:
            ref = stoptime['strategy0']
            treat = stoptime[list(filter(lambda x: x != 'strategy0', i))[0]]
            name = ('strategy0', list(filter(lambda x: x != 'strategy0', i))[0])
        else:
            ref, treat = stoptime[i[0]], stoptime[i[1]]
            name = i
        hz[name] = DPM_analysis_HZ(ref, treat, Simduration)
    return hz, p


def DPM_analysis_KM(data, duration):
    E = [1 if i_val <= duration else 0 for i_val in data]
    epf = ExponentialFitter().fit(data, E)
    hazard = {'mean': epf.hazard_.mean(), 'ci': epf.confidence_interval_hazard_.mean()}

    kmf = KaplanMeierFitter()
    kmf.fit(data, E)
    km_interval = kmf.confidence_interval_survival_function_
    km = kmf.survival_function_

    t = km.index.values
    val = km.iloc[:].values.flatten()
    interval_lower = km_interval.iloc[:, 0].values.flatten()
    interval_upper = km_interval.iloc[:, 1].values.flatten()
    median_survival = kmf.median_survival_time_

    idx = np.where(t <= duration)[0]
    t, val, interval_lower, interval_upper = t[idx], val[idx], interval_lower[idx], interval_upper[idx]
    t, val, interval_lower, interval_upper = np.append(t, duration), \
        np.append(val, val[-1]), \
        np.append(interval_lower, interval_lower[-1]), \
        np.append(interval_upper, interval_upper[-1])

    return {'t': t, 'val': val, 'median_survival': median_survival, 'interval_lower': interval_lower, 'interval_upper': interval_upper,
            'hazard': hazard}


def DPM_analysis_HZ(data_ref, data, Simduration):
    if (len(data_ref) != 0) & (len(data) != 0):
        treat = np.concatenate((np.zeros(len(data_ref)), np.ones(len(data))))
        val = data_ref + data
        E = [1 if i_val <= Simduration else 0 for i_val in data_ref]
        E.extend([1 if i_val <= Simduration else 0 for i_val in data])
        E = np.array(E)
        d = {'val': val, 'E': E, 'treat': treat}
        df = pd.DataFrame(data=d)
        cph = CoxPHFitter()
        cph.fit(df, duration_col='val', event_col='E')
        hz_ratio = cph.hazard_ratios_.values[0]
    else:
        hz_ratio = None
    return hz_ratio


def DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration):
    G, T, E = tuple(), tuple(), tuple()
    flag_empty = False
    for i_strategy in Strategy_name:
        i_stop = stoptime[i_strategy]
        if not i_stop:
            flag_empty = True
            break
        G = G + tuple(len(i_stop) * [str(i_strategy)])
        E = E + tuple([1 if i_val <= Simduration else 0 for i_val in i_stop])
        T = T + tuple(i_stop)
    p = statistics.pairwise_logrank_test(T, G, E) if not flag_empty else None
    p_out = dict(zip(p.name, p.p_value)) if not flag_empty else None
    return p_out


def DPM_analysis_dose(dose, strategyname, inddrug=1):
    firstuse, max_num_change, num_change, average_move_num = [], [], [], []
    use_atbegin = 0
    for i, i_dose in enumerate(dose):
        i_dose = i_dose.split(';')
        if '-1' in i_dose:
            i_dose = i_dose[:i_dose.index('-1')]
        i_firstuse, i_average_move_num, drugovermove, drugtotal, i_num_change,  i_current = None, None, 0, 0, 0, i_dose[0]
        if i_current == '(0.0,1.0)':
            use_atbegin += 1

        for j, j_step in enumerate(i_dose):
            i_val = [float(i_val) for i_val in j_step[1:-2].split(',')]
            drugtotal = drugtotal + i_val[inddrug]
            drugovermove = drugovermove + (j+1) * i_val[inddrug]
            if strategyname == 'strategy0' and i_current not in ['(0.0,1.0)', '(0.0,0.0)', '(1.0,0.0)']:
                assert strategyname == 'strategy0' and i_current not in ['(0.0,1.0)', '(0.0,0.0)', '(1.0,0.0)']
            if i_firstuse is None and i_val[inddrug] != 0:
                i_firstuse = j+1
            if j_step != i_current and j_step != '(0.0,0.0)':
                i_num_change += 1
                i_current = j_step

        i_average_move_num = drugovermove/drugtotal if drugtotal != 0 else None
        i_firstuse = len(i_dose) + 1 if i_firstuse is None else i_firstuse
        i_average_move_num = len(i_dose) + 1 if i_average_move_num is None else i_average_move_num
        firstuse.append(i_firstuse)
        average_move_num.append(i_average_move_num)
        num_change.append(i_num_change/len(i_dose))
        max_num_change.append(i_num_change)

        # if strategyname == 'strategy0' and i_num_change > 1:
        #     print(i_num_change)
        #     print(i_dose)

        # firstuse = np.array(firstuse)
        # num_change = np.array(num_change)
        # a = firstuse[num_change==0]

    return dict(zip(['first drug 2', 'average move number', 'drug changes', 'max drug changes'],
                    [np.mean(firstuse), np.mean(average_move_num), np.mean(num_change), np.max(max_num_change)]))


def DPM_analysis_sigbetter(stoptime, Strategy):
    paramID, stoptime_ref, stoptime_test = stoptime['paramID'], stoptime[Strategy[0]], stoptime[Strategy[1]]
    ind = np.logical_and(np.array(stoptime_test) > np.array(stoptime_ref) + 30*2, np.array(stoptime_test) > 1.25 * np.array(stoptime_ref))
    return list(itertools.compress(paramID, ind))
