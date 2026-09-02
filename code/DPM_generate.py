from DPM_lib import itertools, sys, rv_continuous, loguniform, simpson, plt, deepcopy, math
from DPM_constant import *
""" This script checks arguments used in the DPM (Dynamic Precision Medicine) model. """


# Generate the default discrete dose combination based on drug number.
def DPM_generate_discrete_dose_combination(num_drug):
    dose_list = [float(0)]
    for i in range(int(num_drug)):
        dose_list.append(float(1 / (i + 1)))
    dose_list.sort()
    dose_combination = itertools.product(dose_list, repeat=int(num_drug))
    dose_combination = [x for x in dose_combination if sum(x) == 1]
    return dose_combination


# Generate the continuous dose combination based on drug number and dose interval.
def DPM_generate_continuous_dose_combination(num_drug, dose_interval):
    dose_list = np.arange(0, 1+dose_interval, dose_interval, dtype=float)
    dose_list.sort()
    dose_combination = itertools.product(dose_list, repeat=int(num_drug))
    dose_combination = [x for x in dose_combination if sum(x) == 1]
    return dose_combination


# Generate default parameters in 2 drug case according to PNAS(2012).
def DPM_generate_default_par_2drug():
    # g0: basel proliferation rate for 4 kinds for cells: S, R1, R2, R12.
    g0 = [0.001, 0.0026, 0.007, 0.0184, 0.0487, 0.1287, 0.34]
    # R1 to X0 ratio.
    ratioR1toX0 = [0, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 9e-1]
    # R2 to X0 ratio.
    ratioR2toX0 = [0, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 9e-1]
    # Sa:4×2 matrix of drug sensitivities.
    # Sa(S,D1)/g0.
    Sa_ratio_S_D1tog0 = [5.6e-4, 0.0054, 0.0517, 0.4964, 4.7683, 45.8045, 440]
    # Sa(S,D2)/Sa(S,D1).
    Sa_ratio_S_D2toS_D1 = [4e-4, 0.0015, 0.0054, 0.02, 0.0737, 0.2714, 1e0]
    # Sa(R1,D1)/Sa(S,D1).
    Sa_ratio_R1_D1toS_D1 = [0, 1e-5, 9.5635e-5, 9.1461e-4, 0.0087, 0.0837, 0.8]
    # Sa(R2,D2)/Sa(S,D2).
    Sa_ratio_R2_D2toS_D2 = [0, 1e-5, 9.5635e-5, 9.1461e-4, 0.0087, 0.0837, 0.8]
    # T:4×4 transition rate matrix, R1->R12=S->R2; R2->R12=S->R1.
    T_StoR1 = [1e-11, 2.154e-10, 4.642e-9, 1e-7, 2.154e-6, 4.642e-5, 1e-3]
    T_StoR2 = [1e-11, 2.154e-10, 4.642e-9, 1e-7, 2.154e-6, 4.642e-5, 1e-3]
    return g0, ratioR1toX0, ratioR2toX0, Sa_ratio_S_D1tog0, Sa_ratio_S_D2toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, T_StoR1, T_StoR2


# Generate default parameter criterisa index based on drug number.
def DPM_generate_PAR_criterisa_list(Num_drug):
    if Num_drug == 2:
        PAR_criterisa_list = [0, 1, 2, 3, 4, 5, 6]
    else:
        PAR_criterisa_list = None
    return PAR_criterisa_list


# Generate heading for .csv saving files.
def DPM_generate_heading_csv(Num_drug, Simduration, Stepsize):
    assert Num_drug == 2
    # Heading of param csv file.
    Heading_param_csv: list[str] = HEADING_2DRUG_PARAM_CSV
    # Heading of stopt csv file.
    Heading_stopt_csv = HEADING_STOPT_CSV
    # Heading of dosage csv file.
    Heading_dosage_csv = HEADING_DOSAGE_CSV
    Heading_dosage_csv_str = 'Drug1 dosage,Drug2 dosage'
    Heading_dosage_csv.extend([f'({Heading_dosage_csv_str}) at t={t_i}' for t_i in np.arange(0, Simduration, Stepsize)])
    # Heading of pop csv file.
    Heading_pop_csv = HEADING_POP_CSV
    Heading_pop_csv_str = ','.join(ALL_POSSIBLE_CELLTYPE_2DRUG)
    Heading_pop_csv.extend(f'({Heading_pop_csv_str}) at t={t_i}' for t_i in np.arange(Stepsize, Simduration + Stepsize, Stepsize))
    # Heading of eachtimepoint csv file.
    Heading_eachtimepoint_csv = HEADING_EACHTIMEPOINT_CSV
    Heading_eachtimepoint_csv_str = [f'({Heading_pop_csv_str}) at t={t_i}, ({Heading_dosage_csv_str}) at t={t_i}'.split(', ')
                                     for t_i in np.arange(0, Simduration + SIMTIMESTEP_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL)]
    Heading_eachtimepoint_csv.extend(list(itertools.chain.from_iterable(Heading_eachtimepoint_csv_str)))
    return Heading_param_csv, Heading_stopt_csv, Heading_dosage_csv, Heading_pop_csv, Heading_eachtimepoint_csv


# Generate initial cell number X0 based on input parameter.
def DPM_generate_X0(PAR):
    assert PAR['Num_drug'] == 2
    return DPM_generate_X0_2drug(PAR)


# Generate initial cell number X0 based on input parameter, 2 drug case.
def DPM_generate_X0_2drug(PAR):
    return np.array([PAR['Spop'], PAR['R1pop'], PAR['R2pop'], PAR['R12pop']], dtype=float)


# Generate g0 based on input parameter.
def DPM_generate_g0(PAR):
    assert PAR['Num_drug'] == 2
    return DPM_generate_g0_2drug(PAR)


# Generate g0 based on input parameter, 2 drug case.
def DPM_generate_g0_2drug(PAR):
    return np.array([PAR['g0_S'], PAR['g0_R1'], PAR['g0_R2'], PAR['g0_R12']], dtype=float)


# Generate Sa based on input parameter.
def DPM_generate_Sa(PAR):
    assert PAR['Num_drug'] == 2
    return DPM_generate_Sa_2drug(PAR)


# Generate Sa based on input parameter, 2 drug case.
def DPM_generate_Sa_2drug(PAR):
    Sa = np.zeros((PAR['Num_cell_type'], PAR['Num_drug']), dtype=float)

    # Sensitivity of S cell on drug 1 and drug 2.
    Sa[0, :] = np.array([PAR['Sa.S.D1.'], PAR['Sa.S.D2.']])
    # Sensitivity of R1 cell on drug 1 and drug 2.
    Sa[1, :] = np.array([PAR['Sa.R1.D1.'], PAR['Sa.R1.D2.']])
    # Sensitivity of R2 cell on drug 1 and drug 2.
    Sa[2, :] = np.array([PAR['Sa.R2.D1.'], PAR['Sa.R2.D2.']])
    # Sensitivity of R12 cell on drug 1 and drug 2.
    Sa[3, :] = np.array([PAR['Sa.R12.D1.'], PAR['Sa.R12.D2.']])

    return Sa


# Generate T based on input parameter.
def DPM_generate_T(PAR):
    assert PAR['Num_drug'] == 2
    T = DPM_generate_T_2drug(PAR)
    return T


# Generate T based on input parameter, 2 drug case.
def DPM_generate_T_2drug(PAR):
    T = np.zeros((PAR['Num_cell_type'], PAR['Num_cell_type']), dtype=float)
    # Transition to S.
    T[0, :] = np.array([PAR['T.S..S.'], PAR['T.S..R1.'], PAR['T.S..R2.'], PAR['T.S..R12.']])
    # Transition to R1.
    T[1, :] = np.array([PAR['T.R1..S.'], PAR['T.R1..R1.'], PAR['T.R1..R2.'], PAR['T.R1..R12.']])
    # Transition to R2.
    T[2, :] = np.array([PAR['T.R2..S.'], PAR['T.R2..R1.'], PAR['T.R2..R2.'], PAR['T.R2..R12.']])
    # Transition to R12.
    T[3, :] = np.array([PAR['T.R12..S.'], PAR['T.R12..R1.'], PAR['T.R12..R2.'], PAR['T.R12..R12.']])
    return T


def DPM_generate_pdf(x, k, LOD):
    return (k/x**2)/(1-k/LOD+k)


def DPM_generate_pho_cumulative(p_lower, p_upper, k, pho_min, pho_max):
    return k/p_lower-k/p_upper+2*k*(1/p_lower-1/p_upper)/(1/pho_min-1/pho_max)


def DPM_generate_misspecification_subclone(par, subclone_LOD, misspecification_LOD, mutation_rate, celltype, mis_specfiy_pop, Strategy_name,
                                           norm_pdf=1, Xtotal=None):

    # def p(x, k, LOD):
    #     return (k/x**2)/(1-k/LOD+k)

    class MutationFraction(rv_continuous):
        def _pdf(self, x, k, LOD, const):
            return (1.0/const)*DPM_generate_pdf(x, k, LOD)

    # x = np.linspace(mutation_rate, subclone_LOD, int(1e7))
    # norm_pdf = simps(DPM_generate_pdf(x, mutation_rate, subclone_LOD), x)
    mutationFraction_distribution = MutationFraction(name='mutationFraction_distribution', a=mutation_rate, b=subclone_LOD)
    # pdf = mutationFraction_distribution.pdf(k=mutation_rate, LOD=subclone_LOD, x=x, const=norm_pdf)
    # cdf = mutationFraction_distribution.cdf(k=mutation_rate, LOD=subclone_LOD, x=x, const=norm_pdf)
    # plt.plot(x, cdf)
    # plt.xscale("log")
    # plt.yscale("log")

    for i, i_par in enumerate(par):
        # print(i)
        i_par_ori = deepcopy(i_par)
        flag_mis = False
        # Some virtual patients have a total of 5,000,000,000 cells, while others have 5,000,050,000 cells.
        i_X = i_par['Spop'] + i_par['R1pop'] + i_par['R2pop'] + i_par['R12pop']

        # if Xtotal is None:
        #     i_X = i_par['Spop'] + i_par['R1pop'] + i_par['R2pop'] + i_par['R12pop']
        # else:
        #     i_X = Xtotal
        for i_sub in celltype:
            if i_par[i_sub]/i_X < 2*subclone_LOD:
                flag_mis = True
                mis_specfiy_pop[Strategy_name][i_sub]['total'].append(i_X)
                mis_specfiy_pop[Strategy_name][i_sub]['percent'].append(i_par[i_sub]/i_X)
                if misspecification_LOD == 0:
                    i_par[i_sub] = 0.0
                    mis_specfiy_pop[Strategy_name][i_sub]['est'].append(i_par[i_sub]/i_X)
                elif misspecification_LOD == 'max':
                    if i_sub in ['R1pop', 'R2pop']:
                        val = i_X * (2 * subclone_LOD - 1e-9)
                        i_par[i_sub] = np.ceil(val).item() if val >= 1 else 0
                    else:
                        i_par[i_sub] = 0.0
                    mis_specfiy_pop[Strategy_name][i_sub]['est'].append(i_par[i_sub]/i_X)
                elif misspecification_LOD == 'pdf':
                    if i_sub in ['R1pop', 'R2pop']:
                        i_rvs = mutationFraction_distribution.rvs(k=mutation_rate, LOD=subclone_LOD, const=norm_pdf, size=1)
                        val = i_X * 2 * i_rvs
                        i_par[i_sub] = val[0]  # np.ceil(val).item() if val >= 1 else 0
                    else:
                        i_par[i_sub] = 0.0
                    est = i_par[i_sub]/i_X
                    assert est <= 2 * subclone_LOD
                    mis_specfiy_pop[Strategy_name][i_sub]['est'].append(est)
                elif misspecification_LOD == 'loguni':
                    if i_sub in ['R1pop', 'R2pop']:
                        i_rvs = loguniform.rvs(mutation_rate, subclone_LOD, size=1)
                        val = i_X * 2 * i_rvs
                        i_par[i_sub] = val[0]  # np.ceil(val).item() if val >= 1 else 0
                    else:
                        i_par[i_sub] = 0.0
                    est = i_par[i_sub]/i_X
                    assert est <= 2 * subclone_LOD
                    mis_specfiy_pop[Strategy_name][i_sub]['est'].append(est)
        if misspecification_LOD == 'pdf':
            if i_par['R1pop'] > 0:
                assert i_par['R1pop'] >= 2 * mutation_rate * i_X
            if i_par['R2pop'] > 0:
                assert i_par['R2pop'] >= 2 * mutation_rate * i_X
        if flag_mis:
            i_par['Spop'] = i_X - i_par['R1pop'] - i_par['R2pop'] - i_par['R12pop']
            assert math.isclose(i_par['Spop'] + i_par['R1pop'] + i_par['R2pop'] + i_par['R12pop'],
                                i_par_ori['Spop'] + i_par_ori['R1pop'] + i_par_ori['R2pop'] + i_par_ori['R12pop'],
                                rel_tol=1e-6)
        else:
            assert i_par == i_par_ori

        i_par = {k: float(v) if isinstance(v, np.floating) else v for k, v in i_par.items()}
        par[i] = i_par
    return par, mis_specfiy_pop
