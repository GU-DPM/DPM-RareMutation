from DPM_assign_check import *


def DPM_strategy_pro(X0, g0, Sa, T, Num_drug, tnow, X_Strategy_sim, X_Strategy_true, t_Strategy_sim, d_Strategy_sim, Stepsize, maxthreshold,
                     Simtimestep, Limit_mortality, Limit_radiologicdetection,  subclone_LOD, misspecification_LOD, mutation_rate, celltype,
                     mis_specfiy_pop, Strategy_name):
    flag_treat, flag_break = False, False
    # no treat, let proliferate
    treat_i = np.zeros([Num_drug, 1], dtype=int)
    t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

    X_true_i, t_true_i, _ = DPM_sim_model(X0, T, g0, Sa, d_i, t_i, maxthreshold, True)
    X_true_total_i = np.sum(X_true_i, axis=0)

    # True dead, stop sim.
    if np.any(X_true_total_i >= Limit_mortality):
        flag_break = True

    # True cell number of all types are smaller than 1 (cured), should not happen, error.
    assert not np.all(X_true_i[:, -1] < 1)

    X_sim = np.tile(X_Strategy_sim[:, [-1]], (1, X_true_i.shape[1]))
    X_Strategy = np.append(X_Strategy_sim, X_sim, 1) if tnow == 0 else np.append(X_Strategy_sim, X_sim[:, 1:], 1)
    t_Strategy = np.append(t_Strategy_sim, t_i + tnow) if tnow == 0 else np.append(t_Strategy_sim, t_true_i[1:] + tnow)
    d_Strategy = np.append(d_Strategy_sim, d_i, 1)

    X_Strategy_true = np.append(X_Strategy_true, X_true_i, 1) if tnow == 0 else np.append(X_Strategy_true, X_true_i[:, 1:], 1)

    # true total cell population reemerges from radiologic detection, do subclone detection and treat again
    if X_true_total_i[-1] >= Limit_radiologicdetection:
        flag_treat = True
        X0_LOD, mis_specfiy_pop = DPM_strategy_subclone_detection(X_true_i[:, -1],
                                                                  subclone_LOD,
                                                                  Limit_radiologicdetection,
                                                                  misspecification_LOD,
                                                                  mutation_rate,
                                                                  celltype,
                                                                  mis_specfiy_pop,
                                                                  Strategy_name)
        # did subclone detection, assign X0_LOD to X_Strategy
        X_Strategy[:, -1] = X0_LOD

    X0_true = X_true_i[:, -1]
    return X0_true, X_Strategy, X_Strategy_true, t_Strategy, d_Strategy, flag_treat, flag_break, mis_specfiy_pop


def DPM_strategy_LOD(X_sim, t_sim, d_sim, X_true, t_true, tnow, X_Strategy_sim, X_Strategy_true, t_Strategy_sim, d_Strategy_sim,
                     Stepsize, Limit_mortality, Limit_radiologicdetection, subclone_LOD, misspecification_LOD, mutation_rate, celltype,
                     mis_specfiy_pop, Strategy_name):
    flag_break, subclone_detection, flag_treat = False, False, True
    X_sim_total = np.sum(X_sim, axis=0)
    X_true_total = np.sum(X_true, axis=0)

    # If total true cell population is >= Limit_mortality (death) or
    # the true cell number of all types are smaller than 1 (cured), stop.
    if np.any(X_true_total >= Limit_mortality) or np.all(X_true[:, -1] < 1):
        X_sim = X_sim[:, :X_true.shape[1]]
        t_sim = t_sim[:X_true.shape[1]]
        flag_break = True

    X_Strategy = np.append(X_Strategy_sim, X_sim, 1) if tnow == 0 else np.append(X_Strategy_sim, X_sim[:, 1:], 1)
    t_Strategy = np.append(t_Strategy_sim, t_sim + tnow) if tnow == 0 else np.append(t_Strategy_sim, t_sim[1:] + tnow)
    d_Strategy = np.append(d_Strategy_sim, d_sim, 1)

    X_Strategy_true = np.append(X_Strategy_true, X_true, 1) if tnow == 0 else np.append(X_Strategy_true, X_true[:, 1:], 1)

    # If the simulated cell numbers for all types are less than 1
    if np.all(X_sim[:, -1] < 1):
        # If true cell population is >= radiologicdetection
        if X_true_total[-1] >= Limit_radiologicdetection:
            X_true_total_atStepsize = np.sum(X_Strategy_true[:, range(0, X_Strategy_true.shape[1], Stepsize)][:, :-1], axis=0)
            replapse = np.any(X_true_total_atStepsize < Limit_radiologicdetection)
            # If it is a relapse (true total cell number reemerge from radiologic LOD), do subclone detection
            # If it is a relapse(the true total cell population re-emerges above the radiologic LOD), perform subclone detection.
            if replapse:
                subclone_detection = True
                X_LOD, mis_specfiy_pop = DPM_strategy_subclone_detection(X_true[:, -1],
                                                                         subclone_LOD,
                                                                         Limit_radiologicdetection,
                                                                         misspecification_LOD,
                                                                         mutation_rate,
                                                                         celltype,
                                                                         mis_specfiy_pop,
                                                                         Strategy_name)
                X_sim[:, -1] = X_LOD
        # If true cell population is < radiologicdetection stop treating because doctors cannot not detect the tumor.
        elif X_true_total[-1] < Limit_radiologicdetection:
            flag_treat = False
    # Subclone detection will be performed if the simulated cell population is >= Limit_mortality threshold,
    # while the actual cell population remains below this threshold.
    if np.any(X_sim_total >= Limit_mortality) and np.all(X_true_total < Limit_mortality):
        subclone_detection = True
        X_LOD, mis_specfiy_pop = DPM_strategy_subclone_detection(X_true[:, -1],
                                                                 subclone_LOD,
                                                                 Limit_radiologicdetection,
                                                                 misspecification_LOD,
                                                                 mutation_rate,
                                                                 celltype,
                                                                 mis_specfiy_pop,
                                                                 Strategy_name)
        X_sim[:, -1] = X_LOD
    return X_Strategy, X_Strategy_true, t_Strategy, d_Strategy, X_sim, flag_break, subclone_detection, flag_treat, mis_specfiy_pop


def DPM_strategy_subclone_detection(x, subclone_LOD, Limit_radiologicdetection, misspecification_LOD, mutation_rate, celltype,
                                    mis_specfiy_pop, Strategy_name):
    # If true cell population is bigger or radiologic detection, can do a tissue biopsy otherwise a liquid biopsy, the LOD will be
    # higher in liquid biopsy.
    liquid_biopsy_ratio = LIQUID_BIOPSY_LOD_RATIO_DEFAULT_VAl  # else 1
    subclone_LOD_i = subclone_LOD * liquid_biopsy_ratio if np.sum(x) < Limit_radiologicdetection else subclone_LOD
    if subclone_LOD_i > 1e-1:
        subclone_LOD_i = 1e-1
    par = dict(zip(celltype, x))
    X_LOD, mis_specfiy_pop = DPM_generate_misspecification_subclone([par],
                                                                    subclone_LOD_i,
                                                                    misspecification_LOD,
                                                                    mutation_rate,
                                                                    celltype[1:],
                                                                    mis_specfiy_pop,
                                                                    Strategy_name)
    X_LOD = np.fromiter(X_LOD[0].values(), dtype=float)
    return X_LOD, mis_specfiy_pop


def DPM_strategy_select_drug_majority(x, Num_drug, Sa):
    treat = np.zeros([Num_drug, 1], dtype=int)
    # Index of the majority cell type.
    ind_drug = np.argmax(Sa[np.argmax(x), :])
    treat[ind_drug] = 1
    return treat, ind_drug


# Define Strategy 0.
def DPM_strategy_0(PAR, Simduration, Stepsize, Simtimestep, Limit_mortality, Limit_radiologicdetection, LSsim, misspecification_ofdecision,
                   misspecification_oftrue, mis_PAR):
    # Strategy 0: Current personalized medicine:
    maxthreshold = Limit_mortality
    if misspecification_oftrue:
        X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(mis_PAR), DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)
    else:
        X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_Sa_i = DPM_generate_Sa(mis_PAR) if misspecification_ofdecision else None

    nadir = X0_i.sum()
    tnow = 0.0
    t_Strategy0_i = np.zeros([0], dtype=float)
    X_Strategy0_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy0_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    drug_used = np.zeros([PAR['Num_drug']], dtype=int)

    if misspecification_ofdecision:
        treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_i, PAR['Num_drug'], mis_Sa_i)
    else:
        treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_i, PAR['Num_drug'], Sa_i)
    drug_used[ind_drug] = 1

    while tnow < Simduration:
        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

        X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

        X_Strategy0_i = np.append(X_Strategy0_i, X_i, 1) if tnow == 0 else np.append(X_Strategy0_i, X_i[:, 1:], 1)
        t_Strategy0_i = np.append(t_Strategy0_i, t_i + tnow) if tnow == 0 else np.append(t_Strategy0_i, t_i[1:] + tnow)
        d_Strategy0_i = np.append(d_Strategy0_i, d_i, 1)
        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        # or if the cell number of all types are smaller than 1, cured. Stop.
        if (X_i[:, -1].sum() >= Limit_mortality) or np.all(X_i[:, -1] < 1):
            break

        # (1) If total cell population is bigger than 2 times of nadir, total cell population is bigger than Limit_radiologicdetection
        # (tumor can be detected), there are drugs have not been used (each drug is used only once) or
        # (2) If the total cell population reemerges from a level below the detection,
        # and there are drugs have not been used (each drug is used only once), switch to the other unused drugs.
        if ((X_i[:, -1].sum() >= 2 * nadir and X_i[:, -1].sum() >= Limit_radiologicdetection) or
                (X_i[:, 0].sum() < Limit_radiologicdetection <= X_i[:, -1].sum())) and np.any(drug_used == 0):
            treat_i = np.zeros([PAR['Num_drug'], 1], dtype=int)
            index_drug_not_used = [i for i, x in enumerate(drug_used.tolist()) if x == 0]
            treat_i[index_drug_not_used] = 1
            drug_used[index_drug_not_used] = 1
            nadir = X_i[:, -1].sum()

        # If the total cell population is smaller than the current nadir and bigger than the Limit_radiologicdetection (tumor can be detected),
        # update the nadir.
        if Limit_radiologicdetection <= X_i[:, -1].sum() < nadir:
            nadir = X_i[:, -1].sum()
        # Update X0_i and tnow.
        X0_i = X_i[:, -1]
        tnow += Stepsize
    return t_Strategy0_i, X_Strategy0_i, d_Strategy0_i


# Define Strategy 0, LOD
def DPM_strategy_LOD_0(PAR, Simduration, Stepsize, Simtimestep, Limit_mortality, Limit_radiologicdetection, PAR_LOD, subclone_LOD,
                       misspecification_LOD, mutation_rate, mis_specfiy_pop, Strategy_name):
    # Strategy 0: Current personalized medicine:
    maxthreshold = Limit_mortality
    X0_true_i = DPM_generate_X0(PAR)
    X0_LOD_i = DPM_generate_X0(PAR_LOD)
    g0_i = DPM_generate_g0(PAR)
    Sa_i = DPM_generate_Sa(PAR)
    T_i = DPM_generate_T(PAR)

    celltype = ['Spop', 'R1pop', 'R2pop', 'R12pop']
    nadir = X0_true_i.sum()
    if nadir < Limit_radiologicdetection:
        raise Exception('Initial tumor cell number is smaller than limit of radiologic detction.')
    tnow = 0.0
    t_Strategy0_i = np.zeros([0], dtype=float)
    X_Strategy0_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    X_Strategy0_true_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy0_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    drug_used = np.zeros([PAR['Num_drug']], dtype=int)

    treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_LOD_i, PAR['Num_drug'], Sa_i)
    drug_used[ind_drug] = 1

    flag_treat, subclone_detection = True, False
    while tnow < Simduration:
        if flag_treat:
            t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

            X_i, t_i, d_i = DPM_sim_model(X0_LOD_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, False)
            X_true_i, t_true_i, _ = DPM_sim_model(X0_true_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, True)
            X_true_total_i = np.sum(X_true_i, axis=0)

            X_Strategy0_i, X_Strategy0_true_i, t_Strategy0_i, d_Strategy0_i, X_i, flag_break, subclone_detection, flag_treat, mis_specfiy_pop =\
                DPM_strategy_LOD(X_i, t_i, d_i, X_true_i, t_true_i, tnow, X_Strategy0_i, X_Strategy0_true_i, t_Strategy0_i, d_Strategy0_i,
                                 Stepsize, Limit_mortality, Limit_radiologicdetection, subclone_LOD, misspecification_LOD, mutation_rate,
                                 celltype, mis_specfiy_pop, Strategy_name)
            if flag_break:
                break

            if subclone_detection:
                treat_i, ind_drug = DPM_strategy_select_drug_majority(X_i[:, -1], PAR['Num_drug'], Sa_i)
                drug_used[ind_drug] = 1

            # (1) If the true total cell population is bigger than 2 times of nadir and total cell population is bigger than
            # Limit_radiologicdetection (tumor can be detected), there are drugs have not been used (each drug is used only once) or
            # (2) If the true total cell population reemerges from radiologic detection,
            # and there are drugs have not been used (each drug is used only once), switch to the other unused drugs.
            if ((X_true_total_i[-1] >= 2 * nadir and X_true_total_i[-1] >= Limit_radiologicdetection) or
                    (X_true_total_i[0] < Limit_radiologicdetection <= X_true_total_i[-1])) \
                    and np.any(drug_used == 0) and (not subclone_detection):
                treat_i = np.zeros([PAR['Num_drug'], 1], dtype=int)
                index_drug_not_used = [i for i, x in enumerate(drug_used.tolist()) if x == 0]
                treat_i[index_drug_not_used], drug_used[index_drug_not_used] = 1, 1
                nadir = X_true_total_i[-1]
            # If the total cell population is smaller than the current nadir and bigger than the Limit_radiologicdetection
            # (tumor can be detected), update the nadir.
            if X_true_total_i[-1] < nadir:
                nadir = X_true_i[:, -1].sum() if X_true_i[:, -1].sum() >= Limit_radiologicdetection else 0
            # Update X0_LOD_i, X0_true_i and tnow.
            X0_LOD_i = X_i[:, -1]
            X0_true_i = X_true_i[:, -1]
        else:
            # no treat, let proliferate
            X0_true_i, X_Strategy0_i, X_Strategy0_true_i, t_Strategy0_i, d_Strategy0_i, flag_treat, flag_break, mis_specfiy_pop = \
                DPM_strategy_pro(X0_true_i, g0_i, Sa_i, T_i, PAR['Num_drug'], tnow, X_Strategy0_i, X_Strategy0_true_i, t_Strategy0_i,
                                 d_Strategy0_i, Stepsize, maxthreshold, Simtimestep, Limit_mortality, Limit_radiologicdetection,
                                 subclone_LOD, misspecification_LOD, mutation_rate, celltype, mis_specfiy_pop, Strategy_name)
            if flag_break:
                break

            if flag_treat:
                X0_LOD_i = X_Strategy0_i[:, -1]
                treat_i, ind_drug = DPM_strategy_select_drug_majority(X0_true_i, PAR['Num_drug'], Sa_i)
                drug_used[ind_drug] = 1

        tnow += Stepsize

    assert t_Strategy0_i.shape[0] == X_Strategy0_true_i.shape[1] == X_Strategy0_i.shape[1]
    diffpts = DPM_miscellaneous_treatment_change_time(d_Strategy0_i)/Stepsize
    # Drug changes occur at integer multiples of the step size.
    assert np.all(np.equal(np.mod(diffpts, 1), 0))
    return t_Strategy0_i, X_Strategy0_true_i, d_Strategy0_i, X_Strategy0_i, mis_specfiy_pop


# Define Strategy 1.
def DPM_strategy_1(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # Strategy 1: Minimize the total cell population.
    # In each Stepsize, select the d_i that minimizes the total cell population.
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy1_i = np.zeros([0], dtype=float)
    X_Strategy1_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy1_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

    while tnow < Simduration:
        X_i_total, t_i_total, d_i_total, X_i_end_total, t_sim = [], [], [], [], []
        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        for i_dose_combination in dose_combination:
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
            d_i_total.append(d_i)
            if mis_specification:
                mis_X_i, mis_t_i, _ = DPM_sim_model(X0_i, mis_T_i, mis_g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(mis_X_i)
                X_i_end_total.append(mis_X_i[:, -1].sum())
                t_sim.append(mis_t_i[-1])
            else:
                X_i, t__i, d__i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(X_i)
                t_i_total.append(t__i)
                d_i_total[-1] = d__i
                X_i_end_total.append(X_i[:, -1].sum())
                t_sim.append(t__i[-1])

        # Find the d_i that minimizes total cell population.
        if np.argmin(X_i_end_total) < Limit_mortality:
            index_select = int(np.argmin(X_i_end_total))
        else:
            # This indicates that the simulation for all dose combinations reaches mortality.
            # The combination with the longest survival time is then selected.
            index_select = np.argmax(t_sim)
        d_i_minimum = d_i_total[index_select]

        if mis_specification:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[index_select], t_i_total[index_select]

        X_Strategy1_i = np.append(X_Strategy1_i, X_i_minimum, 1) if tnow == 0 else np.append(X_Strategy1_i, X_i_minimum[:, 1:], 1)
        t_Strategy1_i = np.append(t_Strategy1_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_Strategy1_i, t_i_minimum[1:] + tnow)
        d_Strategy1_i = np.append(d_Strategy1_i, d_i_minimum, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update X0_i and tnow.
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_Strategy1_i, X_Strategy1_i, d_Strategy1_i


# Define Strategy 2.
def DPM_strategy_2(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, Strategy2threshold,
                   misspecification_ofdecision, misspecification_oftrue, mis_PAR):
    # Strategy 2: Minimize the risk of incurable cells developing unless there is an immediate threat of mortality.
    # Strategy 2.1: threshold is 1e9
    # Strategy 2.2: threshold is 1e11
    maxthreshold = Limit_mortality
    if misspecification_oftrue:
        X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(mis_PAR), DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)
        g0_true_i, Sa_true_i, T_true_i = DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    else:
        X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
        g0_true_i, Sa_true_i, T_true_i = None, None, None

    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        misspecification_ofdecision else (None, None, None)

    tnow = 0.0
    t_Strategy2_i = np.zeros([0], dtype=float)
    X_Strategy2_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy2_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    while tnow < Simduration:
        X_i_total, t_i_total, d_i_total, X_i_end_total, X_i_end_multi_resis, t_sim = [], [], [], [], [], []
        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        for i_dose_combination in dose_combination:
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
            d_i_total.append(d_i)
            if misspecification_ofdecision:
                mis_X_i, mis_t_i, _ = DPM_sim_model(X0_i, mis_T_i, mis_g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(mis_X_i)
                X_i_end_total.append(mis_X_i[:, -1].sum())
                X_i_end_multi_resis.append(mis_X_i[-1, -1])
                t_sim.append(mis_t_i[-1])
            elif misspecification_oftrue:
                mis_X_i, mis_t_i, _ = DPM_sim_model(X0_i, T_true_i, g0_true_i, Sa_true_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(mis_X_i)
                X_i_end_total.append(mis_X_i[:, -1].sum())
                X_i_end_multi_resis.append(mis_X_i[-1, -1])
                t_sim.append(mis_t_i[-1])
            else:
                X_i, t__i, d__i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, True)
                X_i_total.append(X_i)
                t_i_total.append(t__i)
                d_i_total[-1] = d__i
                X_i_end_total.append(X_i[:, -1].sum())
                X_i_end_multi_resis.append(X_i[-1, -1])
                t_sim.append(t__i[-1])
        # If one or more doses cause mortality while others do not, exclude the dose(s) that result in mortality.
        assert max(t_sim) <= Stepsize
        if min(X_i_end_total) < Limit_mortality <= max(X_i_end_total):
            ind = [i for i, x in enumerate(X_i_end_total) if x < Limit_mortality]
            X_i_total = [X_i_total[i] for i in ind]
            t_i_total = [t_i_total[i] for i in ind]
            d_i_total = [d_i_total[i] for i in ind]
            X_i_end_total = [X_i_end_total[i] for i in ind]
            X_i_end_multi_resis = [X_i_end_multi_resis[i] for i in ind]
            t_sim = [t_sim[i] for i in ind]
            assert min(X_i_end_total) < Limit_mortality

        if min(X_i_end_total) < Limit_mortality:
            # If the current total cell population does not exceed the threshold, minimize the multiply-resistant population.
            if X0_i.sum() <= Strategy2threshold:
                index_minimum_end_multi_resis = int(np.argmin(X_i_end_multi_resis))
                index_same_minimum_end_multi_resis = [i for i, x in enumerate(X_i_end_multi_resis) if x ==
                                                      X_i_end_multi_resis[index_minimum_end_multi_resis]]
                # If there is only one d_i that gives the minimum multiply-resistant population.
                if len(index_same_minimum_end_multi_resis) == 1:
                    index_select = index_minimum_end_multi_resis
                # If there are multiple d_i that give the same minimum multiply-resistant population.
                else:
                    X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                    # Find the d_i gives the minimum total cell population in all the d_i giving the same minimum multiply-resistant population.
                    index_minimum_end_total_same_minimum_end_multi_resis = int(np.argmin(X_i_end_total_same_minimum_end_multi_resis))
                    index_select = index_same_minimum_end_multi_resis[index_minimum_end_total_same_minimum_end_multi_resis]
            # If the current total cell population exceeds the threshold, minimize the total population.
            else:
                index_select = np.argmin(X_i_end_total)
        # If all dose options result in mortality, select the dose that extends survival the longest.
        else:
            index_select = np.argmax(t_sim)

        d_i_minimum = d_i_total[index_select]
        if misspecification_ofdecision or misspecification_oftrue:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[index_select], t_i_total[index_select]

        X_Strategy2_i = np.append(X_Strategy2_i, X_i_minimum, 1) if tnow == 0 else np.append(X_Strategy2_i, X_i_minimum[:, 1:], 1)
        t_Strategy2_i = np.append(t_Strategy2_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_Strategy2_i, t_i_minimum[1:] + tnow)
        d_Strategy2_i = np.append(d_Strategy2_i, d_i_minimum, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update X0_i and tnow
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_Strategy2_i, X_Strategy2_i, d_Strategy2_i


# Define Strategy 2_LOD.
def DPM_strategy_LOD_2(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, Limit_radiologicdetection,
                       Strategy2threshold, PAR_LOD, subclone_LOD, misspecification_LOD, mutation_rate, mis_specfiy_pop, Strategy_name):
    # Strategy 2: Minimize the risk of incurable cells developing unless there is an immediate threat of mortality.
    # Strategy 2.1: threshold is 1e9
    # Strategy 2.2: threshold is 1e11
    maxthreshold = Limit_mortality
    X0_true_i = DPM_generate_X0(PAR)
    X0_LOD_i = DPM_generate_X0(PAR_LOD)
    g0_i = DPM_generate_g0(PAR)
    Sa_i = DPM_generate_Sa(PAR)
    T_i = DPM_generate_T(PAR)
    assert X0_true_i.sum() > Limit_radiologicdetection

    celltype = ['Spop', 'R1pop', 'R2pop', 'R12pop']
    tnow = 0.0
    t_Strategy2_i = np.zeros([0], dtype=float)
    X_Strategy2_i, X_Strategy2_true_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float), np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy2_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

    flag_treat, subclone_detection = True, False
    while tnow < Simduration:
        if flag_treat:
            X_i_total, t_i_total, d_i_total, X_i_end_total, X_i_end_multi_resis, t_sim = [], [], [], [], [], []
            for i_dose_combination in dose_combination:
                t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
                treat_i = np.array([i_dose_combination], dtype=float).T
                d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
                d_i_total.append(d_i)
                X_i, t_i, d_i = DPM_sim_model(X0_LOD_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, False)   # True
                X_i_total.append(X_i)
                t_i_total.append(t_i)
                if not d_i.any():
                    d_i = np.array(treat_i)
                d_i_total[-1] = d_i
                X_i_end_total.append(X_i[:, -1].sum())
                X_i_end_multi_resis.append(X_i[-1, -1])
                t_sim.append(t_i[-1])
            # If the current total cell population does not exceed the threshold, minimize the multiply-resistant population.
            if X0_true_i.sum() <= Strategy2threshold:
                index_minimum_end_multi_resis = int(np.argmin(X_i_end_multi_resis))
                index_same_minimum_end_multi_resis = [i for i, x in enumerate(X_i_end_multi_resis) if x ==
                                                      X_i_end_multi_resis[index_minimum_end_multi_resis]]
                # If there is only one d_i that gives the minimum multiply-resistant population.
                if len(index_same_minimum_end_multi_resis) == 1:
                    index_select = index_minimum_end_multi_resis
                # If there are multiple d_i that give the same minimum multiply-resistant population.
                else:
                    X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                    # Find the d_i gives the minimum total cell population in all the d_i giving the same minimum multiply-resistant population.
                    index_minimum_end_total_same_minimum_end_multi_resis = int(np.argmin(X_i_end_total_same_minimum_end_multi_resis))
                    index_select = index_same_minimum_end_multi_resis[index_minimum_end_total_same_minimum_end_multi_resis]
            # If the current total cell population exceeds the threshold, minimize the total population.
            else:
                index_select = np.argmin(X_i_end_total)

            t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
            X_i_minimum, t_i_minimum,  d_i_minimum = X_i_total[index_select], t_i_total[index_select], d_i_total[index_select]
            X_true_i, t_true_i, _ = DPM_sim_model(X0_true_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, True)

            X_Strategy2_i, X_Strategy2_true_i, t_Strategy2_i, d_Strategy2_i, X_i_minimum, flag_break, subclone_detection, flag_treat, \
                mis_specfiy_pop = DPM_strategy_LOD(X_i_minimum, t_i_minimum,  d_i_minimum, X_true_i, t_true_i, tnow, X_Strategy2_i,
                                                   X_Strategy2_true_i, t_Strategy2_i, d_Strategy2_i, Stepsize, Limit_mortality,
                                                   Limit_radiologicdetection, subclone_LOD, misspecification_LOD, mutation_rate,
                                                   celltype, mis_specfiy_pop, Strategy_name)
            if flag_break:
                break

            # Update X0_LOD_i, X0_true_i and tnow
            X0_LOD_i = X_i_minimum[:, -1]
            X0_true_i = X_true_i[:, -1]
        else:
            # no treat, let proliferate
            X0_true_i, X_Strategy2_i, X_Strategy2_true_i, t_Strategy2_i, d_Strategy2_i, flag_treat, flag_break, mis_specfiy_pop = \
                DPM_strategy_pro(X0_true_i, g0_i, Sa_i, T_i, PAR['Num_drug'], tnow, X_Strategy2_i, X_Strategy2_true_i, t_Strategy2_i,
                                 d_Strategy2_i, Stepsize, maxthreshold, Simtimestep, Limit_mortality, Limit_radiologicdetection, subclone_LOD,
                                 misspecification_LOD, mutation_rate, celltype, mis_specfiy_pop, Strategy_name)
            if flag_break:
                break

            if flag_treat:
                X0_LOD_i = X_Strategy2_i[:, -1]

        tnow += Stepsize

    assert t_Strategy2_i.shape[0] == X_Strategy2_true_i.shape[1] == X_Strategy2_i.shape[1]
    diffpts = DPM_miscellaneous_treatment_change_time(d_Strategy2_i)/Stepsize
    assert np.all(np.equal(np.mod(diffpts, 1), 0))
    return t_Strategy2_i, X_Strategy2_true_i, d_Strategy2_i, X_Strategy2_i, mis_specfiy_pop


# Define Strategy 3.
def DPM_strategy_3(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # Strategy 3: Minimize the predicted total cell population unless the first multiply-resistant cell will arise by the selection of the d_i
    # which gives the minimum total cell population.

    # At each Stepsize:
    # If the predicted multiply-resistant population < 1, or it is curable, select d_i to minimize the total cell population.
    # If the selected d_i rises the first multiply-resistant cell, re-select d_i to minimize the multiply-resistant population.
    # else if the current multiply-resistant >= 1 and multiply-resistant is not curable, minimize the total cell population.
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy3_i = np.zeros([0], dtype=float)
    X_Strategy3_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy3_i = np.zeros([PAR['Num_drug'], 0], dtype=float)

    t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    while tnow < Simduration:
        X_i_total, t_i_total, d_i_total, X_i_end_total, X_i_end_multi_resis, t_sim = [], [], [], [], [], []
        for i_dose_combination in dose_combination:
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
            d_i_total.append(d_i)
            if mis_specification:
                mis_X_i, mis_t_i, _ = DPM_sim_model(X0_i, mis_T_i, mis_g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(mis_X_i)
                t_i_total.append(mis_t_i)
                X_i_end_total.append(mis_X_i[:, -1].sum())
                X_i_end_multi_resis.append(mis_X_i[-1, -1])
                t_sim.append(mis_t_i[-1])
            else:
                X_i, t__i, d__i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_i_total.append(X_i)
                t_i_total.append(t__i)
                d_i_total[-1] = d__i
                X_i_end_total.append(X_i[:, -1].sum())
                X_i_end_multi_resis.append(X_i[-1, -1])
                t_sim.append(t__i[-1])

        # If the current multiply-resistant < 1 or the multiply-resistant is curable (curable means any Sa[-1, :] >= g0),
        # minimize the total cell population.
        index_minimum_end_total = int(np.argmin(X_i_end_total))
        index_select = index_minimum_end_total
        if (X0_i[-1] < 1 or (np.any(np.greater_equal(mis_Sa_i[-1, :], mis_g0_i[-1]))) if mis_specification
                else np.any(np.greater_equal(Sa_i[-1, :], g0_i[-1]))):
            X_i_minimum = X_i_total[index_minimum_end_total]
            # If the first multiply-resistant cell will arise under the selected d_i, minimize the multiply-resistant population.
            if X0_i[-1] < 1 <= X_i_minimum[-1, -1]:
                index_minimum_end_multi_resis = int(np.argmin(X_i_end_multi_resis))
                index_same_minimum_end_multi_resis = [i for i, x in enumerate(X_i_end_multi_resis) if x ==
                                                      X_i_end_multi_resis[index_minimum_end_multi_resis]]
                # If there is only 1 d_i that gives the minimum multiply-resistant population.
                if len(index_same_minimum_end_multi_resis) == 1:
                    index_select = index_minimum_end_multi_resis
                # If there are multiple d_i that give the same minimum multiply-resistant populaiton, minimize the total cell population in all
                # the d_i giving the same minimum multiply-resistant populaiton.
                else:
                    X_i_end_total_same_minimum_end_multi_resis = [X_i_end_total[i] for i in index_same_minimum_end_multi_resis]
                    # Find the d_i gives the minimum total cell population in all the d_i giving the same minimum multiply-resistant population.
                    index_minimum_end_total_same_minimum_end_multi_resis = int(np.argmin(X_i_end_total_same_minimum_end_multi_resis))
                    index_select = index_same_minimum_end_multi_resis[index_minimum_end_total_same_minimum_end_multi_resis]
                # If minimizing the multiply-resistant population still results in mortality, minimize the total cell population instead.
                if X_i_end_total[index_select] >= Limit_mortality:
                    index_select = index_minimum_end_total
        # If current multiply-resistant >= 1 and multiply-resistant is not curable, minimize the total cell population.
        else:
            pass
        # If all dose combinations result in mortality, select the one with the longest survival time.
        if X_i_end_total[index_select] >= Limit_mortality:
            index_select = np.argmax(t_sim)
        d_i_minimum = d_i_total[index_select]
        if mis_specification:
            X_i_minimum, t_i_minimum, d_i_minimum = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i_minimum, t_i, maxthreshold, LSsim)
        else:
            X_i_minimum, t_i_minimum = X_i_total[index_select], t_i_total[index_select]

        X_Strategy3_i = np.append(X_Strategy3_i, X_i_minimum, 1) if tnow == 0 else np.append(X_Strategy3_i, X_i_minimum[:, 1:], 1)
        t_Strategy3_i = np.append(t_Strategy3_i, t_i_minimum + tnow) if tnow == 0 else np.append(t_Strategy3_i, t_i_minimum[1:] + tnow)
        d_Strategy3_i = np.append(d_Strategy3_i, d_i_minimum, 1)

        # If total cell populatin is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i_minimum[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i_minimum[:, -1] < 1):
            break
        # Update X0_i and tnow.
        X0_i = X_i_minimum[:, -1]
        tnow += Stepsize
    return t_Strategy3_i, X_Strategy3_i, d_Strategy3_i


# Define Strategy 4.
def DPM_strategy_4(PAR, dose_combination, Simduration, Stepsize, Simtimestep, Limit_mortality, LSsim, mis_specification, mis_PAR):
    # Strategy 4: Estimate the time to either incurability or death, and react to the most proximal threat as long as there is a chance of cure.
    # In each Stepsize, evaluate the predicted durations toward incurability (multiply-resistant population >= 1) and
    # mortality (population >= 1e13) dictated by the growth of S, R1, R2, ... and multiply-resistant populations.
    # For each dosage combination d_i, define τ_inc(d) as the predicted time to incurability (multiply-resistant >= 1), given the currently
    # observed population and d_i fixed.
    # Define τ_S(d) as the predicted time to S causing mortality (S > 1e13), given the currently observed population and d_i fixed.
    # τ_R1(d), R1 causing mortality (R1 > 1e13).
    # τ_R2(d), R2 causing mortality (R2 > 1e13).
    # ...
    # τ_X(d), X type causing mortality (X > 1e13).
    # τ_multi_resis(d), multiply-resistant population causing mortality (multiply_resistant > 1e13).
    # If the current multiply-resistant population < 1 or multiply-resistant population is curable, i.e., there exists some d_i such that each
    # component of diag(Sa*d) > g0, vary d to maximize min(τ_inc,τ_s,τ_R1,τ_R2,...,τ_multi_resis) with the constraint that
    # min(τ_S,τ_R1,τ_R2,...,τ_multi_resis) > Stepsize. If such a dosage combination does not exist, maximize min(τ_S,τ_R1,τ_R2,..,τ_multi_resis).
    # If the current R_multi_resis >= 1, and multiply-resistant population is not curable, maximize min(τ_S,τ_R1,τ_R2,...,τ_multi_resis).
    maxthreshold = Limit_mortality
    X0_i, g0_i, Sa_i, T_i = DPM_generate_X0(PAR), DPM_generate_g0(PAR), DPM_generate_Sa(PAR), DPM_generate_T(PAR)
    mis_g0_i, mis_Sa_i, mis_T_i = (DPM_generate_g0(mis_PAR), DPM_generate_Sa(mis_PAR), DPM_generate_T(mis_PAR)) if \
        mis_specification else (None, None, None)

    tnow = 0.0
    t_Strategy4_i = np.zeros([0], dtype=float)
    X_Strategy4_i = np.zeros([PAR['Num_cell_type'], 0], dtype=float)
    d_Strategy4_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
    curable_multi_resis = False

    # If the current multiply-resistant populations < 1 or multiply-resistant populations is curable (curable means any Sa[-1, :] >= g0).
    if np.any(np.greater_equal(mis_Sa_i[-1, :], mis_g0_i[-1])) if mis_specification else np.any(np.greater_equal(Sa_i[-1, :], g0_i[-1])):
        curable_multi_resis = True

    while tnow < Simduration:
        τ_total = np.zeros([PAR['Num_cell_type'] + 1, 0], dtype=float)
        d_τ_i = np.zeros([PAR['Num_drug'], 0], dtype=float)
        timeleft = Simduration - tnow
        t_i = np.arange(0, timeleft + Simtimestep, Simtimestep)
        for i_dose_combination in dose_combination:
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))
            τ_i_drug = np.zeros(PAR['Num_cell_type'] + 1, dtype=float)
            if mis_specification:
                mis_X_τ_i, _, _ = DPM_sim_model(X0_i, mis_T_i, mis_g0_i, mis_Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_τ_i_use = mis_X_τ_i
            else:
                X_τ_i, _, _ = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)
                X_τ_i_use = X_τ_i

            # Calculate the τ_inc,τ_s,τ_R1,τ_R2,...,τ_multi_resis.
            index_celltype = [-1]
            index_celltype.extend(list(range(PAR['Num_cell_type'])))
            limit_celltype = [1]
            limit_celltype.extend([Limit_mortality] * PAR['Num_cell_type'])
            for i in range(len(index_celltype)):
                occurrences_τ, = np.where(X_τ_i_use[index_celltype[i], :] >= limit_celltype[i])
                if not occurrences_τ.size == 0:
                    τ_i_drug[i] = occurrences_τ[0]
                else:
                    τ_i_drug[i] = timeleft

            τ_total = np.append(τ_total, np.reshape(τ_i_drug, (τ_i_drug.shape[0], 1)), 1)
            d_τ_i = np.append(d_τ_i, treat_i, 1)

        τ_S_to_multi_resis = τ_total[1:, :]
        min_τ_S_to_multi_resis = np.amin(τ_S_to_multi_resis, 0)
        index_min_τ_S_to_multi_resis_biggerthanStepsize, = np.where(min_τ_S_to_multi_resis > Stepsize)
        # If the current multiply-resistant population < 1 or it is curable.
        if X0_i[-1] < 1 or curable_multi_resis:
            # If there exists min(τ_S, τ_R1, τ_R2,...,τ_multi_resis) > Stepsize, then maximize min(τ_inc, τ_S, τ_R1, τ_R2, τ_R12) among the d_i
            # that meet the criteria min(τ_S, τ_R1, τ_R2,...,τ_multi_resis) > Stepsize.
            if index_min_τ_S_to_multi_resis_biggerthanStepsize.size:
                τ_total = τ_total[:, index_min_τ_S_to_multi_resis_biggerthanStepsize]
                d_τ_i = d_τ_i[:, index_min_τ_S_to_multi_resis_biggerthanStepsize]
                min_τ_inc_S_to_multi_resis = np.amin(τ_total, 0)
                index_maxmin_τ_inc_S_to_multi_resis = np.argmax(min_τ_inc_S_to_multi_resis)
                d_τ = np.array([d_τ_i[:, index_maxmin_τ_inc_S_to_multi_resis]]).T
            # If there doesn't exist exists min(τ_S, τ_R1, τ_R2,...,τ_multi_resis) > Stepsize,
            # then maximize min(τ_S, τ_R1, τ_R2,...,τ_multi_resis).
            else:
                index_maxmin_τ_S_to_multi_resis = np.argmax(min_τ_S_to_multi_resis)
                d_τ = np.array([d_τ_i[:, index_maxmin_τ_S_to_multi_resis]]).T
        # If the current multiply-resistant population >= 1 and it is not curable, maximize min(τ_S, τ_R1, τ_R2,...,τ_multi_resis).
        elif X0_i[-1] >= 1 and (not curable_multi_resis):
            index_maxmin_τ_S_to_multi_resis = np.argmax(min_τ_S_to_multi_resis)
            d_τ = np.array([d_τ_i[:, index_maxmin_τ_S_to_multi_resis]]).T
        else:
            raise ValueError('Wrong situation')

        t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
        d_i = np.tile(d_τ, (1, t_i.shape[0] - 1))
        X_i, t_i, d_i = DPM_sim_model(X0_i, T_i, g0_i, Sa_i, d_i, t_i, maxthreshold, LSsim)

        X_Strategy4_i = np.append(X_Strategy4_i, X_i, 1) if tnow == 0 else np.append(X_Strategy4_i, X_i[:, 1:], 1)
        t_Strategy4_i = np.append(t_Strategy4_i, t_i + tnow) if tnow == 0 else np.append(t_Strategy4_i, t_i[1:] + tnow)
        d_Strategy4_i = np.append(d_Strategy4_i, d_i, 1)

        # If total cell population is bigger or equal than Limit_mortality, stop. Mortality happens.
        if X_i[:, -1].sum() >= Limit_mortality:
            break
        # If the cell number of all types are smaller than 1, stop. Cured.
        elif all(X_i[:, -1] < 1):
            break
        # Update X0_i and tnow.
        X0_i = X_i[:, -1]
        tnow += Stepsize
    return t_Strategy4_i, X_Strategy4_i, d_Strategy4_i
