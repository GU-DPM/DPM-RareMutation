from DPM_lib import np, groupby, re, deepcopy
from DPM_constant import PARFORMAT_LIST_2DRUG
""" This script contains the miscellaneous functions used in the DPM (Dynamic Precision Medicine) model. """


def DPM_miscellaneous_treatment_change_time(d: np.ndarray):
    d_diff = np.diff(d)
    d_diff_boolen = d_diff.astype(bool)
    if d_diff_boolen.ndim != 1:
        pos_changedose = np.sum(d_diff_boolen, axis=0).astype(bool)
    else:
        pos_changedose = d_diff_boolen
    pos_changedose, = np.where(pos_changedose)
    # Position where drug dose changes.
    pos_changedose += 1
    # At the beginning, index=0, drug dose always changes. It can be seen as the start of applying treatment, previous drug doses are all 0.
    pos_changedose = np.insert(pos_changedose, 0, 0)
    return pos_changedose


def DPM_miscellaneous_allequal(data):
    g = groupby(data)
    return next(g, True) and not next(g, False)


def DPM_miscellaneous_slicedosage(dosage, num_dosage, strategy):
    dosage_sliced = deepcopy(dosage)
    for i_strategy in strategy:
        i_dosage_sliced = dosage_sliced[i_strategy]
        for i, i_dose in enumerate(i_dosage_sliced):
            ind_object = re.finditer(pattern=';', string=i_dose)
            ind = [index.start() for index in ind_object]
            dosage_sliced[i_strategy][i] = i_dose[:ind[num_dosage - 1]]
    return dosage_sliced


def DPM_miscellaneous_fillful(par):
    keys = par.keys()
    if 'g0_R1' not in keys:
        par['g0_R1'] = par['g0_S']
    if 'g0_R2' not in keys:
        par['g0_R2'] = par['g0_S']
    if 'g0_R12' not in keys:
        par['g0_R12'] = par['g0_S']
    if 'Sa.R1.D2.' not in keys:
        par['Sa.R1.D2.'] = par['Sa.S.D2.']
    if 'Sa.R2.D1.' not in keys:
        par['Sa.R2.D1.'] = par['Sa.S.D1.']
    if 'Sa.R12.D1.' not in keys:
        par['Sa.R12.D1.'] = par['Sa.R1.D1.']
    if 'Sa.R12.D2.' not in keys:
        par['Sa.R12.D2.'] = par['Sa.R2.D2.']
    if 'T.S..S.' not in keys:
        par['T.S..S.'] = 0
    if 'T.S..R1.' not in keys:
        par['T.S..R1.'] = 0
    if 'T.S..R2.' not in keys:
        par['T.S..R2.'] = 0
    if 'T.S..R12.' not in keys:
        par['T.S..R12.'] = 0
    if 'T.R1..R1.' not in keys:
        par['T.R1..R1.'] = 0
    if 'T.R1..R2.' not in keys:
        par['T.R1..R2.'] = 0
    if 'T.R1..R12.' not in keys:
        par['T.R1..R12.'] = 0
    if 'T.R2..R1.' not in keys:
        par['T.R2..R1.'] = 0
    if 'T.R2..R2.' not in keys:
        par['T.R2..R2.'] = 0
    if 'T.R2..R12.' not in keys:
        par['T.R2..R12.'] = 0
    if 'T.R12..S.' not in keys:
        par['T.R12..S.'] = 0
    if 'T.R12..R1.' not in keys:
        par['T.R12..R1.'] = par['T.R2..S.']
    if 'T.R12..R2.' not in keys:
        par['T.R12..R2.'] = par['T.R1..S.']
    if 'T.R12..R12.' not in keys:
        par['T.R12..R12.'] = 0
    return par
