from DPM_lib import datetime, os
from DPM_generate import *
from DPM_model import *
from DPM_strategy import DPM_strategy_0
""" This script assigns and checks arguments used in the DPM (Dynamic Precision Medicine) model. """

# Check simulation parameter.
# (i) Total cell number (X0total) should be smaller than Limit_mortality at t = 0, otherwise mortality happens at t = 0.
# (ii) Simduration should be larger than stepsize, otherwise only one step simulation.
# Return the information whether the parameter fullfill the conditions.


def DPM_check_par(X0total, Limit_mortality, Simduration, Stepsize):
    Par_check_index = []
    # (1) X0total >= Limit_mortality. Mortality happens at t = 0.
    if X0total >= Limit_mortality:
        Par_check_index.append('X0total >= Limit_mortality')
    # (2) Stepsize > Simuduration. Single stepsize is larger than the whole simulation period and cannot finish one stepsize.
    if Stepsize > Simduration:
        Par_check_index.append('Stepsize > Simuduration')
    return Par_check_index


# Check whether the type of one par is dictionary.
def DPM_check_dict_par(par):
    par_text = Colored.BOLD + Colored.PURPLE + 'par' + Colored.END
    print_quotation = True
    flag_inlist = True
    beginning_char = '\n'
    if type(par) is not dict:
        print('The type of the inputted ' + par_text + ' should be dictionary.')
        DPM_print_errorandclose()
        return
    keys = par.keys()
    keys_missing = set()
    if 'Num_drug' not in keys:
        Num_drug_text = Colored.BOLD + Colored.PURPLE + 'Num_drug' + Colored.END
        print(Num_drug_text + ' attribute is missing in the input ' + par_text + '.')
        DPM_print_errorandclose()
        return
    Num_drug = DPM_check_one_numerical_par(par['Num_drug'], 'Num_drug', '{:g}'.format, NUM_DRUG_UPPERTHRES, NUM_DRUG_LOWERTHRES, flag_inlist,
                                           NUM_DRUG_LIST, print_quotation, beginning_char)
    if Num_drug is None:
        return
    elif Num_drug == 2:
        keys_missing = set(PARFORMAT_LIST_2DRUG[1:2] + PARFORMAT_LIST_2DRUG[3:]).difference(set(keys))

    if not keys_missing == set():
        if len(keys_missing) == 1:
            print('The following attribute is missing in the input ' + par_text + ':')
        else:
            print('The following attributes are missing in the input ' + par_text + ':')
        DPM_print_list(list(keys_missing), print_quotation)
        DPM_print_errorandclose()
        return
    else:
        return par


# Check whether one numerical parameter is in its range or in a list if it should be. Return back the parameter if it is, otherwise return None.
def DPM_check_one_numerical_par(par, par_name, par_display, upper_limit, lower_limit, flag_inlist, lis=None, print_quotation=None, beginning_char=''):
    if type(par) not in [int, float]:
        DPM_print_notnumerical(par_name)
        DPM_print_errorandclose()
        return
    if flag_inlist:
        par = DPM_check_in_list_par([par], par_name, lis, print_quotation, beginning_char)
        if par is None:
            return
        else:
            return par[0]
    else:
        if par > upper_limit:
            print('The inputted ' + Colored.BOLD + Colored.PURPLE + par_name + Colored.END + ' = ' + Colored.BOLD + Colored.PURPLE +
                  par_display(par) + Colored.END + ' is larger than the upper limit ' + Colored.BOLD + Colored.PURPLE +
                  par_display(upper_limit) + Colored.END + '.')
            DPM_print_errorandclose()
            return
        elif par < lower_limit:
            print('The inputted ' + Colored.BOLD + Colored.PURPLE + par_name + Colored.END + ' = ' + Colored.BOLD + Colored.PURPLE +
                  par_display(par) + Colored.END + ' is lower than the lower limit ' + Colored.BOLD + Colored.PURPLE +
                  par_display(lower_limit) + Colored.END + '.')
            DPM_print_errorandclose()
            return
        else:
            return par


# Check whether the parameter is a string or in a list if it should be. Return back the parameter if it is, otherwise return None.
def DPM_check_one_string_par(par, par_name, flag_inlist, lis, print_quotation, beginning_char):
    if type(par) is not str:
        DPM_print_notstring(par_name)
        DPM_print_errorandclose()
        return
    if flag_inlist:
        par = DPM_check_in_list_par([par], par_name, lis, print_quotation, beginning_char)
        if par is None:
            return
        else:
            return par[0]
    else:
        return par


# Check whether the parameter is a boolean. Return back the parameter if it is, otherwise return None.
def DPM_check_one_boolean_par(par, par_name):
    if type(par) is not bool:
        DPM_print_notboolean(par_name)
        DPM_print_errorandclose()
        return
    else:
        return par


# Check whether the input parameter is in the input list. Return back the parameter if it is in, otherwise return None.
def DPM_check_in_list_par(par, par_name, lis, print_quotation, beginning_char=''):
    if not (all(elem in lis for elem in par)):
        print('The allowed input of ' + Colored.BOLD + Colored.PURPLE + par_name + Colored.END + ': ', end='')
        DPM_print_list(lis, print_quotation)
        element_not_in_list = [elem for elem in par if elem not in lis]
        print('\n', end='')
        print('The following input of ' + Colored.BOLD + Colored.PURPLE + par_name + Colored.END + ' not allowed: ', end='')
        DPM_print_list(element_not_in_list, print_quotation)
        DPM_print_errorandclose(beginning_char)
        return
    else:
        return par


# Check a list of numerical tuple, used for dose combination input.
# (i) Number of elements in each tuple it equal to the required number.
# (ii) Each element in the tuple can be convert to a number.
# Return the list of valid inputs if existed, otherwise return None.
def DPM_check_list_of_tuple_of_numerical_par(par, par_name, num_element_in_tuple, par_element_lowerthres, par_element_upperthres, par_sum_lowerthres,
                                             par_sum_upperthres, print_quotation):
    text_input = Colored.BOLD + Colored.PURPLE + par_name + Colored.END
    valid_input, invalid_input = [], []
    if par:
        for par_tuple in par:
            # Number of element in tuple is equal to the required number.
            if type(par_tuple) is tuple and len(par_tuple) == num_element_in_tuple:
                try:
                    par_tuple = tuple(float(i) for i in par_tuple)
                    if all(par_element_lowerthres <= i <= par_element_upperthres for i in par_tuple) \
                            and par_sum_lowerthres <= sum(par_tuple) <= par_sum_upperthres:
                        valid_input.append(par_tuple)
                    else:
                        invalid_input.append(par_tuple)
                except ValueError:
                    invalid_input.append(par_tuple)
            else:
                invalid_input.append(par_tuple)

        if valid_input:
            if len(valid_input) > 1:
                print('The following ' + text_input + ' are valid inputs and accepted:')
            else:
                print('The following ' + text_input + ' is valid input and accepted:')
            DPM_print_list(valid_input, print_quotation)

        if invalid_input:
            if len(invalid_input) > 1:
                print('\nThe following ' + text_input + ' are invalid inputs and ignored:')
            else:
                print('\nThe following ' + text_input + ' is invalid input and ignored:')
            DPM_print_list(invalid_input, print_quotation)
            print('\n')
        return valid_input
    else:
        return


# Check the parameters in the PARNAME_LIST.
def DPM_check_inputpar(kwargs, parname_list, caller=None):
    parname_list = [i_name.lower() for i_name in parname_list]

    print_quotation, flag_inlist, beginning_char = True, True, '\n'

    # Check parameter Num_drug should be numerical and in its range.
    if 'Num_drug'.lower() in parname_list:
        if 'Num_drug'.lower() in kwargs.keys():
            Num_drug = kwargs['Num_drug'.lower()]
            Num_drug = DPM_check_one_numerical_par(Num_drug, 'Num_drug', '{:g}'.format, NUM_DRUG_UPPERTHRES, NUM_DRUG_LOWERTHRES, flag_inlist,
                                                   NUM_DRUG_LIST, print_quotation, beginning_char)
            if Num_drug is None:
                sys.exit()
        else:
            Num_drug = NUM_DRUG_DEFAULT_VAL
    else:
        Num_drug = ''

    # Check parameter par is inputted or not.
    if 'par'.lower() in parname_list:
        if 'par'.lower() in kwargs.keys():
            par = kwargs['par']
            if 'Num_drug' not in par:
                par['Num_drug'] = Num_drug
            par = DPM_check_dict_par(par)
            if par is None:
                sys.exit()
            Num_drug = par['Num_drug']
        else:
            par_text = Colored.BOLD + Colored.PURPLE + 'par' + Colored.END
            print('No ' + par_text + ' keyword argument input. Need input ' + par_text + ' to run the function.')
            DPM_print_errorandclose()
            sys.exit()
    else:
        par = ''

    # Check parameter dose_method should be a string, and be either continous or discrete.
    if 'dose_method'.lower() in parname_list:
        if 'dose_method'.lower() in kwargs.keys():
            dose_method = kwargs['dose_method'.lower()]
            dose_method = DPM_check_one_string_par(dose_method, 'dose_method', flag_inlist, DOSE_METHOD_LIST, print_quotation, beginning_char)
            if dose_method is None:
                sys.exit()
            elif dose_method in ['continous', 'c']:
                dose_method = 'c'
            else:
                dose_method = 'd'
        else:
            dose_method = DOSE_METHOD_DEFAULT_VAL
    else:
        dose_method = ''

    # Check parameter dose_interval should be numerical and in its range.
    if 'dose_interval'.lower() in parname_list:
        if dose_method == 'c':
            if 'dose_interval'.lower() in kwargs.keys():
                dose_interval = kwargs['dose_interval'.lower()]
                dose_interval = DPM_check_one_numerical_par(dose_interval, 'dose_interval', '{:g}'.format, DOSE_INTERVAL_UPPERTHRES,
                                                            DOSE_INTERVAL_LOWERTHRES, flag_inlist, DOSE_INTERVAL_LIST, print_quotation,
                                                            beginning_char)
                if dose_interval is None:
                    sys.exit()
            else:
                dose_interval = DOSE_INTERVAL_DEFAULT_VAL
        else:
            dose_interval = ''
    else:
        dose_interval = ''

    # Generate default parameter dose_combination.
    if dose_method == 'd' and Num_drug != '':
        dose_combination = DPM_generate_discrete_dose_combination(Num_drug)
    elif dose_method == 'c':
        dose_combination = DPM_generate_continuous_dose_combination(Num_drug, dose_interval)
    else:
        dose_combination = []

    # Check parameter use_input_dose_combination_only should be a boolean.
    if 'use_input_dose_combination_only'.lower() in parname_list:
        if 'use_input_dose_combination_only'.lower() in kwargs.keys():
            use_input_dose_combination_only = kwargs['use_input_dose_combination_only'.lower()]
            use_input_dose_combination_only = DPM_check_one_boolean_par(use_input_dose_combination_only, 'use_input_dose_combination_only')
            if use_input_dose_combination_only is None:
                sys.exit()
        else:
            use_input_dose_combination_only = False
    else:
        use_input_dose_combination_only = ''

    # Check parameter dose_combination should be a list of tuple and in its range.
    if 'dose_combination'.lower() in parname_list:
        if 'dose_combination'.lower() in kwargs.keys():
            dose_combination_input = kwargs['dose_combination'.lower()]
            if type(dose_combination_input) is not list:
                dose_combination_input = [dose_combination_input]
            dose_combination_input = DPM_check_list_of_tuple_of_numerical_par(
                dose_combination_input, 'dose_combination', Num_drug, DOSE_COMBINATION_ELEMENT_LOWERTHRES, DOSE_COMBINATION_ELEMENT_UPPERTHRES,
                DOSE_COMBINATION_SUM_LOWERTHRES, DOSE_COMBINATION_SUM_UPPERTHRES, print_quotation)
            if dose_combination_input:
                if not use_input_dose_combination_only:
                    dose_combination_input = [x for x in dose_combination_input if x not in dose_combination]
                    dose_combination.extend(dose_combination_input)
                elif use_input_dose_combination_only:
                    dose_combination = dose_combination_input
    else:
        dose_combination = ''

    flag_inlist = False

    # Check parameter Simduration should be numerical and in its range.
    if 'Simduration'.lower() in parname_list:
        if 'Simduration'.lower() in kwargs.keys():
            Simduration = kwargs['Simduration'.lower()]
            Simduration = DPM_check_one_numerical_par(Simduration, 'Simduration', '{:g}'.format, SIMDURATION_UPPERTHRES,
                                                      SIMDURATION_LOWERTHRES, flag_inlist)
            if Simduration is None:
                sys.exit()
        else:
            Simduration = SIMDURATION_DEFAULT_VAL
    else:
        Simduration = ''

    # Check parameter Limit_mortality should be numerical and in its range.
    if 'Limit_mortality'.lower() in parname_list:
        if 'Limit_mortality'.lower() in kwargs.keys():
            Limit_mortality = kwargs['Limit_mortality'.lower()]
            Limit_mortality = DPM_check_one_numerical_par(Limit_mortality, 'Limit_mortality', '{:.2e}'.format,
                                                          LIMIT_MORTALITY_UPPERTHRES, LIMIT_MORTALITY_LOWERTHRES, flag_inlist)
            if Limit_mortality is None:
                sys.exit()
        else:
            Limit_mortality = LIMIT_MORTALITY_DEFAULT_VAL
    else:
        Limit_mortality = ''

    # Check parameter Stepsize should be numerical and in its range.
    if 'Stepsize'.lower() in parname_list:
        if 'Stepsize'.lower() in kwargs.keys():
            Stepsize = kwargs['Stepsize'.lower()]
            Stepsize = DPM_check_one_numerical_par(Stepsize, 'Stepsize', '{:g}'.format, STEPSIZE_UPPERTHRES,
                                                   STEPSIZE_LOWERTHRES, flag_inlist)
            if Stepsize is None:
                sys.exit()
        else:
            Stepsize = STEPSIZE_DEFAULT_VAL
    else:
        Stepsize = ''

    # Check parameter Limit_radiologicdetection should be numerical and in its range.
    if 'Limit_radiologicdetection'.lower() in parname_list:
        if 'Limit_radiologicdetection'.lower() in kwargs.keys():
            Limit_radiologicdetection = kwargs['Limit_radiologicdetection'.lower()]
            Limit_radiologicdetection = DPM_check_one_numerical_par(Limit_radiologicdetection, 'Limit_radiologicdetection',
                                                                    '{:.2e}'.format, LIMIT_RADIOLOGICDETECTION_UPPERTHRES,
                                                                    LIMIT_RADIOLOGICDETECTION_LOWERTHRES, flag_inlist)
            if Limit_radiologicdetection is None:
                sys.exit()
        else:
            Limit_radiologicdetection = LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL
    else:
        Limit_radiologicdetection = ''

    # Check parameter Limit_moleculardetection should be numerical and in its range.
    if 'Limit_moleculardetection'.lower() in parname_list:
        if 'Limit_moleculardetection'.lower() in kwargs.keys():
            Limit_moleculardetection = kwargs['Limit_moleculardetection'.lower()]
            Limit_moleculardetection = DPM_check_one_numerical_par(Limit_moleculardetection, 'Limit_moleculardetection',
                                                                   '{:.2e}'.format, LIMIT_MOLECULARDETECTION_UPPERTHRES,
                                                                   LIMIT_MOLECULARDETECTION_LOWERTHRES, flag_inlist)
            if Limit_moleculardetection is None:
                sys.exit()
        else:
            Limit_moleculardetection = LIMIT_MOLECULARDETECTION_DEFAULT_VAL
    else:
        Limit_moleculardetection = ''

    # Check parameter run_sim should be a boolean.
    if 'run_sim'.lower() in parname_list:
        if 'run_sim'.lower() in kwargs.keys():
            run_sim = kwargs['run_sim'.lower()]
            run_sim = DPM_check_one_boolean_par(run_sim, 'run_sim')
            if run_sim is None:
                sys.exit()
        else:
            run_sim = True
    else:
        run_sim = ''

    # Check parameter use_parallel should be a boolean.
    if 'use_parallel'.lower() in parname_list:
        if 'use_parallel'.lower() in kwargs.keys():
            use_parallel = kwargs['use_parallel'.lower()]
            use_parallel = DPM_check_one_boolean_par(use_parallel, 'use_parallel')
            if use_parallel is None:
                sys.exit()
        else:
            use_parallel = False
    else:
        use_parallel = ''

    # Check parameter erase_preresult should be a boolean.
    if 'erase_preresult'.lower() in parname_list:
        if 'erase_preresult'.lower() in kwargs.keys():
            erase_preresult = kwargs['erase_preresult'.lower()]
            erase_preresult = DPM_check_one_boolean_par(erase_preresult, 'erase_preresult')
            if erase_preresult is None:
                sys.exit()
        else:
            erase_preresult = True
    else:
        erase_preresult = ''

    # Generate parameter criterisa based on drug number.
    if Num_drug in [2]:
        PAR_criterisa_list = DPM_generate_PAR_criterisa_list(Num_drug)
    else:
        PAR_criterisa_list = []

    # Check parameter misspecification_sim_only should be a boolean.
    if 'misspecification_sim_only'.lower() in parname_list:
        if 'misspecification_sim_only'.lower() in kwargs.keys():
            misspecification_sim_only = kwargs['misspecification_sim_only'.lower()]
            misspecification_sim_only = DPM_check_one_boolean_par(misspecification_sim_only, 'misspecification_sim_only')
            if misspecification_sim_only is None:
                sys.exit()
        else:
            misspecification_sim_only = False
    else:
        misspecification_sim_only = ''

    # Check the elements of parameter PAR_criterisa are in the list.
    if 'PAR_criterisa'.lower() in parname_list:
        if 'PAR_criterisa'.lower() in kwargs.keys():
            PAR_criterisa = kwargs['PAR_criterisa'.lower()]
            if type(PAR_criterisa) is not list:
                PAR_criterisa = [PAR_criterisa]
            PAR_criterisa = DPM_check_in_list_par(PAR_criterisa, 'PAR_criterisa', PAR_criterisa_list, print_quotation, beginning_char)
            if PAR_criterisa is None:
                sys.exit()
            else:
                # Make PAR_criterisa unique.
                PAR_criterisa = list(set(PAR_criterisa))
        else:
            PAR_criterisa = []
    else:
        PAR_criterisa = ''

    # Check parameter Strategy_name should be in the list.
    if 'Strategy_name'.lower() in parname_list:
        if 'Strategy_name'.lower() in kwargs.keys():
            Strategy_name = kwargs['Strategy_name'.lower()]
            if type(Strategy_name) is not list:
                Strategy_name = list(Strategy_name)
            Strategy_name = DPM_check_in_list_par(Strategy_name, 'Strategy_name', STRATEGY_LIST, print_quotation, beginning_char)
            if Strategy_name is None:
                sys.exit()
        else:
            Strategy_name = ['strategy0', 'strategy1', 'strategy2.1', 'strategy2.2', 'strategy3', 'strategy4']
        if 'all' in Strategy_name:
            Strategy_name = STRATEGY_LIST
    else:
        Strategy_name = ''

    # Check parameter filename_csv is inputted or not.
    if 'filename_csv'.lower() in parname_list:
        if 'filename_csv'.lower() in kwargs.keys():
            filename_csv = kwargs['filename_csv']
        else:
            filename_csv_text: str = Colored.BOLD + Colored.PURPLE + 'filename_csv' + Colored.END
            print('No ' + filename_csv_text + ' keyword argument input. Need at least input the name of the ' +
                  Colored.BOLD + Colored.PURPLE + '.csv' + Colored.END + ' file containing the parameters.')
            DPM_print_errorandclose()
            sys.exit()
    else:
        filename_csv = ''

    # Check parameter misspecification_filename_csv should has the same type and length as filename_csv.
    if 'misspecification_filename_csv'.lower() in parname_list:
        if 'misspecification_filename_csv' in kwargs.keys():
            misspecification_filename_csv = kwargs['misspecification_filename_csv']
            if type(misspecification_filename_csv) != type(filename_csv):
                print('The type of the input keyword argument ' + Colored.BOLD + Colored.PURPLE + 'filename_csv' + Colored.END +
                      ' is different from the type of the input keyword argument ' + Colored.BOLD + Colored.PURPLE +
                      'misspecification_filename_csv' + Colored.END + '.')
                DPM_print_errorandclose()
                sys.exit()
            elif type(misspecification_filename_csv) is list and len(misspecification_filename_csv) != len(filename_csv):
                print('The length of the input keyword argument ' + Colored.BOLD + Colored.PURPLE + 'filename_csv' + Colored.END +
                      ' is different from the length of the input keyword argument ' + Colored.BOLD + Colored.PURPLE +
                      'misspecification_filename_csv' + Colored.END + '.')
                DPM_print_errorandclose()
                sys.exit()
        else:
            misspecification_filename_csv = ''
    else:
        misspecification_filename_csv = ''

    # Check parameter fullinput should be a boolean.
    if 'fullinput'.lower() in parname_list:
        if 'fullinput'.lower() in kwargs.keys():
            fullinput = kwargs['fullinput'.lower()]
            fullinput = DPM_check_one_boolean_par(fullinput, 'fullinput')
            if fullinput is None:
                sys.exit()
        else:
            fullinput = True
    else:
        fullinput = ''

    # Check parameter par_save_block_size should be numerical and in its range.
    if 'par_save_block_size'.lower() in parname_list:
        flag_inlist = False
        if 'par_save_block_size'.lower() in kwargs.keys():
            par_save_block_size = kwargs['par_save_block_size'.lower()]
            par_save_block_size = DPM_check_one_numerical_par(par_save_block_size, 'par_save_block_size', '{:.2e}'.format,
                                                              PAR_SAVE_BLOCK_SIZE_UPPERTHRES, PAR_SAVE_BLOCK_SIZE_LOWERTHRES,
                                                              flag_inlist)
            if par_save_block_size is None:
                sys.exit()
        else:
            par_save_block_size = PAR_SAVE_BLOCK_SIZE_DEFAULT_VAL
    else:
        par_save_block_size = ''

    # Check parameter save_filename_param should be a boolean.
    if 'save_filename_param'.lower() in parname_list:
        if 'save_filename_param'.lower() in kwargs.keys():
            save_filename_param = kwargs['save_filename_param'.lower()]
            save_filename_param = DPM_check_one_boolean_par(save_filename_param, 'save_filename_param')
            if save_filename_param is None:
                sys.exit()
        else:
            save_filename_param = True
    else:
        save_filename_param = ''

    # Check parameter save_filename_stopt should be a boolean.
    if 'save_filename_stopt'.lower() in parname_list:
        if 'save_filename_stopt'.lower() in kwargs.keys():
            save_filename_stopt = kwargs['save_filename_stopt'.lower()]
            save_filename_stopt = DPM_check_one_boolean_par(save_filename_stopt, 'save_filename_stopt')
            if save_filename_stopt is None:
                sys.exit()
        else:
            save_filename_stopt = True
    else:
        save_filename_stopt = ''

    # Check parameter save_filename_pop should be a boolean.
    if 'save_filename_pop'.lower() in parname_list:
        if 'save_filename_pop'.lower() in kwargs.keys():
            save_filename_pop = kwargs['save_filename_pop'.lower()]
            save_filename_pop = DPM_check_one_boolean_par(save_filename_pop, 'save_filename_pop')
            if save_filename_pop is None:
                sys.exit()
        else:
            save_filename_pop = True
    else:
        save_filename_pop = ''

    # Check parameter save_filename_dosage should be a boolean.
    if 'save_filename_dosage'.lower() in parname_list:
        if 'save_filename_dosage'.lower() in kwargs.keys():
            save_filename_dosage = kwargs['save_filename_dosage'.lower()]
            save_filename_dosage = DPM_check_one_boolean_par(save_filename_dosage, 'save_filename_dosage')
            if save_filename_dosage is None:
                sys.exit()
        else:
            save_filename_dosage = True
    else:
        save_filename_dosage = ''

    # Check parameter save_filename_eachtimepoint should be a boolean.
    if 'save_filename_eachtimepoint'.lower() in parname_list:
        if 'save_filename_eachtimepoint'.lower() in kwargs.keys():
            save_filename_eachtimepoint = kwargs['save_filename_eachtimepoint'.lower()]
            save_filename_eachtimepoint = DPM_check_one_boolean_par(save_filename_eachtimepoint, 'save_filename_eachtimepoint')
            if save_filename_eachtimepoint is None:
                sys.exit()
        else:
            save_filename_eachtimepoint = False
    else:
        save_filename_eachtimepoint = ''

    flag_inlist = False
    path_current = os.getcwd()
    # Check parameter pathsave should be a vaild directory. If not existed, try to create the directory otherwise use the current path.
    if 'pathsave'.lower() in parname_list:
        if 'pathsave'.lower() in kwargs.keys():
            pathsave = kwargs['pathsave'.lower()]
            pathsave = DPM_check_one_string_par(pathsave, 'pathsave', flag_inlist, [], print_quotation, beginning_char)
            if pathsave is None:
                sys.exit()
            # If the path to save is not exist.
            pathsave = os.path.abspath(pathsave)
            if not os.path.exists(pathsave):
                try:
                    os.makedirs(pathsave)
                except (OSError, FileExistsError):
                    pathsave_text = Colored.BOLD + Colored.PURPLE + pathsave + Colored.END
                    print('The inputted ' + Colored.BOLD + Colored.PURPLE + 'pathsave' + Colored.END + ':')
                    print(pathsave_text + ' is not exist and cannot be created. The current working directory is used to save.')
                    pathsave = path_current
        else:
            pathsave = path_current
        print('The directory for saving results is set to:')
        print(Colored.BOLD + Colored.PURPLE + pathsave + Colored.END)
    else:
        pathsave = ''

    path_current = os.getcwd()
    # Check parameter pathload should be a vaild directory. If not existed, try load data from the current path.
    if 'pathload'.lower() in parname_list:
        if 'pathload'.lower() in kwargs.keys():
            pathload = kwargs['pathload'.lower()]
            pathload = DPM_check_one_string_par(pathload, 'pathload', flag_inlist, [], print_quotation, beginning_char)
            if pathload is None:
                sys.exit()
            # If the path to load data is not exist.
            pathload = os.path.abspath(pathload)
            if not os.path.exists(pathload):
                pathload_text = Colored.BOLD + Colored.PURPLE + pathload + Colored.END
                print('The inputted ' + Colored.BOLD + Colored.PURPLE + 'pathload' + Colored.END + ':')
                print(pathload_text + ' is not exist. The current working directory is used to save.')
                pathload = path_current
        else:
            pathload = path_current
        print('The directory for loading files is set to:')
        print(Colored.BOLD + Colored.PURPLE + pathload + Colored.END)
    else:
        pathload = ''

    # Check parameter pathload should be a vaild directory. If not existed, try load data from the current path.
    if 'misspecification_fileload'.lower() in parname_list:
        if 'misspecification_fileload'.lower() in kwargs.keys():
            misspecification_fileload = kwargs['misspecification_fileload'.lower()]
            misspecification_fileload = DPM_check_one_string_par(misspecification_fileload, 'misspecification_fileload',
                                                                 flag_inlist, [], print_quotation, beginning_char)
            if misspecification_fileload is None:
                sys.exit()
            # If the path to load data is not exist.
            misspecification_fileload = os.path.abspath(misspecification_fileload)
            if not os.path.exists(misspecification_fileload):
                misspecification_fileload_text = Colored.BOLD + Colored.PURPLE + misspecification_fileload + Colored.END
                print('The inputted ' + Colored.BOLD + Colored.PURPLE + 'misspecification_fileload' + Colored.END + ':')
                print(misspecification_fileload_text+ ' is not exist.')
                DPM_print_errorandclose()
                sys.exit()
        else:
            misspecification_fileload = ''
    else:
        misspecification_fileload = ''

    # Check parameter pathloadmis should be a vaild directory. If not existed, try load data from the current path.
    if 'pathloadmis'.lower() in parname_list:
        if 'pathloadmis'.lower() in kwargs.keys():
            pathloadmis = kwargs['pathloadmis'.lower()]
            pathloadmis = DPM_check_one_string_par(pathloadmis, 'pathloadmis', False, [], print_quotation, beginning_char)
            if pathloadmis is None:
                sys.exit()
            # If the path to load data is not exist.
            pathloadmis = os.path.abspath(pathloadmis) if pathloadmis != '' else ''
            if (not os.path.exists(pathloadmis)) and pathloadmis != '':
                pathloadmis_text = Colored.BOLD + Colored.PURPLE + pathloadmis + Colored.END
                print('The inputted ' + Colored.BOLD + Colored.PURPLE + 'pathloadmis' + Colored.END + ':')
                print(pathloadmis_text + ' is not exist')
                pathloadmis = ''
        else:
            pathloadmis = ''
    else:
        pathloadmis = ''

    # Check parameter par_ind should be numerical or list of numerical.
    if 'par_ind'.lower() in parname_list:
        if 'par_ind'.lower() in kwargs.keys():
            par_ind_accept = []
            par_ind = kwargs['par_ind'.lower()]
            # Change par_ind type to list if not.
            par_ind = [par_ind] if type(par_ind) is not list else par_ind
            for i_par_ind in par_ind:
                par_ind = DPM_check_one_numerical_par(i_par_ind, 'par_ind', '{:g}'.format, np.inf, -np.inf, flag_inlist, [],
                                                      print_quotation, beginning_char)
                if par_ind is not None and par_ind % 1 == 0:
                    par_ind_accept.append(i_par_ind)
            if len(par_ind_accept) == 0:
                print('No index avaiable in the inputted ' + Colored.BOLD + Colored.PURPLE + 'par_ind' + Colored.END)
                DPM_print_errorandclose()
                sys.exit()
            used: set = set()
            par_ind = [int(x) for x in par_ind_accept if int(x) not in used and (used.add(int(x)) or True)]
        else:
            par_ind = ''
    else:
        par_ind = ''

    # Check parameter savename should be a string.
    if 'savename'.lower() in parname_list:
        if 'savename'.lower() in kwargs.keys():
            savename = kwargs['savename'.lower()]
            savename = DPM_check_one_string_par(savename, 'savename', flag_inlist, [], print_quotation, beginning_char)
            if savename is None:
                sys.exit()
        else:
            # Use local time as savename, up to milliseconds.
            savename = datetime.now().strftime('%Y-%m-%d-%H-%M-%S-%f')
    else:
        savename = ''

    # Check parameter save_filename_eachtimepoint should be a boolean.
    if 'plot'.lower() in parname_list:
        if 'plot'.lower() in kwargs.keys():
            plot = kwargs['plot'.lower()]
            plot = DPM_check_one_boolean_par(plot, 'plot')
            if plot is None:
                sys.exit()
        else:
            plot = True if caller == 'DPM_run_plot_1PAR' else False
    else:
        plot = ''

    # Check parameter Num_stepdiff should be numerical and in its range.
    if 'Num_stepdiff'.lower() in parname_list:
        if 'Num_stepdiff'.lower() in kwargs.keys():
            Num_stepdiff = kwargs['Num_stepdiff'.lower()]
            Num_stepdiff = DPM_check_one_numerical_par(Num_stepdiff, 'Num_stepdiff', '{:g}'.format, NUM_STEPDIFF_UPPERTHRES,
                                                       NUM_STEPDIFF_LOWERTHRES, flag_inlist)
            if Num_stepdiff is None:
                sys.exit()
        else:
            Num_stepdiff = NUM_STEPDIFF_DEFAULT_VAL
    else:
        Num_stepdiff = ''

    # Check parameter filename_pattern should be string.
    if 'filename_pattern'.lower() in parname_list:
        if 'filename_pattern'.lower() in kwargs.keys():
            filename_pattern = kwargs['filename_pattern'.lower()]
            filename_pattern = DPM_check_one_string_par(filename_pattern, 'filename_pattern', flag_inlist, [], print_quotation, beginning_char)
            if filename_pattern is None:
                sys.exit()
        else:
            filename_pattern = ''
    else:
        filename_pattern = ''

    # Check parameter misspecification_ofdecision should be a boolean.
    if 'misspecification_ofdecision'.lower() in parname_list:
        if 'misspecification_ofdecision'.lower() in kwargs.keys():
            misspecification_ofdecision = kwargs['misspecification_ofdecision'.lower()]
            misspecification_ofdecision = DPM_check_one_boolean_par(misspecification_ofdecision, 'misspecification_ofdesiction')
            if misspecification_ofdecision is None:
                sys.exit()
        else:
            misspecification_ofdecision = False
    else:
        misspecification_ofdecision = ''

    # Check parameter misspecification_atsim should be a boolean.
    if 'misspecification_atsim'.lower() in parname_list:
        if 'misspecification_atsim'.lower() in kwargs.keys():
            misspecification_atsim = kwargs['misspecification_atsim'.lower()]
            misspecification_atsim = DPM_check_one_boolean_par(misspecification_atsim, 'misspecification_atsim')
            if misspecification_atsim is None:
                sys.exit()
        else:
            misspecification_atsim = False
    else:
        misspecification_atsim = ''

    # Check parameter subclone limit of detection should be numerical and in its range.
    if 'subclone_LOD'.lower() in parname_list:
        flag_inlist = False
        if 'subclone_LOD'.lower() in kwargs.keys():
            subclone_LOD = kwargs['subclone_LOD'.lower()]
            if subclone_LOD:
                subclone_LOD = DPM_check_one_numerical_par(subclone_LOD, 'subclone_LOD', '{:.2e}'.format,
                                                           SUBCOLONE_LOD_UPPERTHRES, SUBCOLONE_LOD_LOWERTHRES, flag_inlist)
                if subclone_LOD is None:
                    sys.exit()
            else:
                subclone_LOD = ''
        else:
            subclone_LOD = ''
    else:
        subclone_LOD = ''

    # Check mutation rate should be numerical and in its range.
    if 'mutation_rate'.lower() in parname_list:
        flag_inlist = False
        if 'mutation_rate'.lower() in kwargs.keys():
            mutation_rate = kwargs['mutation_rate'.lower()]
            mutation_rate = DPM_check_one_numerical_par(mutation_rate, 'mutation_rate', '{:.2e}'.format,
                                                        MUTATION_RATE_UPPERTHRES, MUTATION_RATE_LOWERTHRES, flag_inlist)
            if mutation_rate is None:
                sys.exit()
        else:
            mutation_rate = MUTATION_RATE_DEFAULT_VAL
    else:
        mutation_rate = ''

    # Check parameter misspecification_LOD for subclone under limit of detection
    if 'misspecification_LOD'.lower() in parname_list:
        if 'misspecification_LOD'.lower() in kwargs.keys() and subclone_LOD:
            misspecification_LOD = kwargs['misspecification_LOD'.lower()]
            if type(misspecification_LOD) is not list:
                misspecification_LOD = [misspecification_LOD]
            misspecification_LOD = DPM_check_in_list_par(misspecification_LOD,
                                                         'misspecification_LOD',
                                                         MISSPECIFICATION_LOD,
                                                         print_quotation,
                                                         beginning_char)
            if misspecification_LOD is None:
                sys.exit()
            misspecification_LOD = misspecification_LOD[0]
        else:
            misspecification_LOD = ''
    else:
        misspecification_LOD = ''

    # Check parameter LOD.
    if 'LOD'.lower() in parname_list:
        if 'LOD'.lower() in kwargs.keys():
            LOD = kwargs['LOD'.lower()]
            if type(LOD) is not list:
                LOD = list(LOD)
            LOD = DPM_check_in_list_par(LOD, 'LOD', LOD_LIST, print_quotation, beginning_char)
            if LOD is None:
                sys.exit()
        else:
            LOD = LOD_LIST
    else:
        LOD = ''

    PAR_dict = {'Num_drug': Num_drug, 'dose_method': dose_method, 'dose_interval': dose_interval, 'dose_combination': dose_combination,
                'Simduration': Simduration, 'use_input_dose_combination_only': use_input_dose_combination_only, 'Stepsize': Stepsize,
                'Limit_mortality': Limit_mortality, 'Limit_radiologicdetection': Limit_radiologicdetection,  'run_sim': run_sim,
                'erase_preresult': erase_preresult, 'PAR_criterisa': PAR_criterisa, 'Strategy_name': Strategy_name, 'fullinput': fullinput,
                'par_save_block_size': par_save_block_size, 'save_filename_param': save_filename_param,
                'save_filename_stopt': save_filename_stopt, 'save_filename_pop': save_filename_pop,
                'save_filename_dosage': save_filename_dosage, 'pathsave': pathsave, 'pathload': pathload,
                'save_filename_eachtimepoint': save_filename_eachtimepoint, 'use_parallel': use_parallel, 'filename_csv': filename_csv,
                'misspecification_filename_csv': misspecification_filename_csv, 'Limit_moleculardetection': Limit_moleculardetection,
                'misspecification_sim_only': misspecification_sim_only, 'savename': savename, 'par': par, 'plot': plot, 'par_ind': par_ind,
                'Num_stepdiff': Num_stepdiff, 'filename_pattern': filename_pattern, 'misspecification_ofdecision': misspecification_ofdecision,
                'misspecification_atsim': misspecification_atsim, 'subclone_LOD': subclone_LOD, 'mutation_rate': mutation_rate,
                'misspecification_LOD': misspecification_LOD, 'pathloadmis': pathloadmis, 'LOD': LOD,
                'misspecification_fileload': misspecification_fileload}

    return PAR_dict


# Assign parameter from .csv files.
def DPM_assign_par_csvread(par_csv, Num_drug, PAR_criterisa, Simduration, Limit_mortality, Stepsize):
    PARtotal = []
    for i_par in par_csv:
        if type(i_par) is dict:
            assert Num_drug == 2
            Num_cell_type = 4
            i_par.update({'Num_drug': Num_drug, 'Num_cell_type': Num_cell_type})
            criterisa_i = DPM_check_criterisa_2drug(PAR_criterisa, i_par, Simduration, Limit_mortality)
            par_check_index = DPM_check_par(sum([i_par['Spop'], i_par['R1pop'], i_par['R2pop'], i_par['R12pop']]), Limit_mortality, Simduration,
                                            Stepsize)

            if all(criterisa_i) and not par_check_index:
                PARtotal.append(i_par)
    return PARtotal


# Assign parameters into 2 drug case fomart.
def DPM_assign_par_default_2drug(par, X0total):
    par_out = dict()

    g0_S, ratioR1toX0, ratioR2toX0, Sa_ratio_S_D1tog0_S, Sa_ratio_S_D2toS_D1, Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2, T_StoR1, T_StoR2 = par
    ratioStoX0, ratioR12toX0 = 1-ratioR1toX0-ratioR2toX0, 0

    g0_R1 = g0_R2 = g0_R12 = g0_S
    Sa_ratio_R1_D2toS_D2 = Sa_ratio_R2_D1toS_D1 = 1.0
    Sa_ratio_R12_D1toS_D1, Sa_ratio_R12_D2toS_D2 = Sa_ratio_R1_D1toS_D1, Sa_ratio_R2_D2toS_D2

    T_R1toR12, T_R2toR12 = T_StoR2, T_StoR1
    # X0
    par_out['Spop'], par_out['R1pop'], par_out['R2pop'], par_out['R12pop'] = \
        X0total * np.array([ratioStoX0, ratioR1toX0,  ratioR2toX0, ratioR12toX0]).ravel()
    # g0
    par_out['g0_S'], par_out['g0_R1'], par_out['g0_R2'], par_out['g0_R12'] = np.array([g0_S, g0_R1, g0_R2, g0_R12]).ravel()
    # Sa
    par_out['Sa.S.D1.'] = Sa_ratio_S_D1tog0_S * g0_S
    par_out['Sa.S.D2.'] = Sa_ratio_S_D2toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R1.D1.'] = Sa_ratio_R1_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R1.D2.'] = Sa_ratio_R1_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R2.D1.'] = Sa_ratio_R2_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R2.D2.'] = Sa_ratio_R2_D2toS_D2 * par_out['Sa.S.D2.']
    par_out['Sa.R12.D1.'] = Sa_ratio_R12_D1toS_D1 * par_out['Sa.S.D1.']
    par_out['Sa.R12.D2.'] = Sa_ratio_R12_D2toS_D2 * par_out['Sa.S.D2.']
    # T
    par_out['T.R1..S.'] = T_StoR1
    par_out['T.R2..S.'] = T_StoR2
    par_out['T.R12..R1.'] = T_R1toR12
    par_out['T.R12..R2.'] = T_R2toR12

    T = [i for i in HEADING_2DRUG_INPUT_PARAM_CSV if i.find('T.') != -1]
    for i_T in T:
        if i_T not in par_out.keys():
            par_out[i_T] = 0.0

    return par_out


# Assign strategy results for save into .csv files.
def DPM_assign_saveresult_csv(Strategy, Strategy_text, Simduration, Stepsize, cured_survival, non_report, para_PAR_index, filename_stopt,
                              filename_dosage, filename_pop, filename_eachtimepoint, Limit_mortality):
    if Strategy is not None:
        t, X, d = Strategy
        t_atstepsize = list(np.arange(0, Simduration+Stepsize, Stepsize))

        if filename_dosage:
            # find the index of t
            pos = [int(np.where(t == i_step)[0]) if np.isin(i_step, t, assume_unique=True) and i_step < max(t) else -1
                   for i_step in t_atstepsize[0:-1]]
            dose_Strategy = [tuple(d[:, x].T) if x != -1 else -1 for x in pos]
            # round to 2 digit for the tuple of dose value
            dose_Strategy = [tuple(round(y, 2) for y in x) if type(x) is tuple else -1 for x in dose_Strategy]
            dose_Strategy = para_PAR_index + Strategy_text + ', '.join(str(x).replace(' ', '') for x in dose_Strategy)
        else:
            dose_Strategy = ''

        if filename_pop:
            # if t[-1] not in t_atstepsize:
            #     t_atstepsize.append(t[-1])
            #     t_atstepsize.sort()
            pos = [int(np.where(t == i_step)[0]) if np.isin(i_step, t, assume_unique=True) else -1 for i_step in t_atstepsize[1:]]
            pop_Strategy = [tuple(X[:, x].T) if x != -1 else -1 for x in pos]
            pop_Strategy = [tuple(round(y, 2) for y in x) if type(x) is tuple else -1 for x in pop_Strategy]
            pop_Strategy = para_PAR_index + Strategy_text + ', '.join('(' + ','.join('{:.2e}'.format(y) for y in x) + ')'
                                                                      if type(x) is tuple else str(x) for x in pop_Strategy)
        else:
            pop_Strategy = ''

        if filename_stopt:
            if all(X[:, -1] < 1):
                survival_Strategy = cured_survival
            else:
                survival_Strategy = min(t[-1], Simduration)
                if survival_Strategy == Simduration and sum(X[:, -1]) < Limit_mortality:
                    survival_Strategy += 1
        else:
            survival_Strategy = ''

        if filename_eachtimepoint:
            X = map(tuple, X.T)
            X = tuple(X)
            X = [tuple(round(y, 2) for y in x) if type(x) is tuple else -1 for x in X]
            d = map(tuple, d.T)
            d = tuple(d)
            d = d + (d[-1],)
            a = tuple(zip(d, X))
            eachtimepoint_Strategy = para_PAR_index + Strategy_text + \
                ', '.join(', '.join(str(y).replace(' ', '') for y in x) if type(x) is tuple else str(x) for x in a)
        else:
            eachtimepoint_Strategy = ''
    else:
        survival_Strategy = non_report
        dose_Strategy = para_PAR_index + Strategy_text + non_report
        pop_Strategy = para_PAR_index + Strategy_text + non_report
        eachtimepoint_Strategy = para_PAR_index + Strategy_text + non_report

    return {'survival': survival_Strategy, 'dose': dose_Strategy, 'pop': pop_Strategy, 'eachtimepoint':  eachtimepoint_Strategy}


# Check whether the parameter fullfile the criterisa in 2 drug case.
def DPM_check_criterisa_2drug(PAR_criterisa, PAR, Simduration, Limit_mortality):
    Num_criterisa = 9
    criterisa = [True] * Num_criterisa

    Num_drug, g0_S, g0_R1, g0_R2, g0_R12 = [PAR.get(key) for key in ['Num_drug', 'g0_S', 'g0_R1', 'g0_R2', 'g0_R12']]

    # Criterisa(0): The potency of drugs on S cell follows: drug 1 >= drug 2.
    criterisa_0 = (0 not in PAR_criterisa) or (PAR['Sa.S.D1.'] >= PAR['Sa.S.D2.'])
    criterisa[0] = criterisa_0
    if not criterisa_0:
        return criterisa

    # Criterisa PNAS(1): S, R1, R2, R12 cells >= 0.
    criterisa_1 = (1 not in PAR_criterisa) or all(i >= 0 for i in [PAR['Spop'], PAR['R1pop'], PAR['R2pop'], PAR['R12pop']])
    criterisa[1] = criterisa_1
    if not criterisa_1:
        return criterisa

    # Criterisa PNAS(2): Sa(R1,D1) <= g0_R1 and Sa(R2,D2) <= g0_R2.
    criterisa_2 = (2 not in PAR_criterisa) or (PAR['Sa.R1.D1.'] <= PAR['g0_R1'] and PAR['Sa.R2.D2.'] <= PAR['g0_R2'])
    criterisa[2] = criterisa_2
    if not criterisa_2:
        return criterisa

    # Criterisa PNAS(3): Sa(S,D1) > 0 and Sa(S,D2) > 0.
    criterisa_3 = (3 not in PAR_criterisa) or (PAR['Sa.S.D1.'] > 0 and PAR['Sa.S.D2.'] > 0)
    criterisa[3] = criterisa_3
    if not criterisa_3:
        return criterisa

    # Criterisa Biology Direct(4): Sa(S,D1) > g0_S and Sa(S,D2) > g0_S, both inequalities hold simultaneously.
    criterisa_4 = (4 not in PAR_criterisa) or (PAR['Sa.S.D1.'] > PAR['g0_S'] and PAR['Sa.S.D2.'] > PAR['g0_S'])
    criterisa[4] = criterisa_4
    if not criterisa_4:
        return criterisa

    # Criterisa Biology Direct(5): Sa(R12,D1) < g0_R12 and Sa(R12,D2) < g0_R12, both inequalities hold simultaneously.
    criterisa_5 = (5 not in PAR_criterisa) or (PAR['Sa.R12.D1.'] < PAR['g0_R12'] and PAR['Sa.R12.D2.'] < PAR['g0_R12'])
    criterisa[5] = criterisa_5
    if not criterisa_5:
        return criterisa

    # Criterisa Biology Direct(6): The patient can be cured by simultaneous full dosages of all drugs (an invalid option due to toxicity)
    # but cannot be cured by any valid static treatment.
    if 6 in PAR_criterisa:
        LSsim = True
        X0, g0, Sa, T = DPM_generate_X0_2drug(PAR), DPM_generate_g0_2drug(PAR), DPM_generate_Sa_2drug(PAR), DPM_generate_T_2drug(PAR)
        t = np.arange(0, Simduration+SIMTIMESTEP_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL)
        d_fulldose = np.full([Num_drug, len(t)-1], 1, dtype=float)
        X_fulldose, t_fulldose, _ = DPM_sim_model2012(X0, T, g0, Sa, d_fulldose, t, Limit_mortality, LSsim)
        if not np.all(X_fulldose[:, -1] < 1):
            criterisa[6] = False
            return criterisa
        else:
            d_drug1fulldose = np.tile(np.array([[1, 0]]).T, (1, len(t)-1))
            X_drug1fulldose, t_drug1fulldose, _ = DPM_sim_model2012(X0, T, g0, Sa, d_drug1fulldose, t, Limit_mortality, LSsim)
            if not np.any(X_drug1fulldose[:, -1] > 1):
                criterisa[6] = False
                return criterisa
            else:
                d_drug2fulldose = np.tile(np.array([[0, 1]]).T, (1, len(t)-1))
                X_drug2fulldose, t_drug2fulldose, _ = DPM_sim_model2012(X0, T, g0, Sa, d_drug2fulldose, t, Limit_mortality, LSsim)
                if not np.any(X_drug2fulldose[:, -1] > 1):
                    criterisa[6] = False
                    return criterisa
                else:
                    d_halfdose = np.tile(np.array([[1/2, 1/2]]).T, (1, len(t)-1))
                    X_halfdose, t_halfdose, _ = DPM_sim_model2012(X0, T, g0, Sa, d_halfdose, t, Limit_mortality, LSsim)
                    if not np.any(X_halfdose[:, -1] > 1):
                        criterisa[6] = False
                        return criterisa

    # Criterisa PNAS(7): full dose of single drug1 and drug2 be sufficient to cause a 25% increase in the progression-free survival compared to
    # no treatment.
    if 7 in PAR_criterisa:
        LSsim = True
        X0, g0, Sa, T = DPM_generate_X0_2drug(PAR), DPM_generate_g0_2drug(PAR), DPM_generate_Sa_2drug(PAR), DPM_generate_T_2drug(PAR)
        t = np.arange(0, Simduration+SIMTIMESTEP_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL)
        # No treatment
        d_zerodose = np.zeros([Num_drug, len(t)-1])
        _, t_notreat, _ = DPM_sim_model(X0, T, g0, Sa, d_zerodose, t, Limit_mortality, LSsim)
        if t_notreat[-1] > Simduration/1.25:
            criterisa[7] = False
            return criterisa
        else:
            # Full dose of drug 1
            d_fulldosedrug1 = np.tile(np.array([[1, 0]]).T, (1, t.shape[0] - 1))
            _, t_fulldosedrug1, _ = DPM_sim_model(X0, T, g0, Sa, d_fulldosedrug1, t, Limit_mortality, LSsim)
            if t_fulldosedrug1[-1] < 1.25 * t_notreat[-1]:
                criterisa[7] = False
                return criterisa
            else:
                # Full dose of drug 2
                d_fulldosedrug2 = np.tile(np.array([[0, 1]]).T, (1, t.shape[0] - 1))
                _, t_fulldosedrug2, _ = DPM_sim_model(X0, T, g0, Sa, d_fulldosedrug2, t, Limit_mortality, LSsim)
                if t_fulldosedrug2[-1] < 1.25 * t_notreat[-1]:
                    criterisa[7] = False
                    return criterisa

    # Criterisa PNAS(8): survival time of strategy 0 should be less than 75% of the simulation time to test whether there is
    # 25% relative improvement by other strategies.
    if 8 in PAR_criterisa:
        LSsim = True
        t_Strategy0, X_Strategy0, d_Strategy0 = DPM_strategy_0(PAR, Simduration, STEPSIZE_DEFAULT_VAL, SIMTIMESTEP_DEFAULT_VAL,
                                                               Limit_mortality, LIMIT_RADIOLOGICDETECTION_DEFAULT_VAL, LSsim, False, None)
        if (t_Strategy0[-1] > Simduration/1.25) or (not all(X_Strategy0[:, -1] < 1)):
            criterisa[8] = False

    return criterisa
