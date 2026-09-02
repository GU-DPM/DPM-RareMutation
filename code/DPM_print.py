import sys
from DPM_constant import Colored
""" This script defines the print functions used in DPM (Dynamic Precision Medicine) model."""


# Print message when error happens and will exit the program.
def DPM_print_errorandclose(beginchar=''):
    print(beginchar + Colored.BOLD + Colored.RED + 'ERROR! ' + Colored.END + 'Program terminated.')
    return


# Print error when input positional arguments but only keyword arguments are allowed in the function.
def DPM_print_keywordargument_only(argv, keys, print_quotation):
    if len(argv) > 2:
        print('You entered non-keyword arguments:')
    else:
        print('You entered non-keyword argument:')
    for arg_i in argv:
        if arg_i is not argv[-1]:
            DPM_print_stringandnum(arg_i, ', ', print_quotation)
        else:
            DPM_print_stringandnum(arg_i, '. ', print_quotation)
    print('The function only allows case insensitive keyword arguments or without argument.')
    print('The allowed case insensitive keywords are: ')
    DPM_print_list(keys, False)
    DPM_print_errorandclose()
    return


# Print error when the keyword arguments are not existed in the function.
def DPM_print_keyworderror(kwargs, parname_list, print_quotation, beginning_char):
    input_parname_lowercase = [name_i.lower() for name_i in parname_list]
    falseinput = [key for key, value in kwargs.items() if key.lower() not in input_parname_lowercase]
    if falseinput:
        if len(falseinput) > 1:
            print('The entered keyword arguments', end=' ')
        else:
            print('The entered keyword argument', end=' ')
        for key_i in falseinput:
            if key_i is not falseinput[-1]:
                print(Colored.BOLD + Colored.PURPLE + key_i + Colored.END, end=', ')
            else:
                print(Colored.BOLD + Colored.PURPLE + key_i + Colored.END, end=' ')
        print('are not in the list.\nThe allowed case insensitive keywords are:')
        DPM_print_list(parname_list, print_quotation)
        DPM_print_errorandclose(beginning_char)
        sys.exit()
    else:
        return


# Print a str or number.
def DPM_print_stringandnum(val, end, print_quotation):
    if print_quotation:
        if type(val) is str:
            print(Colored.BOLD + Colored.PURPLE + '\'' + val + '\'' + Colored.END, end=end)
        else:
            print(Colored.BOLD + Colored.PURPLE + str(val) + Colored.END, end=end)
    else:
        print(Colored.BOLD + Colored.PURPLE + str(val) + Colored.END, end=end)
    return


# Print a tuple.
def DPM_print_tuple(tupl, end, print_quotation):
    print(Colored.BOLD + Colored.PURPLE + '(' + Colored.END, end='')
    count = 0
    for i in tupl:
        if count != len(tupl)-1:
            DPM_print_stringandnum(i, '', print_quotation)
            print(Colored.BOLD + Colored.PURPLE + ', ' + Colored.END, end='')
        else:
            DPM_print_stringandnum(i, '', print_quotation)
        count += 1
    print(Colored.BOLD + Colored.PURPLE + ')' + Colored.END, end=end)
    return


# Print a list.
def DPM_print_list(input_list, print_quotation):
    for i in input_list:
        if i is not input_list[-1]:
            if type(i) is not tuple:
                DPM_print_stringandnum(i, ', ', print_quotation)
            else:
                DPM_print_tuple(i, ', ', print_quotation)
        else:
            if type(i) is not tuple:
                DPM_print_stringandnum(i, '. ', print_quotation)
            else:
                DPM_print_tuple(i, '. ', print_quotation)
    return


# Print the message not a number.
def DPM_print_notnumerical(name):
    print('The inputted ' + Colored.BOLD + Colored.PURPLE + name + Colored.END + ' should be a single numerical value.')
    return


# Print the message not a integer.
def DPM_print_notinteger(name):
    print('The inputted ' + Colored.BOLD + Colored.PURPLE + name + Colored.END + ' should be an integer.')
    return


# Print the message not a string.
def DPM_print_notstring(name):
    print('The inputted ' + Colored.BOLD + Colored.PURPLE + name + Colored.END + ' should be a string.')
    return


# Print the message not a boolean.
def DPM_print_notboolean(name):
    print('The inputted ' + Colored.BOLD + Colored.PURPLE + name + Colored.END + ' should be a "True" or "False" boolean value.')
    return


# Print the type of the input.
def DPM_print_val_type(var, name):
    print('The type of ' + Colored.BOLD + Colored.PURPLE + name + Colored.END + ' is ' +
          Colored.BOLD + Colored.PURPLE + type(var).__name__ + Colored.END + '.')
    return
