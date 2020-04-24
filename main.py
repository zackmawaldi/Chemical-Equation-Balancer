try:
    import math
except ImportError:
    print("This script requires 'math' package, and it wasn't found")
    quit()

try:
    import re
except ImportError:
    print("This script requires 're' package, and it wasn't found")
    quit()

try:
    import sympy
except ImportError:
    print("This script requires 'sympy' package, and it wasn't found")
    quit()

try:
    from fractions import Fraction
except ImportError:
    print("This script requires 'Fraction' from fractions package, and it wasn't found")
    quit()


def lcm(lst):  # finds least common multiple in lst
    in_lcm = lst[0]
    for i in lst[1:]:
        in_lcm = abs(in_lcm * i) // math.gcd(in_lcm, i)
    return in_lcm


def full_eq_split(full_eq):  # splits inputted string eq to a list of reactants, '-->', and products.
    global r_compounds, p_compounds, r_compounds_wpai, p_compounds_wpai
    initial_list = full_eq.split()
    initial_list = list(filter(lambda a: a != '+', initial_list))
    n = initial_list.index('-->')
    r_compounds_wpai = initial_list[:n]  # wpai = with poly atomic ion. For the final printed string.
    p_compounds_wpai = initial_list[n + 1:]
    for w, item in enumerate(initial_list):  # this loop multiply through PAI Ca3(PO4)2 --> Ca3P2O8
        if '(' in item:  # finds PAIs
            split_pai = re.split('[(-)]', item)
            first_element = split_pai.pop(0)
            if split_pai[-1] == '':  # split_pai[-1] is cof outside ()
                split_pai[-1] = '1'
            mod = [s for s in re.split("([A-Z][^A-Z]*)", split_pai[0]) if s]  # As1O2Na3 --> [As1, O2, Na3]
            mod_elem = [re.sub('[1-9]', '', x) for x in mod]
            mod_elem_num = [re.sub('[a-zA-Z]', '', x) for x in mod]
            for n, element in enumerate(mod_elem_num):
                try:
                    mod_elem_num[n] = str(int(element) * int(split_pai[-1]))
                except ValueError:  # arises when element == ''
                    mod_elem_num[n] = split_pai[-1]
            final_lst = [mod_elem[x] + mod_elem_num[x] for x in range(len(mod_elem))]
            final_str = ''.join(final_lst)
            initial_list[w] = first_element + final_str
    n = initial_list.index('-->')
    r_compounds = initial_list[:n]
    p_compounds = initial_list[n + 1:]


global r_compounds, p_compounds, r_compounds_wpai, p_compounds_wpai  # this line made PyCharm happy


def dict_convert(lst_compounds):  # converts compounds into dict. Ex. ['NO2'] --> {N: 1, O: 2}
    joined_str = ''.join(lst_compounds)
    mod = [s for s in re.split("([A-Z][^A-Z]*)", joined_str) if s]
    comp_num_lst = []
    comp_num_lst.extend(mod)
    for item in comp_num_lst:
        t = comp_num_lst.index(item)
        check = any(c.isdigit() for c in item)
        if check is False:
            comp_num_lst[t] = '{}1'.format(comp_num_lst[t])
    length = len(comp_num_lst)
    for item in comp_num_lst:
        match = re.match(r"([a-z]+)([0-9]+)", item, re.I)
        if match:
            items = list(match.groups())
            for element in items:
                comp_num_lst.append(element)
    del comp_num_lst[0:length]
    res_dct = {}
    for i in comp_num_lst:
        x = comp_num_lst.index(i)
        entry = comp_num_lst.pop(x)
        if entry in res_dct:
            res_dct[entry] += int(comp_num_lst[x])
        else:
            res_dct[entry] = int(comp_num_lst[x])
    return res_dct


def get_element_list(r_or_p):  # converts list of compounds to list of the elements in the reaction, w/o cof.
    joined_str = ''.join(r_or_p)
    mod = [s for s in re.split("([A-Z][^A-Z]*)", joined_str) if s]
    comp_num_lst = []
    comp_num_lst.extend(mod)
    for item in comp_num_lst:
        t = comp_num_lst.index(item)
        check = any(c.isdigit() for c in item)
        if check is False:
            comp_num_lst[t] = '{}1'.format(comp_num_lst[t])
    for element in comp_num_lst:
        n = comp_num_lst.index(element)
        mod = ''.join([i for i in element if not i.isdigit()])
        comp_num_lst.remove(element)
        comp_num_lst.insert(n, mod)
    comp_num_lst.sort()
    for element in comp_num_lst:  # removes duplicates
        n = comp_num_lst.index(element)
        try:
            if comp_num_lst[n] == comp_num_lst[n+1]:
                comp_num_lst.remove(element)
        except IndexError:
            pass
    return comp_num_lst


def welcome():
    while True:
        eq_input = input('Your equation: \n')
        if '-->' and '+' in eq_input:
            break
        print("Invalid equation; follow example.")
    return eq_input


def main_script():
    full_equation = welcome()

    full_eq_split(full_equation)

    reaction_elements = get_element_list(r_compounds)

    n = len(r_compounds) + len(p_compounds)
    m = len(get_element_list(r_compounds))

    initial_array = [[0] * n for i in range(m)]  # this is the array where RREF will be preformed before values appended
    final_array = [[0] * n for i in range(m)]  # saved for the end

    for item in r_compounds:  # converts compounds in list to dict, then appends cof into array to prepare for RREF
        n = r_compounds.index(item)
        compound_dict = dict_convert(item)
        for key in compound_dict:
            if key in reaction_elements:
                m = reaction_elements.index(key)
                initial_array[m][n] = compound_dict['{}'.format(key)]
            else:
                m = reaction_elements.index(key)
                initial_array[m][n] = 0

    for item in p_compounds:  # same as above but for products' list
        n = p_compounds.index(item) + len(r_compounds)
        compound_dict = dict_convert(item)
        for key in compound_dict:
            if key in reaction_elements:
                m = reaction_elements.index(key)
                initial_array[m][n] = -1 * compound_dict['{}'.format(key)]  # for RREF to work, products have to be * -1
            else:
                m = reaction_elements.index(key)
                initial_array[m][n] = 0

    initial_matrix = sympy.Matrix(initial_array)  # convert array to matrix then RREF
    reduced_matrix = initial_matrix.rref()
    reduced_matrix = list(reduced_matrix[0])

    y = 0  # converts from reduced matrix to a final array
    while True:
        if y == len(reduced_matrix):
            break
        c = y // (n + 1)
        r = y - (c * (n + 1))
        final_array[c][r] = float(reduced_matrix[y])
        y += 1

    n = len(r_compounds) + len(p_compounds)
    m = len(get_element_list(r_compounds))

# At this point, our matrix's results are fractions. The answer works, but fractions aren't expected.
# So, we get lcm of the denominators of those fractions and multiply through by result.
# This gives lowest possible whole number answers.

    denominators = []
    fraction_values = []
    x = 0
    while True:  # extracts denom and fraction values of results into their separate lists
        if x > m - 1:
            break
        insert_den = Fraction(final_array[x][n - 1]).limit_denominator().denominator
        if insert_den != '0':
            denominators.append(insert_den)
        insert_float = Fraction(final_array[x][n - 1]).limit_denominator()
        fraction_values.append(insert_float)
        x += 1

    den_lcm = lcm(denominators)

    coefficient_values = [int(abs(x * den_lcm)) for x in fraction_values]
    for n, item in enumerate(coefficient_values):
        if item == 0:  # if it's 0, means it was a free variable in the RREF matrix
            coefficient_values[n] = den_lcm
    mod_eq = r_compounds_wpai + p_compounds_wpai
    for n, item in enumerate(mod_eq):  # adds cof behind compound to get it ready for ''.join
        try:
            mod_eq[n] = '{}'.format(coefficient_values[n]) + ' ' + item
        except IndexError:
            mod_eq[n] = '{}'.format(den_lcm) + ' ' + item

    for i in range((len(r_compounds_wpai) + len(p_compounds_wpai)) * 2):  # adds the '+' back
        if i % 2 != 0:
            mod_eq.insert(i, ' + ')
    if mod_eq[-1] == ' + ':
        del mod_eq[-1]

    mod_eq.insert(len(r_compounds_wpai) * 2, ' --> ')  # adds the arrow
    del mod_eq[mod_eq.index(' --> ') - 1]

    final_eq = ''.join(mod_eq)
    print('\n' + final_eq)


print("Welcome to Zack's Chemical Equation Balancer!\n"
      "Ex. C7H5(NO2)3 --> N2 + O2 + H2 + CO2\n"
      "*The script is case and capitalization sensitive!\n"
      "*Double salts, and complexes have to be manually multiplied through \n"
      "Ex. K[PtCl3(C2H4)] --> KPtCl3C2H4 and (NH4)3(PO4) --> N3H12PO4\n"
      "Good rule of thumb: multiply polyatomic cat/ion through whenever you're getting funky answers\n")

while True:
    main_script()
    ask_try_again = input('\nWant to try another equation? Type 1 or 2.\n'
                          '1) Yes\n'
                          '2) No\n')
    if ask_try_again == '2':
        break
