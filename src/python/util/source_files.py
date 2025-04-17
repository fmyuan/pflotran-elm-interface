def get_source_file_roots():
    source_file_roots = []
    # Obtain list of source files
    for line in open('pflotran_object_files.txt','r'):
        # find .o file
        # could use re.split() here, but too complicated.
        w = line.split('}')
        if len(w) == 2:
            w2 = w[1].split('.o')
            source_file_roots.append(w2[0])
    source_file_roots.append('pflotran')
    source_file_roots.append('pflotran_rxn')
    source_file_roots.append('pflotran_derivative')

    # Alphabetize
    source_file_roots.sort()
    #print(source_file_roots)
    return source_file_roots

def get_clm_pflotran_source_file_roots():
    source_file_roots = []
    # clm-pflotran files
    source_file_roots.append('../clm-pflotran/pflotran_model')
    source_file_roots.append('../clm-pflotran/pflotran_interface_main')
    source_file_roots.append('../clm-pflotran/mapping')
    source_file_roots.append('../clm-pflotran/clm_pflotran_interface_data')

    # Alphabetize
    source_file_roots.sort()
    #print(source_file_roots)
    return source_file_roots

def get_unit_test_files():
    source_file_roots = []
    # unit test files
    source_file_roots.append('unittests/test_characteristic_curves.pf')
    source_file_roots.append('unittests/test_characteristic_curves_thermal.pf')
    source_file_roots.append('unittests/test_eos_gas.pf')
    source_file_roots.append('unittests/test_eos_water.pf')
    source_file_roots.append('unittests/test_geometry.pf')
    source_file_roots.append('unittests/test_material.pf')
    source_file_roots.append('unittests/test_saturation_function.pf')
    source_file_roots.append('unittests/test_string.pf')
    source_file_roots.append('unittests/test_utility.pf')

    # Alphabetize
    source_file_roots.sort()
    #print(source_file_roots)
    return source_file_roots

