#
# Siril Tests
#

# Criterion based unit-tests
if criterion_dep.found()
     sorting_exec = executable('sorting',
                               'sorting.c',
                               dependencies : [siril_dep, criterion_dep],
                               link_args : siril_link_arg,
                               c_args : siril_c_flag,
                               cpp_args : siril_cpp_flag)

     test('sorting', sorting_exec)
endif


# Performance testing
sorting_perf_exec = executable('sorting_perf',
                               'measure_sorting.c',
                               dependencies : siril_dep,
                               link_args : siril_link_arg,
                               c_args : siril_c_flag,
                               cpp_args : siril_cpp_flag)

test('sorting_perf', sorting_perf_exec, suite: 'perfs', protocol: 'exitcode', is_parallel : false)
