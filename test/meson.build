project('FFT', 'cpp', default_options : ['cpp_std=c++20'])

gtest_subproj = subproject('gtest')
gtest_dep = gtest_subproj.get_variable('gtest_main_dep')

FFT_dir = include_directories('../include')

e = executable(
    'FFTtest',
    'FFTtest.cpp',
    include_directories : FFT_dir,
    dependencies : gtest_dep
)
test('FFTtest', e)