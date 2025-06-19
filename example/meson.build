project('FFT', 'cpp', default_options : ['cpp_std=c++20'])

FFT_dir = include_directories('../include')

cpp_flags = ['-std=c++20', '-O3', '-march=native']

executable('main', 'main.cpp', 
    cpp_args: cpp_flags,
    include_directories : FFT_dir
    )