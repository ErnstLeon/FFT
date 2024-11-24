#include <array>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <FFT.hpp>

int main(int argc, char ** argv){

    constexpr double PI = 3.14159265358979323846264338327950288;

    std::ofstream file;
    file.open("fft.dat");

    file << "t " << "signal " << "freq " << "component" << std::endl;

    constexpr int N = 512;
    constexpr double delta_t = 0.025;
    constexpr double t_0 = 0;

    std::array<double, N> time = {0.};
    std::array<double, N> freq = {0.};
    std::array<double, N * 2> signal = {0.};
    std::array<double, N * 2> component = {0.};

    for(int i = 0; i < N; i++)
    {
        time.at(i) = t_0 + delta_t * i;
        signal.at(i * 2) = sin(t_0 + delta_t * i * (2 * PI) * 2);
    }

    for(int i = 0; i < N / 2; i++)
    {
        freq.at(i) = 1 / (N * delta_t) * i;
        freq.at(N / 2 + i) = - 1 / (2 * delta_t) + 1 / (N * delta_t) * i;
    }
    
    component = signal;

    fft(component);

    for(int i = 0; i < N; i++)
    {
        file << time.at(i) << " " << signal.at(i * 2) << " " << freq.at(i) << " " << std::abs(component.at(i * 2))<< std::endl;
    }
    
    file.close();
}