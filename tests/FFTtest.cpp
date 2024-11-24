#include <array>
#include <cstdlib>
#include <FFT.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(FFTtest, FFTandInverseFFT)
{
    std::array<double, 16> data_orig;

    std::srand(10);

    for(auto &i : data_orig) i = static_cast<double>(std::rand() % 10);

    auto data_copy = data_orig;

    fft(data_copy);
    inverse_fft(data_copy);

    for(int i = 0; i < 16; ++i)
    {
        EXPECT_NEAR(data_orig.at(i), data_copy.at(i), 1e-8)
            << i
            << "th value in the array does not equal to the original data after FFT and inverse FFT "
            << data_orig.at(i) << " " << data_copy.at(i);
    }

}

TEST(FFTtest, FFTzeroPadding)
{
    std::array<double, 16> data_orig;

    std::srand(10);
    for(auto &i : data_orig) i = static_cast<double>(std::rand() % 10);
    data_orig.at(12) = static_cast<double>(0);
    data_orig.at(13) = static_cast<double>(0);
    data_orig.at(14) = static_cast<double>(0);
    data_orig.at(15) = static_cast<double>(0);

    fft(data_orig);

    std::vector<double> data_orig_v;

    std::srand(10);
    for(int i = 0; i < 12; ++i) data_orig_v.push_back(static_cast<double>(std::rand() % 10));
    fft(data_orig_v);

    for(int i = 0; i < 16; ++i)
    {
        EXPECT_NEAR(data_orig.at(i), data_orig_v.at(i), 1e-8)
            << i
            << "th value of the FFT with automiatc zero padding differs from an FFT with manual zero padding"
            << data_orig.at(i) << " " << data_orig.at(i);
    }

}

TEST(FFTtest, InverseFFTzeroPadding)
{
    std::array<double, 16> data_orig;

    std::srand(10);
    for(auto &i : data_orig) i = static_cast<double>(std::rand() % 10);
    data_orig.at(12) = static_cast<double>(0);
    data_orig.at(13) = static_cast<double>(0);
    data_orig.at(14) = static_cast<double>(0);
    data_orig.at(15) = static_cast<double>(0);

    inverse_fft(data_orig);

    std::vector<double> data_orig_v;

    std::srand(10);
    for(int i = 0; i < 12; ++i) data_orig_v.push_back(static_cast<double>(std::rand() % 10));
    inverse_fft(data_orig_v);

    for(int i = 0; i < 16; ++i)
    {
        EXPECT_NEAR(data_orig.at(i), data_orig_v.at(i), 1e-8)
            << i
            << "th value of the FFT with automiatc zero padding differs from an FFT with manual zero padding"
            << data_orig.at(i) << " " << data_orig.at(i);
    }

}
