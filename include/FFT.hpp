#ifndef FFT_H
#define FFT_H

#include <array>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

/*
  fft using Danielson-Lanczos Lemma

  data contains the real and imaginary part of (size / 2) = N discrete
  data points in the time domain f_n

  data[n] = Re(f_n)
  data[n + 1] = Im(f_n)

  data is overwritten with the result, which are the real and imaginary part of N discrete
  data points in the frequency domain F_n

  when delta_t is the time difference of two points in the time domain (i.e. sampling frequency),
  F_0 corresponds to v = 0, F_1 to v = 1/(delta_t * N), ..., F_N to v = (N - 1)/(delta_t * N)
*/
template<typename T, typename S, S size>
  requires (size > 2 && !(size&(size - 1))) //size must be a power of 2 and larger than 2
static inline int fft(std::array<T, size> &data)
{
  constexpr double PI = 3.14159265358979323846264338327950288;
  
  S j = 1;  
  for(S i = 1; i < size; i += 2)
  {
    if ( j > i )
    {
      std::swap(data[j], data[i]);
      std::swap(data[j - 1], data[i - 1]);
    }

    S mid = size / 2;        

    while( mid >= 2 && j > mid)
    { 
      j -= mid;
      mid /= 2; 
    }   
    j += mid;
  }

  //recombine according to Danielson-Lanczos lemma
  S l = 1;
  
  while(l < size / 2)
  { 
    // go through data with blocks of size (2 * l) * 2 (for real and im part)        
  	S block_start = 1;

  	while(block_start < size)
    {                     
      // merge F_k(...e) and F_k(...o) for k = 0 ... (l - 1) 
      // to F_k(...) and F_(k + l)(...)
      // F_k(...e) are replaced by F_k(...) places
      // and F_k(...o) are replaced by F_(k + l)(...)
      // this step works on (2 * l)* 2 (for real and im part) consecutive elements data[]
      S pos = block_start;
      for(S k = 0; k < l; ++k)
      { 
        // data[pos + (2 * l - 1)] contains the real part of F_k(...e)
        // data[pos + (2 * l)] contains the imaginary part of F_k(...e)

        // data[pos] contains the real part of F_k(...o)
        // data[pos -1] contains the imaginary part of F_k(...o)

        // WF_re contains the real part of (W_2l)^k * F_k(...e)
        // WF_im contains the imaginary part of (W_2l)^k * F_k(...e)
        const auto WF_re = data[pos + (2 * l - 1)] * cos((PI * k) / l) - data[pos + (2 * l)] * sin((PI * k) / l);
        const auto WF_im = data[pos + (2 * l - 1)] * sin((PI * k) / l) + data[pos + (2 * l)] * cos((PI * k )/ l);

        // Im(F_k(...o)) - Im((W_2l)^k * F_k(...e)) = Im(F_(k + l)(...))
        // Re(F_k(...o)) - Re((W_2l)^k * F_k(...e)) = Re(F_(k + l)(...))
        data[pos + (2 * l)] = data[pos] - WF_im;   
        data[pos + (2 * l - 1)] = data[pos - 1] - WF_re;  

        // Im(F_k(...o)) + Im((W_2l)^k * F_k(...e)) = Im(F_k(...))
        // Re(F_k(...o)) + Re((W_2l)^k * F_k(...e)) = Re(F_k(...))
        data[pos] += WF_im;  
        data[pos - 1] += WF_re;   

        pos += 2;     
        
      }   
      block_start += l * 4;
    }
    l *= 2;
  }
  return 0;
}

/*
  fft using Danielson-Lanczos Lemma

  data contains the real and imaginary part of (size / 2) = N discrete
  data points in the time domain f_n

  data[n] = Re(f_n)
  data[n + 1] = Im(f_n)

  data is overwritten with the result, which are the real and imaginary part of N discrete
  data points in the frequency domain F_n

  when delta_t is the time difference of two points in the time domain (i.e. sampling frequency),
  F_0 corresponds to v = 0, F_1 to v = 1/(delta_t * N), ..., F_N to v = (N - 1)/(delta_t * N)

  N must be a power of 2, so if this is not the case, 0s are added to the time signal upto the next power in 2
*/
template<typename T>
static inline int fft(std::vector<T> &data)
{
  if(data.size() % 2 != 0) 
  {
    std::cout << "the data vector is not a multiple of 2, "
      "which it must be to contain real and im part" << std::endl;
    return 1;
  } 

  //if the size of data is not a power of 2 or not larger than 2 (N = 1) add zeros
  while((data.size()&(data.size() - 1)) || data.size() <= 2)
  {
    data.push_back(static_cast<T>(0));
  }

  const auto size = data.size();

  constexpr double PI = 3.14159265358979323846264338327950288;
  
  long unsigned int j = 1;  
  for(long unsigned int i = 1; i < size; i += 2)
  {
    if ( j > i )
    {
      std::swap(data[j], data[i]);
      std::swap(data[j - 1], data[i - 1]);
    }

    long unsigned int mid = size / 2;        

    while( mid >= 2 && j > mid)
    { 
      j -= mid;
      mid /= 2; 
    }   
    j += mid;
  }

  //recombine according to Danielson-Lanczos lemma
  long unsigned int l = 1;
  
  while(l < size / 2)
  { 
    // go through data with blocks of size (2 * l) * 2 (for real and im part)        
  	long unsigned int block_start = 1;

  	while(block_start < size)
    {                     
      // merge F_k(...e) and F_k(...o) for k = 0 ... (l - 1) 
      // to F_k(...) and F_(k + l)(...)
      // F_k(...e) are replaced by F_k(...) places
      // and F_k(...o) are replaced by F_(k + l)(...)
      // this step works on (2 * l)* 2 (for real and im part) consecutive elements data[]
      long unsigned int pos = block_start;
      for(long unsigned int k = 0; k < l; ++k)
      { 
        // data[pos + (2 * l - 1)] contains the real part of F_k(...e)
        // data[pos + (2 * l)] contains the imaginary part of F_k(...e)

        // data[pos] contains the real part of F_k(...o)
        // data[pos -1] contains the imaginary part of F_k(...o)

        // WF_re contains the real part of (W_2l)^k * F_k(...e)
        // WF_im contains the imaginary part of (W_2l)^k * F_k(...e)
        const auto WF_re = data[pos + (2 * l - 1)] * cos((PI * k) / l) - data[pos + (2 * l)] * sin((PI * k) / l);
        const auto WF_im = data[pos + (2 * l - 1)] * sin((PI * k) / l) + data[pos + (2 * l)] * cos((PI * k )/ l);

        // Im(F_k(...o)) - Im((W_2l)^k * F_k(...e)) = Im(F_(k + l)(...))
        // Re(F_k(...o)) - Re((W_2l)^k * F_k(...e)) = Re(F_(k + l)(...))
        data[pos + (2 * l)] = data[pos] - WF_im;   
        data[pos + (2 * l - 1)] = data[pos - 1] - WF_re;  

        // Im(F_k(...o)) + Im((W_2l)^k * F_k(...e)) = Im(F_k(...))
        // Re(F_k(...o)) + Re((W_2l)^k * F_k(...e)) = Re(F_k(...))
        data[pos] += WF_im;  
        data[pos - 1] += WF_re;   

        pos += 2;     
        
      }   
      block_start += l * 4;
    }
    l *= 2;
  }
  return 0;
}

/*
  inverse fft using Danielson-Lanczos Lemma

  data contains the real and imaginary part of (size / 2) = N discrete
  data points in the frequency domain F_n

  data[n] = Re(F_n)
  data[n + 1] = Im(F_n)

  data is overwritten with the result, which are the real and imaginary part of N discrete
  data points in the time domain f_n

  when delta_v is the frequency difference of two points in the frequency domain,
  f_0 corresponds to t = 0, f_1 to t = delta_t, ..., f_N to f = (N - 1) * delta_t
*/
template<typename T, typename S, S size>
  requires (size > 2 && !(size&(size - 1))) //size must be a power of 2 and larger than 2
static inline int inverse_fft(std::array<T, size> &data)
{
  constexpr double PI = 3.14159265358979323846264338327950288;
  
  S j = 1;  
  for(S i = 1; i < size; i += 2)
  {
    if ( j > i )
    {
      std::swap(data[j], data[i]);
      std::swap(data[j - 1], data[i - 1]);
    }

    S mid = size / 2;        

    while( mid >= 2 && j > mid)
    { 
      j -= mid;
      mid /= 2; 
    }   
    j += mid;
  }

  //recombine according to Danielson-Lanczos lemma
  S l = 1;
  
  while(l < size / 2)
  { 
    // go through data with blocks of size (2 * l) * 2 (for real and im part)        
  	S block_start = 1;

  	while(block_start < size)
    {                     
      // merge F_k(...e) and F_k(...o) for k = 0 ... (l - 1) 
      // to F_k(...) and F_(k + l)(...)
      // F_k(...e) are replaced by F_k(...) places
      // and F_k(...o) are replaced by F_(k + l)(...)
      // this step works on (2 * l)* 2 (for real and im part) consecutive elements data[]
      S pos = block_start;
      for(S k = 0; k < l; ++k)
      { 
        // data[pos + (2 * l - 1)] contains the real part of F_k(...e)
        // data[pos + (2 * l)] contains the imaginary part of F_k(...e)

        // data[pos] contains the real part of F_k(...o)
        // data[pos -1] contains the imaginary part of F_k(...o)

        // WF_re contains the real part of (W_2l)^k * F_k(...e)
        // WF_im contains the imaginary part of (W_2l)^k * F_k(...e)
        const auto WF_re = data[pos + (2 * l - 1)] * cos((PI * k) / l) + data[pos + (2 * l)] * sin((PI * k) / l);
        const auto WF_im = - data[pos + (2 * l - 1)] * sin((PI * k) / l) + data[pos + (2 * l)] * cos((PI * k )/ l);

        // Im(F_k(...o)) - Im((W_2l)^k * F_k(...e)) = Im(F_(k + l)(...))
        // Re(F_k(...o)) - Re((W_2l)^k * F_k(...e)) = Re(F_(k + l)(...))
        data[pos + (2 * l)] = data[pos] - WF_im;   
        data[pos + (2 * l - 1)] = data[pos - 1] - WF_re;   

        // Im(F_k(...o)) + Im((W_2l)^k * F_k(...e)) = Im(F_k(...))
        // Re(F_k(...o)) + Re((W_2l)^k * F_k(...e)) = Re(F_k(...))
        data[pos] += WF_im;  
        data[pos - 1] += WF_re;  

        pos += 2;        
        
      }   
      block_start += l * 4;
    }
    l *= 2;
  }

  for(auto &i : data) i /= (size / 2);

  return 0;
}

/*
  inverse fft using Danielson-Lanczos Lemma

  data contains the real and imaginary part of (size / 2) = N discrete
  data points in the frequency domain F_n

  data[n] = Re(F_n)
  data[n + 1] = Im(F_n)

  data is overwritten with the result, which are the real and imaginary part of N discrete
  data points in the time domain f_n

  when delta_v is the frequency difference of two points in the frequency domain,
  f_0 corresponds to t = 0, f_1 to t = delta_t, ..., f_N to f = (N - 1) * delta_t

  N must be a power of 2, so if this is not the case, 0s are added to the time signal upto the next power in 2
*/
template<typename T>
static inline int inverse_fft(std::vector<T> &data)
{
  if(data.size() % 2 != 0) 
  {
    std::cout << "the data vector is not a multiple of 2, "
      "which it must be to contain real and im part" << std::endl;
    return 1;
  } 

  //if the size of data is not a power of 2 or not larger than 2 (N = 1) add zeros
  while((data.size()&(data.size() - 1)) || data.size() <= 2)
  {
    data.push_back(static_cast<T>(0));
  }

  const auto size = data.size();

  constexpr double PI = 3.14159265358979323846264338327950288;
  
  long unsigned int j = 1;  
  for(long unsigned int i = 1; i < size; i += 2)
  {
    if ( j > i )
    {
      std::swap(data[j], data[i]);
      std::swap(data[j - 1], data[i - 1]);
    }

    long unsigned int mid = size / 2;        

    while( mid >= 2 && j > mid)
    { 
      j -= mid;
      mid /= 2; 
    }   
    j += mid;
  }

  //recombine according to Danielson-Lanczos lemma
  long unsigned int l = 1;
  
  while(l < size / 2)
  { 
    // go through data with blocks of size (2 * l) * 2 (for real and im part)        
  	long unsigned int block_start = 1;

  	while(block_start < size)
    {                     
      // merge F_k(...e) and F_k(...o) for k = 0 ... (l - 1) 
      // to F_k(...) and F_(k + l)(...)
      // F_k(...e) are replaced by F_k(...) places
      // and F_k(...o) are replaced by F_(k + l)(...)
      // this step works on (2 * l)* 2 (for real and im part) consecutive elements data[]
      long unsigned int pos = block_start;
      for(long unsigned int k = 0; k < l; ++k)
      { 
        // data[pos + (2 * l - 1)] contains the real part of F_k(...e)
        // data[pos + (2 * l)] contains the imaginary part of F_k(...e)

        // data[pos] contains the real part of F_k(...o)
        // data[pos -1] contains the imaginary part of F_k(...o)

        // WF_re contains the real part of (W_2l)^k * F_k(...e)
        // WF_im contains the imaginary part of (W_2l)^k * F_k(...e)
        const auto WF_re = data[pos + (2 * l - 1)] * cos((PI * k) / l) + data[pos + (2 * l)] * sin((PI * k) / l);
        const auto WF_im = - data[pos + (2 * l - 1)] * sin((PI * k) / l) + data[pos + (2 * l)] * cos((PI * k )/ l);

        // Im(F_k(...o)) - Im((W_2l)^k * F_k(...e)) = Im(F_(k + l)(...))
        // Re(F_k(...o)) - Re((W_2l)^k * F_k(...e)) = Re(F_(k + l)(...))
        data[pos + (2 * l)] = data[pos] - WF_im;   
        data[pos + (2 * l - 1)] = data[pos - 1] - WF_re;  

        // Im(F_k(...o)) + Im((W_2l)^k * F_k(...e)) = Im(F_k(...))
        // Re(F_k(...o)) + Re((W_2l)^k * F_k(...e)) = Re(F_k(...))
        data[pos] += WF_im;  
        data[pos - 1] += WF_re;   

        pos += 2;     
        
      }   
      block_start += l * 4;
    }
    l *= 2;
  }

  for(auto &i : data) i /= (size / 2);

  return 0;
}

/*
  fft using Danielson-Lanczos Lemma

  data contains only N real, discrete
  data points in the time domain f_n

  data[n] = f_n = Re(f_n)

  data is overwritten with the result, which are the real and imaginary part of N/2 discrete
  data points in the frequency domain F_n

  when delta_t is the time difference of two points in the time domain (i.e. sampling frequency),
  F_0 corresponds to v = 0, F_1 to v = 1/(delta_t * N), ..., F_N/2 to v = (N/2 - 1)/(delta_t * N)
*/
template<typename T, typename S, S size>
  requires (size > 2 && !(size&(size - 1))) //size must be a power of 2 and larger than 2
static inline int real_fft(std::array<T, size> &data)
{
  constexpr double PI = 3.14159265358979323846264338327950288;

  // copy the array and fourier transform it using a FFT that expects
  // complex data points, i.e. data[n] = Re(f_n), data[n+1] = Im(f_n). 
  // This is equivalent of computing the FFT of (f_(2n) + i f_(2n+1)). 
  // The resulting FFT is then the combination of a FFT over the even data points
  // and a FFT over the odd data points, F_n(even) + i F_n(odd). This can be combined
  // to F_n using the danielson lanczos lemma F_n = F_n(even) + F_n(odd) * exp(i 2Pi n/N)
  std::array<T, size> tmp_data (std::move(data));
  fft(tmp_data);
  
  constexpr S half_size = size >> 1;

  for(S i = 1; i < half_size; ++i)
  {
    T re_H_n = tmp_data[2 * i];
    T im_H_n = tmp_data[(2 * i) + 1];

    T re_H_N = tmp_data[2 * (half_size - i)];
    T im_H_N = tmp_data[(2 * (half_size - i)) + 1];

    T cos_val = cos((2 * PI * i) / size);
    T sin_val = sin((2 * PI * i) / size);

    data[2 * i] =  (re_H_n + re_H_N) + (im_H_n + im_H_N) * cos_val + 
                      (re_H_n - re_H_N) * sin_val;
    data[2 * i] /= 2;

    data[(2 * i) + 1] =  (im_H_n - im_H_N) + (im_H_n + im_H_N) * sin_val + 
                            (- re_H_n + re_H_N) * cos_val;
    data[(2 * i) + 1] /= 2;
  }

  data[0] = tmp_data[0] + tmp_data[1];
  data[1] = 0;

  return 0;
}


#endif
