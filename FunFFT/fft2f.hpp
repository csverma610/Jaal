#pragma once

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <array>
#include <fftw3.h>
#include <iomanip>
#include <fstream>
#include <boost/timer/timer.hpp>
#include <cassert>
#include <algorithm>

struct FFT2F
{
public:
    void readImage( const std::string &filename);
    void setSpatialData( std::vector<float> &data,  std::array<int,2> &dim);
    void setFrequencyData( std::vector<float> &data, std::array<int,2> &dim);

    void setNumThreads(int n) { nthreads = n; }

    void fft( bool newplan = 0);
    void ifft( bool newplan = 0);

    int  getWidth()  const;
    int  getHeight() const;

    std::vector<float> getFFTData();
    std::vector<float> getIFFTData();

private:
    bool inplace  = 1;
    int  nthreads = 8;
    std::array<int,2> dim;

    fftwf_complex   *inData    = nullptr;
    fftwf_complex   *fftData   = nullptr;
    fftwf_complex   *ifftData  = nullptr;

    fftwf_plan  forwardPlan, backwardPlan;
};

void  fft2f(const std::string &infile, std::vector<float> &data, std::array<int,2> &dim);
void  fft2f(const std::string &infile, const std::string &outfile);

void  ifft2f(const std::string &infile, std::vector<float> &data, std::array<int,2> &dim);
void  ifft2f(const std::string &infile, const std::string &outfile);

void  conv2f(const std::string &signal, const std::string &mask, const std::string &outfile, bool expand);
void  corr2f(const std::string &signal, const std::string &mask, const std::string &outfile);

void  genRandom( int Nx, int Ny, const std::string &filename);
