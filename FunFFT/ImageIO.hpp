#pragma once

#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <cassert>
#include <iostream>

double fftdiff2(const std::string &signal, const std::string &mask);

template<typename T>
std::ostream& write_binary(std::ostream& stream, const T& value) {
    return stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T>
std::ostream& write_binary(std::ostream& stream, const std::vector<T>& value) {
    size_t nSize = value.size();
    return stream.write(reinterpret_cast<const char*>(&value[0]), nSize*sizeof(T));
}

template<typename T>
std::istream & read_binary(std::istream& stream, T& value) {
    return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template<typename T>
std::istream & read_binary(std::istream& stream, std::vector<T>& value) {
    size_t nSize = value.size();
    return stream.read(reinterpret_cast<char*>(&value[0]), nSize*sizeof(T));
}

template<class T>
std::vector<T> circShift(const std::vector<T> &in, std::array<int,2> &dim)
{
    int width   = dim[0];
    int height  = dim[1];
    int xshift  = width/2;
    int yshift  = height/2;

    assert(in.size() == width*height);
    std::vector<T> out(width*height);

    for (int j =0; j < height; j++) {
        int jj = (j+ yshift)%height;
        for (int i =0; i < width; i++) {
            int ii = (i+ xshift)%width;
            out[jj*width + ii] = in[j*width+i];
        }
    }
    return out;
}

template<class T>
std::vector<T> reflect(const std::vector<T> &in, std::array<int,2> &dim)
{
    int width   = dim[0];
    int height  = dim[1];

    assert(in.size() == width*height);
    std::vector<T> out(width*height);

    for (int j =0; j < height; j++) {
        int jj = (-j+ height-1)%height;
        for (int i =0; i < width; i++) {
            int ii = (-i+ width-1)%width;
            out[jj*width + ii] = in[j*width+i];
        }
    }
    return out;
}

int  savePGM(const std::string &filename, const std::vector<unsigned char> &image, const std::array<int,2> &dim);
int  readPGM(const std::string &filename, std::vector<unsigned char> &image, std::array<int,2> &dim);
std::array<int,2> getImageSize( std::string &filename);

void padImage( const std::string &f, const std::array<int,2> &p, const std::string &outfile);
void fftImage( const std::string &sigfile, const std::string &kerfile, int shift = 0);
void ifftImage( const std::string &sigfile, const std::string &outfile);

void Image2Float(const std::string &infile,  const std::string &outfile);
