#include "fft2f.hpp"
#include "ImageIO.hpp"

#include <iostream>
#include <fstream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    string inimage = argv[1];

    Image2Float( inimage, "signal.dat");

    fft2f(  "signal.dat", "fft.dat");

    fftImage( "fft.dat", "fft.pgm");
    fftImage( "fft.dat", "cfft.pgm", 1);

    ifft2f(  "fft.dat", "ifft.dat");
    ifftImage( "ifft.dat", "ifft.pgm");

    cout << "MaxDiff " << fftdiff2( "signal.dat", "ifft.dat") << endl;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

