#include "ImageIO.hpp"
#include <iostream>

using namespace std;

//
//////////////////////////////////////////////////////////////////////////////////////////////
//
int savePGM(const string &filename, const vector<unsigned char> &image, const array<int,2> &dim)
{
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) return 1;

    int width  = dim[0];
    int height = dim[1];

    ofile << "P5 \n";
    ofile << width << " " << height << "\n";

    auto val = *max_element( image.begin(), image.end());

    ofile << (int)val << "\n";

    write_binary(ofile, image);
    return 0;
}
//
//////////////////////////////////////////////////////////////////////////////////////////////
//
std::array<int,2> getImageSize( const string &filename)
{
    std::array<int,2> dim = {0,0};
    ifstream ifile(filename.c_str(), ios::out);
    if( ifile.fail() ) return dim;

    string str;
    ifile >> str;
    assert( str == "P5");

    ifile >> dim[0] >> dim[1];
    return dim;
}

//////////////////////////////////////////////////////////////////////////////////////////////

int readPGM(const string &filename, vector<unsigned char> &image, array<int,2> &dim)
{
    ifstream ifile(filename.c_str(), ios::out);
    if( ifile.fail() ) return 1;

    string str;
    ifile >> str;
    assert( str == "P5");

    int width, height;
    ifile >> width >> height;

    int maxVal;
    ifile >> maxVal;

    image.resize(width*height);
    read_binary(ifile, image);

    dim[0] = width;
    dim[1] = height;

    if( maxVal == 1) {
        for( size_t i = 0; i < image.size(); i++)
            image[i] *= 255;
    }

    return 0;
}
//
//////////////////////////////////////////////////////////////////////////////////////////////
//
void padImage(const string &infile, const array<int,2> &paddim, const string &outfile)
{
    vector<unsigned char> inbuf;
    array<int,2>    indim;

    readPGM(infile, inbuf, indim);

    int width  = indim[0];
    int height = indim[1];

    assert( inbuf.size()  == width*height);

    int pwidth  = width  + 2*paddim[0];
    int pheight = height + 2*paddim[1];

    vector<unsigned char> padbuf(pwidth*pheight, 0);

    size_t index = 0;
    for( size_t j = 0; j <  height; j++) {
        int jj = j + paddim[1];
        for( size_t i = 0; i <  width;  i++) {
            int ii = i + paddim[0];
            size_t offset = jj*pwidth + ii;
            padbuf[offset] = inbuf[index++];
        }
    }

    std::array<int,2> newdim = {pwidth, pheight};
    savePGM(outfile, padbuf, newdim);
}
//
//////////////////////////////////////////////////////////////////////////////////////////////
//
void fftImage( const string &datafile, const string &imgfile, int shift)
{
    ifstream ifile(datafile.c_str(), ios::binary);
    if( ifile.fail() ) {
        cout << "Error: Can not open file " << datafile << endl;
        return;
    }

    int width, height;
    read_binary( ifile, width);
    read_binary( ifile, height);

    size_t numPixels = width*height;

    vector<float> buf(2*numPixels);
    read_binary( ifile, buf);

    vector<float> data(numPixels);

    for( size_t i = 0; i < numPixels; i++) {
        double x = buf[2*i];
        double y = buf[2*i+1];
        data[i] = log(sqrt(x*x + y*y) + 1.0E-01);
    }

    std::array<int,2> dim = {width, height};

    if(shift) data = circShift(data, dim);

    float minval = *min_element( data.begin(), data.end());
    float maxval = *max_element( data.begin(), data.end());

    vector<unsigned char> imgdata(numPixels);
    for( size_t i = 0; i < numPixels; i++) {
        double val = 255*(data[i] - minval)/(maxval-minval);
        if( val < 0.0) val = 0;
        if( val > 255) val = 255.0;
        imgdata[i] = val;
    }

    savePGM( imgfile, imgdata, dim);
}
//
//////////////////////////////////////////////////////////////////////////////////////////////
//
void ifftImage( const string &datafile, const string &imgfile)
{
    ifstream ifile(datafile.c_str(), ios::binary);
    if( ifile.fail() ) {
        cout << "Error: Can not open file " << datafile << endl;
        return;
    }

    int width, height;
    read_binary(ifile, width);
    read_binary(ifile, height);

    size_t numPixels = width*height;
    vector<float> data(numPixels);
    read_binary( ifile, data);

    float minval = *min_element( data.begin(), data.end());
    float maxval = *max_element( data.begin(), data.end());

    vector<unsigned char> imgdata(numPixels);
    for( size_t i = 0; i < numPixels; i++)
        imgdata[i] = 255*(data[i] - minval)/(maxval-minval);

    std::array<int,2> dim = {width, height};

    savePGM( imgfile, imgdata, dim);
}
//
//////////////////////////////////////////////////////////////////////////////////////////////
//
void Image2Float(const string &infile, const string &outfile)
{
    vector<unsigned char> cdata;
    array<int,2> dim;
    readPGM( infile, cdata, dim);

    ofstream ofile(outfile.c_str(), ios::out);
    if( ofile.fail() ) return;

    write_binary( ofile, dim[0]);
    write_binary( ofile, dim[1]);

    vector<float> fdata( dim[0]*dim[1] );
    for( size_t i = 0; i < cdata.size(); i++)
        fdata[i] = cdata[i]/255.0;

    write_binary( ofile, fdata);
}
//
//////////////////////////////////////////////////////////////////////////////////////////////
//
double fftdiff2( const string &afile, const string &bfile)
{
    ifstream ifile1(afile.c_str(), ios::binary);
    if( ifile1.fail() ) {
        cout << "Error: Can not open file " << afile << endl;
        return 0.0;
    }

    ifstream ifile2(bfile.c_str(), ios::binary);
    if( ifile2.fail() ) {
        cout << "Error: Can not open file " << bfile << endl;
        return 0.0;
    }
    int Nx, Ny;
    read_binary( ifile1, Nx);
    read_binary( ifile1, Ny);

    int Nx1, Ny1;
    read_binary( ifile2, Nx1);
    read_binary( ifile2, Ny1);

    assert( Nx == Nx1);
    assert( Ny == Ny1);

    vector<float> buf1(Nx*Ny);
    read_binary( ifile1, buf1);

    vector<float> buf2(Nx*Ny);
    read_binary( ifile2, buf2);

    float maxDiff = 0.0;
    for( size_t i = 0; i < Nx*Ny; i++)
        maxDiff = max(maxDiff, fabs(buf1[i]-buf2[i] ));

    return maxDiff;
}
