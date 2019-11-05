#include "fft2f.hpp"
#include "ImageIO.hpp"

using namespace std;

bool inplace  = 1;
int  nthreads = 8;

//////////////////////////////////////////////////////////////////////////////////////////////

void genRandom( int Nx, int Ny, const string &filename)
{
    ofstream ofile(filename.c_str(), ios::binary);

    write_binary(ofile, Nx);
    write_binary(ofile, Ny);

    vector<float> buf(Nx*Ny);

    for( int i = 0; i < Nx*Ny; i++) buf[i] = 0;

    int index = 0;
    for( int i = 0; i < Nx; i++) {
        for( int j = 0; j  < Ny; j++) {
            int dx = abs(i - 0.5*Nx);
            int dy = abs(j - 0.5*Ny);
            float val  = drand48();
            index++;
        }
    }
    write_binary( ofile, buf);
    ofile.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////

void fft2f( const string &infilename, const string &outfilename)
{
    boost::timer::cpu_timer fft_watch;
    fftwf_init_threads();
    fftwf_plan_with_nthreads(nthreads);

    // Input and Plan Creatiion ...
    ifstream ifile(infilename.c_str(), ios::binary);
    int width, height;
    read_binary(ifile, width);
    read_binary(ifile, height);

    size_t nbytes = size_t(width)*size_t(height)*sizeof(fftwf_complex);
    auto in  = (fftwf_complex*) fftwf_malloc( nbytes);
    auto out = in;

    if( !inplace ) {
        out = (fftwf_complex*) fftwf_malloc( nbytes);
    }

    boost::timer::cpu_timer plan_watch;

    auto plan = fftwf_plan_dft_2d( width, height, &in[0], &out[0], FFTW_FORWARD, FFTW_MEASURE);

    cout << "Time: Plan creation  : " << plan_watch.elapsed().wall*1.0E-09 << " seconds " << endl;

    boost::timer::cpu_timer read_watch;

    vector<float> buf(height*width);
    read_binary(ifile, buf);

    size_t index = 0;
    for( int i = 0; i < height; i++) {
        for( int j = 0; j < width; j++) {
            in[index][0] = buf[index];
            in[index][1] = 0.0;
            index++;
        }
    }
    ifile.close();
    cout << "Time: Read Data      : " << read_watch.elapsed().wall*1.0E-09 << " seconds " << endl;

    // Core execution ...
    boost::timer::cpu_timer exe_watch;
    fftwf_execute ( plan );

    cout << "Time: Execution      : " << exe_watch.elapsed().wall*1.0E-09 << " seconds " << endl;

    fftwf_destroy_plan( plan );

    // Output ....
    boost::timer::cpu_timer write_watch;

    buf.resize(2*height*width);

    index = 0;
    for( int i = 0; i < height; i++) {
        for( int j = 0; j < width; j++) {
            buf[2*index]   =  out[index][0];
            buf[2*index+1] =  out[index][1];
            index++;
        }
    }

    ofstream ofile(outfilename.c_str(), ios::binary);
    write_binary( ofile, width);
    write_binary( ofile, height);
    write_binary(ofile, buf);
    ofile.close();

    cout << "Time: Write Data     : " << write_watch.elapsed().wall*1.0E-09 << " seconds " << endl;

    fftwf_cleanup_threads();
    cout << "Total Time: FFT      : " << fft_watch.elapsed().wall*1.0E-09 << " seconds " << endl;
    fftwf_free(in);

    if( !inplace ) fftwf_free(out);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void ifft2f( const string &infilename, const string &outfilename)
{
    cout << "IFFT ... " << endl;
    boost::timer::cpu_timer fft_watch;
    fftwf_init_threads();
    fftwf_plan_with_nthreads(nthreads);

    // Input and Plan Creatiion ...
    ifstream ifile(infilename.c_str(), ios::binary);

    int width, height;
    read_binary( ifile, width);
    read_binary( ifile, height);

    size_t nbytes = size_t(width)*size_t(height)*sizeof(fftwf_complex);

    auto in  = (fftwf_complex*) fftwf_malloc( nbytes);
    auto out = in;

    if( !inplace ) {
        out = (fftwf_complex*) fftwf_malloc( nbytes);
    }

    boost::timer::cpu_timer plan_watch;

    auto plan = fftwf_plan_dft_2d( width, height, &in[0], &out[0], FFTW_BACKWARD, FFTW_MEASURE);

    cout << "Time: Plan creation  : " << plan_watch.elapsed().wall*1.0E-09 << " seconds " << endl;

    boost::timer::cpu_timer read_watch;

    vector<float> buf(2*width*height);
    read_binary(ifile, buf);

    int index = 0;
    for( int i = 0; i < height; i++) {
        for( int j = 0; j < width; j++) {
            in[index][0] = buf[2*index];
            in[index][1] = buf[2*index+1];
            index++;
        }
    }
    cout << "Time: Read Data      : " << read_watch.elapsed().wall*1.0E-09 << " seconds " << endl;

    boost::timer::cpu_timer exe_watch;
    fftwf_execute ( plan );
    cout << "Time: Execution      : " << exe_watch.elapsed().wall*1.0E-09 << " seconds " << endl;

    fftwf_destroy_plan( plan );

    // Output ....
    boost::timer::cpu_timer write_watch;

    buf.resize(width*height);

    index = 0;
    for( int i = 0; i < height; i++) {
        for( int j = 0; j < width; j++) {
            buf[index++] = out[index][0]/(double)width/(double)height;
        }
    }

    ofstream ofile(outfilename.c_str(), ios::binary);
    write_binary( ofile, width);
    write_binary( ofile, height);
    write_binary( ofile, buf);
    ofile.close();

    cout << "Time: Write Data     : " << write_watch.elapsed().wall*1.0E-09 << " seconds " << endl;
    fftwf_cleanup_threads();
    fftwf_free(in);
    cout << "Total Time: IFFT     : " << fft_watch.elapsed().wall*1.0E-09 << " seconds " << endl;

    if( !inplace ) fftwf_free(out);
}

////////////////////////////////////////////////////////////////////////////////

void conv2( const string &afile, const string &bfile, const string &outfilename, bool expand)
{

	/*
   std::array<int,2> adim = getImageSize(afile);
   std::array<int,2> bdim = getImageSize(bfile);

   string kernel = bfile;
   if( adim[0]*adim[1] > bdim[0]*bdim[1] ) {
       enlargeImage( bfile, adim, "kernel.pgm");
       kernel = "kernel2.pgm"l
   }
   bdim = getImageSize(kernel);

   assert( adim[0] == bdim[0] );
   assert( adim[1] == bdim[1] );

   fft2f( afile,  "afft.dat");
   fft2f( kernel, "bfft.dat");
*/
}

////////////////////////////////////////////////////////////////////////////////
