CPPFLAGS = -O3 -fopenmp
CPPFLAGS += -I$(FFTW_DIR)/include 
CPPFLAGS += -I$(BOOST_DIR)/include 

LIBS1   = -L$(FFTW_DIR)/lib -lfftw3f -lfftw3f_omp
LIBS1  += -L$(BOOST_DIR)/lib -lboost_timer
LIBS1  += -lgomp

OBJS2 = main.o fft2f.o ImageIO.o

fft2:$(OBJS2)
	g++ -o fft2 $(OBJS2) $(LIBS1)

LIBS2   = -L$(FFTW_DIR)/lib -lfftw3 -lfftw3_omp
LIBS2  += -L$(BOOST_DIR)/lib -lboost_timer
LIBS2  += -lgomp

.o:.cpp
	g++ $(CPPFLAGS) $<

clean:
	\rm -rf *.o fft2
