# g++ gw_fin.cc -std=c++23 -I/opt/eigen-3.4.0/  -L/mnt/hdd_barracuda/opt//myhdfstuff/hdf5-2.1.0/lib/ -lhdf5
CXX=g++
CXXFLAGS=-Wall -Wextra -std=c++23 -O3 -fopenmp
INCLUDES=-I"/mnt/hdd_barracuda/opt/eigen/" -I"/mnt/hdd_barracuda/opt/highfive/build/install/include/" -I"/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
# LDFLAGS=-L"/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/" -L"/usr/lib/" -L"/usr/local/lib/" -lhdf5 -lgsl -lgslcblas -lm # -lhdf5_cpp
LDFLAGS=-L"/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/" -lhdf5 # -lhdf5_cpp

# g++ dummy.cc -I/mnt/hdd_barracuda/opt/highfive/build/install/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/ -L/mnt/hdd_barracuda/opt//myhdfstuff/hdf5-2.1.0/lib/ -lhdf5
# INCLUDES=-I"/Users/nobuyukimatsumoto/opt/highFive/include/" -I"/Users/nobuyukimatsumoto/opt/hdf5-v1.14.0/include/" -I"/Users/nobuyukimatsumoto/opt/eigen/" -I"/usr/local/include/"
# # INCLUDES=-I"/Users/nobuyukimatsumoto/opt/hdf5-v1.14.0/include/" -I"/Users/nobuyukimatsumoto/opt/eigen/"
# LDFLAGS=-L"/Users/nobuyukimatsumoto/opt/hdf5-v1.14.0/lib/" -lhdf5
# INCLUDES=-I"/Users/nobuyukimatsumoto/opt/eigen/"
# LDFLAGS=

SRCS := $(wildcard *.cc)
OBJS := $(SRCS:%.cc=%.o)

all: $(OBJS)

%.d: %.cc
	$(CXX) $< -o $@ $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -M

include $(SRCS:.cc=.d)

%.o: %.cc
	$(CXX) $< -o $@ $(CXXFLAGS) $(INCLUDES) $(LDFLAGS)

# all:
# 	g++ hist_spline.cc ${INCLUDES} ${LDFLAGS}
# 	g++ get_potential.cc ${INCLUDES} ${LDFLAGS}
