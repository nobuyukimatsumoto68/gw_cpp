CXX=g++-13
CXXFLAGS=-Wall -Wextra -std=c++23 -O3 -fopenmp
# INCLUDES=-I"/mnt/hdd_barracuda/opt/eigen/" -I"/mnt/hdd_barracuda/opt/HighFive/include/" -I"/mnt/hdd_barracuda/opt/hdf5-v1.14.0/include/" -I"~/opt/eigen/"
# LDFLAGS=-L"/mnt/hdd_barracuda/opt/hdf5-v1.14.0/lib/" -L"/usr/lib/" -L"/usr/local/lib/" -lhdf5 -lgsl -lgslcblas -lm # -lhdf5_cpp
INCLUDES=-I"/Users/nobuyukimatsumoto/opt/highFive/include/" -I"/Users/nobuyukimatsumoto/opt/hdf5-v1.14.0/include/" -I"/Users/nobuyukimatsumoto/opt/eigen/" -I"/usr/local/include/"
# INCLUDES=-I"/Users/nobuyukimatsumoto/opt/hdf5-v1.14.0/include/" -I"/Users/nobuyukimatsumoto/opt/eigen/"
LDFLAGS=-L"/Users/nobuyukimatsumoto/opt/hdf5-v1.14.0/lib/" -lhdf5
# INCLUDES=-I"/Users/nobuyukimatsumoto/opt/eigen/"
# LDFLAGS=

SRCS := $(wildcard *.cc)
OBJS := $(SRCS:%.cc=%.o)

all: $(OBJS)

%.d: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -M $< -o $@

include $(SRCS:.cc=.d)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $< -o $@ 

# all:
# 	g++ hist_spline.cc ${INCLUDES} ${LDFLAGS}
# 	g++ get_potential.cc ${INCLUDES} ${LDFLAGS}
