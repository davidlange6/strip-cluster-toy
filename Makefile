CXXFLAGS += -std=c++14 -g #-O3
LDFLAGS += -std=c++14
CC = c++
all : strip-cluster strip-cluster2
strip-cluster : strip-cluster.o Clusterizer.o FEDChannel.o
strip-cluster.o: strip-cluster.cc Clusterizer.h FEDChannel.h FEDZSChannelUnpacker.h
strip-cluster2 : strip-cluster2.o Clusterizer.o FEDChannel.o
strip-cluster2.o: strip-cluster2.cc Clusterizer.h FEDChannel.h FEDZSClusterUnpacker.h
Clusterizer.o: Clusterizer.cc Clusterizer.h
FEDChannel.o : FEDChannel.cc FEDChannel.h Clusterizer.h
