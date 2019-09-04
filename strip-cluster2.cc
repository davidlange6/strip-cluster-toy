#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>

#include "Clusterizer.h"
#include "FEDChannel.h"
#include "FEDZSClusterUnpacker.h"

const unsigned int MAX_CLUSTERS=500000;
const unsigned int MAX_SIZE=255;



FEDSet fillFeds()
{
  std::ifstream fedfile("stripdata.bin", std::ios::in | std::ios::binary);

  FEDSet feds;
  detId_t detid;  

  while (fedfile.read((char*)&detid, sizeof(detid)).gcount() == sizeof(detid)) {
    FEDChannel fed(fedfile);
    feds[detid].push_back(std::move(fed));
  }
  return feds;
}

unsigned int
fillStructures(detId_t idet, Clusterizer& clusterizer, const std::vector<FEDChannel>& channels, 
	       unsigned int c,
	       std::vector<uint16_t> &firstStrip,
	       std::vector<uint16_t> &lastStrip,
	       std::vector<std::array<uint8_t,MAX_SIZE> > &adcs,
	       std::vector<detId_t> &detids,
	       std::ofstream& digidata_out)
{
  for (auto const& chan : channels) {
    auto unpacker = FEDZSClusterUnpacker::zeroSuppressedModeUnpacker(chan);
    while ( unpacker.readNext() ) {
      firstStrip[c] = unpacker.firstSampleNumber() + chan.iPair()*256;
      lastStrip[c] = firstStrip[c] + unpacker.nSamples() -1;
      for ( uint8_t i=0; i< unpacker.nSamples(); i++ ) adcs[c][i]=unpacker.adc(i);
      detids[c]=idet;
      c++;
    }
  }
  std::cout << "nClusters " << c << std::endl;
  return c;
}

int main()
{
  //std::cout << "start initializing clusterizer" << std::endl;
  Clusterizer clusterizer;
  //  Clusterizer::State state;

  std::ofstream digidata_out("digidata.bin", std::ofstream::out | std::ios::binary);
  std::ifstream digidata_in("digidata.bin", std::ofstream::in | std::ios::binary);

  std::vector<uint16_t> firstStrip(MAX_CLUSTERS),lastStrip(MAX_CLUSTERS);
  std::vector<std::array<uint8_t,MAX_SIZE> > adcs(MAX_CLUSTERS);
  std::vector<detId_t> detids(MAX_CLUSTERS);

  FEDSet feds(fillFeds());
  unsigned int c=0;
  for (auto idet : clusterizer.allDetIds()) {
    if (feds.find(idet) != feds.end()) {
      c=fillStructures(idet, clusterizer, feds[idet], c, firstStrip, lastStrip, adcs, detids, digidata_out);
    }
  }

  for ( unsigned int i=0; i< 100; i++ ) {
    std::cout << i << " " << detids[i] << " " << firstStrip[i] << " " << lastStrip[i] << " " << (unsigned int)adcs[i][0] << std::endl;
  }


  return 0;

}
