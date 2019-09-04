#pragma once

class FEDZSClusterUnpacker
{
public:

  static FEDZSClusterUnpacker zeroSuppressedModeUnpacker(const FEDChannel& channel);
  FEDZSClusterUnpacker();
  uint8_t firstSampleNumber() const;
  uint8_t adc(uint8_t n) const;
  uint8_t nSamples() const {return valuesInCluster_;}
  bool readNext();

private:
  //pointer to beginning of FED or FE data, offset of start of channel payload in data and length of channel payload
  FEDZSClusterUnpacker(const uint8_t* payload, const uint16_t channelPayloadOffset, const int16_t channelPayloadLength, const uint16_t offsetIncrement=1);
  bool hasData(uint16_t extra = 0) const;
  void readNewClusterInfo();
  const uint8_t* data_;
  uint16_t currentOffset_;
  uint16_t offsetIncrement_;
  uint8_t currentStrip_;
  uint8_t valuesInCluster_;
  uint16_t channelPayloadOffset_;
  uint16_t channelPayloadLength_;
};

inline void FEDZSClusterUnpacker::readNewClusterInfo()
{
  currentStrip_ = data_[(currentOffset_++)^7];
  valuesInCluster_ = data_[(currentOffset_++)^7];
}


inline
FEDZSClusterUnpacker::FEDZSClusterUnpacker(const uint8_t* payload, const uint16_t channelPayloadOffset, const int16_t channelPayloadLength, const uint16_t offsetIncrement)
: data_(payload),
  currentOffset_(channelPayloadOffset),
  offsetIncrement_(offsetIncrement),
  currentStrip_(0),
  valuesInCluster_(0),
  channelPayloadOffset_(channelPayloadOffset),
  channelPayloadLength_(channelPayloadLength)
{

}

inline FEDZSClusterUnpacker FEDZSClusterUnpacker::zeroSuppressedModeUnpacker(const FEDChannel& channel)
{
  uint16_t length = channel.length();
  FEDZSClusterUnpacker result(channel.data()-channel.offset(),channel.offset()+7,length-7);
  return result;
}

inline uint8_t FEDZSClusterUnpacker::firstSampleNumber() const
{
  return currentStrip_;
}

inline uint8_t FEDZSClusterUnpacker::adc(uint8_t n) const
{
  return data_[(currentOffset_+n)^7];
}

inline bool FEDZSClusterUnpacker::hasData(uint16_t extra) const
{
  return (currentOffset_+valuesInCluster_*offsetIncrement_+extra < channelPayloadOffset_+channelPayloadLength_);
}

inline bool FEDZSClusterUnpacker::readNext() {
  if ( valuesInCluster_ > 0 ) currentOffset_ += valuesInCluster_*offsetIncrement_;
  if (hasData(2)) {
    readNewClusterInfo();
    return true;
  }
  return false;
}


