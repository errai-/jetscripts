#include <Hannouris/QCDPFJet.h>
#include <Hannouris/QCDJet.h>
#include <Hannouris/LorentzVector.h>
#include <vector>
#include <iostream>

class Auxiliary{
  private:
    std::vector<QCDPFJet> vals;

  public:

    Auxiliary(){ }

    ~Auxiliary(){ }

    void clear() { vals.clear(); }

    void addValue( double px, double py, double pz, double e, double chf, 
      double nhf, double phf, double elf, double muf, double hfhf, double hfphf,
      double tightID){
      QCDPFJet tmpJet;
      LorentzVector tmpVector;
      tmpVector.setX( px ); tmpVector.setY( py ); 
      tmpVector.setZ( pz ); tmpVector.setT( e );

      tmpJet.setP4( tmpVector );
      tmpJet.setFrac(chf, nhf, phf, elf, muf);
      tmpJet.setHFFrac( hfhf, hfphf );
      tmpJet.setTightID( tightID );

      vals.push_back( tmpJet );
    }

    std::vector<QCDPFJet> values(){return vals; } 
};
