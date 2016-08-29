#include "../include/PythiaSaver.h"

PythiaSaver::PythiaSaver(int _nEvent, char *_fileName):nEvent(_nEvent){
  if (weightedPt){
    pythia.setUserHooksPtr( &ptGenReweight );
  }
  
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  pythia.readString("particleDecays:limitTau0=on");
  pythia.readString("particleDecays:tauMax=10.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("PartonLevel:ISR=on");
  //pythia.particleData.listAll();

  pythia.readString("Beams:eCM = 14000.");
  pythia.init();
  pythia.settings.listChanged();
  
  output.open(_fileName, std::ofstream::out | std::ofstream::app);

  timer.set_params(nEvent,100);

  output << std::fixed << std::setprecision(60);

  timer.start_timing();  
}

void PythiaSaver::EventLoop(){
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    // event.list();
    
    if (iEvent!=0&&iEvent%100==0){
      timer.print_time();
    }
    
    ParticleLoop();
    
    tmpEvent.Write(&output);
    tmpEvent.Nullify();
  }
}
  
void PythiaSaver::ParticleLoop(){
  for (int i = 0; i != event.size(); ++i) {
    double status = abs( event[i].status() );
    if (event[i].isFinal() && event[i].isVisible()) {
      int tmpId = event[i].id();
      tmpId -= ( ((tmpId == 22) && gammaChecker( event, i )) ? 2 : 0 ); // Indicate pi0 photons
      tmpEvent.SetVals(event[i].px(),event[i].py(),event[i].pz(),event[i].e(), tmpId );
    }
  }
}