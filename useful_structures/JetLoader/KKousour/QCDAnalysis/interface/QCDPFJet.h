#ifndef QCDPFJet_h
#define QCDPFJet_h
#include "KKousour/QCDAnalysis/interface/QCDJet.h"

class QCDPFJet : public QCDJet {
   public:
     //------------ Constructor ------------------------------
     QCDPFJet() {chf_=0;nhf_=0;phf_=0;elf_=0;muf_=0;chm_=0;nhm_=0;phm_=0;elm_=0;mum_=0;}
     //------------ Destructor -------------------------------
     ~QCDPFJet() {}
     //------------ Set methods ------------------------------
     void setFrac(float fchf, float fnhf, float fphf, float felf, float fmuf)  {chf_ = fchf; nhf_ = fnhf; phf_ = fphf; elf_ = felf; muf_ = fmuf;}
     void setMulti(int fncand, int fchm, int fnhm, int fphm, int felm, int fmum) {ncand_ = fncand; chm_ = fchm; nhm_ = fnhm; phm_ = fphm; elm_ = felm; mum_ = fmum;}
     void setBeta(float fbeta) {beta_ = fbeta;}
     void setBetaStar(float fbetaStar) {betaStar_ = fbetaStar;}
     void setHFFrac(float fhf_hf, float fhf_phf) {hf_hf_ = fhf_hf; hf_phf_ = fhf_phf;}
     void setHFMulti(int fhf_hm, int fhf_phm) {hf_hm_ = fhf_hm; hf_phm_ = fhf_phm;}
     //------------ Get methods ------------------------------ 
     float beta()     const {return beta_;}                
     float betaStar() const {return betaStar_;}
     float chf()      const {return chf_;} 
     float nhf()      const {return nhf_;}
     float phf()      const {return phf_;} 
     float elf()      const {return elf_;}
     float muf()      const {return muf_;}
     float hf_hf()    const {return hf_hf_;}
     float hf_phf()   const {return hf_phf_;}
     int chm()        const {return chm_;}
     int nhm()        const {return nhm_;}
     int phm()        const {return phm_;}
     int elm()        const {return elm_;}
     int mum()        const {return mum_;}
     int hf_hm()      const {return hf_hm_;}
     int hf_phm()     const {return hf_phm_;}
     int ncand()      const {return ncand_;}
   private:
     //---- charged hadron energy fraction ----
     float chf_;
     //---- neutral hadron energy fraction ----
     float nhf_;
     //---- photon energy fraction ------------
     float phf_;
     //---- electron energy fraction ----------
     float elf_;
     //---- muon energy fraction --------------
     float muf_;
     //-----HF Hadron Fraction ---------------
     float hf_hf_;
     //-----HF Photon Fraction ------------
     float hf_phf_;
     //-----HF Hadron Multiplicity---------
     int hf_hm_;
     //-----HF Photon Multiplicity --------
     int hf_phm_;
     //---- charged hadron multiplicity -------
     int chm_;
     //---- neutral hadron multiplicity -------
     int nhm_;
     //---- photon multiplicity ---------------
     int phm_;
     //---- electron multiplicity -------------
     int elm_;
     //---- muon multiplicity -----------------
     int mum_;
     //---- number of PF candidates -----------
     int ncand_;
     //---- fraction of track pt coming from the signal vertex ---
     float beta_;
     //---- fraction of track pt NOT coming from the signal vertex ---
     float betaStar_;    
    };
#endif    
