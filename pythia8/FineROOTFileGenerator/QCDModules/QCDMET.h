#ifndef QCDMET_h
#define QCDMET_h

class QCDMET 
{
   public:
     //------------ Constructor ------------------------------
     QCDMET() {et_=0; sumEt_=0;}
     //------------ Destructor -------------------------------
     ~QCDMET() {}
     //------- Set method ------------------------------------
     void setVar(float fEt, float fSumEt, float fPhi) {et_ = fEt; sumEt_ = fSumEt; phi_ = fPhi;} 
     //------- Get methods -----------------------------------
     float met()         const {return et_;}
     float phi()         const {return phi_;}
     float sumet()       const {return sumEt_;}
     float met_o_sumet() const {return et_/sumEt_;}
     
   private:
     //---- size of MET vector ----------
     float et_;
     //---- sumET -----------------------
     float sumEt_;
     //---- phi of MET vector -----------
     float phi_;
};
#endif
