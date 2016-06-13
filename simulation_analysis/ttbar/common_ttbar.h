#include "tdrstyle_mod14.C"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::cerr;

int zeros = 0;
const int ptBins = 45.;//29.;//61.;
const double ptRange[]=
    //{18, 21, 24,·
    {28, 32, 37, 43, 49,
     56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
     2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450,
     3637, 3832, 4037};

static double pz_calc(TLorentzVector& lepton, TLorentzVector& MET, TLorentzVector& neutrino) {
    double l = -pow(lepton.M(),2)+pow(80.13,2)+2*( MET.Px()*lepton.Px()+MET.Py()*lepton.Py() );
    double m = 2*( pow(lepton.E(),2)-pow(lepton.Pz(),2) );

    double t1 = pow(l,2)-2*m*(MET.Px()*MET.Px()+MET.Py()*MET.Py());
    double t2 = l*lepton.Pz()/m;
    if (t1<0) {
        //cerr << "Less-than-zero" << ++zeros << endl;
        return t2;
    }
    double t3 = fabs(lepton.E()*sqrt(t1)/m);
    double pz1 = t2+t3;
    double pz2 = t2-t3;
    double pz = (fabs(pz1-lepton.Pz())<fabs(pz2-lepton.Pz())) ? pz1 : pz2;

    return pz;
}

static inline bool compatibility(double mass_sum, double mass_diff) {
    return (mass_sum < 400 && mass_sum > 300 && mass_diff < 50);
}

static inline int mass_study(double m1, double m2, double n1, double n2, bool noisy, int& id, double& diff) {
    double sum_1 = m1+n2, sum_2 = m2+n1;
    double diff_1 = fabs(m1-n2), diff_2 = fabs(m2-n1);
    unsigned mult = 0;

    if (compatibility(sum_1,diff_1)) {
        if (noisy)
            cout << "Lepton t " << m1 << "Jet t " << n2 << endl;
        id = 0;
        diff = diff_1;
        ++mult;
    }
    if (compatibility(sum_2,diff_2)) {
        if (noisy)
            cout << "Lepton t " << m2 << " Jet t " << n1 << endl;
        id = 1;
        diff = diff_2;
        ++mult;
        if (mult>1) {
            if (diff_1<diff_2) {
                id = 0;
                diff = diff_1;
            }
        }
    }

    return mult;
}


double upperW = 110.0;
double lowerW = 60.0;
unsigned noW = 180;
double upperT = 220.0;
double lowerT = 140.0;
unsigned noT = 200;

unsigned goods = 0;
unsigned bads = 0;
