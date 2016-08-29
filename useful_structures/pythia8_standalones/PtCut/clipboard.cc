

feh->SetParameters(0.5, 0.5,
       1.4160, -1.0004, -0.1800, -2.19, // a
       0.9656, -0.5652, -0.1800, -0.18);// b

double a = max(0.33, (pa[0]+pa[1]*pow(T+pa[3],pa[2])) * max(0., T+pa[3])/T);
double b = max(0.33, (pb[0]+pb[1]*pow(T+pb[3],pb[2])) * max(0., T+pb[3])/T);


// H
TF1 *fhr = new TF1("fhr","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",4,200);
fhr->SetParameters(1.0806,-0.5,-0.15,-0.3);
fhr->FixParameter(0, 1.0806); // constrain to old VLE+ME sample
fhr->SetParLimits(2,-0.18,-0.13);
cout << "*** Fitting fhr ***" << endl;
ghr->Fit(fhr,"RN");
fhr->DrawClone("SAME");
fhr->SetLineStyle(kDashed);
fhr->SetRange(-fhr->GetParameter(3),200);
fhr->DrawClone("SAME");

// HCAL
TF1 *frh = new TF1("frh","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
frh->SetParameters(0.9*fhr->GetParameter(0), fhr->GetParameter(1),
       fhr->GetParameter(2), fhr->GetParameter(3));
frh->FixParameter(1,fhr->GetParameter(1));
frh->SetParLimits(2,-0.18,-0.13);
frh->SetLineColor(kMagenta+1);//Red);
cout << "*** Fitting frh ***" << endl;
grh->Fit(frh,"RN");
frh->DrawClone("SAME");
frh->SetLineStyle(kDashed);
frh->SetRange(-frh->GetParameter(3),200);
frh->DrawClone("SAME");

// ECAL
TF1 *fre = new TF1("fre","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",9,200);
fre->SetParameters(0.95*fr->GetParameter(0), fr->GetParameter(1),
       fr->GetParameter(2), fr->GetParameter(3));
fre->SetParLimits(2,-0.18,-0.13);
fre->SetLineColor(kCyan+1);
cout << "*** Fitting fre ***" << endl;
gre->Fit(fre,"RN");
fre->DrawClone("SAME");
fre->SetLineStyle(kDashed);
fre->SetRange(-fre->GetParameter(3),200);
fre->DrawClone("SAME");

// E
TF1 *fer = new TF1("fer","([0]+[1]*pow(x+[3],[2])) * (x+[3])/x",3,200);
fer->SetParameters(0.98*fre->GetParameter(0), fre->GetParameter(1),
       fre->GetParameter(2), fre->GetParameter(3));
fre->SetParLimits(2,-0.18,-0.13);
fer->SetLineColor(kBlue);
cout << "*** Fitting fer ***" << endl;
ger->Fit(fer,"RN");
fer->DrawClone("SAME");
fer->SetLineStyle(kDashed);
fer->SetRange(-fer->GetParameter(3),200);
fer->DrawClone("SAME");
