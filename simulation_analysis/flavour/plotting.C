
{

TFile *f1 = new TFile("sim_dir/cteq6l1/herwig_Zjet_physics_1000000.root");
JetTree->Draw("fDR>>h1(150,0,1.5)");
TFile *f2 = new TFile("sim_dir/cteq6l1/pythia6_Zjet_physics_1000000.root");
JetTree->Draw("fDR>>h2(150,0,1.5)","","same");
TFile *f3 = new TFile("sim_dir/cteq6l1/pythia8_Zjet_physics_1000000.root");
JetTree->Draw("fDR>>h3(150,0,1.5)","","same");

gPad->SetLogy();
h1->SetLineColor(6);
h2->SetLineColor(12);

cout << "p8: " << h3->Integral(30,150)/h3->Integral() << endl;
cout << "p6: " << h2->Integral(30,150)/h2->Integral() << endl;
cout << "hpp: " << h1->Integral(30,150)/h1->Integral() << endl;
}