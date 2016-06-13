
{

TFile *f1 = new TFile("sim_dir/cteq6l1/pythia8_dijet_physics_1000000.root");
JetTree->Draw("fPartonPT>>h1(100,0,2000)");
TFile *f2 = new TFile("sim_dir/cteq6l1/pythia6_dijet_physics_1000000.root");
JetTree->Draw("fPartonPT>>h2(100,0,2000)","","same");
TFile *f3 = new TFile("sim_dir/cteq6l1/herwig_dijet_physics_1000000.root");
JetTree->Draw("fPartonPT>>h3(100,0,2000)","","same");

gPad->SetLogy();
h1->SetLineColor(6);
h2->SetLineColor(12);

cout << "p8: " << h1->Integral(2,31)/h1->Integral() << endl;
cout << "p6: " << h2->Integral(2,31)/h2->Integral() << endl;
cout << "hpp: " << h3->Integral(2,31)/h3->Integral() << endl;
}
