//PU_Plots.C

void PU_Plots(){
  TCanvas *c1 = new TCanvas("c1","",700,600);

   // auto pu1 = new TF1("f2","ROOT::Math::crystalball_function(-x, 1, 3, 1, 1.7)",-3,4);

   TLine* Line1 = new TLine(-1.5, 0, -1.5, 1.2);
   Line1->SetLineWidth(2);
   TLine* Line2 = new TLine(-0.5, 0, -0.5, 1.2);
   Line2->SetLineWidth(2);
   TLine* Line3 = new TLine(0.5, 0, 0.5, 1.2);
   Line3->SetLineWidth(2);
   TLine* Line4 = new TLine(1.5, 0, 1.5, 1.2);
   Line4->SetLineWidth(2);

   TF1 *pu1 = new TF1("pu1","TMath::Landau(x,[0],[1],1)",-2.6,3.7);
   pu1->SetParameters(-1,0.15);
   pu1->SetFillColor(kRed);
   pu1->SetLineColor(kRed);
   pu1->SetFillColorAlpha(kRed, 0.3);
   pu1->SetFillStyle(3001);
   pu1->Draw("c");


   TF1 *pu2 = new TF1("pu2","TMath::Landau(x,[0],[1],1)",-2.6,3.7);
   pu2->SetParameters(0,0.15);
   pu2->SetLineColor(kBlue);
   pu2->SetFillColor(kBlue);
   pu2->SetFillColorAlpha(kBlue, 0.3);
   pu2->SetFillStyle(3001);
   pu2->Draw("c same");

   TF1 *pu3 = new TF1("pu3","TMath::Landau(x,[0],[1],1)",-2.6,3.7);
   pu3->SetParameters(1,0.15);
   pu3->SetFillColor(kGreen);
   pu3->SetLineColor(kGreen);
   pu3->SetFillColorAlpha(kGreen, 0.3);
   pu3->SetFillStyle(3001);
   pu3->Draw("c same");

   Line1->Draw("SAME");
   Line2->Draw("SAME");
   Line3->Draw("SAME");
   Line4->Draw("SAME");

}
