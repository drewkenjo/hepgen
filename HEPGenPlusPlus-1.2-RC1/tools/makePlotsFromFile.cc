#include <iostream>
#include <fstream>
#include <string>

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TText.h>
#include <TAxis.h>

using namespace std;

int main(int argc, char** argv) {
if (argc < 2){
  printf("Usage: %s [dataTableFromXSecPlotter]\n",argv[0]);
  return -1;
}
ifstream inFile;
inFile.open(argv[1]);
if (!inFile.good()){
  printf("Could not open file: %s \n",argv[1]);
  return -1;
}

double qsq,phi;
string dummy;
inFile >> dummy;
inFile >> qsq;
inFile >> dummy;
inFile >> phi;

double t[300],moutarde[300],mosse[300];
for (int i = 0; i < 250;i++){
  inFile >> t[i] >> mosse[i] >> moutarde[i];
  t[i] *= -1;
  if (moutarde[i] == -1)
    moutarde[i] = 0;
  printf("%f %f %f \n",t[i],mosse[i],moutarde[i]);
}
TCanvas* myCanvas = new TCanvas("myCanvas","Plots",1280,720);
myCanvas->SetGrid();
TGraph* mosseGraph = new TGraph(250,t,mosse);
mosseGraph->SetLineColor(1);
mosseGraph->SetLineWidth(3);

mosseGraph->SetTitle("L. Mosse");
TGraph* moutardeGraph = new TGraph(250,t,moutarde);
moutardeGraph->SetLineColor(2);
moutardeGraph->SetLineWidth(3);
moutardeGraph->SetTitle("H. Moutarde");
TMultiGraph* myMG = new TMultiGraph;

myMG->Add(mosseGraph);
myMG->Add(moutardeGraph);
myMG->Draw("APL");
gPad->SetLogy(true);
myMG->SetMinimum(0.00001);
myMG->SetMaximum(10000000);
myMG->GetYaxis()->SetTitle("#sigma [nbarn]");
myMG->GetXaxis()->SetTitle("-t [GeV^2]");
char title[200];
sprintf(title,"q2: %f phi: %f",qsq,phi);
TText* myText = new TText(1.5,30000000,title);
myText->Draw();

myCanvas->SetLogy(true);

// myCanvas->BuildLegend();

myCanvas->Update();
myCanvas->Draw();

char fileName[200];
sprintf(fileName,"%s.png",argv[1]);

myCanvas->SaveAs(fileName);


return 0;



}
