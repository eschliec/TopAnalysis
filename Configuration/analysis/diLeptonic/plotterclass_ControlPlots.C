#include "plotterclass_ControlPlots.h"



Plotter::Plotter(TString name_, TString XAxis_,TString YAxis_, double rangemin_, double rangemax_)
{
  name=name_;
  rangemin=rangemin_;
  rangemax=rangemax_;
  XAxis=XAxis_;
  YAxis=YAxis_;
  initialized=false;
}

Plotter::Plotter()
{
  name="defaultName";
  specialComment="Standard";
  rangemin=0;
  rangemax=3;
  YAxis="N_{events}";
  initialized=false;
  datafiles = 0;
  
  // DAVID  
  outpath = "";
  outpathPlots = "Plots";
  subfolderChannel = "";
  subfolderSpecial = "";

}


Plotter::~Plotter()
{
}


//IVAN
void Plotter::SetDataLumi(double Lumi){
   lumi=Lumi;
}


// DAVID
void Plotter::SetOutpath(TString path)
{
  outpath = path;
}



void Plotter::DYScaleFactor(){
    doDYScale = false;
    if (!doDYScale){
        //Don't do the DY SF calculation: Set all DY SFs to 1 and exit
        DYScale[0]=1.;
        DYScale[1]=1.;
        DYScale[2]=1.;
        DYScale[3]=1.;
        return;
    }
    
    TH1::AddDirectory(kFALSE);
    ifstream FileList("FileLists/HistoFileList_Nominal_combined.txt");
    TString filename;
    
    double NoutEEDYMC=0, NinEEDYMC=0, NoutMuMuDYMC=0, NinMuMuDYMC=0;//Number of events in/out of z-veto region for the DY MC
    double NinEE=0, NinMuMu=0, NinEMu=0;//Number of events in z-veto region for data
    double NinEEloose=0, NinMuMuloose=0;//Number of data events in Z-Veto region with MET cut
    double NinEEMC=0, NinMuMuMC=0;//All other MC events

    while(!FileList.eof()){
        FileList>>filename;
        if(filename!=""){
            TFile *ftemp = TFile::Open(filename);
            if(filename.Contains("ee")){
                if(filename.Contains("run")){
                    TH1D *htemp = (TH1D*)ftemp->Get("Zh1");      NinEE+=htemp->Integral();
                    TH1D *htemp1 = (TH1D*)ftemp->Get("Looseh1"); NinEEloose+=htemp1->Integral();	  
                }
                else if(filename.Contains("dy")){
                    if(filename.Contains("50inf")){
                        TH1D *htemp = (TH1D*)ftemp->Get("Zh1");   NinEEDYMC+=htemp->Integral();
                        TH1D *htemp1 = (TH1D*)ftemp->Get("TTh1"); NoutEEDYMC+=htemp1->Integral();
                    }
                    else{TH1D *htemp = (TH1D*)ftemp->Get("TTh1"); NoutEEDYMC+=htemp->Integral();}
                }
                else{TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinEEMC+=htemp->Integral();}
            }

            if(filename.Contains("emu") && filename.Contains("run")){TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinEMu+=htemp->Integral();}

            if(filename.Contains("mumu")){
                if(filename.Contains("run")){
                    TH1D *htemp = (TH1D*)ftemp->Get("Zh1");      NinMuMu+=htemp->Integral();
                    TH1D *htemp1 = (TH1D*)ftemp->Get("Looseh1"); NinMuMuloose+=htemp1->Integral();
                }
                else if(filename.Contains("dy")){
                    if(filename.Contains("50inf")){
                        TH1D *htemp = (TH1D*)ftemp->Get("Zh1");     NinMuMuDYMC+=htemp->Integral();
                        TH1D *htemp1 = (TH1D*)ftemp->Get("TTh1");   NoutMuMuDYMC+=htemp1->Integral();
                    }
                    else{TH1D *htemp = (TH1D*)ftemp->Get("TTh1"); NoutMuMuDYMC+=htemp->Integral();}
                }
                else{TH1D *htemp = (TH1D*)ftemp->Get("Zh1"); NinMuMuMC+=htemp->Integral();}
            }
            delete ftemp;
        }
    }
    double NoutMCEE = (NoutEEDYMC/NinEEDYMC)*(NinEE - 0.5*NinEMu*sqrt(NinEEloose/NinMuMuloose));
    double NoutMCMuMu = (NoutMuMuDYMC/NinMuMuDYMC)*(NinMuMu - 0.5*NinEMu*sqrt(NinMuMuloose/NinEEloose));
    
    double DYSFEE = NoutMCEE/NoutEEDYMC;
    double DYSFMuMu = NoutMCMuMu/NoutMuMuDYMC;

    cout << endl;
    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    cout << "Calculation of DY Scale Factors for '" << name << "'  " << endl;
    cout<<"DYSFEE:                 "<<DYSFEE<<endl;
    cout<<"DYSFMuMu:               "<<DYSFMuMu<<endl;
    cout<<"NinEEloose:             "<<NinEEloose<<endl;
    cout<<"NinMMloose:             "<<NinMuMuloose<<endl;
    cout<<"kee:                    "<<sqrt(NinEEloose/NinMuMuloose)<<" +- "<<0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<endl;
    cout<<"kmumu:                  "<<sqrt(NinMuMuloose/NinEEloose)<<" +- "<<0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<endl;
    cout<<"Rout/Rin Mumu:          "<<(NoutMuMuDYMC/NinMuMuDYMC)<<endl;
    cout<<"Rout/Rin ee:            "<<(NoutEEDYMC/NinEEDYMC)<<endl;
    cout<<"Est. From Data(ee):     "<<NoutMCEE<<endl;
    cout<<"Est. From Data(mumu):   "<<NoutMCMuMu<<endl;
    cout<<"Est. From MC(ee):       "<<NoutEEDYMC<<endl;
    cout<<"Est. From MC(mumu):     "<<NoutMuMuDYMC<<endl;
    cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
    cout << endl;

    DYScale[0]=DYSFEE;
    DYScale[1]=DYSFMuMu;
    DYScale[2]=1.;
    DYScale[3]=(DYSFEE+DYSFMuMu)/2;

}


void Plotter::setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_, int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_, double rangemax_, int bins_, std::vector<double>XAxisbins_, std::vector<double>XAxisbinCenters_)
{
    XAxisbins.clear();
    XAxisbinCenters.clear();
    XAxisbinCenters = XAxisbinCenters_;
    XAxisbins = XAxisbins_;
    rebin=rebin_;
    name=name_;
    specialComment=specialComment_;
    bins=bins_;
    doDYScale = doDYScale_;
    logX = logX_;
    logY = logY_;
    ymin = ymin_;
    ymax = ymax_;
    rangemin=rangemin_;
    rangemax=rangemax_;
    YAxis=YAxis_;
    XAxis=XAxis_;
    if(XAxis.Contains("band#bar{b}")){//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("band#bar{b}",11,"b and #bar{b}",13);
    }
    if(XAxis.Contains("tand#bar{t}")){//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("tand#bar{t}",11,"t and #bar{t}",13);
    }
    if(XAxis.Contains("l^{+}andl^{-}")){//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("l^{+}andl^{-}",13,"l^{+} and l^{-}",15);
    }
    if(YAxis.Contains("toppairs")){
        YAxis.ReplaceAll("toppairs",8,"top-quark pairs",15);
    }
    if(YAxis.Contains("topquarks")){
        YAxis.ReplaceAll("topquarks",9, "top quarks",10);
    }

    DYScale.insert(DYScale.begin(), 4, 1.);
}


void Plotter::setDataSet(std::vector<TString> dataset_, std::vector<double> scales_, std::vector<TString> legends_, std::vector<int> colors_, TString DYEntry_)
{
  dataset.clear();
  scales.clear();
  legends.clear();
  colors.clear();
  dataset=dataset_;
  scales=scales_;
  legends=legends_;
  colors=colors_;
  DYEntry=DYEntry_;
  
}


void Plotter::setDataSet(TString mode)
{
    initialized=false;

    channelLabel.insert(channelLabel.begin(), 4, "");
    channel=mode;
    if(channel =="ee")      {channelType=0; channelLabel[0]="ee";}
    if(channel =="mumu")    {channelType=1; channelLabel[1]="#mu#mu";}
    if(channel =="emu")     {channelType=2; channelLabel[2]="e#mu";}
    if(channel =="combined"){channelType=3; channelLabel[3]="Dilepton Combined";}

    // Set dataset specific subfolders
    outpathPlots = "./Plots";
    subfolderChannel = channel; 
    subfolderChannel.Prepend("/");
    subfolderSpecial = "";
    if ( specialComment.CompareTo("Standard") != 0 ) {subfolderSpecial = specialComment.Prepend("/");}

    DYEntry = "Z / #gamma* #rightarrow ee/#mu#mu";

    if(channel!="ee" && channel!="emu" && channel!="mumu" && channel!="combined"){
        cout<<"ERROR: wrong channel. You want"<<channel<<" but 'ee', 'emu', 'mumu' or 'combined' only accepted"<<endl;
        return;
    }
    ifstream FileList("FileLists/HistoFileList_Nominal_"+channel+".txt");
    TString filename;

    datafiles=0;
    dataset.clear();
    legends.clear();
    colors.clear();

    while(!FileList.eof()){
        FileList>>filename;
        if(filename!=""){
            dataset.push_back(filename);
            if(filename.Contains("run")){             legends.push_back("Data"); colors.push_back(kBlack);datafiles++;}
            else if(filename.Contains("ttbarsignal")){legends.push_back("t#bar{t} Signal"); colors.push_back(kRed+1);}
            else if(filename.Contains("ttbarbg")){    legends.push_back("t#bar{t} Other"); colors.push_back(kRed-7);}
            else if(filename.Contains("single")){     legends.push_back("Single Top"); colors.push_back(kMagenta);}
            else if(filename.Contains("wtolnu")){     legends.push_back("W+Jets"); colors.push_back(kGreen-3);}
            else if(filename.Contains("qcd")){        legends.push_back("QCD Multijet"); colors.push_back(kYellow);}
            else if(filename.Contains("dytautau")){   legends.push_back("Z / #gamma* #rightarrow #tau#tau"); colors.push_back(kAzure+8);}
            else if(filename.Contains("dymumu") || filename.Contains("dyee")){legends.push_back("DYEntry"); colors.push_back(kAzure-2);}
            else if(filename.Contains("ww") || filename.Contains("wz") || filename.Contains("zz")){legends.push_back("Diboson"); colors.push_back(10);}
        }
    }
}


void Plotter::fillHisto()
{
    if (initialized) { return; }
    TH1::AddDirectory(kFALSE);
    
    hists.clear();
    for(unsigned int i=0; i<dataset.size(); i++){
        TFile *ftemp = TFile::Open(dataset[i]);
        TH1D *hist = (TH1D*)ftemp->Get(name)->Clone();
        if(name.Contains("Lepton") || name.Contains("BJet") || name.Contains("Top")){
            TString stemp = name;

            if(name.Contains("Lepton"))    {stemp.ReplaceAll("Lepton",6,"AntiLepton",10);}
            else if(name.Contains("BTag")) {stemp.ReplaceAll("BJet",4,"AntiBJet",8);}
            else if(name.Contains("Top"))  {stemp.ReplaceAll("Top",3,"AntiTop",7);}

            TH1D *hist2 = (TH1D*)ftemp->Get(stemp)->Clone();
            hist->Add(hist2);
        }

        //Rescaling to the data luminosity
        double LumiWeight = CalcLumiWeight(dataset[i]);
        cout << "LumiWeight for " << dataset[i] << " = " << LumiWeight << endl;
        ApplyFlatWeights(hist, LumiWeight);

//         //Apply any other flat weight (still to do: define the weights somewhere else)
//         ApplyFlatWeights(hist, triggSF);
//         ApplyFlatWeights(hist, lepSF);
//         ApplyFlatWeights(hist, KinRecSF);

        setHHStyle(*gStyle);

        hists.push_back(*hist);
        delete ftemp;
    }
    initialized=true;
    
}



void Plotter::write() // do scaling, stacking, legending, and write in file 
{
  if(initialized){

  TCanvas *c = new TCanvas("","");

  THStack *stack=  new THStack("def", "def");
  TLegend *leg =  new TLegend(0.70,0.55,0.98,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
  leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() - 0.05 - leg->GetNRows()*0.04);
  leg->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
  leg->SetY2NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength());
  TH1D *drawhists[hists.size()];

  int legchange=0;
  leg->Clear();
  c->Clear();
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  c->SetName("");
  c->SetTitle("");
  for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
    drawhists[i]=(TH1D*) hists[i].Clone();
    if(rebin>1) drawhists[i]->Rebin(rebin);
    setStyle(*drawhists[i], i);
    if(legends[i] != "Data"){
        if(legends[i] == legends[0]){drawhists[0]->Add(drawhists[i]);}
        if(legends[i] == "t#bar{t} Signal"){signalHist = i;}
        if((legends[i] == DYEntry) && channelType!=2 ){drawhists[i]->Scale(DYScale[channelType]);}
        if(i > 1){
            if(legends[i] != legends[i-1]){
                legchange = i; 
                if((legends[i] == DYEntry)&& DYScale[channelType] != 1) leg->AddEntry(drawhists[i], legends[i],"f");
                else leg->AddEntry(drawhists[i], legends[i] ,"f");
            }
            else{drawhists[legchange]->Add(drawhists[i]);}
        }

        if(i!=(hists.size()-1)){
            if(legends[i]!=legends[i+1]){drawhists[i]->SetLineColor(1);}
        }
        else{drawhists[i]->SetLineColor(1);}

        if(legends.at(i) != legends.at(i-1)){
            drawhists[i]->SetLineColor(1);
            stack->Add(drawhists[i]); 
        }
    }
    else{
        if(i==0) leg->AddEntry(drawhists[i], legends[i] ,"pe");
        if(i>0){
            if(legends[i] != legends[i-1]){leg->AddEntry(drawhists[i], legends[i] ,"pe");}
            if(legends[i] == legends[0])  {drawhists[0]->Add(drawhists[i]);}
        }
    }
  } 
  
  leg = ControlLegend(hists.size(), drawhists, legends, leg);

    //What is this for??
    TFile *f0 = new TFile("SigBackground.root","UPDATE");

    TList* l = stack->GetHists(); 
    TH1D* stacksum = (TH1D*) l->At(0)->Clone();

    TString aaa = "a";
    for (int i = 1; i < l->GetEntries(); ++i) {
    aaa=aaa+"a"; 
    stacksum->Add((TH1D*)l->At(i));
    if(legends[datafiles+i] == "t#bar{t} Other") {
        stacksum->Write(name+"_"+channel+aaa+"_Background");
    }
    if(legends[datafiles+i] == "t#bar{t} Signal") {
        stacksum->Write(name+"_"+channel+aaa+"_Signal");
    }

    } 
    f0->Close();

  //stat uncertainty::make a function 
  TH1* syshist = (TH1*)stacksum->Clone();
  double lumierr = 0.045;
  //Kidonakis
  double topxsec = SampleXSection(TString("ttbar"));
  double topxsecErr2 = 2.2*2.2 + 4.4*4.4 + 5.5*5.5;

  for(Int_t i=0; i<=syshist->GetNbinsX(); ++i){
    
    Double_t binc = stacksum->GetBinContent(i);
    // calculate uncertainty: lumi uncertainty
    Double_t binerr2 = binc*binc*lumierr*lumierr;
    Double_t topunc = 0; // uncertainty on top xsec
    
    double topRelUnc =  TMath::Sqrt(topxsecErr2)/topxsec;
    topunc += drawhists[signalHist]->GetBinContent(i)*topRelUnc;
    binerr2 += (topunc*topunc);
    syshist->SetLineColor(1);
    syshist->SetBinError(i, TMath::Sqrt(binerr2));
  }

  if(logY)c->SetLogy();
  syshist->SetFillStyle(3004);
  syshist->SetFillColor(kBlack);
  syshist->SetMarkerStyle(0);
  leg->AddEntry( syshist, "Uncertainty", "f" );
  
  // Set histogram ranges in X- and Y-axis
  if(rangemin!=0 || rangemax!=0) drawhists[0]->SetAxisRange(rangemin, rangemax, "X");
  drawhists[0]->SetMinimum(ymin);
  if(ymax==0){
    if(logY){ymax = 18*drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin());}
    else{    ymax = 1.5*drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin());}
  }
  drawhists[0]->SetMaximum(ymax);
  drawhists[0]->GetXaxis()->SetNoExponent(kTRUE);
  
  TGaxis::SetMaxDigits(2);

  //Removal of extra ticks in JetMult plots
  if(name.Contains("jet") && name.Contains("Multi")){ drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(),0,0, 1);}
  
  drawhists[0]->Draw("e1"); //############## 
  
  stack->Draw("same HIST");
  gPad->RedrawAxis();
  TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
  setex1->Draw();
  
  syshist->Draw("same,E2");
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
  setex2->Draw();
  drawhists[0]->Draw("same,e1"); //#############
  
  DrawCMSLabels(false, lumi);
  
  DrawDecayChLabel(channelLabel[channelType]);    
  leg->Draw("SAME");  
  //drawRatio(drawhists[0], stacksum, 0.5, 1.9, *gStyle);

  // Create Directory for Output Plots 
  gSystem->MakeDirectory(outpathPlots);
  gSystem->MakeDirectory(outpathPlots+"/"+subfolderChannel+"/"+subfolderSpecial);  
  c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/"+name+".eps");  
  c->Print(outpathPlots+subfolderChannel+subfolderSpecial+"/"+name+".C");  
  c->Clear();
  leg->Clear();
  delete c;
  delete leg;
  
  }
  else std::cout << "Histogram " << name << " not filled during the process." << std::endl;
}

void Plotter::setStyle(TH1D &hist, unsigned int i)
{
  hist.SetFillColor(colors[i]);
  hist.SetLineColor(colors[i]);
  
  if(legends[i] == "Data"){
    hist.SetMarkerStyle(20); 
    hist.SetMarkerSize(1.);
    hist.SetLineWidth(1);
    hist.GetXaxis()->SetLabelFont(42);
    hist.GetYaxis()->SetLabelFont(42);
    hist.GetXaxis()->SetTitleSize(0.04);
    hist.GetYaxis()->SetTitleSize(0.04);
    hist.GetXaxis()->SetTitleFont(42);
    hist.GetYaxis()->SetTitleFont(42);
    hist.GetYaxis()->SetTitle(YAxis);
    hist.GetYaxis()->SetTitleOffset(1.7);
    hist.GetXaxis()->SetTitleOffset(1.25);
    
    //set XAxis label
    TString h_name= TString(hist.GetName());
    if(h_name.Contains("pT") || h_name.Contains("Mass")){hist.GetXaxis()->SetTitle(XAxis+" #left[GeV#right]");}
    else {hist.GetXaxis()->SetTitle(XAxis);};
    
    //set YAxis label with binwidth
    double binwidth = hist.GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width<<binwidth;
    TString ytitle = TString(hist.GetYaxis()->GetTitle()).Copy();
    ytitle.Append(" / ").Append(width.str());
    if(h_name.Contains("pT") || h_name.Contains("Mass") || h_name.Contains("mass") || h_name.Contains("MET") || h_name.Contains("HT")){
        ytitle.Append(" / ").Append(width.str()).Append(" GeV");
    };
    hist.GetYaxis()->SetTitle(ytitle);
  }
}


TLegend* Plotter::ControlLegend(int HistsSize, TH1* drawhists[], std::vector<TString> Legends, TLegend *leg){
    //hardcoded ControlPlot legend
    std::vector<TString> OrderedLegends;    
    OrderedLegends.push_back("Data");
    OrderedLegends.push_back("t#bar{t} Signal");
    OrderedLegends.push_back("t#bar{t} Other");
    OrderedLegends.push_back("Single Top");
    OrderedLegends.push_back("W+Jets");
    OrderedLegends.push_back("Z / #gamma* #rightarrow ee/#mu#mu");
    OrderedLegends.push_back("Z / #gamma* #rightarrow #tau#tau");
    OrderedLegends.push_back("Diboson");
    OrderedLegends.push_back("QCD Multijet");
      
    leg->Clear();
    leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
    leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()  -gStyle->GetTickLength()-0.30);
    leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()  -gStyle->GetTickLength());
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);
    for(int i=0; i<(int)OrderedLegends.size(); ++i){
        for(int j=0; j<HistsSize; ++j){
            if (OrderedLegends[i] == Legends[j]){
                if( OrderedLegends[i] == "Data"){
                    leg->AddEntry(drawhists[j], OrderedLegends[i], "pe");
                    break;
                }
                else{
                    leg->AddEntry(drawhists[j], OrderedLegends[i], "f");
                    break;
                }
            }
        }
    }
    return leg;
}

TLegend* Plotter::ControlLegend(int HistsSize, TH1D* drawhists[], std::vector<TString> Legends, TLegend *leg){
    
    //hardcoded ControlPlot legend
    std::vector<TString> OrderedLegends;    
    OrderedLegends.push_back("Data");
    OrderedLegends.push_back("t#bar{t} Signal");
    OrderedLegends.push_back("t#bar{t} Other");
    OrderedLegends.push_back("Single Top");
    OrderedLegends.push_back("W+Jets");
    OrderedLegends.push_back(DYEntry);
    OrderedLegends.push_back("Z / #gamma* #rightarrow #tau#tau");
    OrderedLegends.push_back("Diboson");
    OrderedLegends.push_back("QCD Multijet");
    
    leg->Clear();
    leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
    leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.30);
    leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);
    for(int i=0; i<(int)OrderedLegends.size(); ++i){
        for(int j=0; j<HistsSize; ++j){
            if (OrderedLegends[i] == Legends[j]){
                if( OrderedLegends[i] == "Data"){
                    leg->AddEntry(drawhists[j], OrderedLegends[i], "pe");
                    break;
                }
                else{
                    leg->AddEntry(drawhists[j], OrderedLegends[i], "f");
                    break;
                }
            }
        }
    }
    return leg;
}


void Plotter::ApplyFlatWeights(TH1* varhists, const double weight){

    if(weight == 0) {cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<endl;}
    if(weight >=1e3){cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<endl;}
    varhists->Scale(weight);
}


void Plotter::ApplyFlatWeights(TH1* varhists[], const double weight){

    if(weight == 0) {cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<endl;}
    if(weight >=1e3){cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<endl;}
    for(size_t i=0; i<hists.size(); i++){
        varhists[i]->Scale(weight);
    }
}

double Plotter::CalcLumiWeight(TString WhichSample){
    if (WhichSample.Contains("run")) return 1;
    
    double lumiWeight=0;
    if(WhichSample!=""){
        double XSection = SampleXSection(WhichSample);
        if(XSection <= 0.){
            cout<<"Sample XSection is <0. Can't calculate luminosity weight!! returning"<<endl;
            return 0;
        };
        
        //From 'filename' get the number of weighted (MC weights) event processed.
        TFile *f = TFile::Open(WhichSample);
        if(!f || f->IsZombie()){
            cout<<"The file you requested is not valid"<<endl;
            return 0;
        };
        TH1 *h_NrOfEvts = dynamic_cast<TH1*>(f->Get("weightedEvents"));
        if(!h_NrOfEvts){
            cout<<"The histogram you requested is not valid"<<endl;
            return 0;
        };
        double NrOfEvts = h_NrOfEvts->GetBinContent(1);
        //double XSection = SampleXSection(WhichSample);

        lumiWeight = lumi*XSection/NrOfEvts;

        f->Close();
        h_NrOfEvts->Delete();
    }
    
    if (lumiWeight == 0) {
        cout << WhichSample << " has lumi weight 0\n";
    }

    return lumiWeight;
}


double Plotter::SampleXSection(TString filename){
    
    //MC cross sections taken from:
    //  https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    //  AN-12/194    AN-12/228
    
    double XSec=-1.;
    if(filename.Contains("run"))              {XSec = 1;}
    else if(filename.Contains("ttbar"))       {XSec = 225.197;}
    else if(filename.Contains("single"))      {XSec = 22.2;}
    else if(filename.Contains("ww"))          {XSec = 54.838;}
    else if(filename.Contains("wz"))          {XSec = 33.21;}
    else if(filename.Contains("zz"))          {XSec = 17.654;}
    else if(filename.Contains("1050"))        {XSec = 860.5;}
    else if(filename.Contains("50inf"))       {XSec = 3532.8;}
    else if(filename.Contains("wtolnu"))      {XSec = 36257.2;}
    else if(filename.Contains("qcdmu15"))     {XSec = 3.640E8*3.7E-4;}
    else if(filename.Contains("qcdmu2030"))   {XSec = 2.870E8*6.500E-3;}
    else if(filename.Contains("qcdmu3050"))   {XSec = 6.609E7*12.20E-3;}
    else if(filename.Contains("qcdmu5080"))   {XSec = 8.802E6*21.80E-3;}
    else if(filename.Contains("qcdmu80120"))  {XSec = 1.024E6*39.50E-3;}
    else if(filename.Contains("qcdmu120170")) {XSec = 1.578E5*47.30E-3;}
    else if(filename.Contains("qcdem2030"))   {XSec = 2.886E8*10.10E-3;}
    else if(filename.Contains("qcdem3080"))   {XSec = 7.433E7*62.10E-3;}
    else if(filename.Contains("qcdem80170"))  {XSec = 1.191E6*153.9E-3;}
    else if(filename.Contains("qcdbcem2030")) {XSec = 2.886E8*5.800E-4;}
    else if(filename.Contains("qcdbcem3080")) {XSec = 7.424E7*2.250E-3;}
    else if(filename.Contains("qcdbcem80170")){XSec = 1.191E6*10.90E-3;}
    
    return XSec;
}



void Plotter::DrawDecayChLabel(TString decaychannel, double textSize){
    // Draw label for Decay Channel in upper left corner of plot

    TPaveText *decch = new TPaveText();

    decch -> AddText(decaychannel);

    decch -> SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch -> SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch -> SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch -> SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch -> SetFillStyle(0);
    decch -> SetBorderSize(0);
    if(textSize!=0) decch->SetTextSize(textSize);
    decch -> SetTextAlign(12);
    decch -> Draw("same");
}
