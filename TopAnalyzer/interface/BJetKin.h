#ifndef JetKin_h
#define JetKin_h

using namespace std; 

class BJetKin{

 public:

  explicit BJetKin();
  ~BJetKin();

  void book(ofstream&, char*);
  void fill(const pat::Jet&);
  
 private:

  TH1F *e_;
  TH1F *et_;
  TH1F *p_;
  TH1F *pt_;
  TH1F *eta_;
  TH1F *theta_;
  TH1F *rapidity_;
  TH1F *phi_;
  TH1F *mass_;
  TH1F *mt_;
  TH1F *n60_;
  TH1F *n90_;
  TH1F *nConst_;
};

#endif
