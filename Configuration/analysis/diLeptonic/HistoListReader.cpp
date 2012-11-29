#include "HistoListReader.h"
#include <TH1.h>
#include <iostream>

HistoListReader::HistoListReader(const char* filename) : 
    filename_ ( filename ),
    isZombie_ ( false )
{
    std::ifstream controlHistStream(filename, std::ifstream::in);
    if (!controlHistStream.good()) {
        isZombie_ = true;
        std::cerr << "HistoListReader: cannot read " << filename << std::endl;
        return;
    }
    plots_.clear();
    
    while(!controlHistStream.eof()){
        // Read HistoList-File
        PlotProperties m;
//        # Name, Extra, axis labels (y,x), rebin, do_dyscale, logx, logy, ymin, ymax, xmin, xmax, nbins, xbins, bcs
        controlHistStream 
            >> m.name
            >> m.specialComment
            >> m.ytitle
            >> m.xtitle
            >> m.rebin
            >> m.do_dyscale
            >> m.logX
            >> m.logY
            >> m.ymin
            >> m.ymax
            >> m.xmin
            >> m.xmax
            >> m.bins;

        // Avoid running over empty lines in 'HistoList'-File
        if (m.name == "") continue;

        m.xbinbounds.clear();
        m.bincenters.clear();

        for(int i = 0; i <= m.bins; ++i){
            double temp;
            controlHistStream>>temp;
            m.xbinbounds.push_back(temp);
        }
        for(int i = 0; i < m.bins; i++){//only until bincenter code is finalized
            double temp;
            controlHistStream>>temp;
            m.bincenters.push_back(temp);
        }
        plots_[m.name] = m;
    }
}

bool HistoListReader::IsZombie() const
{
    return isZombie_;
}

/////////////////////////////////////////

PlotProperties& HistoListReader::getPlotProperties(TString name)
{
    std::map <TString, PlotProperties >::iterator it = plots_.find(name);
    if (it == plots_.end()) { std::cerr << "no such plot in " << filename_ << ": ``" << name << "''" <<std::endl; exit(671); }
    return it->second;
}

PlotProperties::PlotProperties() :
    histo_ (0)
{
}

TH1* PlotProperties::getHistogram()
{
    if (!histo_) MakeHisto();
    return histo_;
}

void PlotProperties::MakeHisto()
{
    bool old = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
    histo_ = new TH1D(name, name, bins, &xbinbounds[0]);
    TH1::AddDirectory(old);
    histo_->GetXaxis()->SetTitle(xtitle);
    histo_->GetYaxis()->SetTitle(ytitle);
}

PlotProperties::~PlotProperties()
{
    delete histo_;
}
