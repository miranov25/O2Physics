// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Both tasks, ATask and BTask create two histograms. But whereas in
///        the first case (ATask) the histograms are not saved to file, this
///        happens automatically if OutputObj<TH1F> is used to create a
///        histogram. By default the histogram is saved to file
///        AnalysisResults.root. HistogramRegistry is yet an other possibility
///        to deal with histograms. See tutorial example histogramRegistery.cxx
///        for details.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "CommonUtils/TreeStream.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "TFile.h"
#include "TRandom.h"
 
using namespace o2;
using namespace o2::framework;


/// Tsalis/Hagedorn function describing charged pt spectra (m s = 62.4 GeV to 13 TeV) as in https://iopscience.iop.org/article/10.1088/2399-6528/aab00f/pdf
/// https://github.com/alisw/AliPhysics/blob/523f2dc8b45d913e9b7fda9b27e746819cbe5b09/PWGPP/AliAnalysisTaskFilteredTree.h#L145
/// \param pt     - transverse momentum
/// \param mass   - mass of particle
/// \param sqrts  -
/// \return       - invariant yields of the charged particle *pt
///    n(sqrts)= a + b/sqrt(s)                             - formula 6
///    T(sqrts)= c + d/sqrt(s)                             - formula 7
///    a = 6.81 ± 0.06       and b = 59.24 ± 3.53 GeV      - for charged particles page 3
///    c = 0.082 ± 0.002 GeV and d = 0.151 ± 0.048 (GeV)   - for charged particles page 4
Double_t TsalisCharged(Double_t pt, Double_t mass, Double_t sqrts){
  const Double_t a=6.81,   b=59.24;
  const Double_t c=0.082,  d=0.151;
  Double_t mt=TMath::Sqrt(mass*mass+pt*pt);
  Double_t n=a+b/sqrts;
  Double_t T=c+d/sqrts;
  Double_t p0 = n*T;
  Double_t result=TMath::Power((1.+mt/p0),-n);
  result*=pt;
  return result;
}

/// Random downsampling trigger function using Tsalis/Hagedorn spectra fit (sqrt(s) = 62.4 GeV to 13 TeV) as in https://iopscience.iop.org/article/10.1088/2399-6528/aab00f/pdf
/// \param pt
/// \param mass
/// \param sqrts
/// \param factorPt
/// \param factor1Pt
/// \return trigger bitmask
///         bit 1 - flat pt   trigger
///         bit 2 - flat q/pt trigger
///         bit 3 - MB trigger
Int_t  DownsampleTsalisCharged(Double_t pt, Double_t factorPt, Double_t factor1Pt, Double_t sqrts, Double_t mass, Double_t *weight){
  Double_t prob=TsalisCharged(pt,mass,sqrts);
  Double_t probNorm=TsalisCharged(1.,mass,sqrts);
  Int_t triggerMask=0;
  (*weight)=prob/probNorm;
  if (gRandom->Rndm()*prob/probNorm<factorPt) triggerMask|=1;
  if ((gRandom->Rndm()*((prob/probNorm)*pt*pt))<factor1Pt) triggerMask|=2;
  if (gRandom->Rndm()<factorPt) triggerMask|=4;
  return triggerMask;
}



struct OutputObjects {
  // histogram created with OutputObj<TH1F>
  o2::utils::TreeStreamRedirector *pcstream = nullptr;
  int counter=0;
  void init(o2::framework::InitContext& ic)
  {
    auto finishFunction = [this]() {
      if (pcstream) {
        pcstream->Close();
        delete pcstream;
        pcstream = nullptr;
      }
      LOGP(info, "OutputObject finishFunction");
    };
    ic.services().get<CallbackService>().set(CallbackService::Id::Stop, finishFunction);
  }

  void process(aod::Tracks const& tracks)
  {
    Float_t mass=0.139;
    Float_t sqrts=14400;
    float factor1Pt=0.01;
    float factorPt=0.01;
    //
    if (pcstream== nullptr) pcstream = new o2::utils::TreeStreamRedirector("tpcqcskimming.root", "recreate");
	  pcstream->GetFile()->cd();
    for (auto& track : tracks) {
      float phi=track.phi();
      float pt=track.pt();
      pcstream->GetFile()->cd();
      auto track2=track;
      //if (track.pt()<2) continue;
      double weight = 0;
      Int_t triggerMask = DownsampleTsalisCharged(pt, factorPt, factor1Pt, sqrts, mass, &weight);
      if (weight==0) continue;
      if (counter%100==0) LOGP(info, "11- Track {} has pT = {}", track.index(), track.pt());
      (*pcstream) << "tracks" <<
        "triggerMask="<<triggerMask<<
        "weight="<<weight<<
        "phi=" << phi <<
        "pt="<<pt<<
         "\n";
      counter++;
    }
    //pcstream->Close();
  }
  ~OutputObjects(){
    LOGP(info, "OutputObject destructor {}",(void*)pcstream);
    //if (pcstream) pcstream->Close();
    //delete pcstream;
  }
};


struct OutputTracks {
  // histogram created with OutputObj<TH1F>
  o2::utils::TreeStreamRedirector *pcstream = nullptr;
  int counter=0;
  void init(o2::framework::InitContext& ic)
  {
    auto finishFunction = [this]() {
      if (pcstream) {
        pcstream->Close();
        delete pcstream;
        pcstream = nullptr;
      }
      LOGP(info, "OutputObject finishFunction");
    };
    ic.services().get<CallbackService>().set(CallbackService::Id::Stop, finishFunction);
  }

  void process(soa::Join<aod::Tracks, aod::TracksExtra> const& tracks)
  {
    Float_t mass=0.139;
    Float_t sqrts=14400;
    float factor1Pt=0.01;
    float factorPt=0.01;
    //
    if (pcstream== nullptr) pcstream = new o2::utils::TreeStreamRedirector("tpcqcskimmingTracks.root", "recreate");
	  pcstream->GetFile()->cd();
    for (auto& track : tracks) {
      float phi=track.phi();
      float pt=track.pt();
      pcstream->GetFile()->cd();
      //if (track.pt()<2) continue;
      double weight = 0;
      Int_t triggerMask = DownsampleTsalisCharged(pt, factorPt, factor1Pt, sqrts, mass, &weight);
      if (triggerMask==0) continue;
      if (weight==0) continue;
      if (counter%100==0) LOGP(info, "11- Track {} has pT = {}", track.index(), track.pt());
      // Extra - table varaibles as in the https://github.com/AliceO2Group/AliceO2/blob/cab4b330a261046eba9aa628781754da4e849205/Framework/Core/include/Framework/AnalysisDataModel.h#L202
       uint32_t flags= track.flags();
       uint8_t itsClusterMap=track.itsClusterMap();
        uint8_t tpcNClsFindable=track.tpcNClsFindable();
        short  tpcNClsFindableMinusFound=track.tpcNClsFindableMinusFound();             //! TPC Clusters: Findable - Found
        short tpcNClsFindableMinusCrossedRows=track.tpcNClsFindableMinusCrossedRows(); //! TPC Clusters: Findable - crossed rows
        uint8_t tpcNClsShared =track.tpcNClsShared();                                    //! Number of shared TPC clusters
        uint8_t trdPattern= track.trdPattern();                                       //! Contributor to the track on TRD layer in bits 0-5, starting from the innermost
        float itsChi2NCl=track.itsChi2NCl();                                          //! Chi2 / cluster for the ITS track segment
        float tpcChi2NCl=track.tpcChi2NCl();                                          //! Chi2 / cluster for the TPC track segment
        float trdChi2=track.trdChi2();                                                //! Chi2 for the TRD track segment
        float tofChi2=track.tofChi2();                                                //! Chi2 for the TOF track segment
        float tpcSignal=track.tpcSignal();                                            //! dE/dx signal in the TPC
        float trdSignal=track.trdSignal();                                            //! dE/dx signal in the TRD
        float length=track.length();                                                  //! Track length
        float tofExpMom=track.tofExpMom();                                            //! TOF expected momentum obtained in tracking, used to compute the expected times
        float trackEtaEmcal=track.trackEtaEmcal();                                    //!
        float trackPhiEmcal=track.trackPhiEmcal();                                    //!
        float trackTime=track.trackTime();                                            //! Estimated time of the track in ns wrt collision().bc() or ambiguoustrack.bcSlice()[0]
        float trackTimeRes=track.trackTimeRes();                                      //! Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)
      /*
       DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);                                      //! Momentum at inner wall of the TPC
DECLARE_SOA_COLUMN(Flags, flags, uint32_t);                                                   //! Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, uint8_t);                                    //! ITS cluster map, one bit per a layer, starting from the innermost
DECLARE_SOA_COLUMN(TPCNClsFindable, tpcNClsFindable, uint8_t);                                //! Findable TPC clusters for this track geometry
DECLARE_SOA_COLUMN(TPCNClsFindableMinusFound, tpcNClsFindableMinusFound, int8_t);             //! TPC Clusters: Findable - Found
DECLARE_SOA_COLUMN(TPCNClsFindableMinusCrossedRows, tpcNClsFindableMinusCrossedRows, int8_t); //! TPC Clusters: Findable - crossed rows
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, uint8_t);                                    //! Number of shared TPC clusters
DECLARE_SOA_COLUMN(TRDPattern, trdPattern, uint8_t);                                          //! Contributor to the track on TRD layer in bits 0-5, starting from the innermost
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);                                            //! Chi2 / cluster for the ITS track segment
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);                                            //! Chi2 / cluster for the TPC track segment
DECLARE_SOA_COLUMN(TRDChi2, trdChi2, float);                                                  //! Chi2 for the TRD track segment
DECLARE_SOA_COLUMN(TOFChi2, tofChi2, float);                                                  //! Chi2 for the TOF track segment
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);                                              //! dE/dx signal in the TPC
DECLARE_SOA_COLUMN(TRDSignal, trdSignal, float);                                              //! dE/dx signal in the TRD
DECLARE_SOA_COLUMN(Length, length, float);                                                    //! Track length
DECLARE_SOA_COLUMN(TOFExpMom, tofExpMom, float);                                              //! TOF expected momentum obtained in tracking, used to compute the expected times
DECLARE_SOA_COLUMN(TrackEtaEMCAL, trackEtaEmcal, float);                                      //!
DECLARE_SOA_COLUMN(TrackPhiEMCAL, trackPhiEmcal, float);                                      //!
DECLARE_SOA_COLUMN(TrackTime, trackTime, float);                                              //! Estimated time of the track in ns wrt collision().bc() or ambiguoustrack.bcSlice()[0]
DECLARE_SOA_COLUMN(TrackTimeRes, trackTimeRes, float);                                        //! Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)
       */
      (*pcstream) << "tracks" <<
        "triggerMask="<<triggerMask<<
        "weight="<<weight<<
        "phi=" << phi <<
        "pt="<<pt<<
         // extra
         "flags="<<flags<<
         "itsClusterMap="<<itsClusterMap<<
          "tpcNClsFindable="<<tpcNClsFindable<<
        "tpcNClsFindableMinusFound="<<tpcNClsFindableMinusFound<<
            "tpcNClsFindableMinusCrossedRows="<<tpcNClsFindableMinusCrossedRows<<
            "tpcNClsShared=" <<tpcNClsShared<<
          "trdPattern="<<trdPattern<<
        "itsChi2NCl="<<itsChi2NCl<<
        "tpcChi2NCl="<<tpcChi2NCl<<
        "trdChi2="<<trdChi2<<
        "tofChi2="<<tofChi2<<
        "tpcSignal="<<tpcSignal<<
        "trdSignal="<<trdSignal<<
        "length="<<length<<
        "tofExpMom="<<tofExpMom<<
        "trackEtaEmcal="<<trackEtaEmcal<<
        "trackPhiEmcal="<<trackPhiEmcal<<
        "trackTime="<<trackTime<<
        "trackTimeRes="<<trackTimeRes<<
         "\n";
      counter++;
    }
    //pcstream->Close();
  }
  ~OutputTracks(){
    LOGP(info, "OutputTracks destructor {}",(void*)pcstream);
    //if (pcstream) pcstream->Close();
    //delete pcstream;
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    // adaptAnalysisTask<OutputObjects>(cfgc),
       adaptAnalysisTask<OutputTracks>(cfgc),
  };
}
