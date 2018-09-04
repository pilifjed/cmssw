// system include files
#include <memory>

// user include files
#include <map>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CTPPSDigi/interface/TotemTimingDigi.h"
#include "DataFormats/CTPPSDetId/interface/TotemTimingDetId.h"

#include "DataFormats/CTPPSReco/interface/TotemTimingRecHit.h"

#include "RecoCTPPS/TotemRPLocal/interface/TotemTimingConversions.h"

#include "TH1.h"
#include "TGraph.h"

namespace pt = boost::property_tree;

class SimpleAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SimpleAnalyzer(const edm::ParameterSet&);
      ~SimpleAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
     struct calFileId{
       int db;
       int sampic;
       int channel;

       calFileId()
       : db(-1), sampic(-1), channel(-1){}

       calFileId(int db, int sampic, int channel)
       : db(db), sampic(sampic), channel(channel){}
     };

      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void initHistograms( const TotemTimingDetId&);

      // ---------- constants ---------------------------
      static const double SAMPIC_MAX_NUMBER_OF_SAMPLES;
      static const double SAMPIC_ADC_V;
      static const double SAMPIC_SAMPLING_PERIOD_NS;

      // ---------- objects to retrieve ---------------------------
      edm::EDGetTokenT< edm::DetSetVector<TotemTimingDigi> > tokenDigi_;
      edm::EDGetTokenT< edm::DetSetVector<TotemTimingRecHit> > tokenRecHit_;


      //----------- dependencies----------------------------
      std::map<TotemTimingDetId, std::vector<TotemTimingDetId>> substr_map_;

      // ---------- directories ---------------------------
      std::map< TotemTimingDetId, TFileDirectory > maindir_map_;

      // ---------- histograms ---------------------------
      std::map< TotemTimingDetId, TH1F*> yHisto_map_;
      std::map< TotemTimingDetId, TH1F*> tHisto_map_;
      std::map< std::pair<TotemTimingDetId, TotemTimingDetId>, TH1F*> tDiffHisto_map_;

      TH1F* tHisto_tot_;

      // ---------- graphs ---------------------------
      std::map< TotemTimingDetId, TGraph*> samples_graph_map_;
      std::map< TotemTimingDetId, unsigned int> samples_graphCtr_map_;
      //TotemTimingConversions conv, conv2;
      int ct;
      std::map<TotemTimingDetId, calFileId> detId_calFileId_map_;
};

template <typename Tree>
Tree query(Tree& pt, typename Tree::path_type path) {
    if (path.empty())
        return pt;

    auto const head = path.reduce();

    auto subscript = head.find('[');
    auto name      = head.substr(0, subscript);
    auto index     = std::string::npos != subscript && head.back() == ']'
        ? std::stoul(head.substr(subscript+1))
        : 0u;

    auto matches = pt.equal_range(name);
    if (matches.first==matches.second)
        throw std::out_of_range("name:" + name);

    for (; matches.first != matches.second && index; --index)
        ++matches.first;

    return query(matches.first->second, path);
}

//
// constants, enums and typedefs
//
const double    SimpleAnalyzer::SAMPIC_MAX_NUMBER_OF_SAMPLES = 64;
const double    SimpleAnalyzer::SAMPIC_ADC_V = 1./256;
const double    SimpleAnalyzer::SAMPIC_SAMPLING_PERIOD_NS = 1./7.8;

//
// static data member definitions
//

//
// constructors and destructor
//
SimpleAnalyzer::SimpleAnalyzer(const edm::ParameterSet& iConfig)
 :
 tokenDigi_            ( consumes< edm::DetSetVector<TotemTimingDigi> >      ( iConfig.getParameter<edm::InputTag>( "tagDigi" ) ) ),
 tokenRecHit_          ( consumes< edm::DetSetVector<TotemTimingRecHit> >      ( iConfig.getParameter<edm::InputTag>( "tagRecHit" ) ) )
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tHisto_tot_ = fs->make<TH1F>("t total", "t total", 100, -20, 20 );

}


SimpleAnalyzer::~SimpleAnalyzer()
{}


//
// member functions
//

void
SimpleAnalyzer::initHistograms(const TotemTimingDetId& detId)
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  if (maindir_map_.find(detId) == maindir_map_.end())
  {
    std::string dirName;
    detId.rpName(dirName, TotemTimingDetId::nPath);
    std::string chName;
    detId.channelName(chName, TotemTimingDetId::nFull);

    // create directory for the detector, if not already done
    maindir_map_[detId] = fs->mkdir( dirName);

    // create all histograms
    // std::string yHisto_name(chName);
    // yHisto_name.insert(0, "yDistribution_");
    // yHisto_map_[ detId ] = maindir_map_[ detId ].make<TH1F>(yHisto_name.c_str(), yHisto_name.c_str(), 120, -60, 60 );

    std::string tHisto_name(chName);
    tHisto_name.insert(0, "tDistribution_");
    tHisto_map_[ detId ] = maindir_map_[ detId ].make<TH1F>(tHisto_name.c_str(), tHisto_name.c_str(), 100, -4, 4); //-4 4


    // std::string plName;
    // detId.planeName(plName, TotemTimingDetId::nFull);
    // std::string tDiffHisto_name(plName);
    // tDiffHisto_name.insert(0, "tDiffWithPl+1_");
    // tDiffHisto_map_[ detId ] = maindir_map_[ detId ].make<TH1F>(tDiffHisto_name.c_str(), tDiffHisto_name.c_str(), 100, -20, 20 ); //-2 2

  }

}

// ------------ method called for each event  ------------
void
SimpleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  edm::Service<TFileService> fs;
  edm::Handle< edm::DetSetVector<TotemTimingDigi> > timingDigi;
  edm::Handle< edm::DetSetVector<TotemTimingRecHit> > timingRecHit;
  iEvent.getByToken( tokenDigi_, timingDigi );
  iEvent.getByToken( tokenRecHit_, timingRecHit );


  for (const auto& digis : *timingDigi)
  {
    // const TotemTimingDetId detId( digis.detId() );
    // if (maindir_map_.find(detId) == maindir_map_.end())
    //   initHistograms( detId );
    // std::string chName;
    // detId.channelName(chName, TotemTimingDetId::nFull);
    // std::string samples_graph_name(chName);
    // samples_graph_name.insert(0, "_cal_");
    // samples_graph_name.insert(0, std::to_string(ct));
    // samples_graph_name.insert(0, "samples_graph_");
    // TGraph *graph_hndl = maindir_map_[ detId ].make<TGraph>();
    // graph_hndl->SetName(samples_graph_name.c_str());
    // graph_hndl->SetTitle(samples_graph_name.c_str());
    // std::vector<float> voltage;
    //
    // std::string samples_graph_name2(chName);
    // samples_graph_name.insert(0, "_noncal_");
    // samples_graph_name.insert(0, std::to_string(ct));
    // samples_graph_name.insert(0, "samples_graph_");
    // TGraph *graph_hndl2 = maindir_map_[ detId ].make<TGraph>();
    // graph_hndl2->SetName(samples_graph_name.c_str());
    // graph_hndl2->SetTitle(samples_graph_name.c_str());
    // std::vector<float> voltage2;
    // ct++;
    const TotemTimingDetId detId(digis.detId());
    for (const auto& digi : digis )
    {
      detId_calFileId_map_[detId] = calFileId(digi.getHardwareBoardId(), digi.getHardwareSampicId(), digi.getHardwareChannelId());
      // // Do stuff on DIGIs
      // voltage = conv.getVoltSamples(digi);
      // voltage2 = conv2.getVoltSamples(digi);
      // int i=0;
      // for (auto sampIt = voltage.begin(); sampIt != voltage.end(); ++sampIt){
      //   // std::cout << *sampIt << ", ";
      //   graph_hndl->SetPoint(i, SAMPIC_SAMPLING_PERIOD_NS*(i), *sampIt  );
      //   i++;
      // }
      // // std::cout << "\n";
      //
      // i=0;
      // for (auto sampIt = voltage2.begin(); sampIt != voltage2.end(); ++sampIt){
      //   graph_hndl2->SetPoint(i, SAMPIC_SAMPLING_PERIOD_NS*(i), *sampIt  );
      //   i++;
      // }
    }
  }

  std::map< TotemTimingDetId, bool> onehit;


  for (const auto& recHits : *timingRecHit){
    const TotemTimingDetId planeId(recHits.detId());
    if(onehit.find(planeId.getPlaneId()) == onehit.end())
      onehit[planeId.getPlaneId()]=true;
    else
      onehit[planeId.getPlaneId()]=false;
  }
  for (const auto& recHits : *timingRecHit){
    for (const auto& recHit : recHits ){
      const TotemTimingDetId detId1(recHits.detId());
      if(onehit[detId1.getPlaneId()]){
        for (const auto& recHits2 : *timingRecHit){
          const TotemTimingDetId detId2( recHits2.detId());
          if  (onehit[detId2.getPlaneId()] &&
              (detId1.arm() == detId2.arm()) &&
              (detId1.station() == detId2.station()) &&
              (detId1.rp() == detId2.rp()) &&
              (detId1.plane()!= detId2.plane())){
            std::string pl1Name, pl2Name;
            detId1.planeName(pl1Name, TotemTimingDetId::nFull);
            detId2.planeName(pl2Name, TotemTimingDetId::nFull);
            std::string tDiffHisto_name( pl1Name + "-" + pl2Name);
            if(maindir_map_.find(detId2.getPlaneId()) == maindir_map_.end())
              maindir_map_[detId2.getPlaneId()] = fs->mkdir("-"+ pl2Name);
            if(tDiffHisto_map_.find(std::make_pair(detId1.getPlaneId(), detId2.getPlaneId())) == tDiffHisto_map_.end()){
              tDiffHisto_map_[std::make_pair(detId1.getPlaneId(), detId2.getPlaneId())] = maindir_map_[detId2.getPlaneId()].make<TH1F>(tDiffHisto_name.c_str(), tDiffHisto_name.c_str(), 100, -20, 20 );
              if(substr_map_.find(detId2.getPlaneId()) == substr_map_.end()){
                std::vector<TotemTimingDetId> a;
                substr_map_[detId2.getPlaneId()] = a;
              }
              substr_map_[detId2.getPlaneId()].push_back(detId1.getPlaneId());
            }
            for (const auto& recHit2 : recHits2 ){
              tDiffHisto_map_[std::make_pair(detId1.getPlaneId(), detId2.getPlaneId())]->Fill(recHit.getT() - recHit2.getT());
            }
          }
        }
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void
SimpleAnalyzer::beginJob()
{
  // ct=0;
  // conv = TotemTimingConversions();
  // conv2 = TotemTimingConversions();
  // conv.openCalibrationFile("/afs/cern.ch/user/f/fdej/sampic_cal/correct_offsets.cal.json");
}

// ------------ method called once each job just after ending the event loop  ------------
void
SimpleAnalyzer::endJob(){
  for(const auto& substrahend : substr_map_){
    std::cout << substrahend.first;
    std::string samples_graph_name;
    substrahend.first.planeName(samples_graph_name, TotemTimingDetId::nFull);
    TGraph *graph_hndl = maindir_map_[substrahend.first].make<TGraph>();
    graph_hndl->SetName(samples_graph_name.c_str());
    graph_hndl->SetTitle(samples_graph_name.c_str());
    //create scatter plot in maindir_map_
    for(const auto& minuend : substrahend.second){
      graph_hndl->SetPoint(minuend.plane(), minuend.plane(), tDiffHisto_map_[std::make_pair(minuend, substrahend.first)]->GetRMS());
      std::cout  << " - "<< minuend;
    }
    std::cout << "\n";
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleAnalyzer);
