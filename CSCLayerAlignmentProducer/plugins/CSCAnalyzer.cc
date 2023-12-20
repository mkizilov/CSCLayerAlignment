#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h>
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"



//Refitter

#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h" //from MK for Refit
 #include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h" //TA: chck if needed for refit
 #include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h" //TA: check if needed for refit
 #include "TrackingTools/KalmanUpdators/interface/KFUpdator.h" //TA: check if needed for refit
 #include "TrackingTools/Records/interface/TransientRecHitRecord.h" //TA: check if needed for refit
 #include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h" //TA: check if needed for refit
 #include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h" //TA: check if needed for refit
 #include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h" //TA: check if needed for refit
 #include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmBase.h" //TA: check if needed for refit:ConstTrajTrackPairs type


using namespace std;
using namespace edm;


struct CSCData
{
  void init();
  TTree* book(TTree *t);
  //Muon Info//////////////////////////////////////////////////////
  int muon_charge; float muon_pt; float muon_eta; float muon_momentum;
  unsigned long long  evtNum; unsigned long long  lumiBlock; int muonIdx;
  int runNum;
  //Propagation Info//////////////////////////////////////////////////////
  float prop_GP[3]; float prop_LP[3]; float prop_startingPoint_GP[3];
  float prop_localphi_rad; float prop_localphi_deg;
  float prop_globalphi_rad;
  bool has_prop; bool has_fidcut;
  int prop_location[5];
  //Track Info//////////////////////////////////////////////////////
  float track_chi2; float track_ndof; int n_ME11_segment; int which_track;
  int hasME11; int hasME11RecHit; int hasME11A; int hasME11ARecHit;
  int nCSCSeg; int nDTSeg; int nME11RecHits; float ME11_BunchX; int ME11_strip;
  int ME11_location[5];
  //Rechit Info//////////////////////////////////////////////////////
  float rechit_GP[3]; float rechit_LP[3];
  float rechit_localphi_rad; float rechit_localphi_deg;
  bool has_rechit;
  float RdPhi; int rechit_detId;
  int rechit_location[5];
  //Sim info for MC
  float sim_GP[3]; float sim_LP[3];
  float simDy; float sim_yroll; int nSim;
};

void CSCData::init()
{
  //Muon Info//////////////////////////////////////////////////////
  muon_charge = 9999; muon_pt = 9999; muon_eta = 9999; muon_momentum = 9999;
  evtNum = 99999999; lumiBlock = 99999999; muonIdx = 99999999; runNum = 99999999;
  //Propagation Info//////////////////////////////////////////////////////
  for(int i=0; i<3; ++i){
    prop_GP[i] = 99999; prop_LP[i] = 99999; prop_startingPoint_GP[i] = 99999;
  }
  prop_localphi_rad = 99999; prop_localphi_deg = 99999;
  prop_globalphi_rad = 99999;
  has_prop = false; has_fidcut = false;
  for(int i=0; i<5; ++i){
    prop_location[i] = 99999;
  }
  //Track Info//////////////////////////////////////////////////////
  track_chi2 = 999999; track_ndof = 999999; n_ME11_segment = 999999; which_track = 999999;
  hasME11 = 0; hasME11RecHit = 0; hasME11A = 0; hasME11ARecHit = 0;
  nCSCSeg = 999999; nDTSeg = 999999; nME11RecHits = 999999; ME11_BunchX = 999999; ME11_strip = 999999;
  for(int i=0; i<5; ++i){
    ME11_location[i] = 999999;
  }
  //Rechit Info//////////////////////////////////////////////////////
  for(int i=0; i<3; ++i){
    rechit_GP[i] = 999999; rechit_LP[i] = 999999;
  }
  rechit_localphi_rad = 999999; rechit_localphi_deg = 999999;
  has_rechit = false;
  RdPhi = 999999; rechit_detId = 999999;
  for(int i=0; i<5; ++i){
    rechit_location[i] = 999999;
  }
  //Sim info for MC
  for(int i=0; i<3; ++i){
    sim_GP[i] = 9999999; sim_LP[i] = 9999999;
  }
  simDy = 9999999; sim_yroll = 9999999; nSim = 9999999;
}

TTree* CSCData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("Inner_Prop", "Inner_Prop");
  //Muon Info//////////////////////////////////////////////////////
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt);
  t->Branch("muon_eta", &muon_eta); t->Branch("muon_momentum", &muon_momentum);
  t->Branch("evtNum", &evtNum); t->Branch("lumiBlock", &lumiBlock); t->Branch("muonIdx", &muonIdx);
  t->Branch("runNum", &runNum);
  //Propagation Info//////////////////////////////////////////////////////
  t->Branch("prop_GP", &prop_GP, "prop_GP[3] (x,y,z)/F");
  t->Branch("prop_LP", &prop_LP, "prop_LP[3] (x,y,z)/F");
  t->Branch("prop_startingPoint_GP", &prop_startingPoint_GP, "prop_startingPoint_GP[3] (x,y,z)/F");
  t->Branch("prop_localphi_rad", &prop_localphi_rad);
  t->Branch("prop_localphi_deg", &prop_localphi_deg);
  t->Branch("prop_globalphi_rad", &prop_globalphi_rad);
  t->Branch("has_prop", &has_prop);
  t->Branch("has_fidcut", &has_fidcut);
  t->Branch("prop_location", &prop_location, "prop_location[5] (reg, sta, ring, cha, lay)/I");
  //Track Info//////////////////////////////////////////////////////
  t->Branch("track_chi2", &track_chi2); t->Branch("track_ndof", &track_ndof);
  t->Branch("n_ME11_segment", &n_ME11_segment); t->Branch("which_track", &which_track);
  t->Branch("hasME11", &hasME11); t->Branch("hasME11RecHit", &hasME11RecHit);
  t->Branch("hasME11A", &hasME11A); t->Branch("hasME11ARecHit", &hasME11ARecHit);
  t->Branch("nCSCSeg", &nCSCSeg); t->Branch("nDTSeg", &nDTSeg);
  t->Branch("nME11RecHits", &nME11RecHits); t->Branch("ME11_BunchX", &ME11_BunchX);
  t->Branch("ME11_strip", &ME11_strip);
  t->Branch("ME11_location", &ME11_location, "ME11_location[5] (end, sta, ring, cha, lay)/I");
  //Rechit Info//////////////////////////////////////////////////////
  t->Branch("rechit_GP", &rechit_GP, "rechit_GP[3] (x,y,z)/F");
  t->Branch("rechit_LP", &rechit_LP, "rechit_LP[3] (x,y,z)/F");
  t->Branch("rechit_localphi_rad", &rechit_localphi_rad);
  t->Branch("rechit_localphi_deg", &rechit_localphi_deg);
  t->Branch("has_rechit", &has_rechit);
  t->Branch("RdPhi", &RdPhi);
  t->Branch("rechit_detId", &rechit_detId);
  t->Branch("rechit_location", &rechit_location, "rechit_location[5] (reg, sta, ring, cha, lay)/I");
  //Sim info for MC
  t->Branch("sim_GP", &sim_GP, "sim_GP[3] (x,y,z)/F");
  t->Branch("sim_LP", &sim_LP, "sim_LP[3] (x,y,z)/F");
  t->Branch("simDy", &simDy);
  t->Branch("sim_yroll", &sim_yroll);
  t->Branch("nSim", &nSim);
  return t;
}


class CSCAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit CSCAnalyzer(const edm::ParameterSet&);
  ~CSCAnalyzer(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  void propagate(const reco::Muon* mu, const edm::Event& iEvent, int i, const Trajectory* traj_of_muon);
  void CSCSegmentCounter(const reco::Muon* mu, CSCData& data_);
  void propagate_to_ME11(const reco::Muon* mu, const CSCLayer* ch, bool &tmp_has_prop, GlobalPoint &pos_GP, CSCData& data_, const Trajectory* traj_of_muon);
  void ME11_rechit_matcher(const CSCLayer* ch, LocalPoint prop_LP, CSCData& data_);
  float RdPhi_func(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, const GEMEtaPartition* ch);
  bool fidcutCheck(float local_y, float local_x, float localphi_deg);

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::Handle<GEMRecHitCollection> gemRecHits;

  edm::EDGetTokenT<vector<PSimHit> > gemSimHits_;
  edm::Handle<vector<PSimHit> > gemSimHits;

  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::Handle<View<reco::Muon> > muons;

  edm::Handle<TrajTrackAssociationCollection> ref_track;
  edm::EDGetTokenT<TrajTrackAssociationCollection> ref_track_;

  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  edm::Handle<CSCSegmentCollection> cscSegments;

  edm::EDGetTokenT<CSCRecHit2DCollection> csc2DRecHits_;
  edm::Handle<CSCRecHit2DCollection> csc2DRecHits;


  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;

  edm::ESHandle<GEMGeometry> GEMGeometry_;
  edm::ESHandle<CSCGeometry> CSCGeometry_;

  bool CSC_prop; bool tracker_prop; bool Segment_prop;
  vector<int> prop_list;
  bool debug;
  bool refitter;
  bool isCosmic;

  CSCData data_;
  TTree* Tracker_tree;
  TH2D* nME11_col_vs_matches = new TH2D("nME11_test", "nME11_test", 5, 0, 5, 5, 0, 5);

  bool isMC;
  const CSCSegment *ME11_segment;

  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;


};


CSCAnalyzer::CSCAnalyzer(const edm::ParameterSet& iConfig)
  : gemGeomToken_(esConsumes()),
    cscGeomToken_(esConsumes()),
    ttkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    geomToken_(esConsumes())
{
  cout << "Begin analyzer" << endl;
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());

  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  gemSimHits_ = consumes<vector<PSimHit> >(iConfig.getParameter<edm::InputTag>("gemSimHits"));
  cscSegments_ = consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"));
  csc2DRecHits_ = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("csc2DRecHits"));
  ref_track_ = consumes<TrajTrackAssociationCollection>(iConfig.getParameter<InputTag>("ref_track"));

  debug = iConfig.getParameter<bool>("debug");
  isCosmic = iConfig.getParameter<bool>("isCosmic");
  refitter = iConfig.getParameter<bool>("refitter");


  Tracker_tree = data_.book(Tracker_tree);

}


void CSCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);
  //iSetup.get<MuonGeometryRecord>().get(CSCGeometry_);
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  //iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  GEMGeometry_ = &iSetup.getData(gemGeomToken_);
  CSCGeometry_ = &iSetup.getData(cscGeomToken_);
  ttrackBuilder_ = &iSetup.getData(ttkToken_);
  theTrackingGeometry = &iSetup.getData(geomToken_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  isMC = false;
  if (! iEvent.eventAuxiliary().isRealData()) isMC = true;
  iEvent.getByToken(gemRecHits_, gemRecHits);
  if (isMC) {
    iEvent.getByToken(gemSimHits_, gemSimHits);
  }
  if (! iEvent.getByToken(muons_, muons)) return;
  if (muons->size() == 0) return;
  
  edm::Handle<TrajTrackAssociationCollection> ref_track;
  iEvent.getByToken(cscSegments_, cscSegments);
  iEvent.getByToken(csc2DRecHits_, csc2DRecHits);

  iEvent.getByToken(ref_track_, ref_track);
  ConstTrajTrackPairs ref_track_pairs;
  for (auto it = ref_track->begin(); it != ref_track->end(); ++it) {
    ref_track_pairs.push_back(ConstTrajTrackPair(&(*(*it).key), &(*(*it).val)));
  } //the loop goes over tracks and saves the key and value of each track as a pair in ref_track_pairs.
  //end of Refit trajectory

 

  if (debug) cout << "New! EvtNumber = " << iEvent.eventAuxiliary().event() << " LumiBlock = " << iEvent.eventAuxiliary().luminosityBlock() << " RunNumber = " << iEvent.run() << endl;

  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();
    if (not mu->standAloneMuon()) continue;
    if (!(mu->passed(reco::Muon::PFIsoTight))) continue; //MK ???
    if (debug) cout << "new standalone" << endl;

    //trajectory muon matching--we want to match muons with the tracks
    const Trajectory* traj_of_muon;
    const Trajectory* traj_of_Track;
    const reco::Track* track_of_Track;

        // propagate(mu, iEvent, i);

    // for (auto it = std::begin(prop_list); it != std::end(prop_list); ++it){
      // if (debug) std::cout << "\tprop " << *it << "about to start propagate" << std::endl;
      // int prop_type = *it;
      for (ConstTrajTrackPairs::const_iterator it = ref_track_pairs.begin(); it != ref_track_pairs.end(); ++it) {
        traj_of_Track = (*it).first;
        track_of_Track = (*it).second;
        if (track_of_Track == mu->track().get()) {
          traj_of_muon = traj_of_Track;
          if (debug) cout << "mu, prop_type, i, traj_of_muon: " << mu << i << "," << traj_of_muon << endl;
          propagate(mu, /*prop_type,*/ iEvent, i, traj_of_muon); //taking the matched trajectory of muon to the propagate function
        }
      }
    }
//}


  if (debug){
    int muon_size = muons->size();
    int CSCSegments_size = cscSegments->size();
    int muon_STA_counter = 0;
    int ME11_counter = 0;
    int Matches_counter_norepeat = 0;
    int ME11_counter_norepeat = 0;
    int tot_muon_ME11_counter = 0;

    //Commented out for github, may be many print statements
    //cout << "Starting loop over all gemRecHits " << gemRecHits->size() << std::endl;
    //for(GEMRecHitCollection::const_iterator gemRecHit = gemRecHits->begin(); gemRecHit != gemRecHits->end(); gemRecHit++){
    //  cout << gemRecHit->gemId() << endl;
    //}

    std::vector<int> list_of_matches = {};
    for (size_t i = 0; i < muons->size(); ++i){
      edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
      const reco::Muon* mu = muRef.get();
      if (mu->isStandAloneMuon()){
        muon_STA_counter ++;
        int muon_ME11_counter_tmp = 0;
        auto matches = mu->matches();
        for (auto MCM : matches){
          if(MCM.detector() != 2) continue;
          for(auto MSM : MCM.segmentMatches){
            auto cscSegRef = MSM.cscSegmentRef;
            auto cscDetID = cscSegRef->cscDetId();
            std::cout << cscDetID.endcap() << cscDetID.station() << cscDetID.ring() << cscDetID.chamber() << std::endl;
            if(cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
              muon_ME11_counter_tmp++;
              int tmp_ID = cscDetID.endcap()*10000+cscDetID.station()*1000+cscDetID.ring()*100+cscDetID.chamber();
              if (std::find(list_of_matches.begin(), list_of_matches.end(), tmp_ID) == list_of_matches.end()){
                list_of_matches.push_back(tmp_ID);
                // std::cout << "New Det Triggered " << tmp_ID << std::endl;
                Matches_counter_norepeat++;
              }
            }
          }
        }
        // std::cout << "Muon had " << muon_ME11_counter_tmp << " segments" << std::endl;
        tot_muon_ME11_counter += muon_ME11_counter_tmp;
      }
    }

    // std::cout << "SEGMENT COLLECTION" << std::endl;
    std::vector<int> list_of_collection = {};
    for(CSCSegmentCollection::const_iterator segmentCSC = cscSegments->begin(); segmentCSC != cscSegments->end(); segmentCSC++){
      auto cscDetID = (segmentCSC)->cscDetId();
      std::cout << cscDetID.endcap() << cscDetID.station() << cscDetID.ring() << cscDetID.chamber() << std::endl;
      if(cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
        ME11_counter++;
        int tmp_ID = cscDetID.endcap()*10000+cscDetID.station()*1000+cscDetID.ring()*100+cscDetID.chamber();
        if (std::find(list_of_collection.begin(), list_of_collection.end(), tmp_ID) == list_of_collection.end()){
          list_of_collection.push_back(tmp_ID);
          std::cout << "New Det Triggered " << tmp_ID << std::endl;
          ME11_counter_norepeat++;
        }
      }
    }
  // std::cout << "Muon size           = " << muon_size << std::endl;
  // std::cout << "Segment Size        = " << CSCSegments_size << std::endl;
  // std::cout << "muon STA  counter   = " << muon_STA_counter << std::endl;
  // std::cout << "Muon Matches ME11   = " << tot_muon_ME11_counter << std::endl;
  // std::cout << "Seg Collection ME11 = " << ME11_counter << std::endl;
  // std::cout << "Muon Matches no rep = " << Matches_counter_norepeat << std::endl;
  // std::cout << "Coll Matches no rep = " << ME11_counter_norepeat << std::endl;
  nME11_col_vs_matches->Fill(Matches_counter_norepeat, ME11_counter_norepeat);
  }
}


float CSCAnalyzer::RdPhi_func(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, const GEMEtaPartition* ch){
  GEMDetId gemid((rechit)->geographicalId());
  const auto& etaPart = GEMGeometry_->etaPartition(gemid);
  const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
  float deltay_roll =  etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp() - etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
  return cos(stripAngle) * (prop_localx - (rechit)->localPosition().x()) - sin(stripAngle) * (prop_localy + deltay_roll);
}
void CSCAnalyzer::CSCSegmentCounter(const reco::Muon* mu, CSCData& data_){
  const reco::Track* Track = mu->outerTrack().get();
  int tmp_CSC_counter = 0; int tmp_DT_counter = 0; int tmp_ME11_counter = 0; int tmp_ME11RecHit_counter = 0; float tmp_ME11_BunchX = 99999; int tmp_ME11_strip = 99999; bool tmp_hasME11A = 0;
  if(isCosmic){
    tmp_CSC_counter = mu->numberOfSegments(1,2) + mu->numberOfSegments(2,2) + mu->numberOfSegments(3,2) + mu->numberOfSegments(4,2);
    tmp_DT_counter = mu->numberOfSegments(1,1) + mu->numberOfSegments(2,1) + mu->numberOfSegments(3,1) + mu->numberOfSegments(4,1);
    auto matches = mu->matches();
    for (auto MCM : matches){
      if(MCM.detector() != 2) continue;
      for(auto MSM : MCM.segmentMatches){
        auto cscSegRef = MSM.cscSegmentRef;
        auto cscDetID = cscSegRef->cscDetId();
        if(cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
          if(cscDetID.ring() == 4){tmp_hasME11A = 1;}
          tmp_ME11_counter++;
          if (debug){std::cout << "isCosmic = True! Getting ME11 Segment" << std::endl;}
          ME11_segment = cscSegRef.get();
          tmp_ME11RecHit_counter = (cscSegRef.get())->nRecHits(); // Find the real function for this. Bad if multiple segments.
          tmp_ME11_BunchX = ME11_segment->time();
          auto cscDetID_FAKE = CSCDetId(cscDetID.endcap(), cscDetID.station(), cscDetID.ring(), cscDetID.chamber(), 3);
          const CSCLayer* tmp_ME11_layer = CSCGeometry_->layer(cscDetID_FAKE);
          const CSCLayerGeometry* tmp_ME11_layer_geo = tmp_ME11_layer->geometry();
          tmp_ME11_strip = tmp_ME11_layer_geo->nearestStrip(ME11_segment->localPosition());
          data_.ME11_location[0] = cscDetID.endcap(); data_.ME11_location[1] = cscDetID.station(); data_.ME11_location[2] = cscDetID.ring(); data_.ME11_location[3] = cscDetID.chamber(); data_.ME11_location[4] = cscDetID.layer();
        }
      }
    }
  }
  else{
    for (size_t RecHit_iter = 0; RecHit_iter != Track->recHitsSize(); RecHit_iter++){
      const TrackingRecHit* RecHit = (Track->recHit(RecHit_iter)).get();
      DetId RecHitId = RecHit->geographicalId();
      uint16_t RecHitDetId = RecHitId.det();
      if (RecHitDetId == DetId::Muon){
        uint16_t RecHitSubDet = RecHitId.subdetId();
        if (RecHitSubDet == (uint16_t)MuonSubdetId::CSC){
          if (CSCDetId(RecHitId).station() == 1 and CSCDetId(RecHitId).ring() == 1 and RecHit->dimension() == 4){
            tmp_ME11_counter++;
            if (debug){std::cout << "isCosmic = False! Getting ME11 Segment" << std::endl;}
            RecSegment* Rec_segment = (RecSegment*)RecHit;
            ME11_segment = (CSCSegment*)Rec_segment;
            tmp_ME11_BunchX = ((CSCRecHit2D*)RecHit)->wgroupsBX();
            auto cscDetID_FAKE = CSCDetId(CSCDetId(RecHitId).endcap(), CSCDetId(RecHitId).station(), CSCDetId(RecHitId).ring(), CSCDetId(RecHitId).chamber(), 3);
            const CSCLayer* tmp_ME11_layer = CSCGeometry_->layer(cscDetID_FAKE);
            const CSCLayerGeometry* tmp_ME11_layer_geo = tmp_ME11_layer->geometry();
            tmp_ME11_strip = tmp_ME11_layer_geo->nearestStrip(ME11_segment->localPosition());
            data_.ME11_location[0] = CSCDetId(RecHitId).endcap(); data_.ME11_location[1] = CSCDetId(RecHitId).station(); data_.ME11_location[2] = CSCDetId(RecHitId).ring(); data_.ME11_location[3] = CSCDetId(RecHitId).chamber(); data_.ME11_location[4] = CSCDetId(RecHitId).layer();
          }
          if (CSCDetId(RecHitId).station() == 1 and CSCDetId(RecHitId).ring() == 1){tmp_ME11RecHit_counter++;}
          if (RecHit->dimension() == 4){tmp_CSC_counter++;}
        }
        if (RecHitSubDet == (uint16_t)MuonSubdetId::DT){
          if (RecHit->dimension() > 1){tmp_DT_counter++;}
        }
      }
    }
  }
  data_.nCSCSeg = tmp_CSC_counter; data_.nDTSeg = tmp_DT_counter;
  data_.n_ME11_segment = tmp_ME11_counter;
  data_.nME11RecHits = tmp_ME11RecHit_counter;
  data_.ME11_BunchX = tmp_ME11_BunchX;
  data_.ME11_strip = tmp_ME11_strip;
  data_.hasME11A =  tmp_hasME11A;
  if(data_.n_ME11_segment >= 1 and data_.n_ME11_segment < 1000){data_.hasME11 = 1;}
}

void CSCAnalyzer::propagate(const reco::Muon* mu, const edm::Event& iEvent, int i, const Trajectory* traj_of_muon){
  const reco::Track* Track;
  reco::TransientTrack ttTrack;
  TTree* tree;
  tree = Tracker_tree;
  // std::cout << "propagate start" << std::endl;
  // std::cout<<"mu->track().isNonnull() = "<<mu->track().isNonnull()<<std::endl;
  // std::cout<<"mu->Track is available"<<mu->track().isAvailable()<<std::endl;
  if(!(mu->track().isNonnull())){return;}
  if (!(mu->isTrackerMuon())) {return;}
  // std::cout << "propagate Track" << std::endl;
  Track = mu->track().get();
  // std::cout << "propagate build Track" << std::endl;
  //Print Track
  std::cout << "Track->pt() = " << Track->pt() << std::endl;

  // std::cout << "Track->pt() = " << Track->pt() << std::endl;
  // std::cout << "Track->eta() = " << Track->eta() << std::endl;
  // std::cout << "Track->phi() = " << Track->phi() << std::endl;
  // std::cout << "Track->charge() = " << Track->charge() << std::endl;


  ttTrack = ttrackBuilder_->build(Track); // Crashes here

  if(!ttTrack.isValid()){std::cout << "BAD EVENT! NO TRACK" << std::endl;}
  data_.init();
  //Muon Info//////////////////////////////////////////////////////
  data_.muon_charge = mu->charge(); data_.muon_pt = mu->pt(); data_.muon_eta = mu->eta(); data_.muon_momentum = mu->momentum().mag2();
  data_.evtNum = iEvent.eventAuxiliary().event(); data_.lumiBlock = iEvent.eventAuxiliary().luminosityBlock(); data_.muonIdx = data_.evtNum*100 + i;
  data_.runNum = iEvent.run();
  //Track Info//////////////////////////////////////////////////////
  data_.track_chi2 = Track->chi2(); data_.track_ndof = Track->ndof();
  //Propagation Info//////////////////////////////////////////////////////
  for (const auto& ch : CSCGeometry_->layers()) {
    if (!(ch->id().station() == 1 and (ch->id().ring() == 1 or ch->id().ring() == 4))) continue;
    if (debug) std::cout << "Looping over CSC layers, at detID " << ch->id() << std::endl;
    GlobalPoint tmp_prop_GP; bool tmp_has_prop = 0;
    propagate_to_ME11(mu, ch, tmp_has_prop, tmp_prop_GP, data_, traj_of_muon);
    if (!tmp_has_prop) continue;
    if (debug) std::cout << "Found a prop!" << std::endl;
    LocalPoint tmp_prop_LP = ch->toLocal(tmp_prop_GP);
    //Rechit Info//////////////////////////////////////////////////////
    ME11_rechit_matcher(ch, tmp_prop_LP, data_);
    tree->Fill();
  }
}


bool CSCAnalyzer::fidcutCheck(float local_y, float local_x, float localphi_deg){
  // Fidcut same as in TBMA https://github.com/cms-sw/cmssw/blob/deaac86743f80cda845f794c9335ba27a4d50417/Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsFitter.h
  const float fidcut_angle = 1.0;
  const float ymin = -80.0;
  const float ymax =  80.0;
  const float xmin = -80.0;
  const float xmax =  80.0;
  const float cut_angle = 5.0 - fidcut_angle;
  if ((abs(localphi_deg) < cut_angle) && (ymin < local_y && local_y < ymax) && (xmin<local_x && local_x < xmax)){return 1;}
  else{return 0;}
}

void CSCAnalyzer::propagate_to_ME11(const reco::Muon* mu, const CSCLayer* ch, bool &tmp_has_prop, GlobalPoint &pos_GP, CSCData& data_, const Trajectory* traj_of_muon){
  const reco::Track* Track;
  reco::TransientTrack track;
  tmp_has_prop = false;
  const BoundPlane& bps(ch->surface());
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  double previous_trackTSOS_globalPositionR = 0.0; 

  TrajectoryStateOnSurface tsos;
  TrajectoryStateOnSurface previous_trackTSOS;
  std::vector<TrajectoryMeasurement> traj_measurement = traj_of_muon->measurements();
  // std::out<<"debug traj measurment size: "<< traj_measurement.size() << std::endl;

  TrajectoryStateOnSurface tsos_ch; TrajectoryStateOnSurface tsos_seg;
  GlobalPoint pos_startingPoint_GP;

  GlobalPoint pos_center_of_layer = ch->centerOfStrip(ch->geometry()->numberOfStrips()/2.0);

  // if (!(mu->isTrackerMuon())) return; //for refitter

  Track = mu->track().get();
  track = ttrackBuilder_->build(Track);

 //For refitter start
      //  for (std::vector<TrajectoryMeasurement>::const_iterator it_traj_measurement = traj_measurement.begin(); it_traj_measurement != traj_measurement.end(); ++it_traj_measurement){
      //   TrajectoryMeasurement iTraj_measurement = *it_traj_measurement;
      //   TrajectoryStateOnSurface tsos = TrajectoryStateCombiner().combine(iTraj_measurement.forwardPredictedState(), iTraj_measurement.backwardPredictedState());
      //   if (tsos.isValid()){
      //     double tsosGlobalPositionR = sqrt(tsos.globalPosition().x() * tsos.globalPosition().x() +
      //                                       tsos.globalPosition().y() * tsos.globalPosition().y());
       
      //     if (tsosGlobalPositionR > previous_trackTSOS_globalPositionR){
      //       previous_trackTSOS = tsos;
      //       previous_trackTSOS_globalPositionR = tsosGlobalPositionR;
      //     }          
      //   }
      // }

      // tsos_ch = propagator->propagate(previous_trackTSOS, ch->surface());

      // if (tsos_ch.isValid()){
      //   const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
      //   const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
      //   if (!(tsos_ch.globalPosition().z() * previous_trackTSOS.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1) {
      //     tmp_has_prop = true;
      //     pos_GP = tsos_ch.globalPosition();
      //     pos_startingPoint_GP = previous_trackTSOS.globalPosition();
      //   }
      // }
  //for refitter end

  if(!refitter){
    float inner_delta = abs(track.innermostMeasurementState().globalPosition().z() - pos_center_of_layer.z());
    float outer_delta = abs(track.outermostMeasurementState().globalPosition().z() - pos_center_of_layer.z());
    float used_delta = 0;

    if (inner_delta < outer_delta){
      tsos_seg = track.innermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface()); used_delta = inner_delta;
      data_.which_track = 0;
    }
    else{
      tsos_seg = track.outermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface()); used_delta = outer_delta;
      data_.which_track = 1;
    }
    if (tsos_ch.isValid()){
      const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
      const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
      if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1){
        tmp_has_prop = true;
        if(debug){std::cout << "Delta to CSC!!! = " << used_delta << std::endl;}
        pos_GP = tsos_ch.globalPosition();
        pos_startingPoint_GP = tsos_seg.globalPosition();
      }
    }
  }
  else{
    for (std::vector<TrajectoryMeasurement>::const_iterator it_traj_measurement = traj_measurement.begin(); it_traj_measurement != traj_measurement.end(); ++it_traj_measurement){
        TrajectoryMeasurement iTraj_measurement = *it_traj_measurement;
        TrajectoryStateOnSurface tsos = TrajectoryStateCombiner().combine(iTraj_measurement.forwardPredictedState(), iTraj_measurement.backwardPredictedState());
        if (tsos.isValid()){
          double tsosGlobalPositionR = sqrt(tsos.globalPosition().x() * tsos.globalPosition().x() +
                                            tsos.globalPosition().y() * tsos.globalPosition().y());
       
          if (tsosGlobalPositionR > previous_trackTSOS_globalPositionR){
            previous_trackTSOS = tsos;
            previous_trackTSOS_globalPositionR = tsosGlobalPositionR;
          }          
        }
      }

      tsos_ch = propagator->propagate(previous_trackTSOS, ch->surface());

      if (tsos_ch.isValid()){
        const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
        const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
        if (!(tsos_ch.globalPosition().z() * previous_trackTSOS.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1) {
          tmp_has_prop = true;
          pos_GP = tsos_ch.globalPosition();
          pos_startingPoint_GP = previous_trackTSOS.globalPosition();
        }
      }
    }
  

  if (tmp_has_prop){
    LocalPoint tmp_prop_LP = ch->toLocal(pos_GP);
    data_.prop_GP[0] = pos_GP.x(); data_.prop_GP[1] = pos_GP.y(); data_.prop_GP[2] = pos_GP.z();
    data_.prop_LP[0] = tmp_prop_LP.x(); data_.prop_LP[1] = tmp_prop_LP.y(); data_.prop_LP[2] = tmp_prop_LP.z();
    data_.prop_startingPoint_GP[0] = pos_startingPoint_GP.x(); data_.prop_startingPoint_GP[1] = pos_startingPoint_GP.y(); data_.prop_startingPoint_GP[2] = pos_startingPoint_GP.z();
    //std::cout << "GP     = " << pos_GP << std::endl;
    //std::cout << "mag    = " << pos_GP.mag() << std::endl;
    //std::cout << "mag2   = " << pos_GP.mag2() << std::endl;
    //std::cout << "perp2  = " << pos_GP.perp2() << std::endl;
    //std::cout << "perp   = " << pos_GP.perp() << std::endl;
    float radius = pos_GP.perp();
    LocalPoint local_to_center(tmp_prop_LP.x(), tmp_prop_LP.y() + radius, 0); //tmp_prop_LP is from center of chamber, not center of endcap
    float local_phi = local_to_center.phi();
    if(debug){
      std::cout << "Local nofix = " << tmp_prop_LP << std::endl;
      std::cout << "Local point = " << local_to_center << std::endl;
      std::cout << "Local phi   = " << ((3.14159265/2.) - local_phi)*(180./3.14159265) << std::endl;
    }

    //int strip = ch->geometry()->nearestStrip(tmp_prop_LP);
    //double stripAngle = ch->geometry()->stripAngle(strip) - M_PI/2.;


    //std::cout << "Strip angle = " << stripAngle << std::endl;
    //std::cout << "Local point = " << local_to_center << std::endl;
    data_.prop_localphi_rad = (3.14159265/2.) - local_phi;
    data_.prop_localphi_deg = ((3.14159265/2.) - local_phi)*(180./3.14159265);
    data_.prop_globalphi_rad = pos_GP.phi();
    data_.has_prop = tmp_has_prop;
    data_.has_fidcut = fidcutCheck(tmp_prop_LP.y(),tmp_prop_LP.x(), ((3.14159265/2.) - local_phi)*(180./3.14159265));
    data_.prop_location[0] = ch->id().zendcap(); data_.prop_location[1] = ch->id().station(); data_.prop_location[2] = ch->id().ring(); data_.prop_location[3] = ch->id().chamber(); data_.prop_location[4] = ch->id().layer();
  }
}



// Need to add     ME11_rechit_matcher(ch, tmp_prop_LP, data_); here!!! Basically uses the GEM_rechit_matcher code
void CSCAnalyzer::ME11_rechit_matcher(const CSCLayer* ch, LocalPoint prop_LP, CSCData& data_){
  float tmp_rechit_GP_x; float tmp_rechit_GP_y; float tmp_rechit_GP_z;
  float tmp_rechit_LP_x; float tmp_rechit_LP_y; float tmp_rechit_LP_z;
  float tmp_rechit_localphi_rad; float tmp_rechit_localphi_deg;
  bool tmp_has_rechit = false;
  float tmp_RdPhi = 9999.; int tmp_rechit_detId;
  int tmp_rechit_region; int tmp_rechit_station; int tmp_rechit_ring; int tmp_rechit_chamber; int tmp_rechit_layer;
  for (auto hit = csc2DRecHits->begin(); hit != csc2DRecHits->end(); hit++){
    CSCDetId cscid((hit)->geographicalId());
    if (!(cscid == ch->id())) continue;
    if (!(ch->id().station() == 1 and ch->id().ring() == 1 and fabs((hit)->localPosition().x() - prop_LP.x()) < 999.0)) continue;

    int strip = ch->geometry()->nearestStrip(hit->localPosition());
    double stripAngle = ch->geometry()->stripAngle(strip) - M_PI/2.;
    if (abs(tmp_RdPhi) < abs(cos(stripAngle) * (prop_LP.x() - (hit)->localPosition().x()) + sin(stripAngle) * (prop_LP.y() - (hit)->localPosition().y()))) continue; //MK added abs for right part to test
    tmp_RdPhi = cos(stripAngle) * (prop_LP.x() - (hit)->localPosition().x()) + sin(stripAngle) * (prop_LP.y() - (hit)->localPosition().y());  // yes, that's +sin() //Gotten from CSC Alignment Code
    tmp_rechit_GP_x = ch->toGlobal((hit)->localPosition()).x(); tmp_rechit_GP_y = ch->toGlobal((hit)->localPosition()).y(); tmp_rechit_GP_z = ch->toGlobal((hit)->localPosition()).z();
    tmp_rechit_LP_x = (hit)->localPosition().x(); tmp_rechit_LP_y = (hit)->localPosition().y(); tmp_rechit_LP_z = (hit)->localPosition().z();
    float local_phi = prop_LP.phi();
    tmp_rechit_localphi_rad = (3.14159265/2.) - local_phi;
    tmp_rechit_localphi_deg = ((3.14159265/2.) - local_phi)*(180./3.14159265);
    tmp_has_rechit = true;
    tmp_rechit_detId = cscid.zendcap()*(cscid.station()*100 + cscid.chamber());
    tmp_rechit_region = cscid.zendcap(); tmp_rechit_station = cscid.station(); tmp_rechit_ring = cscid.ring(); tmp_rechit_chamber = cscid.chamber(); tmp_rechit_layer = cscid.layer();


  }
  if(tmp_has_rechit){
    data_.rechit_GP[0] = tmp_rechit_GP_x; data_.rechit_GP[1] = tmp_rechit_GP_y; data_.rechit_GP[2] = tmp_rechit_GP_z;
    data_.rechit_LP[0] = tmp_rechit_LP_x; data_.rechit_LP[1] = tmp_rechit_LP_y; data_.rechit_LP[2] = tmp_rechit_LP_z;
    data_.rechit_localphi_rad = tmp_rechit_localphi_rad;
    data_.rechit_localphi_deg = tmp_rechit_localphi_deg;
    data_.has_rechit = tmp_has_rechit;
    data_.RdPhi = tmp_RdPhi;
    data_.rechit_detId = tmp_rechit_detId;
    data_.rechit_location[0] = tmp_rechit_region; data_.rechit_location[1] = tmp_rechit_station; data_.rechit_location[2] = tmp_rechit_ring; data_.rechit_location[3] = tmp_rechit_chamber; data_.rechit_location[4] = tmp_rechit_layer;
  }
}




void CSCAnalyzer::beginJob(){}
void CSCAnalyzer::endJob(){
  if(debug){nME11_col_vs_matches->Write();}
  }

DEFINE_FWK_MODULE(CSCAnalyzer);