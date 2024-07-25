#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Tracking and geometry includes
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

// Data formats
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/ESHandle.h"

// ROOT includes
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

// Additional CMSSW includes for track refitting
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmBase.h"

using namespace std;
using namespace edm;

// Structure to hold various CSC data
struct CSCData {
    void init();
    TTree* book(TTree* t);

    // Muon information
    int muon_charge; 
    float muon_pt; 
    float muon_eta; 
    float muon_momentum;
    unsigned long long evtNum; 
    unsigned long long lumiBlock; 
    int muonIdx;
    int runNum;

    // Propagation information
    float prop_GP[3]; 
    float prop_LP[3]; 
    float prop_startingPoint_GP[3];
    float prop_localphi_rad; 
    float prop_localphi_deg;
    float prop_globalphi_rad;
    bool has_prop; 
    bool has_fidcut;
    int prop_location[5];

    // Track information
    float track_chi2; 
    float track_ndof; 
    int n_ME11_segment; 
    int which_track;
    int hasME11; 
    int hasME11RecHit; 
    int hasME11A; 
    int hasME11ARecHit;
    int nCSCSeg; 
    int nDTSeg; 
    int nME11RecHits; 
    float ME11_BunchX; 
    int ME11_strip;
    int ME11_location[5];

    // Rechit information
    float rechit_GP[3]; 
    float rechit_LP[3];
    float rechit_localphi_rad; 
    float rechit_localphi_deg;
    bool has_rechit;
    float RdPhi; 
    int rechit_detId;
    int rechit_location[5];

    // Simulation information for MC
    float sim_GP[3]; 
    float sim_LP[3];
    float simDy; 
    float sim_yroll; 
    int nSim;
};

void CSCData::init() {
    // Initialize all members with default values
    muon_charge = 9999; 
    muon_pt = 9999; 
    muon_eta = 9999; 
    muon_momentum = 9999;
    evtNum = 99999999; 
    lumiBlock = 99999999; 
    muonIdx = 99999999; 
    runNum = 99999999;

    for (int i = 0; i < 3; ++i) {
        prop_GP[i] = 99999; 
        prop_LP[i] = 99999; 
        prop_startingPoint_GP[i] = 99999;
    }

    prop_localphi_rad = 99999; 
    prop_localphi_deg = 99999;
    prop_globalphi_rad = 99999;
    has_prop = false; 
    has_fidcut = false;

    for (int i = 0; i < 5; ++i) {
        prop_location[i] = 99999;
    }

    track_chi2 = 999999; 
    track_ndof = 999999; 
    n_ME11_segment = 999999; 
    which_track = 999999;
    hasME11 = 0; 
    hasME11RecHit = 0; 
    hasME11A = 0; 
    hasME11ARecHit = 0;
    nCSCSeg = 999999; 
    nDTSeg = 999999; 
    nME11RecHits = 999999; 
    ME11_BunchX = 999999; 
    ME11_strip = 999999;

    for (int i = 0; i < 5; ++i) {
        ME11_location[i] = 999999;
    }

    for (int i = 0; i < 3; ++i) {
        rechit_GP[i] = 999999; 
        rechit_LP[i] = 999999;
    }

    rechit_localphi_rad = 999999; 
    rechit_localphi_deg = 999999;
    has_rechit = false;
    RdPhi = 999999; 
    rechit_detId = 999999;

    for (int i = 0; i < 5; ++i) {
        rechit_location[i] = 999999;
    }

    for (int i = 0; i < 3; ++i) {
        sim_GP[i] = 9999999; 
        sim_LP[i] = 9999999;
    }

    simDy = 9999999; 
    sim_yroll = 9999999; 
    nSim = 9999999;
}

TTree* CSCData::book(TTree* t) {
    edm::Service<TFileService> fs_;
    t = fs_->make<TTree>("Inner_Prop", "Inner_Prop");

    // Book branches for the tree
    // Muon Info
    t->Branch("muon_charge", &muon_charge); 
    t->Branch("muon_pt", &muon_pt);
    t->Branch("muon_eta", &muon_eta); 
    t->Branch("muon_momentum", &muon_momentum);
    t->Branch("evtNum", &evtNum); 
    t->Branch("lumiBlock", &lumiBlock); 
    t->Branch("muonIdx", &muonIdx);
    t->Branch("runNum", &runNum);
    // Propagation Info
    t->Branch("prop_GP", &prop_GP, "prop_GP[3] (x,y,z)/F");
    t->Branch("prop_LP", &prop_LP, "prop_LP[3] (x,y,z)/F");
    t->Branch("prop_startingPoint_GP", &prop_startingPoint_GP, "prop_startingPoint_GP[3] (x,y,z)/F");
    t->Branch("prop_localphi_rad", &prop_localphi_rad);
    t->Branch("prop_localphi_deg", &prop_localphi_deg);
    t->Branch("prop_globalphi_rad", &prop_globalphi_rad);
    t->Branch("has_prop", &has_prop);
    t->Branch("has_fidcut", &has_fidcut);
    t->Branch("prop_location", &prop_location, "prop_location[5] (reg, sta, ring, cha, lay)/I");
    // Track Info
    t->Branch("track_chi2", &track_chi2); 
    t->Branch("track_ndof", &track_ndof);
    t->Branch("n_ME11_segment", &n_ME11_segment); 
    t->Branch("which_track", &which_track);
    t->Branch("hasME11", &hasME11); 
    t->Branch("hasME11RecHit", &hasME11RecHit);
    t->Branch("hasME11A", &hasME11A); 
    t->Branch("hasME11ARecHit", &hasME11ARecHit);
    t->Branch("nCSCSeg", &nCSCSeg); 
    t->Branch("nDTSeg", &nDTSeg);
    t->Branch("nME11RecHits", &nME11RecHits); 
    t->Branch("ME11_BunchX", &ME11_BunchX);
    t->Branch("ME11_strip", &ME11_strip);
    t->Branch("ME11_location", &ME11_location, "ME11_location[5] (end, sta, ring, cha, lay)/I");
    // Rechit Info
    t->Branch("rechit_GP", &rechit_GP, "rechit_GP[3] (x,y,z)/F");
    t->Branch("rechit_LP", &rechit_LP, "rechit_LP[3] (x,y,z)/F");
    t->Branch("rechit_localphi_rad", &rechit_localphi_rad);
    t->Branch("rechit_localphi_deg", &rechit_localphi_deg);
    t->Branch("has_rechit", &has_rechit);
    t->Branch("RdPhi", &RdPhi);
    t->Branch("rechit_detId", &rechit_detId);
    t->Branch("rechit_location", &rechit_location, "rechit_location[5] (reg, sta, ring, cha, lay)/I");
    // Sim info for MC
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
    ~CSCAnalyzer() {};

private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();

    void propagate(const reco::Muon* mu, const edm::Event& event, int index, const Trajectory* trajOfMuon);
    void countCSCSegments(const reco::Muon* mu, CSCData& data);
    void propagateToME11(const reco::Muon* mu, const CSCLayer* layer, bool& hasProp, GlobalPoint& globalPos, CSCData& data, const Trajectory* trajOfMuon);
    void matchME11Rechit(const CSCLayer* layer, LocalPoint localPos, CSCData& data);
    float computeRdPhi(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit>>::const_iterator rechit, float localX, float localY, const GEMEtaPartition* partition);
    bool checkFiducialCut(float localY, float localX, float localPhiDeg);

    // Tokens for input collections
    edm::EDGetTokenT<GEMRecHitCollection> gemRecHitsToken_;
    edm::Handle<GEMRecHitCollection> gemRecHits_;
    edm::EDGetTokenT<vector<PSimHit>> gemSimHitsToken_;
    edm::Handle<vector<PSimHit>> gemSimHits_;
    edm::EDGetTokenT<edm::View<reco::Muon>> muonsToken_;
    edm::Handle<View<reco::Muon>> muons_;
    edm::Handle<TrajTrackAssociationCollection> refTrack_;
    edm::EDGetTokenT<TrajTrackAssociationCollection> refTrackToken_;
    edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsToken_;
    edm::Handle<CSCSegmentCollection> cscSegments_;
    edm::EDGetTokenT<CSCRecHit2DCollection> csc2DRecHitsToken_;
    edm::Handle<CSCRecHit2DCollection> csc2DRecHits_;

    edm::Service<TFileService> fs_;

    // Handles for various geometry and tracking services
    MuonServiceProxy* serviceProxy_;
    edm::ESHandle<Propagator> propagator_;
    edm::ESHandle<TransientTrackBuilder> transientTrackBuilder_;
    ESHandle<GlobalTrackingGeometry> trackingGeometry_;
    edm::ESHandle<GEMGeometry> gemGeometry_;
    edm::ESHandle<CSCGeometry> cscGeometry_;

    bool cscPropagation_; 
    bool trackerPropagation_; 
    bool segmentPropagation_;
    vector<int> propagationList_;
    bool debug_;
    bool refitter_;
    bool isCosmic_;

    CSCData data_;
    TTree* trackerTree_;
    TH2D* nME11ColVsMatches_;

    bool isMC_;
    const CSCSegment* me11Segment_;

    // Tokens for ES data
    const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeometryToken_;
    const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeometryToken_;
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;
    const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> trackingGeometryToken_;
};

CSCAnalyzer::CSCAnalyzer(const edm::ParameterSet& config)
    : gemGeometryToken_(esConsumes()),
      cscGeometryToken_(esConsumes()),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      trackingGeometryToken_(esConsumes()) {

    cout << "Begin analyzer" << endl;
    edm::ParameterSet serviceParameters = config.getParameter<edm::ParameterSet>("ServiceParameters");
    serviceProxy_ = new MuonServiceProxy(serviceParameters, consumesCollector());

    muonsToken_ = consumes<View<reco::Muon>>(config.getParameter<InputTag>("muons"));
    gemRecHitsToken_ = consumes<GEMRecHitCollection>(config.getParameter<edm::InputTag>("gemRecHits"));
    gemSimHitsToken_ = consumes<vector<PSimHit>>(config.getParameter<edm::InputTag>("gemSimHits"));
    cscSegmentsToken_ = consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"));
    csc2DRecHitsToken_ = consumes<CSCRecHit2DCollection>(config.getParameter<edm::InputTag>("csc2DRecHits"));
    refTrackToken_ = consumes<TrajTrackAssociationCollection>(config.getParameter<InputTag>("ref_track"));

    debug_ = config.getParameter<bool>("debug");
    isCosmic_ = config.getParameter<bool>("isCosmic");
    refitter_ = config.getParameter<bool>("refitter");

    trackerTree_ = data_.book(trackerTree_);
}

void CSCAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
    // Retrieve geometry and tracking services
    gemGeometry_ = &setup.getData(gemGeometryToken_);
    cscGeometry_ = &setup.getData(cscGeometryToken_);
    transientTrackBuilder_ = &setup.getData(transientTrackBuilderToken_);
    trackingGeometry_ = &setup.getData(trackingGeometryToken_);

    serviceProxy_->update(setup);
    auto propagator = serviceProxy_->propagator("SteppingHelixPropagatorAny");

    // Determine if the event is MC
    isMC_ = !event.eventAuxiliary().isRealData();
    event.getByToken(gemRecHitsToken_, gemRecHits_);
    if (isMC_) {
        event.getByToken(gemSimHitsToken_, gemSimHits_);
    }
    if (!event.getByToken(muonsToken_, muons_)) return;
    if (muons_->size() == 0) return;

    event.getByToken(cscSegmentsToken_, cscSegments_);
    event.getByToken(csc2DRecHitsToken_, csc2DRecHits_);
    event.getByToken(refTrackToken_, refTrack_);
    ConstTrajTrackPairs refTrackPairs;
    for (auto it = refTrack_->begin(); it != refTrack_->end(); ++it) {
        refTrackPairs.push_back(ConstTrajTrackPair(&(*(*it).key), &(*(*it).val)));
    }

    if (debug_) {
        cout << "New! EvtNumber = " << event.eventAuxiliary().event() 
             << " LumiBlock = " << event.eventAuxiliary().luminosityBlock() 
             << " RunNumber = " << event.run() << endl;
    }

    for (size_t i = 0; i < muons_->size(); ++i) {
        edm::RefToBase<reco::Muon> muRef = muons_->refAt(i);
        const reco::Muon* mu = muRef.get();
        if (!mu->standAloneMuon()) continue;
        if (!(mu->passed(reco::Muon::PFIsoTight))) continue; // Filter muons based on PFIsoTight criteria
        if (debug_) cout << "new standalone" << endl;

        const Trajectory* trajOfMuon = nullptr;
        const Trajectory* trajOfTrack = nullptr;
        const reco::Track* trackOfTrack = nullptr;

        // Match muons with tracks
        for (ConstTrajTrackPairs::const_iterator it = refTrackPairs.begin(); it != refTrackPairs.end(); ++it) {
            trajOfTrack = (*it).first;
            trackOfTrack = (*it).second;
            if (trackOfTrack == mu->track().get()) {
                trajOfMuon = trajOfTrack;
                if (debug_) cout << "mu, i, trajOfMuon: " << mu << i << "," << trajOfMuon << endl;
                propagate(mu, event, i, trajOfMuon);
            }
        }
    }

    if (debug_) {
        // int muonCount = muons_->size();
        // int cscSegmentCount = cscSegments_->size();
        int standaloneMuonCount = 0;
        int me11SegmentCount = 0;
        int uniqueMatchCount = 0;
        int uniqueMe11Count = 0;
        int totalMe11Count = 0;

        vector<int> matchList;
        for (size_t i = 0; i < muons_->size(); ++i) {
            edm::RefToBase<reco::Muon> muRef = muons_->refAt(i);
            const reco::Muon* mu = muRef.get();
            if (mu->isStandAloneMuon()) {
                standaloneMuonCount++;
                int me11SegmentCountTmp = 0;
                auto matches = mu->matches();
                for (auto& match : matches) {
                    if (match.detector() != 2) continue;
                    for (auto& segmentMatch : match.segmentMatches) {
                        auto cscSegRef = segmentMatch.cscSegmentRef;
                        auto cscDetID = cscSegRef->cscDetId();
                        cout << cscDetID.endcap() << cscDetID.station() << cscDetID.ring() << cscDetID.chamber() << endl;
                        if (cscDetID.station() == 1 && (cscDetID.ring() == 1 || cscDetID.ring() == 4)) {
                            me11SegmentCountTmp++;
                            int tmpID = cscDetID.endcap() * 10000 + cscDetID.station() * 1000 + cscDetID.ring() * 100 + cscDetID.chamber();
                            if (find(matchList.begin(), matchList.end(), tmpID) == matchList.end()) {
                                matchList.push_back(tmpID);
                                uniqueMatchCount++;
                            }
                        }
                    }
                }
                totalMe11Count += me11SegmentCountTmp;
            }
        }

        vector<int> collectionList;
        for (const auto& segment : *cscSegments_) {
            auto cscDetID = segment.cscDetId();
            cout << cscDetID.endcap() << cscDetID.station() << cscDetID.ring() << cscDetID.chamber() << endl;
            if (cscDetID.station() == 1 && (cscDetID.ring() == 1 || cscDetID.ring() == 4)) {
                me11SegmentCount++;
                int tmpID = cscDetID.endcap() * 10000 + cscDetID.station() * 1000 + cscDetID.ring() * 100 + cscDetID.chamber();
                if (find(collectionList.begin(), collectionList.end(), tmpID) == collectionList.end()) {
                    collectionList.push_back(tmpID);
                    uniqueMe11Count++;
                }
            }
        }

        nME11ColVsMatches_->Fill(uniqueMatchCount, uniqueMe11Count);
    }
}

float CSCAnalyzer::computeRdPhi(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit>>::const_iterator rechit, float localX, float localY, const GEMEtaPartition* partition) {
    GEMDetId gemid((rechit)->geographicalId());
    const auto& etaPart = gemGeometry_->etaPartition(gemid);
    const auto& etaPartPartition = gemGeometry_->etaPartition(partition->id());
    float deltaYRoll = etaPartPartition->toGlobal(etaPartPartition->centreOfStrip(etaPartPartition->nstrips() / 2)).perp() - etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips() / 2)).perp();
    return cos(stripAngle) * (localX - (rechit)->localPosition().x()) - sin(stripAngle) * (localY + deltaYRoll);
}

void CSCAnalyzer::countCSCSegments(const reco::Muon* mu, CSCData& data) {
    const reco::Track* track = mu->outerTrack().get();
    int cscSegmentCount = 0; 
    int dtSegmentCount = 0; 
    int me11SegmentCount = 0; 
    int me11RecHitCount = 0; 
    float me11BunchX = 99999; 
    int me11Strip = 99999; 
    bool hasMe11A = false;

    if (isCosmic_) {
        cscSegmentCount = mu->numberOfSegments(1, 2) + mu->numberOfSegments(2, 2) + mu->numberOfSegments(3, 2) + mu->numberOfSegments(4, 2);
        dtSegmentCount = mu->numberOfSegments(1, 1) + mu->numberOfSegments(2, 1) + mu->numberOfSegments(3, 1) + mu->numberOfSegments(4, 1);
        auto matches = mu->matches();
        for (const auto& match : matches) {
            if (match.detector() != 2) continue;
            for (const auto& segmentMatch : match.segmentMatches) {
                auto cscSegRef = segmentMatch.cscSegmentRef;
                auto cscDetID = cscSegRef->cscDetId();
                if (cscDetID.station() == 1 && (cscDetID.ring() == 1 || cscDetID.ring() == 4)) {
                    if (cscDetID.ring() == 4) { hasMe11A = true; }
                    me11SegmentCount++;
                    if (debug_) { cout << "isCosmic = True! Getting ME11 Segment" << endl; }
                    me11Segment_ = cscSegRef.get();
                    me11RecHitCount = (cscSegRef.get())->nRecHits();
                    me11BunchX = me11Segment_->time();
                    auto cscDetIDFake = CSCDetId(cscDetID.endcap(), cscDetID.station(), cscDetID.ring(), cscDetID.chamber(), 3);
                    const CSCLayer* me11Layer = cscGeometry_->layer(cscDetIDFake);
                    const CSCLayerGeometry* me11LayerGeo = me11Layer->geometry();
                    me11Strip = me11LayerGeo->nearestStrip(me11Segment_->localPosition());
                    data.ME11_location[0] = cscDetID.endcap(); 
                    data.ME11_location[1] = cscDetID.station(); 
                    data.ME11_location[2] = cscDetID.ring(); 
                    data.ME11_location[3] = cscDetID.chamber(); 
                    data.ME11_location[4] = cscDetID.layer();
                }
            }
        }
    } else {
        for (size_t i = 0; i < track->recHitsSize(); ++i) {
            const TrackingRecHit* recHit = track->recHit(i).get();
            DetId recHitId = recHit->geographicalId();
            uint16_t recHitDetId = recHitId.det();
            if (recHitDetId == DetId::Muon) {
                uint16_t recHitSubDet = recHitId.subdetId();
                if (recHitSubDet == MuonSubdetId::CSC) {
                    if (CSCDetId(recHitId).station() == 1 && CSCDetId(recHitId).ring() == 1 && recHit->dimension() == 4) {
                        me11SegmentCount++;
                        if (debug_) { cout << "isCosmic = False! Getting ME11 Segment" << endl; }
                        RecSegment* recSegment = (RecSegment*)recHit;
                        me11Segment_ = (CSCSegment*)recSegment;
                        me11BunchX = ((CSCRecHit2D*)recHit)->wgroupsBX();
                        auto cscDetIDFake = CSCDetId(CSCDetId(recHitId).endcap(), CSCDetId(recHitId).station(), CSCDetId(recHitId).ring(), CSCDetId(recHitId).chamber(), 3);
                        const CSCLayer* me11Layer = cscGeometry_->layer(cscDetIDFake);
                        const CSCLayerGeometry* me11LayerGeo = me11Layer->geometry();
                        me11Strip = me11LayerGeo->nearestStrip(me11Segment_->localPosition());
                        data.ME11_location[0] = CSCDetId(recHitId).endcap(); 
                        data.ME11_location[1] = CSCDetId(recHitId).station(); 
                        data.ME11_location[2] = CSCDetId(recHitId).ring(); 
                        data.ME11_location[3] = CSCDetId(recHitId).chamber(); 
                        data.ME11_location[4] = CSCDetId(recHitId).layer();
                    }
                    if (CSCDetId(recHitId).station() == 1 && CSCDetId(recHitId).ring() == 1) { me11RecHitCount++; }
                    if (recHit->dimension() == 4) { cscSegmentCount++; }
                }
                if (recHitSubDet == MuonSubdetId::DT) {
                    if (recHit->dimension() > 1) { dtSegmentCount++; }
                }
            }
        }
    }
    data.nCSCSeg = cscSegmentCount; 
    data.nDTSeg = dtSegmentCount;
    data.n_ME11_segment = me11SegmentCount;
    data.nME11RecHits = me11RecHitCount;
    data.ME11_BunchX = me11BunchX;
    data.ME11_strip = me11Strip;
    data.hasME11A = hasMe11A;
    if (data.n_ME11_segment >= 1 && data.n_ME11_segment < 1000) { data.hasME11 = true; }
}

void CSCAnalyzer::propagate(const reco::Muon* mu, const edm::Event& event, int index, const Trajectory* trajOfMuon) {
    const reco::Track* track;
    reco::TransientTrack ttTrack;
    TTree* tree = trackerTree_;

    if (!(mu->track().isNonnull())) { return; }
    if (!(mu->isTrackerMuon())) { return; }

    track = mu->track().get();

    // Print Track details for debugging
    cout << "Track->pt() = " << track->pt() << endl;

    ttTrack = transientTrackBuilder_->build(track);

    if (!ttTrack.isValid()) { cout << "BAD EVENT! NO TRACK" << endl; }
    data_.init();

    // Muon Info
    data_.muon_charge = mu->charge(); 
    data_.muon_pt = mu->pt(); 
    data_.muon_eta = mu->eta(); 
    data_.muon_momentum = mu->momentum().mag2();
    data_.evtNum = event.eventAuxiliary().event(); 
    data_.lumiBlock = event.eventAuxiliary().luminosityBlock(); 
    data_.muonIdx = data_.evtNum * 100 + index;
    data_.runNum = event.run();

    // Track Info
    data_.track_chi2 = track->chi2(); 
    data_.track_ndof = track->ndof();

    for (const auto& layer : cscGeometry_->layers()) {
        if (!(layer->id().station() == 1 && (layer->id().ring() == 1 || layer->id().ring() == 4))) continue;
        if (debug_) cout << "Looping over CSC layers, at detID " << layer->id() << endl;
        GlobalPoint globalPos; 
        bool hasProp = false;
        propagateToME11(mu, layer, hasProp, globalPos, data_, trajOfMuon);
        if (!hasProp) continue;
        if (debug_) cout << "Found a prop!" << endl;
        LocalPoint localPos = layer->toLocal(globalPos);
        matchME11Rechit(layer, localPos, data_);
        tree->Fill();
    }
}

bool CSCAnalyzer::checkFiducialCut(float localY, float localX, float localPhiDeg) {
    // Fidcut same as in TBMA https://github.com/cms-sw/cmssw/blob/deaac86743f80cda845f794c9335ba27a4d50417/Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsFitter.h
    const float fiducialCutAngle = 1.0;
    const float ymin = -80.0;
    const float ymax = 80.0;
    const float xmin = -80.0;
    const float xmax = 80.0;
    const float cutAngle = 5.0 - fiducialCutAngle;
    return (abs(localPhiDeg) < cutAngle) && (ymin < localY && localY < ymax) && (xmin < localX && localX < xmax);
}

void CSCAnalyzer::propagateToME11(const reco::Muon* mu, const CSCLayer* layer, bool& hasProp, GlobalPoint& globalPos, CSCData& data, const Trajectory* trajOfMuon) {
    const reco::Track* Track;
    reco::TransientTrack transientTrack;
    hasProp = false;
    const BoundPlane& boundPlane(layer->surface());
    auto propagator = serviceProxy_->propagator("SteppingHelixPropagatorAny");
    double previousTSOSGlobalPositionR = 0.0;

    TrajectoryStateOnSurface tsos;
    TrajectoryStateOnSurface previousTSOS;
    vector<TrajectoryMeasurement> trajMeasurements = trajOfMuon->measurements();

    TrajectoryStateOnSurface tsosLayer; 
    TrajectoryStateOnSurface tsosSegment;
    GlobalPoint startingPoint;

    GlobalPoint centerOfLayer = layer->centerOfStrip(layer->geometry()->numberOfStrips() / 2.0);


    Track = mu->track().get();
    transientTrack = transientTrackBuilder_->build(Track);
    if (!refitter_) {
        float innerDelta = abs(transientTrack.innermostMeasurementState().globalPosition().z() - centerOfLayer.z());
        float outerDelta = abs(transientTrack.outermostMeasurementState().globalPosition().z() - centerOfLayer.z());
        float usedDelta = 0;

        if (innerDelta < outerDelta) {
            tsosSegment = transientTrack.innermostMeasurementState(); 
            tsosLayer = propagator->propagate(tsosSegment, layer->surface()); 
            usedDelta = innerDelta;
            data.which_track = 0;
        } else {
            tsosSegment = transientTrack.outermostMeasurementState(); 
            tsosLayer = propagator->propagate(tsosSegment, layer->surface()); 
            usedDelta = outerDelta;
            data.which_track = 1;
        }
        if (tsosLayer.isValid()) {
            const LocalPoint localPosLayer = layer->toLocal(tsosLayer.globalPosition());
            const LocalPoint local2DPosLayer(localPosLayer.x(), localPosLayer.y(), 0);
            if (!(tsosLayer.globalPosition().z() * tsosSegment.globalPosition().z() < 0) && boundPlane.bounds().inside(local2DPosLayer) && layer->id().station() == 1 && layer->id().ring() == 1) {
                hasProp = true;
                if (debug_) { cout << "Delta to CSC!!! = " << usedDelta << endl; }
                globalPos = tsosLayer.globalPosition();
                startingPoint = tsosSegment.globalPosition();
            }
        }
    } else {
        for (const auto& measurement : trajMeasurements) {
            TrajectoryStateOnSurface tsos = TrajectoryStateCombiner().combine(measurement.forwardPredictedState(), measurement.backwardPredictedState());
            if (tsos.isValid()) {
                double tsosGlobalPositionR = sqrt(tsos.globalPosition().x() * tsos.globalPosition().x() + tsos.globalPosition().y() * tsos.globalPosition().y());
                if (tsosGlobalPositionR > previousTSOSGlobalPositionR) {
                    previousTSOS = tsos;
                    previousTSOSGlobalPositionR = tsosGlobalPositionR;
                }
            }
        }

        tsosLayer = propagator->propagate(previousTSOS, layer->surface());

        if (tsosLayer.isValid()) {
            const LocalPoint localPosLayer = layer->toLocal(tsosLayer.globalPosition());
            const LocalPoint local2DPosLayer(localPosLayer.x(), localPosLayer.y(), 0);
            if (!(tsosLayer.globalPosition().z() * previousTSOS.globalPosition().z() < 0) && boundPlane.bounds().inside(local2DPosLayer) && layer->id().station() == 1 && layer->id().ring() == 1) {
                hasProp = true;
                globalPos = tsosLayer.globalPosition();
                startingPoint = previousTSOS.globalPosition();
            }
        }
    }

    if (hasProp) {
        LocalPoint localPos = layer->toLocal(globalPos);
        data.prop_GP[0] = globalPos.x(); 
        data.prop_GP[1] = globalPos.y(); 
        data.prop_GP[2] = globalPos.z();
        data.prop_LP[0] = localPos.x(); 
        data.prop_LP[1] = localPos.y(); 
        data.prop_LP[2] = localPos.z();
        data.prop_startingPoint_GP[0] = startingPoint.x(); 
        data.prop_startingPoint_GP[1] = startingPoint.y(); 
        data.prop_startingPoint_GP[2] = startingPoint.z();
        float radius = globalPos.perp();
        LocalPoint localToCenter(localPos.x(), localPos.y() + radius, 0);
        float localPhi = localToCenter.phi();
        if (debug_) {
            cout << "Local nofix = " << localPos << endl;
            cout << "Local point = " << localToCenter << endl;
            cout << "Local phi   = " << ((M_PI / 2.) - localPhi) * (180. / M_PI) << endl;
        }

        data.prop_localphi_rad = (M_PI / 2.) - localPhi;
        data.prop_localphi_deg = ((M_PI / 2.) - localPhi) * (180. / M_PI);
        data.prop_globalphi_rad = globalPos.phi();
        data.has_prop = hasProp;
        data.has_fidcut = checkFiducialCut(localPos.y(), localPos.x(), ((M_PI / 2.) - localPhi) * (180. / M_PI));
        data.prop_location[0] = layer->id().zendcap(); 
        data.prop_location[1] = layer->id().station(); 
        data.prop_location[2] = layer->id().ring(); 
        data.prop_location[3] = layer->id().chamber(); 
        data.prop_location[4] = layer->id().layer();
    }
}

void CSCAnalyzer::matchME11Rechit(const CSCLayer* layer, LocalPoint localPos, CSCData& data) {
    float rechitGlobalPosX; 
    float rechitGlobalPosY; 
    float rechitGlobalPosZ;
    float rechitLocalPosX; 
    float rechitLocalPosY; 
    float rechitLocalPosZ;
    float rechitLocalPhiRad; 
    float rechitLocalPhiDeg;
    bool hasRechit = false;
    float rdPhi = 9999.; 
    int rechitDetId;
    int rechitRegion; 
    int rechitStation; 
    int rechitRing; 
    int rechitChamber; 
    int rechitLayer;

    for (const auto& hit : *csc2DRecHits_) {
        CSCDetId cscId(hit.geographicalId());
        if (!(cscId == layer->id())) continue;
        if (!(layer->id().station() == 1 && layer->id().ring() == 1 && fabs(hit.localPosition().x() - localPos.x()) < 999.0)) continue;

        int strip = layer->geometry()->nearestStrip(hit.localPosition());
        double stripAngle = layer->geometry()->stripAngle(strip) - M_PI / 2.;
        if (abs(rdPhi) < abs(cos(stripAngle) * (localPos.x() - hit.localPosition().x()) + sin(stripAngle) * (localPos.y() - hit.localPosition().y()))) continue;
        rdPhi = cos(stripAngle) * (localPos.x() - hit.localPosition().x()) + sin(stripAngle) * (localPos.y() - hit.localPosition().y());
        rechitGlobalPosX = layer->toGlobal(hit.localPosition()).x(); 
        rechitGlobalPosY = layer->toGlobal(hit.localPosition()).y(); 
        rechitGlobalPosZ = layer->toGlobal(hit.localPosition()).z();
        rechitLocalPosX = hit.localPosition().x(); 
        rechitLocalPosY = hit.localPosition().y(); 
        rechitLocalPosZ = hit.localPosition().z();
        float localPhi = localPos.phi();
        rechitLocalPhiRad = (M_PI / 2.) - localPhi;
        rechitLocalPhiDeg = ((M_PI / 2.) - localPhi) * (180. / M_PI);
        hasRechit = true;
        rechitDetId = cscId.zendcap() * (cscId.station() * 100 + cscId.chamber());
        rechitRegion = cscId.zendcap(); 
        rechitStation = cscId.station(); 
        rechitRing = cscId.ring(); 
        rechitChamber = cscId.chamber(); 
        rechitLayer = cscId.layer();
    }
    if (hasRechit) {
        data.rechit_GP[0] = rechitGlobalPosX; 
        data.rechit_GP[1] = rechitGlobalPosY; 
        data.rechit_GP[2] = rechitGlobalPosZ;
        data.rechit_LP[0] = rechitLocalPosX; 
        data.rechit_LP[1] = rechitLocalPosY; 
        data.rechit_LP[2] = rechitLocalPosZ;
        data.rechit_localphi_rad = rechitLocalPhiRad;
        data.rechit_localphi_deg = rechitLocalPhiDeg;
        data.has_rechit = hasRechit;
        data.RdPhi = rdPhi;
        data.rechit_detId = rechitDetId;
        data.rechit_location[0] = rechitRegion; 
        data.rechit_location[1] = rechitStation; 
        data.rechit_location[2] = rechitRing; 
        data.rechit_location[3] = rechitChamber; 
        data.rechit_location[4] = rechitLayer;
    }
}

void CSCAnalyzer::beginJob() {}

void CSCAnalyzer::endJob() {
    if (debug_) { nME11ColVsMatches_->Write(); }
}

DEFINE_FWK_MODULE(CSCAnalyzer);
