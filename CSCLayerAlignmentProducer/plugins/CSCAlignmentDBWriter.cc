#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "Alignment/MuonAlignment/interface/AlignableMuon.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Geometry/CommonTopologies/interface/GeometryAligner.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Alignment/MuonAlignment/interface/AlignableMuon.h"
#include "Alignment/CommonAlignment/interface/AlignableModifier.h"

#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>

class GEMAlDBWriter : public edm::one::EDAnalyzer<> {
public:
  GEMAlDBWriter(const edm::ParameterSet&);
  ~GEMAlDBWriter() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  const std::string chamberFile, endcapFile, ME11ChamberFile, CSCEndcapFile;
  std::string theDTAlignRecordName, theDTErrorRecordName;
  std::string theCSCAlignRecordName, theCSCErrorRecordName;
  std::string theGEMAlignRecordName, theGEMErrorRecordName;
  const bool doChamber, doEndcap, doME11Chamber, doCSCEndcap;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> esTokenDT_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> esTokenCSC_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> esTokenGEM_;
  AlignableMuon* theAlignableMuon;
  AlignableModifier theMuonModifier;
  Alignments* dt_Alignments;
  AlignmentErrorsExtended* dt_AlignmentErrorsExtended;
  Alignments* csc_Alignments;
  AlignmentErrorsExtended* csc_AlignmentErrorsExtended;
  Alignments* gem_Alignments;
  AlignmentErrorsExtended* gem_AlignmentErrorsExtended;
};

GEMAlDBWriter::GEMAlDBWriter(const edm::ParameterSet& p)
  : chamberFile(p.getUntrackedParameter<std::string>("chamberFile")),
    endcapFile(p.getUntrackedParameter<std::string>("endcapFile")),
    ME11ChamberFile(p.getUntrackedParameter<std::string>("ME11ChamberFile")),
    CSCEndcapFile(p.getUntrackedParameter<std::string>("CSCEndcapFile")),
    theDTAlignRecordName("DTAlignmentRcd"),
    theDTErrorRecordName("DTAlignmentErrorExtendedRcd"),
    theCSCAlignRecordName("CSCAlignmentRcd"),
    theCSCErrorRecordName("CSCAlignmentErrorExtendedRcd"),
    theGEMAlignRecordName("GEMAlignmentRcd"),
    theGEMErrorRecordName("GEMAlignmentErrorExtendedRcd"),
    doChamber(p.getUntrackedParameter<bool>("doChamber")),
    doEndcap(p.getUntrackedParameter<bool>("doEndcap")),
    doME11Chamber(p.getUntrackedParameter<bool>("doME11Chamber")),
    doCSCEndcap(p.getUntrackedParameter<bool>("doCSCEndcap")),
    esTokenDT_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))),
    esTokenCSC_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))),
    esTokenGEM_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))) {}

GEMAlDBWriter::~GEMAlDBWriter() {}

void GEMAlDBWriter::analyze(const edm::Event& event, const edm::EventSetup& eventSetup) {
  edm::ESHandle<DTGeometry> theDTGeometry = eventSetup.getHandle(esTokenDT_);
  edm::ESHandle<CSCGeometry> theCSCGeometry = eventSetup.getHandle(esTokenCSC_);
  edm::ESHandle<GEMGeometry> theGEMGeometry = eventSetup.getHandle(esTokenGEM_);
  AlignableMuon* theAlignableMuon = new AlignableMuon(&(*theDTGeometry), &(*theCSCGeometry), &(*theGEMGeometry));
  theAlignableMuon = dynamic_cast<AlignableMuon*>(theAlignableMuon);
  if (!theAlignableMuon)
    throw cms::Exception("TypeMismatch") << "Argument is not an AlignableMuon";
  if (doChamber) {
    const auto& GEMChambers = theAlignableMuon->GEMChambers();
    int detNum, endcap, station;
    std::string line, DetNum, dx, dy, dz, dphix, dphiy, dphiz;
    std::ifstream maptype(chamberFile);
    std::map<GEMDetId, std::vector<float>> alPar;
    while(std::getline(maptype, line)){
      std::cout << line << std::endl;
      std::stringstream ssline(line);
      getline(ssline, DetNum, ',');
      getline(ssline, dx, ',');
      getline(ssline, dy, ',');
      getline(ssline, dz, ',');
      getline(ssline, dphix, ',');
      getline(ssline, dphiy, ',');
      getline(ssline, dphiz, ',');
      detNum = (float)atof(DetNum.c_str());
      float xShift = (float)atof(dx.c_str());
      float yShift = (float)atof(dy.c_str());
      float zShift = (float)atof(dz.c_str());
      float rotX = (float)atof(dphix.c_str());
      float rotY = (float)atof(dphiy.c_str());
      float rotZ = (float)atof(dphiz.c_str());
      endcap = (detNum > 0) ? 1 : -1;
      station = (abs(detNum)/1000)%10;
      std::cout << "endcap is " << endcap << ", station is " << station << std::endl;
      //GEMDetId(int region, int ring, int station, int layer, int chamber, int ieta)
      GEMDetId id = GEMDetId(endcap, 1, station, abs(detNum%10), abs((detNum/10)%100), 0);
      std::vector<float> tmp = {xShift, yShift, zShift, rotX, rotY, rotZ};
      alPar[id.rawId()] = tmp;
      std::cout << "detNum:" << detNum << " rawId: " << id.rawId()<<  " xShift:" << xShift << " yShift:" << yShift << " zShift:" << zShift << " rotX:" << rotX << " rotY:" << rotY << " rotZ:" << rotZ << std::endl;
    }
    for (const auto& chamber : GEMChambers) {
      auto gemId = chamber->id();
      const GEMChamber* gemChamber = theGEMGeometry->chamber(chamber->geomDetId());
      std::cout << "Doing chamber " << gemChamber->id().region() << gemChamber->id().ring() << gemChamber->id().station() << gemChamber->id().chamber() << gemChamber->id().layer() << std::endl;
      //if ((gemChamber->id()).station() == 2){
      //  std::cout << "Station 2 GEM, skip for now, FIX THIS LATER IF YOU USE DATA" << std::endl;
      //  continue;
      //}
      //std::cout << "Testing! new gemId = " << gemId << std::endl;
      if (alPar.count(gemId) < 1){
        std::cout << "Skipping detId " << GEMDetId(gemId) << std::endl;
        continue;
        //throw cms::Exception("NotAvailable") << "can't find detId " << GEMDetId(gemId) ;
      }
      auto par = alPar[gemId];
      std::cout << gemId << ": "<< par.at(0) << ", " << par.at(1) << ", " << par.at(2) << ", " << par.at(3) << ", " << par.at(4) << ", " << par.at(5) << std::endl;
      theMuonModifier.moveAlignableLocal(chamber, false, false, par.at(0), par.at(1), par.at(2));
      theMuonModifier.rotateAlignableLocal(chamber, false, false, par.at(3), par.at(4), par.at(5));
    }
  }
  if (doEndcap){
    const auto& GEMendcaps = theAlignableMuon->GEMEndcaps();
    int r, detNum;
    std::string line, DetNum, dx, dy;
    std::ifstream maptype(endcapFile);
    while(std::getline(maptype, line)){
      std::cout << line << std::endl;
      std::stringstream ssline(line);
      getline(ssline, DetNum, ',');
      getline(ssline, dx, ',');
      getline(ssline, dy, ',');
      detNum = (float)atof(DetNum.c_str());
      float xShift = (float)atof(dx.c_str());
      float yShift = (float)atof(dy.c_str());
      r = (detNum < 0) ? 0 : 1;
      theMuonModifier.moveAlignable(GEMendcaps[r], false, false, xShift, yShift, 0.0);
    }
  }
  if (doME11Chamber) {
    const auto& CSCChambers = theAlignableMuon->CSCLayers();
    int detNum, endcap;
    std::string line, DetNum, dx, dy, dz, dphix, dphiy, dphiz;
    std::ifstream maptype(ME11ChamberFile);
    std::map<CSCDetId, std::vector<float>> alPar;
    while(std::getline(maptype, line)){
      std::cout << line << std::endl;
      std::stringstream ssline(line);
      getline(ssline, DetNum, ',');
      getline(ssline, dx, ',');
      getline(ssline, dy, ',');
      getline(ssline, dz, ',');
      getline(ssline, dphix, ',');
      getline(ssline, dphiy, ',');
      getline(ssline, dphiz, ',');
      detNum = (float)atof(DetNum.c_str());
      float xShift = (float)atof(dx.c_str());
      float yShift = (float)atof(dy.c_str());
      float zShift = (float)atof(dz.c_str());
      float rotX = (float)atof(dphix.c_str());
      float rotY = (float)atof(dphiy.c_str());
      float rotZ = (float)atof(dphiz.c_str());
      endcap = (detNum > 0) ? 1 : 2;
      std::cout << "endcap: " << endcap << " chamber: " << abs(detNum%1000)-abs(detNum%10) <<" layer: "<< abs(detNum%10) << std::endl;
      CSCDetId id = CSCDetId(endcap, 1, 1, (abs(detNum%1000)-abs(detNum%10))/10, abs(detNum%10));
      CSCDetId id2 = CSCDetId(endcap, 1, 4, (abs(detNum%1000)-abs(detNum%10))/10, abs(detNum%10));
      std::vector<float> tmp = {xShift, yShift, zShift, rotX, rotY, rotZ};
      alPar[id.rawId()] = tmp;
      alPar[id2.rawId()] = tmp;
      std::cout << "detNum:" << detNum << " rawId: " << id.rawId()<<  " xShift:" << xShift << " yShift:" << yShift << " zShift:" << zShift << " rotX:" << rotX << " rotY:" << rotY << " rotZ:" << rotZ << std::endl;
    }
    for (const auto& chamber : CSCChambers) {
      auto cscId = chamber->id();
      const CSCLayer* cscChamber = theCSCGeometry->layer(chamber->geomDetId());
      std::cout << "Starting chamber " << cscChamber->id().endcap() << cscChamber->id().station() << cscChamber->id().ring() << cscChamber->id().chamber() << cscChamber->id().layer() << std::endl;
      if (alPar.count(cscId) < 1)
      {
        continue;
      }
      auto par = alPar[cscId];
      std::cout << cscId << ": "<< par.at(0) << ", " << par.at(1) << ", " << par.at(2) << ", " << par.at(3) << ", " << par.at(4) << ", " << par.at(5) << std::endl;
      theMuonModifier.moveAlignableLocal(chamber, false, false, par.at(0), par.at(1), par.at(2));
      theMuonModifier.rotateAlignableLocal(chamber, false, false, par.at(3), par.at(4), par.at(5));
    }
  }
  if (doCSCEndcap){
    const auto& CSCendcaps = theAlignableMuon->CSCEndcaps();
    int r, detNum;
    std::string line, DetNum, dx, dy;
    std::ifstream maptype(CSCEndcapFile);
    while(std::getline(maptype, line)){
      std::cout << line << std::endl;
      std::stringstream ssline(line);
      getline(ssline, DetNum, ',');
      getline(ssline, dx, ',');
      getline(ssline, dy, ',');
      detNum = (float)atof(DetNum.c_str());
      float xShift = (float)atof(dx.c_str());
      float yShift = (float)atof(dy.c_str());
      r = (detNum < 0) ? 1 : 0;
      theMuonModifier.moveAlignable(CSCendcaps[r], false, false, xShift, yShift, 0.0);
    }
  }
  dt_Alignments = theAlignableMuon->dtAlignments();
  dt_AlignmentErrorsExtended = theAlignableMuon->dtAlignmentErrorsExtended();
  csc_Alignments = theAlignableMuon->cscAlignments();
  csc_AlignmentErrorsExtended = theAlignableMuon->cscAlignmentErrorsExtended();
  gem_Alignments = theAlignableMuon->gemAlignments();
  gem_AlignmentErrorsExtended = theAlignableMuon->gemAlignmentErrorsExtended();
  edm::Service<cond::service::PoolDBOutputService> poolDbService;
  if (!poolDbService.isAvailable())  // Die if not available
    throw cms::Exception("NotAvailable") << "PoolDBOutputService not available";
  poolDbService->writeOneIOV<Alignments>((*dt_Alignments), poolDbService->beginOfTime(), theDTAlignRecordName);
  poolDbService->writeOneIOV<AlignmentErrorsExtended>(
      (*dt_AlignmentErrorsExtended), poolDbService->beginOfTime(), theDTErrorRecordName);
  poolDbService->writeOneIOV<Alignments>((*csc_Alignments), poolDbService->beginOfTime(), theCSCAlignRecordName);
  poolDbService->writeOneIOV<AlignmentErrorsExtended>(
      (*csc_AlignmentErrorsExtended), poolDbService->beginOfTime(), theCSCErrorRecordName);
  poolDbService->writeOneIOV<Alignments>((*gem_Alignments), poolDbService->beginOfTime(), theGEMAlignRecordName);
  poolDbService->writeOneIOV<AlignmentErrorsExtended>(
      (*gem_AlignmentErrorsExtended), poolDbService->beginOfTime(), theGEMErrorRecordName);
}

DEFINE_FWK_MODULE(GEMAlDBWriter);
