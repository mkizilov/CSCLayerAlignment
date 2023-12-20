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

class CSCAlignmentDBWriter : public edm::one::EDAnalyzer<> {
public:
  CSCAlignmentDBWriter(const edm::ParameterSet&);
  ~CSCAlignmentDBWriter() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  const std::string chamberFile, endcapFile, CSCLayerAlignmentFile, CSCEndcapFile;
  std::string theDTAlignRecordName, theDTErrorRecordName;
  std::string theCSCAlignRecordName, theCSCErrorRecordName;
  std::string theGEMAlignRecordName, theGEMErrorRecordName;
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

CSCAlignmentDBWriter::CSCAlignmentDBWriter(const edm::ParameterSet& p)
  : CSCLayerAlignmentFile(p.getUntrackedParameter<std::string>("CSCLayerAlignmentFile")),
    theDTAlignRecordName("DTAlignmentRcd"),
    theDTErrorRecordName("DTAlignmentErrorExtendedRcd"),
    theCSCAlignRecordName("CSCAlignmentRcd"),
    theCSCErrorRecordName("CSCAlignmentErrorExtendedRcd"),
    theGEMAlignRecordName("GEMAlignmentRcd"),
    theGEMErrorRecordName("GEMAlignmentErrorExtendedRcd"),
    esTokenDT_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))),
    esTokenCSC_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))),
    esTokenGEM_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))) {}

CSCAlignmentDBWriter::~CSCAlignmentDBWriter() {}

void CSCAlignmentDBWriter::analyze(const edm::Event& event, const edm::EventSetup& eventSetup) {
  edm::ESHandle<DTGeometry> theDTGeometry = eventSetup.getHandle(esTokenDT_);
  edm::ESHandle<CSCGeometry> theCSCGeometry = eventSetup.getHandle(esTokenCSC_);
  edm::ESHandle<GEMGeometry> theGEMGeometry = eventSetup.getHandle(esTokenGEM_);
  AlignableMuon* theAlignableMuon = new AlignableMuon(&(*theDTGeometry), &(*theCSCGeometry), &(*theGEMGeometry));
  theAlignableMuon = dynamic_cast<AlignableMuon*>(theAlignableMuon);
  if (!theAlignableMuon)
    throw cms::Exception("TypeMismatch") << "Argument is not an AlignableMuon";
  //Starting CSC alignment writing
  const auto& CSCChambers = theAlignableMuon->CSCLayers();
  int detNum, endcap;
  std::string line, DetNum, dx, dy, dz, dphix, dphiy, dphiz;
  std::ifstream maptype(CSCLayerAlignmentFile);
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
  //End CSC alignment writing
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

DEFINE_FWK_MODULE(CSCAlignmentDBWriter);
