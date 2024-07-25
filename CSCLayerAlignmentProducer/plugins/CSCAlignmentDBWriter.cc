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
#include <map>

class CSCAlignmentDBWriter : public edm::one::EDAnalyzer<> {
public:
    explicit CSCAlignmentDBWriter(const edm::ParameterSet&);
    ~CSCAlignmentDBWriter() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
    // Configuration parameters
    const std::string cscLayerAlignmentFile_;
    
    // Record names for alignment
    const std::string dtAlignRecordName_;
    const std::string dtErrorRecordName_;
    const std::string cscAlignRecordName_;
    const std::string cscErrorRecordName_;
    const std::string gemAlignRecordName_;
    const std::string gemErrorRecordName_;

    // Tokens for geometry
    edm::ESGetToken<DTGeometry, MuonGeometryRecord> esTokenDT_;
    edm::ESGetToken<CSCGeometry, MuonGeometryRecord> esTokenCSC_;
    edm::ESGetToken<GEMGeometry, MuonGeometryRecord> esTokenGEM_;

    // Alignment data structures
    AlignableMuon* alignableMuon_;
    AlignableModifier muonModifier_;
    Alignments* dtAlignments_;
    AlignmentErrorsExtended* dtAlignmentErrorsExtended_;
    Alignments* cscAlignments_;
    AlignmentErrorsExtended* cscAlignmentErrorsExtended_;
    Alignments* gemAlignments_;
    AlignmentErrorsExtended* gemAlignmentErrorsExtended_;
};

CSCAlignmentDBWriter::CSCAlignmentDBWriter(const edm::ParameterSet& p)
    : cscLayerAlignmentFile_(p.getUntrackedParameter<std::string>("CSCLayerAlignmentFile")),
      dtAlignRecordName_("DTAlignmentRcd"),
      dtErrorRecordName_("DTAlignmentErrorExtendedRcd"),
      cscAlignRecordName_("CSCAlignmentRcd"),
      cscErrorRecordName_("CSCAlignmentErrorExtendedRcd"),
      gemAlignRecordName_("GEMAlignmentRcd"),
      gemErrorRecordName_("GEMAlignmentErrorExtendedRcd"),
      esTokenDT_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))),
      esTokenCSC_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))),
      esTokenGEM_(esConsumes(edm::ESInputTag("", "idealForMuonMisalignedProducer"))) {}

CSCAlignmentDBWriter::~CSCAlignmentDBWriter() {
    delete alignableMuon_;
}

void CSCAlignmentDBWriter::analyze(const edm::Event& event, const edm::EventSetup& eventSetup) {
    // Retrieve geometry handles
    edm::ESHandle<DTGeometry> dtGeometry = eventSetup.getHandle(esTokenDT_);
    edm::ESHandle<CSCGeometry> cscGeometry = eventSetup.getHandle(esTokenCSC_);
    edm::ESHandle<GEMGeometry> gemGeometry = eventSetup.getHandle(esTokenGEM_);

    // Create AlignableMuon from geometry
    alignableMuon_ = new AlignableMuon(&(*dtGeometry), &(*cscGeometry), &(*gemGeometry));

    // Validate AlignableMuon creation
    if (!alignableMuon_) {
        throw cms::Exception("TypeMismatch") << "Failed to create AlignableMuon";
    }

    // Read alignment parameters from file
    std::map<CSCDetId, std::vector<float>> alignmentParameters;
    std::ifstream alignmentFile(cscLayerAlignmentFile_);
    std::string line;

    while (std::getline(alignmentFile, line)) {
        std::stringstream ss(line);
        std::string detNumStr, dxStr, dyStr, dzStr, dphixStr, dphiyStr, dphizStr;
        std::getline(ss, detNumStr, ',');
        std::getline(ss, dxStr, ',');
        std::getline(ss, dyStr, ',');
        std::getline(ss, dzStr, ',');
        std::getline(ss, dphixStr, ',');
        std::getline(ss, dphiyStr, ',');
        std::getline(ss, dphizStr, ',');

        int detNum = std::stoi(detNumStr);
        float xShift = std::stof(dxStr);
        float yShift = std::stof(dyStr);
        float zShift = std::stof(dzStr);
        float rotX = std::stof(dphixStr);
        float rotY = std::stof(dphiyStr);
        float rotZ = std::stof(dphizStr);

        int endcap = (detNum > 0) ? 1 : 2;
        CSCDetId id(endcap, 1, 1, (std::abs(detNum % 1000) - std::abs(detNum % 10)) / 10, std::abs(detNum % 10));
        CSCDetId id2(endcap, 1, 4, (std::abs(detNum % 1000) - std::abs(detNum % 10)) / 10, std::abs(detNum % 10));
        std::vector<float> params = {xShift, yShift, zShift, rotX, rotY, rotZ};

        alignmentParameters[id.rawId()] = params;
        alignmentParameters[id2.rawId()] = params;
    }

    // Apply alignment parameters
    const auto& cscChambers = alignableMuon_->CSCLayers();
    for (const auto& chamber : cscChambers) {
        auto cscId = chamber->id();
        if (alignmentParameters.count(cscId) == 0) {
            continue;
        }

        const auto& params = alignmentParameters[cscId];
        muonModifier_.moveAlignableLocal(chamber, false, false, params[0], params[1], params[2]);
        muonModifier_.rotateAlignableLocal(chamber, false, false, params[3], params[4], params[5]);
    }

    // Retrieve alignment and alignment errors
    dtAlignments_ = alignableMuon_->dtAlignments();
    dtAlignmentErrorsExtended_ = alignableMuon_->dtAlignmentErrorsExtended();
    cscAlignments_ = alignableMuon_->cscAlignments();
    cscAlignmentErrorsExtended_ = alignableMuon_->cscAlignmentErrorsExtended();
    gemAlignments_ = alignableMuon_->gemAlignments();
    gemAlignmentErrorsExtended_ = alignableMuon_->gemAlignmentErrorsExtended();

    // Write alignments to DB
    edm::Service<cond::service::PoolDBOutputService> poolDbService;
    if (!poolDbService.isAvailable()) {
        throw cms::Exception("NotAvailable") << "PoolDBOutputService not available";
    }

    poolDbService->writeOneIOV<Alignments>(*dtAlignments_, poolDbService->beginOfTime(), dtAlignRecordName_);
    poolDbService->writeOneIOV<AlignmentErrorsExtended>(*dtAlignmentErrorsExtended_, poolDbService->beginOfTime(), dtErrorRecordName_);
    poolDbService->writeOneIOV<Alignments>(*cscAlignments_, poolDbService->beginOfTime(), cscAlignRecordName_);
    poolDbService->writeOneIOV<AlignmentErrorsExtended>(*cscAlignmentErrorsExtended_, poolDbService->beginOfTime(), cscErrorRecordName_);
    poolDbService->writeOneIOV<Alignments>(*gemAlignments_, poolDbService->beginOfTime(), gemAlignRecordName_);
    poolDbService->writeOneIOV<AlignmentErrorsExtended>(*gemAlignmentErrorsExtended_, poolDbService->beginOfTime(), gemErrorRecordName_);
}

DEFINE_FWK_MODULE(CSCAlignmentDBWriter);
