import FWCore.ParameterSet.Config as cms


vertexRefitted = cms.EDProducer("TrackReFitter",
      secondaryVertices = cms.InputTag("inclusiveCandidateSecondaryVertices"),
      primaryVertices  = cms.InputTag("offlinePrimaryVertices"),
      beamSpot = cms.InputTag("offlineBeamSpot")
      
)