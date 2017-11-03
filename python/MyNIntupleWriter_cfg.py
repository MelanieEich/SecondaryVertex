import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.nuclearInteractionIdentifier_cfi import *
from RecoBTag.SecondaryVertex.vertexAndTracksCleaner_cfi import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from RecoBTag.SecondaryVertex.trackRefitter_cfi import *


#===========================================
# NI rejection and tracks cleaning procedure
#
#===========================================


processName = "RECODEBUG"
process = cms.Process(processName)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source ("PoolSource",
    fileNames=cms.untracked.vstring(
        'file:////nfs/dust/cms/user/eichm/btag/data/bgd/84B70053-DA98-E711-8525-00259029E87C.root'
    )
)



# NI identifiers

nuclearInteractionIdentifier0 = nuclearInteractionCandIdentifier.clone(
	selection = cms.PSet(
		nuclearInteractionCandIdentifier.selection,
#		position = cms.vdouble(2.65, 3.22, 3.52, 5.11, 6.64, 8.01, 9.53, 10.64),
		position = cms.vdouble(2.17, 2.25, 2.781, 3.014, 6.616, 6.898, 10.729, 11.016, 15.836, 16.125),
		minNctau = cms.double(2.0)
	)
)

# track re-fitting


vertexRefitted0 = vertexRefitted.clone(secondaryVertices = "nuclearInteractionIdentifier0")

nuclearInteractionIdentifierAfterRefit = nuclearInteractionIdentifier0.clone(secondaryVertices = "vertexRefitted0")

# vertex and pfcandidates cleaning steps

vertexAndTracksCleaned0 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifierAfterRefit")

# re-run IVF

inclusiveVertexFinderCleaned0 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned0"))

vertexMergerCleaned0 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned0"))

trackVertexArbitratorCleaned0 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned0"), secondaryVertices = cms.InputTag("vertexMergerCleaned0"))

inclusiveSecondaryVerticesCleaned0 = candidateVertexMerger.clone(
	secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned0"),
	maxFraction = cms.double(0.2),
	minSignificance = cms.double(10.)
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/test.root")
)

process.nuclearInteractionIdentifier0 = cms.EDProducer("NuclearInteractionCandidateIdentifier",
    primaryVertices  = cms.InputTag("offlinePrimaryVertices"),
    secondaryVertices = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    selection = cms.PSet(
        nuclearInteractionCandIdentifier.selection,
#       maxZ = cms.double(100.)				# maximum Z of SV to be identified as NI
        position = cms.vdouble(2.65, 3.22, 3.52, 5.11, 6.64, 8.01, 9.53, 10.64),
        minNctau = cms.double(2.0)
    )
)


process.p = cms.Path(process.nuclearInteractionIdentifier0)
process.e = cms.EndPath(process.out)

# all NI rejection sequences

nuclearInteractionsRemoved0 = cms.Sequence(
	inclusiveCandidateVertexing *
	nuclearInteractionIdentifier0 *
	vertexRefitted0 *
	nuclearInteractionIdentifierAfterRefit *
	vertexAndTracksCleaned0 *
	inclusiveVertexFinderCleaned0 *
	vertexMergerCleaned0 *
	trackVertexArbitratorCleaned0 *
	inclusiveSecondaryVerticesCleaned0
)
