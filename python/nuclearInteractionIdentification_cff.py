import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.nuclearInteractionIdentifier_cfi import *
from RecoBTag.SecondaryVertex.vertexAndTracksCleaner_cfi import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *


#===========================================
# NI rejection and tracks cleaning procedure
#
#===========================================




# IVF run with relaxed cuts to identify more NIs

#inclusiveVertexFinderRelaxed = inclusiveCandidateVertexFinder.clone(
	#maximumLongitudinalImpactParameter = 3.0,
	#vertexMinAngleCosine = 0.7
#)

#vertexMergerRelaxed = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderRelaxed") )

#trackVertexArbitratorRelaxed = candidateVertexArbitrator.clone( secondaryVertices = cms.InputTag("vertexMergerRelaxed"))

#inclusiveSecondaryVerticesRelaxed = candidateVertexMerger.clone(
	#secondaryVertices = cms.InputTag("trackVertexArbitratorRelaxed"),
	#maxFraction = cms.double(0.2),
	#minSignificance = cms.double(10.)
#)

#inclusiveCandidateVertexingRelaxed = cms.Sequence(inclusiveVertexFinderRelaxed*vertexMergerRelaxed*trackVertexArbitratorRelaxed*inclusiveSecondaryVerticesRelaxed)



# NI identifiers

nuclearInteractionIdentifier0 = nuclearInteractionCandIdentifier.clone(
	selection = cms.PSet(
		nuclearInteractionCandIdentifier.selection,
		position = cms.vdouble(2.65, 3.22, 3.52, 5.11, 6.64, 8.01, 9.53, 10.64),
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
