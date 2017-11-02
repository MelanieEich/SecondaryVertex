import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.nuclearInteractionIdentifier_cfi import *
from RecoBTag.SecondaryVertex.vertexAndTracksCleaner_cfi import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from RecoBTag.SecondaryVertex.trackRefitter_cfi import *


#===========================================
# NI rejection and tracks cleaning procedure
#
# short descriptions of the different versions (0 to 7)
# version 0: identify solely based on position
# version 1: identify based on position, mass & ntracks
# version 2: identify based on position & Nctau (=flightDistance2D/(gamma*Bctau) with gamma = pt/mass and Bctau = 0.05 cm
# version 3: identify based on position, mass, ntracks & Nctau
# version 4-7: same as versions 0-3, but first run IVF with relaxed cuts (see below)
#
#===========================================

# IVF run with relaxed cuts to identify more NIs

inclusiveVertexFinderRelaxed = inclusiveCandidateVertexFinder.clone(
	maximumLongitudinalImpactParameter = 25.0,
	vertexMinAngleCosine = 0.
)

vertexMergerRelaxed = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderRelaxed") )

trackVertexArbitratorRelaxed = candidateVertexArbitrator.clone( secondaryVertices = cms.InputTag("vertexMergerRelaxed"))

inclusiveSecondaryVerticesRelaxed = candidateVertexMerger.clone(
	secondaryVertices = cms.InputTag("trackVertexArbitratorRelaxed"),
	maxFraction = cms.double(0.2),
	minSignificance = cms.double(10.)
)

inclusiveCandidateVertexingRelaxed = cms.Sequence(inclusiveVertexFinderRelaxed*vertexMergerRelaxed*trackVertexArbitratorRelaxed*inclusiveSecondaryVerticesRelaxed)

# NI identifiers

Nversion = 2
Nversion2 = 2*Nversion
rhoCuts = ['Rho25', 'Rho9999']

nuclearInteractionIdentifier0 = nuclearInteractionCandIdentifier.clone(
	selection = cms.PSet(
		nuclearInteractionCandIdentifier.selection,
		position = cms.vdouble(2.65, 3.22, 3.52, 5.11, 6.64, 8.01, 9.53, 10.64),
		minNctau = cms.double(2.0)
	)
)

nuclearInteractionIdentifier1 = nuclearInteractionCandIdentifier.clone(
	selection = cms.PSet(
		nuclearInteractionCandIdentifier.selection,
		position = cms.vdouble(2.87, 3.0, 3.65, 3.75, 4.1, 4.25, 4.55, 4.7, 7.05, 7.15, 7.45, 7.6, 9.85, 10., 10.35, 10.45),
		minNctau = cms.double(2.0)
	)
)

#nuclearInteractionIdentifier2 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionCandIdentifier.selection,
		#position = cms.vdouble(2.65, 3.22, 3.52, 5.11, 6.64, 8.01, 9.53, 10.64),
		#minNctau = cms.double(2.0)
	#)
#)

#nuclearInteractionIdentifier3 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionCandIdentifier.selection,
		#position = cms.vdouble(2.87, 3.0, 3.65, 3.75, 4.1, 4.25, 4.55, 4.7, 7.05, 7.15, 7.45, 7.6, 9.85, 10., 10.35, 10.45),
		#minNctau = cms.double(2.0)
	#)
#)

#nuclearInteractionIdentifier4 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionCandIdentifier.selection,
		#position = cms.vdouble(2.65, 3.22, 3.52, 5.11, 6.64, 8.01, 9.53, 10.64),
		#minNctau = cms.double(2.5)
	#)
#)

#nuclearInteractionIdentifier5 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionCandIdentifier.selection,
		#position = cms.vdouble(2.87, 3.0, 3.65, 3.75, 4.1, 4.25, 4.55, 4.7, 7.05, 7.15, 7.45, 7.6, 9.85, 10., 10.35, 10.45),
		#minNctau = cms.double(2.5)
	#)
#)

#nuclearInteractionIdentifier5 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionIdentifier1.selection,
	#)
#)

#nuclearInteractionIdentifier6 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionIdentifier2.selection,
	#)
#)

#nuclearInteractionIdentifier5 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionIdentifier3.selection,
		#distToNI = cms.double(0.3)
	#)
#)

#nuclearInteractionIdentifier6 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionIdentifier1.selection,
		#maxNtracks = cms.int32(4),
		#maxMass = cms.double(2.0)
	#)
#)

#nuclearInteractionIdentifier7 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionIdentifier2.selection,
		#maxNtracks = cms.int32(4),
		#maxMass = cms.double(2.0)
	#)
#)

#nuclearInteractionIdentifier8 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionIdentifier3.selection,
		#maxNtracks = cms.int32(4),
		#maxMass = cms.double(2.0)
	#)
#)

#nuclearInteractionIdentifier9 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionIdentifier4.selection,
		#maxNtracks = cms.int32(4),
		#maxMass = cms.double(2.0)
	#)
#)

#nuclearInteractionIdentifier10 = nuclearInteractionCandIdentifier.clone(
	#selection = cms.PSet(
		#nuclearInteractionIdentifier5.selection,
		#maxNtracks = cms.int32(4),
		#maxMass = cms.double(2.0)
	#)
#)


for i in range(0, Nversion) :
    globals()['vertexAndTracksCleaned%s'%i] = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier%s"%i)

for i in range(Nversion, Nversion2) :
	#globals()['nuclearInteractionIdentifier%s'%i] = globals()['nuclearInteractionIdentifier%s'%(i-Nversion)].clone(secondaryVertices = "inclusiveSecondaryVerticesRelaxed")
	globals()['vertexRefitted%s'%i] = vertexRefitted.clone(secondaryVertices = "nuclearInteractionIdentifier%s"%(i-Nversion))
	globals()['nuclearInteractionIdentifier%s'%i] = globals()['nuclearInteractionIdentifier%s'%(i-Nversion)].clone(secondaryVertices = "vertexRefitted%s"%(i))
	globals()['vertexAndTracksCleaned%s'%i] = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier%s"%i)

# vertex and pfcandidates cleaning steps

for i in range(0, Nversion2) :
	#globals()['vertexAndTracksCleaned%s'%i] = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier%s"%i)
	# re-run IVF
	globals()['inclusiveVertexFinderCleaned%s'%i] = inclusiveCandidateVertexFinder.clone(tracks = "vertexAndTracksCleaned%s"%i)
	globals()['vertexMergerCleaned%s'%i] = candidateVertexMerger.clone(secondaryVertices = "inclusiveVertexFinderCleaned%s"%i)
	globals()['trackVertexArbitratorCleaned%s'%i] = candidateVertexArbitrator.clone(tracks = "vertexAndTracksCleaned%s"%i, secondaryVertices = "vertexMergerCleaned%s"%i)
	globals()['inclusiveSecondaryVerticesCleaned%s'%i] = inclusiveCandidateSecondaryVertices.clone(secondaryVertices = "trackVertexArbitratorCleaned%s"%i)

# all NI rejection sequences

for i in range(0, Nversion) :
	globals()['nuclearInteractionsRemoved%s'%i] = cms.Sequence(
		inclusiveCandidateVertexing *
		globals()['nuclearInteractionIdentifier%s'%i] *
		globals()['vertexAndTracksCleaned%s'%i] *
		globals()['inclusiveVertexFinderCleaned%s'%i] *
		globals()['vertexMergerCleaned%s'%i] *
		globals()['trackVertexArbitratorCleaned%s'%i] *
		globals()['inclusiveSecondaryVerticesCleaned%s'%i]
	)
	
for i in range(Nversion, Nversion2) :
    globals()['nuclearInteractionsRemoved%s'%i] = cms.Sequence(
        inclusiveCandidateVertexing *
        globals()['nuclearInteractionIdentifier%s'%(i-Nversion)] *
        globals()['vertexRefitted%s'%i] *
        globals()['nuclearInteractionIdentifier%s'%i] *
        globals()['vertexAndTracksCleaned%s'%i] *
        globals()['inclusiveVertexFinderCleaned%s'%i] *
        globals()['vertexMergerCleaned%s'%i] *
        globals()['trackVertexArbitratorCleaned%s'%i] *
        globals()['inclusiveSecondaryVerticesCleaned%s'%i]
    )
	
	


#nuclearInteractionIdentifier4 = nuclearInteractionIdentifier0.clone(
	#secondaryVertices = cms.InputTag("inclusiveSecondaryVerticesRelaxed")
#)

#nuclearInteractionIdentifier5 = nuclearInteractionIdentifier1.clone(
	#secondaryVertices = cms.InputTag("inclusiveSecondaryVerticesRelaxed")
#)

#nuclearInteractionIdentifier6 = nuclearInteractionIdentifier2.clone(
	#secondaryVertices = cms.InputTag("inclusiveSecondaryVerticesRelaxed")
#)

#nuclearInteractionIdentifier7 = nuclearInteractionIdentifier3.clone(
	#secondaryVertices = cms.InputTag("inclusiveSecondaryVerticesRelaxed")
#)


#vertexAndTracksCleaned0 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier0")
#vertexAndTracksCleaned1 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier1")
#vertexAndTracksCleaned2 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier2")
#vertexAndTracksCleaned3 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier3")
#vertexAndTracksCleaned4 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier4")
#vertexAndTracksCleaned5 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier5")
#vertexAndTracksCleaned6 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier6")
#vertexAndTracksCleaned7 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier7")

#inclusiveVertexFinderCleaned0 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned0"))
#inclusiveVertexFinderCleaned1 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned1"))
#inclusiveVertexFinderCleaned2 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned2"))
#inclusiveVertexFinderCleaned3 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned3"))
#inclusiveVertexFinderCleaned4 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned4"))
#inclusiveVertexFinderCleaned5 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned5"))
#inclusiveVertexFinderCleaned6 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned6"))
#inclusiveVertexFinderCleaned7 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned7"))
	
#vertexMergerCleaned0 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned0"))
#vertexMergerCleaned1 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned1"))
#vertexMergerCleaned2 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned2"))
#vertexMergerCleaned3 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned3"))
#vertexMergerCleaned4 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned4"))
#vertexMergerCleaned5 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned5"))
#vertexMergerCleaned6 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned6"))
#vertexMergerCleaned7 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned7"))

#trackVertexArbitratorCleaned0 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned0"), secondaryVertices = cms.InputTag("vertexMergerCleaned0"))
#trackVertexArbitratorCleaned1 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned1"), secondaryVertices = cms.InputTag("vertexMergerCleaned1"))
#trackVertexArbitratorCleaned2 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned2"), secondaryVertices = cms.InputTag("vertexMergerCleaned2"))
#trackVertexArbitratorCleaned3 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned3"), secondaryVertices = cms.InputTag("vertexMergerCleaned3"))
#trackVertexArbitratorCleaned4 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned4"), secondaryVertices = cms.InputTag("vertexMergerCleaned4"))
#trackVertexArbitratorCleaned5 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned5"), secondaryVertices = cms.InputTag("vertexMergerCleaned5"))
#trackVertexArbitratorCleaned6 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned6"), secondaryVertices = cms.InputTag("vertexMergerCleaned6"))
#trackVertexArbitratorCleaned7 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned7"), secondaryVertices = cms.InputTag("vertexMergerCleaned7"))

#inclusiveSecondaryVerticesCleaned2 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned2"))
#inclusiveSecondaryVerticesCleaned3 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned3"))
#inclusiveSecondaryVerticesCleaned4 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned4"))
#inclusiveSecondaryVerticesCleaned5 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned5"))
#inclusiveSecondaryVerticesCleaned6 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned6"))
#inclusiveSecondaryVerticesCleaned7 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned7"))


#nuclearInteractionsRemoved0 = cms.Sequence(
	#inclusiveCandidateVertexing *
	#nuclearInteractionIdentifier0 *
	#vertexAndTracksCleaned0 *
	#inclusiveVertexFinderCleaned0 *
	#vertexMergerCleaned0 *
	#trackVertexArbitratorCleaned0 *
	#inclusiveSecondaryVerticesCleaned0
#)

#nuclearInteractionsRemoved1 = cms.Sequence(
	#inclusiveCandidateVertexing *
	#nuclearInteractionIdentifier1 *
	#vertexAndTracksCleaned1 *
	#inclusiveVertexFinderCleaned1 *
	#vertexMergerCleaned1 *
	#trackVertexArbitratorCleaned1 *
	#inclusiveSecondaryVerticesCleaned1
#)

#nuclearInteractionsRemoved2 = cms.Sequence(
	#inclusiveCandidateVertexing *
	#nuclearInteractionIdentifier2 *
	#vertexAndTracksCleaned2 *
	#inclusiveVertexFinderCleaned2 *
	#vertexMergerCleaned2 *
	#trackVertexArbitratorCleaned2 *
	#inclusiveSecondaryVerticesCleaned2
#)

#nuclearInteractionsRemoved3 = cms.Sequence(
	#inclusiveCandidateVertexing *
	#nuclearInteractionIdentifier3 *
	#vertexAndTracksCleaned3 *
	#inclusiveVertexFinderCleaned3 *
	#vertexMergerCleaned3 *
	#trackVertexArbitratorCleaned3 *
	#inclusiveSecondaryVerticesCleaned3
#)

#nuclearInteractionsRemoved4 = cms.Sequence(
	#inclusiveCandidateVertexingRelaxed *
	#nuclearInteractionIdentifier4 *
	#vertexAndTracksCleaned4 *
	#inclusiveVertexFinderCleaned4 *
	#vertexMergerCleaned4 *
	#trackVertexArbitratorCleaned4 *
	#inclusiveSecondaryVerticesCleaned4
#)

#nuclearInteractionsRemoved5 = cms.Sequence(
	#inclusiveCandidateVertexingRelaxed *
	#nuclearInteractionIdentifier5 *
	#vertexAndTracksCleaned5 *
	#inclusiveVertexFinderCleaned5 *
	#vertexMergerCleaned5 *
	#trackVertexArbitratorCleaned5 *
	#inclusiveSecondaryVerticesCleaned5
#)

#nuclearInteractionsRemoved6 = cms.Sequence(
	#inclusiveCandidateVertexingRelaxed *
	#nuclearInteractionIdentifier6 *
	#vertexAndTracksCleaned6 *
	#inclusiveVertexFinderCleaned6 *
	#vertexMergerCleaned6 *
	#trackVertexArbitratorCleaned6 *
	#inclusiveSecondaryVerticesCleaned6
#)

#nuclearInteractionsRemoved7 = cms.Sequence(
	#inclusiveCandidateVertexingRelaxed *
	#nuclearInteractionIdentifier7 *
	#vertexAndTracksCleaned7 *
	#inclusiveVertexFinderCleaned7 *
	#vertexMergerCleaned7 *
	#trackVertexArbitratorCleaned7 *
	#inclusiveSecondaryVerticesCleaned7
#)
