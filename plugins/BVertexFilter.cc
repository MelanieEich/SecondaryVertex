// -*- C++ -*-
//
// Package:    BVertexFilter
// Class:      BVertexFilter
//
/**\class BVertexFilter BVertexFilter.cc DPGAnalysis/BVertexFilter/src/BVertexFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andrea RIZZI
//         Created:  Mon Dec  7 18:02:10 CET 2009
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoBTag/SecondaryVertex/interface/VertexFilter.h"

//
// class declaration
//

template<typename VTX>
class BVertexFilterT : public edm::stream::EDFilter<> {
   public:
      explicit BVertexFilterT(const edm::ParameterSet&);
      ~BVertexFilterT();

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      edm::EDGetTokenT<reco::VertexCollection> token_primaryVertex;
      edm::EDGetTokenT<edm::View<VTX> > token_secondaryVertex;
      reco::VertexFilter                      svFilter;
      bool                                    useVertexKinematicAsJetAxis;
      int                                     minVertices;
};


// BVertexFilter::BVertexFilter(const edm::ParameterSet& params):
//       svFilter(params.getParameter<edm::ParameterSet>("vertexFilter")),
//       useVertexKinematicAsJetAxis(params.getParameter<bool>("useVertexKinematicAsJetAxis")),
//       minVertices(params.getParameter<int>("minVertices"))

// {
//       token_primaryVertex = consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("primaryVertices"));
//       token_secondaryVertex = consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("secondaryVertices"));
//       produces<reco::VertexCollection>();

// }

template<typename VTX>
BVertexFilterT<VTX>::BVertexFilterT(const edm::ParameterSet& params):
      svFilter(params.getParameter<edm::ParameterSet>("vertexFilter")),
      useVertexKinematicAsJetAxis(params.getParameter<bool>("useVertexKinematicAsJetAxis")),
      minVertices(params.getParameter<int>("minVertices"))

{
      token_primaryVertex = consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("primaryVertices"));
      token_secondaryVertex = consumes<edm::View<VTX> >(params.getParameter<edm::InputTag>("secondaryVertices"));
      produces<std::vector<VTX> >();

}

template<typename VTX>
BVertexFilterT<VTX>::~BVertexFilterT()
{
}

template<typename VTX>
bool
BVertexFilterT<VTX>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 int count = 0;
 edm::Handle<reco::VertexCollection> pvHandle;
 iEvent.getByToken(token_primaryVertex, pvHandle);
 edm::Handle<edm::View<VTX> > svHandle;
 iEvent.getByToken(token_secondaryVertex, svHandle);

 // std::auto_ptr<reco::VertexCollection> recoVertices(new reco::VertexCollection);
 auto recoVertices = std::make_unique<std::vector<VTX>>();

 if(pvHandle->size()!=0) {
   const reco::Vertex & primary = (*pvHandle.product())[0];
   const edm::View<VTX>  & vertices = *svHandle.product();


   if(! primary.isFake())
   {
     for(typename edm::View<VTX>::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
      {
            GlobalVector axis(0,0,0);
            if(useVertexKinematicAsJetAxis) axis = GlobalVector(it->p4().X(),it->p4().Y(),it->p4().Z());
            if(svFilter(primary,reco::TemplatedSecondaryVertex<reco::Vertex>(primary,*it,axis,true),axis))  {
                  count++;
                  recoVertices->push_back(*it);
             }
     }
   }
 }
 iEvent.put(std::move(recoVertices));

 return(count >= minVertices);
}

typedef BVertexFilterT<reco::Vertex>                      BVertexFilter;

//define this as a plug-in
DEFINE_FWK_MODULE(BVertexFilter);
