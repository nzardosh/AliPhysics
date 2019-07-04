#ifndef ALIHFJETFINDER_H
#define ALIHFJETFINDER_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFJetFinder
// \helper class to handle jet finding, matching and reclustering
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODMCParticle.h"
#include "AliHFJet.h"
#include "AliFJWrapper.h"
#include "FJ_includes.h"

class AliHFJetFinder : public TObject
{
  public:
  


  AliHFJetFinder();
  AliHFJetFinder(char *name);

  virtual ~AliHFJetFinder();

  
  AliHFJet GetHFMesonJet(TClonesArray *array, AliAODRecoDecayHF *cand);
  AliHFJet GetHFMesonMCJet(TClonesArray *array, AliAODMCParticle *mcpart);
  std::vector<AliHFJet> GetHFMesonJets(TClonesArray *array, AliAODRecoDecayHF *cand);
  std::vector<AliHFJet> GetHFMesonMCJets(TClonesArray *array, AliAODMCParticle *mcpart);
  std::vector<AliHFJet> GetJets(TClonesArray *array);
  std::vector<AliHFJet> GetMCJets(TClonesArray *array);
  void FindJets(TClonesArray *array, AliAODRecoDecayHF *cand=nullptr);
  void FindMCJets(TClonesArray *array, AliAODMCParticle *mcpart=nullptr);
  void SetJetVariables(AliHFJet& HFJet, std::vector<fastjet::PseudoJet> Constituents, fastjet::PseudoJet Jet, Int_t JetID, AliAODRecoDecayHF *cand=nullptr);
  void SetMCJetVariables(AliHFJet& HFJet, std::vector<fastjet::PseudoJet> Constituents, fastjet::PseudoJet Jet, Int_t JetID, AliAODMCParticle *mcpart=nullptr);
  void SetJetSubstructureVariables(AliHFJet& HFJet, std::vector<fastjet::PseudoJet> Constituents);
  Bool_t CheckTrack(AliAODTrack *track);
  Bool_t CheckParticle(AliAODMCParticle *particle);
  Int_t Find_Candidate_Jet();
  Float_t RelativePhi(Float_t Phi1, Float_t Phi2);

  void SetHF(Bool_t b) {fHF=b;}
  void SetMC(Bool_t b) {fMC=b;}
  void SetDoSubstructure(Bool_t b) {fDoSubstructure=b;}
  void SetRecoDecayHFCand(AliAODRecoDecayHF *RecoDecayHFCand) {fRecoDecayHFCand=RecoDecayHFCand;}
  void SetMCParticle(AliAODMCParticle *MCParticle) {fMCParticle=MCParticle;}
  void SetJetRadius (Float_t f) {fJetRadius = f; fFastJetWrapper->SetR(fJetRadius);}
  

  Float_t                  fJetRadius;
  
  Bool_t                   fHF;
  Bool_t                   fMC;
  Bool_t                   fDoSubstructure;
  AliAODRecoDecayHF       *fRecoDecayHFCand;
  AliAODMCParticle        *fMCParticle;
  AliFJWrapper            *fFastJetWrapper;




  /// \cond CLASSIMP
  ClassDef(AliHFJetFinder,1); ///
  /// \endcond
};
#endif
