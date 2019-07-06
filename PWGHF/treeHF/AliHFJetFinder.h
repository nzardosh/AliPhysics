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

  void SetFJWrapper(); 
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

  void SetDoJetSubstructure(Bool_t b)      {fDoJetSubstructure=b;}
  void SetJetRadius(Float_t f)             {fJetRadius = f;}
  void SetJetAlgorithm(Int_t i)            {fJetAlgorithm=i;}
  void SetJetRecombScheme(Int_t i)         {fJetRecombScheme = i;}
  void SetGhostArea(Float_t f)             {fJetGhostArea = f;}
  void SetJetAreaType(Int_t i)             {fJetAreaType = i;}
  

  Float_t                  fJetRadius;
  Int_t                    fJetAlgorithm;
  Int_t                    fJetRecombScheme;
  Float_t                  fJetGhostArea;
  Int_t                    fJetAreaType;
  
  Bool_t                   fDoJetSubstructure;
  AliFJWrapper            *fFastJetWrapper;




  /// \cond CLASSIMP
  ClassDef(AliHFJetFinder,1); ///
  /// \endcond
};
#endif
