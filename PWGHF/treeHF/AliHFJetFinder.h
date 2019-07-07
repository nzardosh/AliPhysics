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
  void SetMinJetPt(Float_t f)              {fMinJetPt = f;}
  void SetJetRadius(Float_t f)             {fJetRadius = f;}
  void SetJetAlgorithm(Int_t i)            {fJetAlgorithm=i;}
  void SetJetRecombScheme(Int_t i)         {fJetRecombScheme = i;}
  void SetGhostArea(Float_t f)             {fJetGhostArea = f;}
  void SetJetAreaType(Int_t i)             {fJetAreaType = i;}
  void SetMinSubJetPt(Float_t f)           {fMinSubJetPt = f;}
  void SetSubJetRadius(Float_t f)          {fSubJetRadius = f;}
  void SetSubJetAlgorithm(Int_t i)         {fSubJetAlgorithm=i;}
  void SetSubJetRecombScheme(Int_t i)      {fSubJetRecombScheme = i;}
  void SetSoftDropZCut(Float_t f)          {fSoftDropZCut = f;}
  void SetSoftDropBeta(Float_t f)          {fSoftDropBeta = f;}
  void SetMinTrackPt(Float_t f)            {fMinTrackPt = f;}
  void SetMaxTrackPt(Float_t f)            {fMaxTrackPt = f;}
  void SetMaxTrackEta(Float_t f)           {fMaxTrackEta = f;}
  void SetMinParticlePt(Float_t f)         {fMinParticlePt = f;}
  void SetMaxParticlePt(Float_t f)         {fMaxParticlePt = f;}
  void SetMaxParticleEta(Float_t f)        {fMaxParticleEta = f;}
  void SetCharged(Int_t i)                 {fCharged = i;}
  

  Float_t                  fMinJetPt;
  Float_t                  fJetRadius;
  Int_t                    fJetAlgorithm;
  Int_t                    fJetRecombScheme;
  Float_t                  fJetGhostArea;
  Int_t                    fJetAreaType;

  Float_t                  fMinSubJetPt;
  Float_t                  fSubJetRadius;
  Int_t                    fSubJetAlgorithm;
  Int_t                    fSubJetRecombScheme;
  Float_t                  fSoftDropZCut;
  Float_t                  fSoftDropBeta;

  Float_t                  fMinTrackPt;
  Float_t                  fMaxTrackPt;
  Float_t                  fMaxTrackEta;

  Float_t                  fMinParticlePt;
  Float_t                  fMaxParticlePt;
  Float_t                  fMaxParticleEta;
  Float_t                  fCharged;
  
  Bool_t                   fDoJetSubstructure;
  AliFJWrapper            *fFastJetWrapper;




  /// \cond CLASSIMP
  ClassDef(AliHFJetFinder,1); ///
  /// \endcond
};
#endif
