/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFJetFinder
// \helper class to handle jet finding, matching and reclustering
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include "AliHFJetFinder.h"
#include "AliAODRecoDecayHF.h"
#include "TMath.h"

/// \cond CLASSIMP
ClassImp(AliHFJetFinder);
/// \endcond

//________________________________________________________________
AliHFJetFinder::AliHFJetFinder():
  fJetRadius(0.4),
  fHF(false),
  fMC(false),
  fRecoDecayHFCand(0x0),
  fMCParticle(0x0),
  fFastJetWrapper(0x0)
{
  //
  // Default constructor
  fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
  fFastJetWrapper->SetAreaType(fastjet::active_area); 
  fFastJetWrapper->SetGhostArea(0.005);  
  fFastJetWrapper->SetR(fJetRadius);
  fFastJetWrapper->SetAlgorithm(fastjet::antikt_algorithm);
  fFastJetWrapper->SetRecombScheme(static_cast<fastjet::RecombinationScheme>(0));
  //
}

//________________________________________________________________
AliHFJetFinder::AliHFJetFinder(char *name):
  fJetRadius(0.4),
  fHF(false),
  fMC(false),
  fRecoDecayHFCand(0x0),
  fMCParticle(0x0),
  fFastJetWrapper(0x0)
{
  fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
  fFastJetWrapper->SetAreaType(fastjet::active_area); 
  fFastJetWrapper->SetGhostArea(0.005);  
  fFastJetWrapper->SetR(fJetRadius);
  fFastJetWrapper->SetAlgorithm(fastjet::antikt_algorithm);
  fFastJetWrapper->SetRecombScheme(static_cast<fastjet::RecombinationScheme>(0));
}

//________________________________________________________________
AliHFJetFinder::~AliHFJetFinder()
{
  //
  // Destructor
  //
  delete fFastJetWrapper;
}



//________________________________________________________________
AliHFJet AliHFJetFinder::GetHFMesonJet(TClonesArray *array, AliAODRecoDecayHF *cand){
  //Jet is clustered with heavy flavour meson and the corresponding variables are set
  AliHFJet HFJet;
  if (!cand) return HFJet;
  fFastJetWrapper->Clear();
  FindJets(array,cand);
  Int_t Jet_Index=Find_Candidate_Jet();
  if (Jet_Index==-1) return HFJet;

 
  std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets();
  fastjet::PseudoJet Jet = Inclusive_Jets[Jet_Index];
  std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(Jet_Index));

  SetJetVariables(HFJet, Constituents, Jet, 0, cand); 

    return HFJet;
}



//________________________________________________________________
AliHFJet AliHFJetFinder::GetHFMesonMCJet(TClonesArray *array, AliAODMCParticle *mcpart){
  //Jet is clustered with heavy flavour meson and the corresponding variables are set
  AliHFJet HFJet;
  if (!mcpart) return HFJet;
  fFastJetWrapper->Clear();
  FindMCJets(array,mcpart);
  Int_t Jet_Index=Find_Candidate_Jet();
  if (Jet_Index==-1) return HFJet;

 
  std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets();
  fastjet::PseudoJet Jet = Inclusive_Jets[Jet_Index];
  std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(Jet_Index));

  SetMCJetVariables(HFJet, Constituents, Jet, 0, mcpart);
 
    return HFJet;
}


//________________________________________________________________
std::vector<AliHFJet> AliHFJetFinder::GetHFMesonJets(TClonesArray *array, AliAODRecoDecayHF *cand) {
  //Jet is clustered with heavy flavour meson and the corresponding variables are set
  std::vector<AliHFJet> HFJets;
  HFJets.clear();
  if (!cand) return HFJets;
  fFastJetWrapper->Clear();
  FindJets(array, cand);
  Int_t Jet_Index=Find_Candidate_Jet();
  if (Jet_Index==-1) return HFJets;

  std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets();
  
  for (Int_t i=0; i<Inclusive_Jets.size(); i++){
    fastjet::PseudoJet Jet = Inclusive_Jets[i];
    std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(i));

    AliHFJet HFJet;

    if (i==Jet_Index) SetJetVariables(HFJet, Constituents, Jet, i, cand); 
    else SetJetVariables(HFJet, Constituents, Jet, i, nullptr); 
    HFJets.push_back(HFJet);
  }
    return HFJets;
}


//________________________________________________________________
std::vector<AliHFJet> AliHFJetFinder::GetHFMesonMCJets(TClonesArray *array, AliAODMCParticle *mcpart) {
  //Jet is clustered with heavy flavour meson and the corresponding variables are set
  std::vector<AliHFJet> HFJets;
  HFJets.clear();
  if (!mcpart) return HFJets;
  fFastJetWrapper->Clear();
  FindMCJets(array, mcpart);
  Int_t Jet_Index=Find_Candidate_Jet();
  if (Jet_Index==-1) return HFJets;

  std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets();
  
  for (Int_t i=0; i<Inclusive_Jets.size(); i++){
    fastjet::PseudoJet Jet = Inclusive_Jets[i];
    std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(i));

    AliHFJet HFJet;
  
    if (i==Jet_Index) SetMCJetVariables(HFJet, Constituents, Jet, i, mcpart); 
    else SetMCJetVariables(HFJet, Constituents, Jet, i, nullptr);
    
    HFJets.push_back(HFJet);
  }
    return HFJets;
}



//________________________________________________________________
std::vector<AliHFJet> AliHFJetFinder::GetJets(TClonesArray *array) {
  //Jet is clustered with heavy flavour meson and the corresponding variables are set
  std::vector<AliHFJet> HFJets;
  HFJets.clear();
  fFastJetWrapper->Clear();
  FindJets(array,nullptr);

  std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets();
  for (Int_t i=0; i<Inclusive_Jets.size(); i++){
    fastjet::PseudoJet Jet = Inclusive_Jets[i];
    std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(i));

    AliHFJet HFJet;
    if (fDoSubstructure) SetJetSubstructureVariables(HFJet,Constituents);

    SetJetVariables(HFJet, Constituents, Jet, i, nullptr);
  
    HFJets.push_back(HFJet);
  }
    return HFJets;
}


//________________________________________________________________
std::vector<AliHFJet> AliHFJetFinder::GetMCJets(TClonesArray *array) {
  //Jet is clustered with heavy flavour meson and the corresponding variables are set
  std::vector<AliHFJet> HFJets;
  HFJets.clear();
  fFastJetWrapper->Clear();
  FindMCJets(array,nullptr);

  std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets();
  for (Int_t i=0; i<Inclusive_Jets.size(); i++){
    fastjet::PseudoJet Jet = Inclusive_Jets[i];
    std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(i));

    AliHFJet HFJet;
    if (fDoSubstructure) SetJetSubstructureVariables(HFJet,Constituents);

    SetMCJetVariables(HFJet, Constituents, Jet, i, nullptr);
  
    HFJets.push_back(HFJet);
  }
    return HFJets;
}



//________________________________________________________________
void AliHFJetFinder::FindJets(TClonesArray *array, AliAODRecoDecayHF *cand) {
  //Performs jet finding. Jets are stored in the fFastJetWrapper object

  std::vector<Int_t> daughters;
  if (cand){
  
    daughters.clear();

    AliVTrack *daughter;
    for (Int_t i = 0; i < cand->GetNDaughters(); i++) {
      daughter = dynamic_cast<AliVTrack *>(cand->GetDaughter(i));   
      if (!daughter) continue;
      daughters.push_back(daughter->GetID());
    }
    fFastJetWrapper->AddInputVector(cand->Px(), cand->Py(), cand->Pz(), cand->E(cand->PdgCode()),0); 
  }

    
  bool IsDaughter;
  AliAODTrack *track=NULL;
  for (Int_t i=0; i<array->GetEntriesFast(); i++) {
 
    track= dynamic_cast<AliAODTrack*>(array->At(i));
    if(!CheckTrack(track)) continue; 
    IsDaughter=false;
    if (cand){
      for (Int_t j=0; j<daughters.size(); j++){
	if (track->GetID()==daughters[j]) IsDaughter=true;
      }
      if(IsDaughter) continue;
    }
    fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(),i+100); 
  }
  fFastJetWrapper->Run();
  //delete track;
}






//________________________________________________________________
void AliHFJetFinder::FindMCJets(TClonesArray *array, AliAODMCParticle *mcpart) {
  //Performs jet finding. Jets are stored in the fFastJetWrapper object

  std::vector<Int_t> daughters;
  if (mcpart){
  
    daughters.clear();

    AliAODMCParticle *daughter;
    for (Int_t i = 0; i < mcpart->GetNDaughters(); i++) {
      daughter = dynamic_cast<AliAODMCParticle *>(array->At(mcpart->GetDaughterLabel(i)));   
      if (!daughter) continue;
      daughters.push_back(daughter->GetLabel());
    }
    fFastJetWrapper->AddInputVector(mcpart->Px(), mcpart->Py(), mcpart->Pz(), mcpart->E(),0); 
  }

    
  bool IsDaughter;
  AliAODMCParticle *particle=NULL;
  for (Int_t i=0; i<array->GetEntriesFast(); i++) {
 
    particle= dynamic_cast<AliAODMCParticle*>(array->At(i));
    if(!CheckParticle(particle)) continue; 
    IsDaughter=false;
    if (mcpart){
      for (Int_t j=0; j<daughters.size(); j++){
	if (particle->GetLabel()==daughters[j]) IsDaughter=true;
      }
      if(IsDaughter) continue;
    }
    fFastJetWrapper->AddInputVector(particle->Px(), particle->Py(), particle->Pz(), particle->E(),i+100); 
  }
  fFastJetWrapper->Run();
  //delete track;
}


//________________________________________________________________
void AliHFJetFinder::SetJetVariables(AliHFJet& HFJet, std::vector<fastjet::PseudoJet> Constituents, fastjet::PseudoJet Jet, Int_t JetID, AliAODRecoDecayHF *cand) {

  HFJet.fID=JetID;
  if (cand)HFJet.fHFMeson=1;
  else HFJet.fHFMeson=0;
  HFJet.fPt=Jet.perp();
  HFJet.fEta=Jet.pseudorapidity();
  HFJet.fPhi=Jet.phi();
  if (cand) HFJet.fDeltaEta=Jet.pseudorapidity()-cand->Eta();
  else HFJet.fDeltaEta=-99;
  if (cand) HFJet.fDeltaPhi=Jet.phi()-cand->Phi();
  else HFJet.fDeltaPhi=-99;
  if (cand) HFJet.fDeltaR=TMath::Sqrt(HFJet.fDeltaEta*HFJet.fDeltaEta + HFJet.fDeltaPhi*HFJet.fDeltaPhi);
  else HFJet.fDeltaR=-99;
  HFJet.fN=Constituents.size();

  if (fDoSubstructure) SetJetSubstructureVariables(HFJet,Constituents);
  
}


//________________________________________________________________
void AliHFJetFinder::SetMCJetVariables(AliHFJet& HFJet, std::vector<fastjet::PseudoJet> Constituents, fastjet::PseudoJet Jet, Int_t JetID, AliAODMCParticle *mcpart) {

  HFJet.fID=JetID;
  if (mcpart)HFJet.fHFMeson=1;
  else HFJet.fHFMeson=0;
  HFJet.fPt=Jet.perp();
  HFJet.fEta=Jet.pseudorapidity();
  HFJet.fPhi=Jet.phi();
  if (mcpart) HFJet.fDeltaEta=Jet.pseudorapidity()-mcpart->Eta();
  else HFJet.fDeltaEta=-99;
  if (mcpart) HFJet.fDeltaPhi=Jet.phi()-mcpart->Phi();
  else HFJet.fDeltaPhi=-99;
  if (mcpart) HFJet.fDeltaR=TMath::Sqrt(HFJet.fDeltaEta*HFJet.fDeltaEta + HFJet.fDeltaPhi*HFJet.fDeltaPhi);
  else HFJet.fDeltaR=-99;
  HFJet.fN=Constituents.size();

  if (fDoSubstructure) SetJetSubstructureVariables(HFJet,Constituents);
  
}




//________________________________________________________________
void AliHFJetFinder::SetJetSubstructureVariables(AliHFJet& HFJet, std::vector<fastjet::PseudoJet> Constituents) {


  Bool_t SoftDropSet=kFALSE;
  Float_t Zg=0;

fastjet::JetDefinition SubJet_Definition(fastjet::cambridge_algorithm, fJetRadius*2.5,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best); 
  try{
    fastjet::ClusterSequence Cluster_Sequence(Constituents, SubJet_Definition);
    std::vector<fastjet::PseudoJet> Reclustered_Jet =  Cluster_Sequence.inclusive_jets(0.0);
    Reclustered_Jet = sorted_by_pt(Reclustered_Jet);
         
    fastjet::PseudoJet Daughter_Jet = Reclustered_Jet[0];
    fastjet::PseudoJet Parent_SubJet_1; 
    fastjet::PseudoJet Parent_SubJet_2;  
	  
    while(Daughter_Jet.has_parents(Parent_SubJet_1,Parent_SubJet_2)){
      if(Parent_SubJet_1.perp() < Parent_SubJet_2.perp()) std::swap(Parent_SubJet_1,Parent_SubJet_2);
      Zg=Parent_SubJet_2.perp()/(Parent_SubJet_1.perp()+Parent_SubJet_2.perp());

      if (Zg >= 0.1 && !SoftDropSet){ //Access to the values can be implemneted as setters once the code structure is finalised
	HFJet.fZg = Zg;
	HFJet.fRg = Parent_SubJet_1.delta_R(Parent_SubJet_2);
	SoftDropSet=kTRUE;
      }
      Daughter_Jet=Parent_SubJet_1;
    }

         
  } catch (fastjet::Error) { /*return -1;*/ }


  
}


//________________________________________________________________
Bool_t AliHFJetFinder::CheckTrack(AliAODTrack *track) { //add all cuts properly later
  if(!track) return false;
  if(track->Pt() > 100.0) return false;
  if(track->Pt() < 1e-6) return false;
  if(TMath::Abs(track->Eta()) > 0.9) return false;
  return true;
}

//________________________________________________________________
Bool_t AliHFJetFinder::CheckParticle(AliAODMCParticle *particle) {
  if(!particle) return false;
  if(!particle->IsPrimary()) return false;
  if(particle->Pt() < 1e-6) return false;
  if (particle->Charge()==0) return false;
  return true;
}



//________________________________________________________________
Int_t AliHFJetFinder::Find_Candidate_Jet() {
  //Finds the label of the jet with the heavy flvaour meson
  std::vector<fastjet::PseudoJet> Inclusive_Jets = fFastJetWrapper->GetInclusiveJets(); 
  for (UInt_t i_Jet=0; i_Jet < Inclusive_Jets.size(); i_Jet++){
    std::vector<fastjet::PseudoJet> Constituents(fFastJetWrapper->GetJetConstituents(i_Jet));
    for (UInt_t i_Constituents = 0; i_Constituents < Constituents.size(); i_Constituents++) { 
      if (Constituents[i_Constituents].user_index() == 0) {
	return i_Jet; 
      }
    }
  }
  return -1;
}



//________________________________________________________________________
Float_t AliHFJetFinder::RelativePhi(Float_t Phi1, Float_t Phi2){

  if(Phi1 < -1*TMath::Pi()) Phi1 += (2*TMath::Pi()); // Turns the range of 0to2Pi into -PitoPi ???????????                                                             
  else if (Phi1 > TMath::Pi()) Phi1 -= (2*TMath::Pi());
  if(Phi2 < -1*TMath::Pi()) Phi2 += (2*TMath::Pi());
  else if (Phi2 > TMath::Pi()) Phi2 -= (2*TMath::Pi());
  Double_t DeltaPhi=Phi2-Phi1;
  if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
  else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  return DeltaPhi;
}
