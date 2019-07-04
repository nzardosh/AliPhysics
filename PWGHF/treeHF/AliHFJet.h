#ifndef ALIHFJET_H
#define ALIHFJET_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFJet
// \helper class to handle jets
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include "TObject.h"
class AliHFJet : public TObject
{
  public:
  


  AliHFJet();
  AliHFJet(const AliHFJet &source);
  virtual ~AliHFJet();
  void Reset();

  Float_t GetID() {return fID;}
  Float_t HFMeson() {return fHFMeson;}
  Float_t GetPt() {return fPt;}
  Float_t GetEta() {return fEta;}
  Float_t GetPhi() {return fPhi;}
  Float_t GetDeltaEta() {return fDeltaEta;}
  Float_t GetDeltaPhi() {return fDeltaPhi;}
  Float_t GetDeltaR() {return fDeltaR;}
  Float_t GetN() {return fN;}
  Float_t GetZg() {return fZg;}
  Float_t GetRg() {return fRg;}

    

  Float_t fID;
  Float_t fHFMeson;
  Float_t fPt;
  Float_t fEta;
  Float_t fPhi;
  Float_t fDeltaEta;
  Float_t fDeltaPhi;
  Float_t fDeltaR;
  Float_t fN;
  Float_t fZg;
  Float_t fRg;



  /// \cond CLASSIMP
  ClassDef(AliHFJet,1); ///
  /// \endcond
};
#endif
