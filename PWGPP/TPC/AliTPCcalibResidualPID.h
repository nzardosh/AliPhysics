/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>             *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// The task:
// stores TPC PID quantities in a THnSparse
//
//  Author:
//  Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>
//

#ifndef ALITPCCALIBRESIDUALPID_H
#define ALITPCCALIBRESIDUALPID_H
#include "AliAnalysisTaskSE.h"

#include <TTreeStream.h>
#include "AliInputEventHandler.h"

class TArrayF;
template <class X>
class THnSparseT;
typedef class THnSparseT<TArrayF> THnSparseF;
class TFile;
class TGraphErrors;
class AliESDEvent;
class AliMCEvent;
class AliESDtrackCuts;
class AliESDpid;
class AliESD;
class AliAnalysisTask;
class AliESDInputHandler;
class AliESDv0KineCuts;
class AliAnalysisManager;
class AliCentrality;
class TTree;
class TSystem;
class TStyle;
class TROOT;
class Riostream;
class TChain;
class TH2;
class TF1;
class TH1;
class TObjArray;


class AliTPCcalibResidualPID : public AliAnalysisTaskSE {
 public:
  enum FitType { kAleph = 0, kLund = 1, kSaturatedLund = 2, kAlephWithAdditionalParam = 3 };
  AliTPCcalibResidualPID();
  AliTPCcalibResidualPID(const char *name);
  virtual ~AliTPCcalibResidualPID();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0);
  virtual void   Terminate(const Option_t*);
  Int_t          CompareFloat(Float_t f1=1, Float_t f2=0) const;
  //setter
  virtual void   SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  virtual void   SetESDtrackCutsV0(AliESDtrackCuts * trackCutsV0){fESDtrackCutsV0 = trackCutsV0;};
  virtual void   SetProduceTPCsignalTHnSparse(Int_t producetpcsignal){fProduceTPCSignalSparse = producetpcsignal;};
  virtual void   SetProducePIDqa(Int_t produceGlobal){fProduceGlobal = produceGlobal;};
  virtual void   SetProduceAllPadsPID(Int_t produceAllpadTypes){fProduceAllPadTypes = produceAllpadTypes;};
  virtual void   SetProduceShortPadsPID(Int_t produceShortpads){fProduceShortPads = produceShortpads;};
  virtual void   SetProduceMediumPadsPID(Int_t produceMediumpads){fProduceMediumPads = produceMediumpads;};
  virtual void   SetProduceLongPadsPID(Int_t produceLongpads){fProduceLongPads = produceLongpads;};
  virtual void   SetProduceOrocPID(Int_t produceOroc){fProduceOroc = produceOroc;};
  
  virtual Bool_t GetVertexIsOk(AliVEvent* event) const;
  
  virtual Bool_t GetUseTPCCutMIGeo() const { return fUseTPCCutMIGeo; };
  virtual void SetUseTPCCutMIGeo(Bool_t newValue) { fUseTPCCutMIGeo = newValue; };
  
  virtual Bool_t GetIsPbpOrpPb() const { return fIsPbpOrpPb; };
  virtual void SetIsPbpOrpPb(Bool_t newValue) { fIsPbpOrpPb = newValue; };
  
  Double_t GetZvtxCutEvent() const { return fZvtxCutEvent; };
  virtual void SetZvtxCutEvent(Double_t newValue) { fZvtxCutEvent = newValue; };
  
  Bool_t GetCorrectdEdxEtaDependence() const { return fCorrectdEdxEtaDependence; };
  virtual void SetCorrectdEdxEtaDependence(Bool_t flag) { fCorrectdEdxEtaDependence = flag; };
  
  Bool_t GetCorrectdEdxMultiplicityDependence() const { return fCorrectdEdxMultiplicityDependence; };
  virtual void SetCorrectdEdxMultiplicityDependence(Bool_t flag) { fCorrectdEdxMultiplicityDependence = flag; };
  
  virtual Char_t GetV0tag(Int_t trackIndex) const;

  Bool_t GetUseMCinfo() const { return fUseMCinfo; };
  virtual void SetUseMCinfo(Bool_t flag) { fUseMCinfo = flag; };
  
  virtual Int_t GetV0motherIndex(Int_t trackIndex) const;
  virtual Int_t GetV0motherPDG(Int_t trackIndex) const;
  
  //
  // static functions for postprocessing
  //
  static Double_t* ExtractResidualPID(THnSparseF * histPidQA,
                                      const Bool_t useV0s = kTRUE,
                                      const Char_t * outFile = "out.root",
                                      const Char_t * type    = "MC",
                                      const Char_t * period  = "LHC10H8",
                                      const Char_t * pass    = "PASS1",
                                      const Char_t * system  = "PBPB",
                                      const Double_t * initialParameters = 0x0,
                                      const Char_t * dedxtype= "",
                                      FitType = kSaturatedLund);
  static  TObjArray * GetResidualGraphs(THnSparseF * histPidQA, const Char_t * system, const Bool_t useV0s);
  static  TObjArray * GetResidualGraphsMC(THnSparseF * histPidQA, const Char_t * system);
  static  TObjArray * GetSeparation(THnSparseF * histPidQA, Int_t kParticle1, Int_t kParticle2);
  static  TObjArray * GetResponseFunctions(TF1* parametrisation, TObjArray* inputGraphs, const Char_t * type, const Char_t * period, const Char_t * pass, const Char_t * system, const Char_t * dedxtype);
  static  TF1*        FitBB(TObjArray* inputGraphs, Bool_t isMC, Bool_t isPPb, const Bool_t useV0s,
                            const Double_t * initialParameters = 0x0, FitType = kSaturatedLund);
  static Int_t MergeGraphErrors(TGraphErrors* mergedGraph, TCollection* li);
  
  static Double_t GetCutGeo() { return fgCutGeo; };
  static Double_t GetCutNcr() { return fgCutNcr; };
  static Double_t GetCutNcl() { return fgCutNcl; };
  
  static void SetCutGeo(Double_t value) { fgCutGeo = value; };
  static void SetCutNcr(Double_t value) { fgCutNcr = value; };
  static void SetCutNcl(Double_t value) { fgCutNcl = value; };
  
  static Bool_t TPCCutMIGeo(const AliVTrack* track, const AliVEvent* evt, TTreeStream* streamer = 0x0);
  static Bool_t TPCCutMIGeo(const AliVTrack* track, const AliInputEventHandler* evtHandler, TTreeStream* streamer = 0x0)
    { if (!evtHandler) return kFALSE; return TPCCutMIGeo(track, evtHandler->GetEvent(), streamer); };

  protected:
  static Double_t fgCutGeo;  // Cut variable for TPCCutMIGeo concerning geometry
  static Double_t fgCutNcr;  // Cut variable for TPCCutMIGeo concerning num crossed rows
  static Double_t fgCutNcl;  // Cut variable for TPCCutMIGeo concerning num clusters
  
  private:
  static Double_t Lund(Double_t* xx, Double_t* par);
  static Double_t SaturatedLund(Double_t* xx, Double_t* par);
  
  void  BinLogAxis(const THnSparseF *h, Int_t axisNumber);
  enum {kElectron=0, kPion, kKaon, kProton} kParticle ;

  static void FitSlicesY(TH2 *hist, Double_t heightFractionForRange, Int_t cutThreshold, TString fitOption, TObjArray *arr);

  void FillV0PIDlist(AliESDEvent* esdEvent = 0x0);
  void ClearV0PIDlist();

  //
  //
  AliESDEvent *fESD;                   //! ESD object
  AliMCEvent  *fMC;                    //! MC object
  TObjArray * fOutputContainer;        //! output data container
  AliESDtrackCuts * fESDtrackCuts;     // basic cut variables for all non-V0 tracks
  AliESDtrackCuts * fESDtrackCutsV0;   // basic cut variables for all V0 tracks
  AliESDpid * fESDpid;                 //! PID handling
  //
  
  Bool_t fUseTPCCutMIGeo;   // Use geometrical cut for TPC 
  
  Bool_t fUseMCinfo;         // Use MC info, if available
  
  Bool_t fIsPbpOrpPb;      // Pbp/pPb collision or something else?
  Double_t fZvtxCutEvent;  // Vertex z cut for the event (cm)
  
  AliESDv0KineCuts *fV0KineCuts;       //! ESD V0 kine cuts
  Int_t fNumTagsStored;     // Number of entries of fV0tags
  Char_t* fV0tags;         //! Pointer to array with tags for identified particles from V0 decays
  Int_t* fV0motherIndex;   //! Pointer to array with index of the mother V0
  Int_t* fV0motherPDG;     //! Pointer to array with pdg of the mother V0

  Bool_t fProduceAllPadTypes, fProduceGlobal, fProduceShortPads, fProduceMediumPads, fProduceLongPads,fProduceOroc;
  THnSparseF * fHistPidQA;             //! histogram for the QA of the PID
  THnSparseF * fHistPidQAshort;        //! histogram for the QA of the PID short pads
  THnSparseF * fHistPidQAmedium;       //! histogram for the QA of the PID med pads
  THnSparseF * fHistPidQAlong;         //! histogram for the QA of the PID long pads
  THnSparseF * fHistPidQAoroc;         //! histogram for the QA of the PID full oroc
  //
  Bool_t fProduceTPCSignalSparse;      //for setter
  Bool_t fCorrectdEdxEtaDependence;    // Correct eta dependence for fHistPidQA (NOTE: Not done for the pad-specific THnSparses)
  Bool_t fCorrectdEdxMultiplicityDependence; // Correct multiplicity dependence for fHistPidQA (NOTE: Not done for the pad-specific THnSparses)
  THnSparseF * fThnspTpc;              //! thnsparse containing the data
  //
  //
  
  // QA histos
  TObjArray* fQAList;           //! Array with QA histos
  TH1F* fhInvMassGamma;         //! Histogram with inv. mass of gamma
  TH1F* fhInvMassK0s;           //! Histogram with inv. mass of K0s
  TH1F* fhInvMassLambda;        //! Histogram with inv. mass of lambda
  TH1F* fhInvMassAntiLambda;    //! Histogram with inv. mass of anti-lambda
  
  TH2F* fhArmenterosAll;        //! Histogram with armenteros plot for all V0s
  TH2F* fhArmenterosGamma;      //! Histogram with armenteros plot for gamma
  TH2F* fhArmenterosK0s;        //! Histogram with armenteros plot for K0s
  TH2F* fhArmenterosLambda;     //! Histogram with armenteros plot for lambda
  TH2F* fhArmenterosAntiLambda; //! Histogram with armenteros plot for anti-lambda
  
  // QA histos for shared clusters
  THnSparseF* fHistSharedClusQAV0Pi;  //! Histogram with shared clusters QA for V0 pi
  THnSparseF* fHistSharedClusQAV0Pr;  //! Histogram with shared clusters QA for V0 pr
  THnSparseF* fHistSharedClusQAV0El;  //! Histogram with shared clusters QA for V0 el
  
  AliTPCcalibResidualPID(const AliTPCcalibResidualPID&); // not implemented
  AliTPCcalibResidualPID& operator=(const AliTPCcalibResidualPID&); // not implemented
  
  ClassDef(AliTPCcalibResidualPID, 4); 
};
#endif
