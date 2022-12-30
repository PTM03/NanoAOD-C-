#include <cmath>
#include <math.h>






void nano9Ana::BookHistograms()
{
  //The histograms are booked here.
  //Binning etc are done here.
  //These histograms are stored in the hst_<process name>.root file in the same order.

  //Example : new TH1F ("hst_name", "hst title", NBins, startVal, EndVal);
  
  h.muprop[0] = new TH1F("nmuons", "Number of Muons", 10, 0, 10);
  h.muprop[1] = new TH1F("muonpt", "Muon pT", 200, 0, 200);
  h.muprop[2] = new TH1F("GM_mu0_pt","Leading muon pT",200,0,200);
  h.muprop[3] = new TH1F("GM_mu0_eta","Leading muon Eta",120,-3.,3.);
  h.muprop[4] = new TH1F("GM_mu0_isol","Leading muon PFRelIso",100,0,2);  
  h.muprop[5] = new TH1F("Di_Mu_mass","di_muon_mass",200,0,200);
  h.muprop[6] = new TH1F("D_Phi","d_phi",100,-4,4);
  h.muprop[7] = new TH1F("D_Eta","d_eta",100,-4,4);
  h.muprop[8] = new TH1F("D_R","d_R",100,0,5);
  h.muprop[9] =  new TH1F("BM_mu0_pt","BM_Leading muon pT",200,0,200);
  h.muprop[10] = new TH1F("BM_mu0_eta","BM_Leading muon Eta",120,-3.,3.);
  h.muprop[11] = new TH1F("BM_mu0_isol","BM_Leading muon PFRelIso",100,0,2);  
  h.muprop[12] = new TH1F("recoMu_Mothers", "Reco_Mu_mother", 2000,-1000, 1000);
  h.muprop[13] = new TH1F("recoMu_goodMother","reco_Mu_goodMother", 10,0,10);

  h.Electronprop[0] = new TH1F("nElectrons","Number of Electrons",10,0,10);
  h.Electronprop[1] = new TH1F("Electronpt","Electron pT",200,0,200);
  // h.Electronprop[2] = new TH1F("Electron_CutBased","Electron_cutBased",100,0,10);
  h.Electronprop[2] =new TH1F ("Electron0_pt","Leading Electron pT",200,0,200);
  h.Electronprop[3] = new TH1F("Electron0_eta","Leading Electron eta",120,-3,3.);
  h.Electronprop[4] = new TH1F("Di-electron_Mass","Dielectron_mass",200,0,200);
  h.Electronprop[5] = new TH1F("Delta_phi_electron","delta_phi_electron",100,-4,4);
  h.Electronprop[6] = new TH1F("Delta-eta","delta_eta_electron",100,0,3);
  h.Electronprop[7] = new TH1F("delta_R_Electrons", "delta_R_Electrons",100,0,6);
  h.Electronprop[8] =new TH1F("delta_R_Electron","delta_R_electron",100,0,6);
  // h.Di_Electron_R = new TH1F("Di-electron_R","di-electron_R",0,50,25);
  // h.Electronprop[2] = new TH1F("Electron0_isol","Leading Electron PFRelIso",100,0,2);  
  // h.Electron_cutBased = new TH1F("Electron_CutBased","electron_cutbased",1000,0,1000);



  //Prachu's additions:
  h.prachu[0] = new TH1F("p_dphi", "dPhi(Muons)", 200, 0, 4);
  
  
  
  //Photon plots
  h.Photonprop[0] = new TH1F("nphotons","Number of photons",10,0,10);
  h.Photonprop[1] = new TH1F("PhotonPt","Photon pT",200,0,200);
  h.Photonprop[2] = new TH1F("photon_pt","Leading Photon pT",200,0,200);
  h.Photonprop[3] = new TH1F("photon_Eta","Leading Photon Eta",120,-3.,3.);
  h.Photonprop[4] = new TH1F("Photon_Mass","Photon_Mass",200,0,200);  
  h.Photonprop[5] = new TH1F("Delta_phi_Photon","delta_phi_photon",100,-4,4);
  h.Photonprop[6] = new TH1F("Delta_R_Photons","Delta_R_photons", 100,0,6);
 //h.photonprop[2] =new TH1F("Photon_isol","Leading Photon PFRelIso", 100,0,2);
  
  
  
  //Jet Plots
  h.Jet[0] = new TH1F("njets","njet",10,0,10);
  h.Jet[1] = new TH1F("JetpT","Jet pT",200,0,200);
  h.Jet[2] = new TH1F("leading_Jet_pT","Jet_pt",200,0,200);
  h.Jet[3] = new TH1F("leading_Jet_Eta","jet_eta",120,-3,3);
  h.Jet[4] = new TH1F("DiJet_mass","di-jet_mass",200,0,300);
  h.Jet[5] = new TH1F("delta_phi_jets","delta_phi_Jets",100,-4,4);
  h.Jet[6] = new TH1F("delta_R_jets","delta_R_Jets",100,0,6);
  h.Jet[7] = new TH1F("DiJet_Mass_RCut_deltaR>1","di_Jet_Mass_Rcut_1",200,0,200);
  h.Jet[8] = new TH1F("DiJet_Mass_RCut_deltaR>2","di_Jet_Mass_Rcut_2",200,0,200);
  
  h.Delta_R_Mu_Ele = new TH1F("Delta_R_Mu_ELe","delta_R_Mu_ELe",100,0,6);  
  
  //GenParticles
  
  
   //Gen Muons
  h.GenPart_StatusID = new TH1F("GenMu_status","genmu_status",10,0,10);
  h.GenMuprop[0] = new TH1F ("nGenMu","Number of GenMuons",10,0,10);
  h.GenMuprop[1] = new TH1F("GenMu_pT","GenMupT",200,0,200);
  h.GenMuprop[2] = new TH1F ("leading_GenMu_pt","Leading_GenMuon_pT",200,0,200);
  h.GenMuprop[3] = new TH1F ("Leading_GenMu_Eta", "Leading_GenMuon_Eta",120,-3,3);
  h.GenMuprop[4] = new TH1F("Di_GenMu_Mass", "di_genmu_mass", 200,0,200);
  h.GenMuprop[5] = new TH1F("GenMother_pdgID_Muons", "Gen_Mother_pdgID_muon",2000,-1000,1000);
  
  
  
  
  //GenElectron
  h.GenElectronprop[0] =new TH1F("nGen_Electrons", "Number_GenElectrons" ,10,0,10);
  h.GenElectronprop[1] = new TH1F("GenElectron_pT","GenElectronpT",200,0,200);
  h.GenElectronprop[2] = new TH1F("Leading_GenElectron_Pt","leading_gen_electron_pt",200,0,200);

  //Resolution
  h.Reso_leading_Pt[0] = new TH1F("Resolution_Pt_leading_Muon_goodMother", "Resolution in pT (leading Muon with good Mother)", 400,-2,2);
  h.Reso_leading_Pt[1] = new TH1F("Resolution_Pt_leading_Muon_badMother","Resolution in pT(leadingMuon with bad Mother)", 400,-2,2);
  h.Reso_leading_Pt[2] = new TH1F("Reso_leading_Pt_Electron", "Reso_leading_Pt_Electron", 400,-2,2); 
  h.Reso_leading_Pt[3] = new TH1F("goodMothers","good_mothers",2000,-1000,1000);
  h.Reso_leading_Pt[4] = new TH1F("badMothers","bad_mothers",2000,-1000,1000);
  h.Min_Delta_R_GenMu_RecoMu = new TH1F("Min_Delta_R_GenMu_RecoMu","min_delta_R_GenMu_RecoMu",400,0,5);
  h.Delta_R_GenMu_RecoMu = new TH1F("Delta_R_GenMu_RecoMu","delta_R_GenMu_RecoMu",100,0,10);

  h.Min_Delta_R_GenElectron_RecoElectron =new TH1F("Min_Delta_R_GenElectron_RecoElectron","min_delta_R_GenElectron_RecoElectron",400,0,5);
  
  h.GenPart_pdgID = new TH1F("GenPart_pdgID","gen_part_pdgID",50,-25,25); 
  
}
  
    
