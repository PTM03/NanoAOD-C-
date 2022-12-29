#define nano9Ana_cxx
// The class definition in nano9Ana.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.



#include "nano9Ana.h"
#include <TH2.h>
#include <TStyle.h>
#include <cmath>
#include <math.h>

void nano9Ana::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
}

void nano9Ana::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  //Initialization of the counters:
  nEvtRan        = 0;
  nEvtTotal      = 0;
  //Other custom counters can be initialized here.

   _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
}

void nano9Ana::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
   _HstFile->Write();
  _HstFile->Close();

  //The following lines are displayed on the root prompt.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total good events = "<<nEvtTotal<<endl;

  //The following lines are written on the sum_<process name>.txt file
  ofstream fout(_SumFileName);
  fout<<"Total events ran = "<<nEvtRan<<endl;
  fout<<"Total good events  = "<<nEvtTotal<<endl;
}

void nano9Ana::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file. 
}

Bool_t nano9Ana::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC.SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);

  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  //The following flags throws away some events based on unwanted properties (such as detector problems)
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  
  nEvtRan++;                             //Total number of events containing everything (including the trash events).
  
  if(GoodEvt){
    nEvtTotal++;                         //Total number of events containing goodEvents
                                         //The analysis is done for these good events.


    //Construction of the arrays:
    
    ////////////////////////////////////////////goodMu array//////////////////////////////////////

    int nmu = 0;                         // This counts the number of muons in each event.
    goodMu.clear();                      // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nMuon); i++){
                                         // This loop runs over all the muon candidates. Some of them will pass our selection criteria.
                                         // These will be stored in the goodMu array.
      Lepton temp;                       // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); //the muon mass in GeV is 0.105
      temp.id = -13*Muon_charge[i];      //pdgID for mu- = 13, pdgID for mu+ = -13  
      temp.ind = i;

      //These are the flags the 'temp' object i.e. the muon candidate has to pass.
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
      //passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1 &&  Muon_pfRelIso04_all[i]<0.25;
      
      if(passCuts){
	goodMu.push_back(temp);
        nmu++;                                  // If 'temp' satisfies all the conditions, it is pushed back into goodMu
      }
    }                                    // This 'for' loop has created a goodMu array.
    
    //Now we sort the goodMu in decreasing pT order
    Sort(0);                           

    //Other arrays, such as RecoEle, GenMu, GenEle can be constructed here.
    



    ///////////////////////////////////////////////goodElectron Array///////////////////////////////
    
    //goodElectron array :
    int nEle = 0;                         // This counts the number of electrons in each event.
    goodElectron.clear();                      // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nElectron); i++){
      // This loop runs over all the muon candidates. Some of them will pass our selection criteria.
      // These will be stored in the goodMu array.
      Lepton temp;                       // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.005);//Electron mass in GeV is 0.000000511
      temp.id = -11*Electron_charge[i];  //pdgID for mu- = 11, pdgID for mu+ = -11  
      temp.ind = i;
      
      
      //  h.Electron_cutBased->Fill(Electron_cutBased[i]);
      
      //These are the flags the 'temp' object i.e. the muon candidate has to pass.
      bool passCuts = temp.v.Pt() > 15 && fabs(temp.v.Eta()) < 2.4  && Electron_cutBased[i]>2 ;

      bool Well_Separated = true;

      for(int j=0; j<(int)goodMu.size(); j++){
	float Electron_Muon_Separation = temp.v.DeltaR(goodMu.at(j).v);
	if(Electron_Muon_Separation < 0.4) Well_Separated = false;
      }
      bool isprompt = false;
      if(fabs(temp.v.Eta())<=1.479){
	if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	  isprompt = true;
      }
      if(fabs(temp.v.Eta())>1.479){
	if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
	isprompt = true;
      }
      
      passCuts = passCuts && fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1;
      passCuts = passCuts && Well_Separated && isprompt;
      
      if(passCuts ){
	goodElectron.push_back(temp);
	nEle ++;
      }  
    } // This 'for' loop has created a goodMu array.
    
      //Now we sort the goodMu in decreasing pT order
    Sort(1);                           
    
  /////////////////////////////////////////goodPhoton array//////////////////////////////////////////
    
    
    
    //Good Photon Array
    int nphoton = 0;                   //counts no. of Photons in each events
    goodPhoton.clear();                //emptying the array
    for(unsigned int i=0; i<(*nPhoton);i++){
      
      
      Lepton temp;
      temp.v.SetPtEtaPhiM(Photon_pt[i],Photon_eta[i],Photon_phi[i],0);
      //   temp.id= 22;                //pdg ID dor photon is 22.....since Photon has no Charge keeping the id constant over the entire array
      temp.ind =i;                                                          
      
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4;                  // && Photon_mediumID[i];
      bool Well_Separated_Mu = true;
      for(int j=0;j<(int)goodMu.size();j++){
	float Photon_Muon_Separation = temp.v.DeltaR(goodMu.at(j).v);
	if(Photon_Muon_Separation > 0.4) Well_Separated_Mu =false;
      }
      
      bool Well_Separated_Ele = true;
      for(int j=0;j<(int)goodElectron.size();j++){
	float Photon_Electron_Separation = temp.v.DeltaR(goodElectron.at(j).v);
	if(Photon_Electron_Separation > 0.4) Well_Separated_Ele =false;
      }
      passCuts = passCuts &&  Well_Separated_Mu &&  Well_Separated_Ele;
      if(passCuts){
	goodPhoton.push_back(temp);
	nphoton++;	}               //If template satisfies the conditions in cuts then photon is pushed back into goodPhoton
    }
    
    // sorting the goodPhotons
    Sort(2);  



    ///////////////////////////////////////////////goodJet Array/////////////////////////////////////


    
    int njet = 0;                   //counts no. of Photons in each events
    goodJet.clear();                //emptying the array
    for(unsigned int i=0; i<(*nJet);i++){
      
      
      Lepton temp,temp2;
      temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],0);
      // temp.id= ;                //pdg ID for Jet is .....since Photon has no Charge keeping the id constant over the entire array
      temp.ind =i;                                                          
      
      bool Separated_Jets = true;
      for(unsigned int j=0;j!=i && j<(*nJet); j++){
	temp2.v.SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j],0);
	temp2.ind =j;
	if(temp.v.DeltaR(temp2.v)<2) Separated_Jets=false;
      }
      
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4;//  && Separated_Jets;
      bool Well_Separated_Mu = true;
      bool Well_Separated_Ele = true;
      bool Well_Separated_Photon =true;

      for(int j=0;j<(int)goodMu.size();j++){
	float Jet_Muon_Separation = temp.v.DeltaR(goodMu.at(j).v);
	if(Jet_Muon_Separation > 0.4) Well_Separated_Mu =false;
      }
      for(int j=0;j<(int)goodElectron.size();j++){
	float Jet_Electron_Separation = temp.v.DeltaR(goodElectron.at(j).v);
	if(Jet_Electron_Separation > 0.4) Well_Separated_Ele =false;
	
      }
      for(int j=0;j<(int)goodPhoton.size();j++){
	float Jet_Photon_Separation = temp.v.DeltaR(goodPhoton.at(j).v);
	if(Jet_Photon_Separation > 0.4) Well_Separated_Photon =false;
      }
      passCuts = passCuts &&  Well_Separated_Mu &&  Well_Separated_Ele &&  Well_Separated_Photon;
      if(passCuts){
	goodJet.push_back(temp);
	njet++;
      }
      
      if(passCuts && Jet_btagDeepB[i]>0.4184){
	good_bJet.push_back(temp);
      }                             //If template satisfies the conditions in cuts then jet is pushed back into goodJet
    }
    
    // sorting the goodJets
    
    Sort(3);
    
      ////////////////////////////////GenMuon array//////////////////////////////////////////////////

    
    //GenMu array :
    int ngenPart = 0;       
    GenMu.clear();          
    for(unsigned int i=0; i<(*nGenPart); i++){
      
      
      Lepton temp,newtemp;          
      temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],0.105);
      
      temp.status = GenPart_status[i];
      
      temp.id = GenPart_pdgId[i];
      temp.ind = i;
      temp.momid = GenPart_pdgId[GenPart_genPartIdxMother[i]];
      
      
      // cout<<temp.momid<<endl;
      h.GenPart_pdgID->Fill(temp.id);
      h.GenMu_StatusID->Fill(temp.status);
      
      int motherid = GenMother(i,GenPart_genPartIdxMother[i]);
      temp.momid = motherid;
      temp.goodMother = 0;
      if( fabs(motherid) <100 && motherid != 22)temp.goodMother = 1;
      
      
      // h.Mu->Fill(temp.momid);
      //These are the flags the 'temp' object i.e. the muon candidate has to pass.
      bool passCuts = fabs(temp.id)==13 && temp.status == 1 && temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4 ;
      bool goodMother = motherid < 100 && motherid != 0;
      
      if(passCuts){
	GenMu.push_back(temp);
	ngenPart++;
	
	  // If 'temp' satisfies all the conditions, it is pushed back into goodMu
	}
    } 
    
    Sort(4);
    // Combine GenParts

    

    //// /////////////////////////////////////GenElectrons///////////////////////////////////////////
    //GenElectron
    int ngenPart_e = 0;
    GenElectron.clear();
    for(unsigned int i=0;i<(*nGenPart); i++){
      
      
      
      
      Lepton temp;//
      temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
      
      temp.status =GenPart_status[i];
      
      temp.id =GenPart_pdgId[i];
      temp.ind =i;
      // These are the flags that the 'temp' object i.e. the Gen Electron candidate should pass
      bool passCuts = fabs(temp.id)==11 && temp.status == 1 && temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4;
      
      if(passCuts){
	GenElectron.push_back(temp);
	ngenPart_e ++;
	// If temp satisfies all the conditions then it is pushed back into the goodElectron Array
      }
    }
    
    Sort(5);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
 
    
    //############################
    // Analysis: goodMuon
    //############################
    
    //Fill the size of goodMu array
    h.muprop[0]->Fill((int)goodMu.size());
    //Plot pT of all muons in same plot
    for(int i=0; i<(int)goodMu.size(); i++){
      h.muprop[1]->Fill(goodMu.at(i).v.Pt());
    }

    
    if((int)goodMu.size()>0 && (int)GenMu.size()>0){
      int Min_d_R = 1000;
      for(int i=0; i<(int)GenMu.size(); i++){
	float dR = goodMu.at(0).v.DeltaR(GenMu.at(i).v);
	if(dR< Min_d_R){
	  Min_d_R =dR;	
	  goodMu.at(0).goodMother = GenMu.at(i).goodMother;
	  goodMU.at(0).momid = GenMu.at(i).momid;
	}
      }
      
      
      if(goodMu.at(0).goodMother==1){
	//Plotting the leading muon pT in each event.
	
	h.muprop[2]->Fill(goodMu.at(0).v.Pt());
	h.muprop[3]->Fill(goodMu.at(0).v.Eta());
	h.muprop[4]->Fill(Muon_pfRelIso04_all[goodMu.at(0).ind]);
	
      }
      else{
	h.muprop[9]->Fill(goodMu.at(0).v.Pt());
	h.muprop[10]->Fill(goodMu.at(0).v.Eta());
	h.muprop[11]->Fill(Muon_pfRelIso04_all[goodMu.at(0).ind]);
      }
    }
      
      // This is code for filling Di-Muon Mass
      if((int)goodMu.size()>1){
	float dimuon_mass = (goodMu.at(0).v +goodMu.at(1).v).M();
	h.muprop[5]->Fill(dimuon_mass);
	
	
	//Filling delta Phi
	float delta_phi =fabs( goodMu.at(0).v.Phi() -goodMu.at(1).v.Phi());
	if(delta_phi>TMath::Pi()){ delta_phi = 2*(TMath::Pi())-delta_phi;}
	h.muprop[6]->Fill(delta_phi);
	
	//Filling Delta_Eta
	
	float delta_eta = goodMu.at(0).v.Eta() - goodMu.at(1).v.Eta();
	h.muprop[7]->Fill(delta_eta);
	// Filling Delta_R
	float d_R= sqrt(pow(delta_phi,2) +pow(delta_eta,2));
	h.muprop[8]->Fill(d_R);
	
	
	// FIlling delta_Phi
	float delta_phi_muons = fabs(goodMu.at(0).v.Phi() - goodMu.at(1).v.Phi());
	if(delta_phi_muons > TMath::Pi()){ delta_phi_muons = 2*TMath::Pi() - delta_phi_muons; }
	h.prachu[0] -> Fill(delta_phi_muons);
	
      } 
      //muon analysis ends here
      
      
   
      
    
    
    
    //############################
    // Analysis: goodElectron
    //############################

    //Fill the size of goodMu array
    h.Electronprop[0]->Fill((int)goodElectron.size());
    
    //Plot pT of all muons in same plot
    for(int i=0;i<(int)goodElectron.size();i++){
      h.Electronprop[1]->Fill(goodElectron.at(i).v.Pt());
    }
    //Plotting the leading muon pT in each event.
    if((int)goodElectron.size()>0){
      h.Electronprop[2]->Fill(goodElectron.at(0).v.Pt());
      h.Electronprop[3]->Fill(goodElectron.at(0).v.Eta());
      
    }
    if((int)goodElectron.size()>1){
      float dielectron_mass = (goodElectron.at(0).v +goodElectron.at(1).v).M();
      h.Electronprop[4]->Fill(dielectron_mass);
      
      float d_phi =fabs(goodElectron.at(0).v.Phi() - goodElectron.at(1).v.Phi());
      if(d_phi>TMath::Pi()){d_phi = 2*(TMath::Pi())-d_phi;}
      h.Electronprop[5]->Fill(d_phi);
      
      float d_eta =fabs(goodElectron.at(0).v.Eta() - goodElectron.at(1).v.Eta());
      h.Electronprop[6]->Fill(d_eta);
      
      float delta_electron_R= sqrt(pow(d_eta,2)+pow(d_phi,2));
      h.Electronprop[7]->Fill(delta_electron_R);
      
      float delta_R= goodElectron.at(0).v.DeltaR(goodElectron.at(1).v);
      h.Electronprop[8]->Fill(delta_R);
      
    }
    
    //############################
    // Analysis: goodPhoton
    //############################
    
    //Filling the size of goodPhoton Array
    h.Photonprop[0]->Fill((int)goodPhoton.size());
    
    //plotting pT of all photons in same plot
    for(int i=0; i<(int)goodPhoton.size(); i++){
      h.Photonprop[1]->Fill(goodPhoton.at(i).v.Pt());
    }
    
    //plotting the leading photon pT in each event.
    if((int)goodPhoton.size()>1){
      h.Photonprop[2]->Fill(goodPhoton.at(0).v.Pt());
      h.Photonprop[3]->Fill(goodPhoton.at(1).v.Eta());
      
      // h.photonprop[2]->Fill(photon_pfrelIso04_all[goodPhoton.at(0).ind]);
      
      float diphoton_mass =( goodPhoton.at(0).v + goodPhoton.at(1).v).M();
      h.Photonprop[4]->Fill(diphoton_mass);
      
      float d_phi= fabs(goodPhoton.at(0).v.Phi() - goodPhoton.at(1).v.Phi());
      if(d_phi>TMath::Pi()){ d_phi = 2*(TMath::Pi()) - d_phi;}
      h.Photonprop[5]->Fill(d_phi);
      
      float d_eta = fabs(goodPhoton.at(0).v.Eta()- goodPhoton.at(1).v.Eta());
      float delta_R = goodPhoton.at(0).v.DeltaR(goodPhoton.at(1).v);
      h.Photonprop[6]->Fill(delta_R);
      
    }
    
    //############################
    // Analysis: goodJet
    //############################

    //Filling the size of goodJet array
    h.Jet[0]->Fill((int)goodJet.size());
    
    //plotting the leading JetPT in the same plot
    for(int i=0; i<(int)goodJet.size(); i++){
      h.Jet[1]->Fill((int)goodJet.at(i).v.Pt());
    }
    
    //Plotting leading Jet Pt in each event.
    if((int)goodJet.size()>1){
      h.Jet[2]->Fill(goodJet.at(0).v.Pt());
      h.Jet[3]->Fill(goodJet.at(1).v.Eta());
      // h.jetprop[2]->Fill(jet_pfrelIso04_all[goodJet.at(0).ind]);
      
      float dijet_mass =( goodJet.at(0).v + goodJet.at(1).v).M();
      h.Jet[4]->Fill(dijet_mass);
      
      
      float d_phi= fabs(goodJet.at(0).v.Phi() - goodJet.at(1).v.Phi());
      if(d_phi>TMath::Pi()){ d_phi = 2*(TMath::Pi()) - d_phi;}
      h.Jet[5]->Fill(d_phi);
      
      float d_eta = fabs(goodJet.at(0).v.Eta()- goodJet.at(1).v.Eta());
      float delta_R = goodJet.at(0).v.DeltaR(goodJet.at(1).v);
      h.Jet[6]->Fill(delta_R); 
      
      
      if(goodJet.at(0).v.DeltaR(goodJet.at(1).v)>1){
	h.Jet[7]->Fill(dijet_mass);
	
      }if(goodJet.at(0).v.DeltaR(goodJet.at(1).v)>2){
	h.Jet[8]->Fill(dijet_mass);
	
      }
      
      
    }
    
    if(goodElectron.size()>0 && goodMu.size()>0 && goodPhoton.size()>0 && goodJet.size()>0){
      float d_R_Mu_Ele = goodElectron.at(0).v.DeltaR(goodMu.at(0).v);
      h.Delta_R_Mu_Ele->Fill(d_R_Mu_Ele);
    }
    


   
    //############################
    // Analysis: GenMuon
    //############################
    h.GenMuprop[0]->Fill((int)GenMu.size());
    for(int i=0;i<(int)GenMu.size();i++){
      h.GenMuprop[1]->Fill(GenMu.at(i).v.Pt());
    }
    
    if(GenMu.size()>1){
      h.GenMuprop[2]->Fill(GenMu.at(0).v.Pt());
      h.GenMuprop[3]->Fill(GenMu.at(0).v.Eta());
    }
    //Di Gen_Mu mass
    if((int)GenMu.size()>1){
      float Di_GenMu_Mass = (GenMu.at(0).v+ GenMu.at(1).v).M();
      h.GenMuprop[4]->Fill(Di_GenMu_Mass);
      h.GenMuprop[5]->Fill(GenMu.at(0).momid);
      }
   
    
    //############################
    // Analysis: GenElectron
    //############################

    h.GenElectronprop[0]->Fill((int)GenElectron.size());
    for(int i=0;i<(int)GenElectron.size();i++){
      h.GenElectronprop[1]->Fill(GenElectron.at(i).v.Pt());
    }
    

    if((int)GenElectron.size()>0){
      h.GenElectronprop[2]->Fill(GenElectron.at(0).v.Pt());
    }


 
    
    //Leading Electron PT Resolution:
    
   
   if((int)goodElectron.size()>0 && (int)GenElectron.size()>0){
     int Min_d_R =1000;
     int Min_index_e =-1;  
     for(int i=0; i<(int)GenElectron.size(); i++){
       float d_R =goodElectron.at(0).v.DeltaR(GenElectron.at(i).v);
       
       if(d_R< Min_d_R){
	 Min_d_R = d_R;
	 Min_index_e= i;
       }
     }
     float Rel_err_e = (goodElectron.at(0).v.Pt()-GenElectron.at(Min_index_e).v.Pt())/goodElectron.at(0).v.Pt();
     h.Reso_leading_Pt[2]->Fill(Rel_err_e);
     h.Min_Delta_R_GenElectron_RecoElectron-> Fill(Min_d_R);  
   }
   
   //##################################
      //       Resolution plots
      //##################################
      
      
      
   if((int)goodMu.size()>0 && (int)GenMu.size()>0){
     int Min_index=-1;      
     int Min_d_R=1000;
     for(int i=0; i<(int)GenMu.size(); i++){
       float dR = goodMu.at(0).v.DeltaR(GenMu.at(i).v);
       
       if(dR< Min_d_R){
	 Min_d_R =dR;	
	 Min_index = i ;
       }
     }
     h.Min_Delta_R_GenMu_RecoMu-> Fill(Min_d_R);  
     
     
     float Rel_err = (goodMu.at(0).v.Pt()-GenMu.at(Min_index).v.Pt())/goodMu.at(0).v.Pt();
     if(GenMu.at(Min_index).goodMother == 1){
       h.Reso_leading_Pt[0]->Fill(Rel_err);
       h.Reso_leading_Pt[3]->Fill(GenMu.at(Min_index).momid)  ;
     }
     else{
       h.Reso_leading_Pt[1]->Fill(Rel_err);
       h.Reso_leading_Pt[4]->Fill(GenMu.at(Min_index).momid)	;
     }
     
   }
 
  
   //##################################
   //Some other Plots
   //##################################
   if(goodMu.size()>0 && GenMu.size()>0){
     float Delta_R_goodMu_GenMu =(goodMu.at(0).v. DeltaR(GenMu.at(0).v));
     h.Delta_R_GenMu_RecoMu-> Fill(Delta_R_goodMu_GenMu);
   }
   
   //  if((int)goodMu.size()>0 && (int)GenMu.size()){
   // int
   
   
   
   
   //########### ANALYSIS ENDS HERE ##############
  }//GoodEvt
  
  
  return kTRUE;
}


//######################################
//        USER DEFINED FUNCTIONS
//######################################
void nano9Ana::Sort(int opt)
{
  //This functions sorts an array in the decreasing order of pT.
  //For goodMu:
  if(opt==0){
    for(int i=0; i<(int)goodMu.size()-1; i++){
      for(int j=i+1; j<(int)goodMu.size(); j++){
	if( goodMu[i].v.Pt() < goodMu[j].v.Pt() ) swap(goodMu.at(i),goodMu.at(j));
      }
    }
  }
  //For goodElectron:
  if(opt==1){
    for(int i=0;i<(int)goodElectron.size()-1;i++){
      for(int j=i+1;j<(int)goodElectron.size();j++){
        if(goodElectron[i].v.Pt()<goodElectron[j].v.Pt()) swap(goodElectron.at(i),goodElectron.at(j));      }
    }
  }
  //For good Photon
  if(opt==2){
    for(int i=0; i<(int)goodPhoton.size()-1; i++){
      for(int j=i+1; j<(int)goodPhoton.size(); j++){
	if(goodPhoton[i].v.Pt() < goodPhoton[j].v.Pt() ) swap(goodPhoton.at(i),goodPhoton.at(j));
      }
    }
  }
  //For goodJet
  if(opt==3){
    for(int i=0;i<(int)goodJet.size()-1;i++){
      for(int j=i+1; j<(int)goodJet.size();j++){
	if(goodJet[i].v.Pt() <goodJet[j].v.Pt()) swap(goodJet.at(i),goodJet.at(j));
      }
    }
  }
  //For GenMu
  if(opt==4){
    for(int i=0;i<(int)GenMu.size()-1;i++){
      for(int j=i+1; j<(int)GenMu.size();j++){
	if(GenMu[i].v.Pt() <GenMu[j].v.Pt()) swap(GenMu.at(i),GenMu.at(j));
      }
    }
  }
    //For GenElectron
  if(opt==5){
    for(int i=0;i<(int)GenElectron.size()-1;i++){
      for(int j=i+1; j<(int)GenElectron.size();j++){
	if(GenElectron[i].v.Pt() < GenElectron[j].v.Pt()) swap(GenElectron.at(i),GenElectron.at(j));      }
    }
  }
   //Repeat this for the other arrays here.
}



int nano9Ana::GenMother(int ind, int mom_ind)
{
  int p_id = GenPart_pdgId[ind];
  int m_id = GenPart_pdgId[mom_ind];
  while(p_id==m_id){
    ind = mom_ind;
    mom_ind = GenPart_genPartIdxMother[ind];
    p_id = GenPart_pdgId[ind];
    m_id = GenPart_pdgId[mom_ind];
  }
  return m_id;
}

float nano9Ana::delta_phi(float phi1, float phi2)
{
  //The correct deltaPhi falls in the interval [0 , pi]
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

float nano9Ana::transv_mass(float E_lep, float MET, float dphi)
{
  //The inputs are the Energy of the lepton, MET and dPhi between the lepton and MET
  float mT = sqrt(2* E_lep * MET *(1-cos(dphi)));
  return mT;
}

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
  h.muprop[6]= new TH1F("D_Phi","d_phi",100,-4,4);
  h.muprop[7]= new TH1F("D_Eta","d_eta",100,-4,4);
  h.muprop[8]= new TH1F("D_R","d_R",100,0,5);
  h.muprop[9]=  new TH1F("BM_mu0_pt","BM_Leading muon pT",200,0,200);
  h.muprop[10]= new TH1F("BM_mu0_eta","BM_Leading muon Eta",120,-3.,3.);
  h.muprop[11]= new TH1F("BM_mu0_isol","BM_Leading muon PFRelIso",100,0,2);  


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
  h.GenMu_StatusID = new TH1F("GenMu_status","genmu_status",10,0,10);
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
  
    
