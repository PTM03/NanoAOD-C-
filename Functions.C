#include<cmath>
#include<math.h>

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

int nano9Ana::GenMatching(Lepton reco, vector <Lepton> gen)
{
  // This is a Snippet for the GenMatrching code.
  int Min_index = -1;
  if((int)gen.size()>0){
    float  Min_d_R =1000;
    
    
    for(int j=0; j<(int)gen.size(); j++){
      float d_R =reco.v.DeltaR(gen.at(j).v);
      
      if(d_R< Min_d_R){
	Min_d_R =d_R;
	Min_index=j;
      }
	//reco.at(i).goodMother = gen.at(Min_index).goodMother;
      //	reco.at(i).momid = gen.at(Min_index).momid;
    }
    
  }
  return Min_index;
}


  
