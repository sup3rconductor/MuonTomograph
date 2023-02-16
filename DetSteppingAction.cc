
#include "DetSteppingAction.hh"
#include "DetDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include <stdio.h> 

extern FILE *st, *stEl, *stMu, *stPr, *stPi;
extern G4int Nphot[8][128], Nfe[8][128], FlMu[8], MatrMu[8][128], MatrAll[8][128];
extern G4int NcopyToSlStr[1024][2], SlStrToNcopy[8][128], Nevent;
extern G4double KoordMu[8][3];

G4double Energy, KvEff, eDeposit;
G4int nEn, NslSA, NstSA, copyNb;

G4String vname, Name;
G4StepPoint* prePoint;
G4ParticleDefinition* particleType;
G4Track * theTrack;

G4double KEF[43] = {0.0019, 0.00523, 0.01028, 0.01784, 0.02788, 0.03755,
					0.04786, 0.06235, 0.07419, 0.09055, 0.10446, 0.11532,
					0.13034, 0.14509, 0.15766, 0.16796, 0.18146, 0.18797,
					0.19339, 0.198,   0.20167, 0.20376, 0.20487, 0.20542,
					0.20542, 0.20542, 0.20542, 0.20542, 0.20542, 0.20542,
					0.20482, 0.20411, 0.20306, 0.20198, 0.19963, 0.19729,
					0.19388, 0.19039, 0.18594, 0.1811,  0.17615, 0.17102,
					0.16589};
G4double  Ef[43] = {1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 
	                2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85,
					2.9, 2.95, 3.0, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35,
					3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85,
					3.9, 3.95, 4.0};
char FElMuName[500];

DetSteppingAction::DetSteppingAction()
{ 
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetSteppingAction::~DetSteppingAction()
{ 
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetSteppingAction::UserSteppingAction(const G4Step* step)
{  
  
  theTrack = step->GetTrack();
  Name = theTrack->GetDefinition()->GetParticleName(); 
  particleType = theTrack->GetDefinition();
  prePoint  = step->GetPreStepPoint();
  vname = prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
  
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition())
  {
	if(vname(0,5)=="FKFEU")
	{
	  copyNb = prePoint->GetTouchableHandle()->GetCopyNumber(); 	  

	  NslSA = NcopyToSlStr[copyNb][0];
	  NstSA = NcopyToSlStr[copyNb][1];
		
	  Nphot[NslSA][NstSA]++;

	  eDeposit = theTrack->GetKineticEnergy();	   
	  Energy = eDeposit/eV;
	   
	  if((Energy-1.9)>1e-20 && (4.0-Energy)>1e-20)
	  {
	    nEn = int((Energy-1.9)/0.05);
	    if(nEn>-1&& nEn<42)
	    {
		  KvEff = ((KEF[nEn+1]-KEF[nEn])*Energy+Ef[nEn+1]*KEF[nEn]-Ef[nEn]*KEF[nEn+1])/0.05;
		  if(KvEff>=(1.0*G4UniformRand())) 
		  {
            Nfe[NslSA][NstSA]++;
		  }
	    }
	  } 	  
	  theTrack->SetTrackStatus(fStopAndKill); 
	}
  }
  if(Name(0,2) == "mu")
  {	
	 if(vname(0,5)=="Strip")
	 {
	   copyNb = prePoint->GetTouchableHandle()->GetCopyNumber(); 	  

	   NslSA = NcopyToSlStr[copyNb][0];
	   NstSA = NcopyToSlStr[copyNb][1];

	   MatrMu[NslSA][NstSA] = 1;

	   FlMu[NslSA] = 1;

	   KoordMu[NslSA][0] = prePoint->GetPosition().x();
	   KoordMu[NslSA][1] = prePoint->GetPosition().y();
	   KoordMu[NslSA][2] = prePoint->GetPosition().z();
	 }
	 if(vname(0,5)=="Strip" || vname(0,5)=="PenSc" || vname(0,5)=="AirSc" || vname(0,5)=="PokSt" || vname(0,5)=="SkScM" || vname(0,5)=="SlScM")
	 {	   
	   if(stMu==NULL)
	   {
	     sprintf(FElMuName, "Ev%07dKMu.dat", Nevent);
         stMu = fopen(FElMuName, "w");  
         if(stMu)
           fprintf(stMu, "Nevent\tFlMu\tx\ty\tz\tenergy\n");   
	   }

	   fprintf(stMu,"%d\t%s\t%g\t%g\t%g\t%g\n", Nevent, Name.c_str(), prePoint->GetPosition().x(), prePoint->GetPosition().y(), prePoint->GetPosition().z(), theTrack->GetKineticEnergy()/MeV);
	   fflush(stMu);
	 }
  }

  if(Name(0,2)=="e-" || Name(0,2)=="e+")
  {	
	 if(vname(0,5)=="Strip" || vname(0,5)=="PenSc" || vname(0,5)=="AirSc" || vname(0,5)=="PokSt" || vname(0,5)=="SkScM" || vname(0,5)=="SlScM")
	 {
	   if(stEl==NULL)   
	   {
	     sprintf(FElMuName, "Ev%07dKEl.dat", Nevent);
         stEl = fopen(FElMuName, "w");  
         if(stEl)
            fprintf(stEl, "Nevent\tFlEl\tx\ty\tz\tenergy\n"); 
	   }

	   fprintf(stEl,"%d\t%s\t%g\t%g\t%g\t%g\n", Nevent, Name.c_str(), prePoint->GetPosition().x(), prePoint->GetPosition().y(), prePoint->GetPosition().z(), theTrack->GetKineticEnergy()/MeV);
	   fflush(stEl);
	 }
  }

  if(Name(0,6)=="proton")
  {	
	 if(vname(0,5)=="Strip" || vname(0,5)=="PenSc" || vname(0,5)=="AirSc" || vname(0,5)=="PokSt" || vname(0,5)=="SkScM" || vname(0,5)=="SlScM")
	 {
	   if(stPr==NULL)   
	   {
	     sprintf(FElMuName, "Ev%07dKPr.dat", Nevent);
         stPr = fopen(FElMuName, "w");  
         if(stPr)
            fprintf(stPr, "Nevent\tFlEl\tx\ty\tz\tenergy\n"); 
	   }

	   fprintf(stPr,"%d\t%s\t%g\t%g\t%g\t%g\n", Nevent, Name.c_str(), prePoint->GetPosition().x(), prePoint->GetPosition().y(), prePoint->GetPosition().z(), theTrack->GetKineticEnergy()/MeV);
	   fflush(stPr);
	 }
  }

  if(Name(0,3)=="pi-" || Name(0,3)=="pi+")
  {	
	 if(vname(0,5)=="Strip" || vname(0,5)=="PenSc" || vname(0,5)=="AirSc" || vname(0,5)=="PokSt" || vname(0,5)=="SkScM" || vname(0,5)=="SlScM")
	 {
	   if(stPi==NULL)   
	   {
	     sprintf(FElMuName, "Ev%07dKPi.dat", Nevent);
         stPi = fopen(FElMuName, "w");  
         if(stPi)
            fprintf(stPi, "Nevent\tFlEl\tx\ty\tz\tenergy\n"); 
	   }

	   fprintf(stPi,"%d\t%s\t%g\t%g\t%g\t%g\n", Nevent, Name.c_str(), prePoint->GetPosition().x(), prePoint->GetPosition().y(), prePoint->GetPosition().z(), theTrack->GetKineticEnergy()/MeV);
	   fflush(stPi);
	 }
  }
  
}
