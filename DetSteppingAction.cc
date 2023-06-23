 
#include "DetSteppingAction.hh"
#include "DetEventAction.hh"
#include "DetDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern G4int Nevent;
extern G4double Z0const;
extern char fpartname[7];
G4int copyNo, strCopy;
long int Nph[1152];
float Tph[1152][200000];
float Energy, Epart, Eph[1152][200000];
extern FILE * CHARGED;
//G4double strx, stry, strz, strux, struy, struz;

DetSteppingAction::DetSteppingAction(DetEventAction* eventAction)
    : G4UserSteppingAction(), fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetSteppingAction::~DetSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* theTrack = step->GetTrack();
  G4String Name = theTrack->GetDefinition()->GetParticleName();
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  G4StepPoint* prePoint = step->GetPreStepPoint();
  G4String vname = prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName(); 



  if (vname == "sipm_l") // SiPM volume
  {
    if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
    {
      Energy = theTrack->GetKineticEnergy()/eV;
      copyNo = prePoint->GetTouchableHandle()->GetVolume()->GetCopyNo();
      if(Energy>=1.9 && Energy<=4.0 && Nph[copyNo]< 200000)
      {
        Eph[copyNo][Nph[copyNo]] = Energy;                             //Energy of Num N photon in Num copyNo of PMT
        Tph[copyNo][Nph[copyNo]] = theTrack->GetGlobalTime()/ns;        //Time of Num N photon in Num copyNo of PMT
        Nph[copyNo]++;                                                  //Plus one photon registered by Num copyNo of PMT
      }
    // collect optical photons 
      theTrack->SetTrackStatus(fStopAndKill);
    }   
  }

  if (vname == "strip_l") // Strip Volume
  {
      strCopy = prePoint->GetTouchableHandle()->GetVolume()->GetCopyNo();

      if (Name == "mu+" || Name == "mu-" || Name == "e+" || Name == "e-" || Name == "proton" || Name == "pi+" || Name == "pi-")
      {
          if (CHARGED) fprintf(CHARGED, "%s\t%d\t%lf\t%lf\t%lf\t%lf\n", Name.c_str(), strCopy, theTrack->GetKineticEnergy() / MeV, prePoint->GetPosition().x(), prePoint->GetPosition().y(), prePoint->GetPosition().z());
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

