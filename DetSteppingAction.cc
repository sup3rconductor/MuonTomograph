 
#include "DetSteppingAction.hh"
#include "DetEventAction.hh"
#include "DetDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern G4int Nevent;
extern G4double Z0const;
extern char fpartname[7];
extern G4int NRows;
G4int copyNo;
G4int Nsipm[NRows];
G4double Tph[NRows][150000];
G4double Energy, Eph[NRows][150000];
G4double scx, scy, scz, scux, scuy, scuz;

DetSteppingAction::DetSteppingAction(DetEventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction), fScoringVolume(0)
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

  if (vname == "strip_l") //Strip volume
  {
    if (Name == fpartname)
    {
      G4ThreeVector worldPosition = theTrack->GetPosition();
      scx = worldPosition.x();
      scy = worldPosition.y();
      scz = worldPosition.z();
      G4ThreeVector direction = theTrack->GetMomentumDirection();
      scux = direction.x();
      scuy = direction.y();
      scuz = direction.z(); 
    }
  }  

  if (vname == "tdlr_l") //Tedlar volume
  {
      if (particleType == G4OpticalPhoton::OpticalPhotonDefinition)
      {
          theTrack->SetStatus(fStopAndKill); //Killing optical photons in Tedlar
      }
  }


  if (vname == "sipm_l") // SiPM volume
  {
    if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
    {
      Energy = theTrack->GetKineticEnergy()/eV;
      copyNo = prePoint->GetTouchableHandle()->GetVolume()->GetCopyNo();
      if(Energy1>=1.9 && Energy1<=4.0 && Nph[copyNo]<150000)
      {
        Eph[copyNo][Nph[copyNo]] = Energy;                             //Energy of Num N photon in Num copyNo of PMT
        Tph[copyNo][Nph[copyNo]] = theTrack->GetGlobalTime()/ns;        //Time of Num N photon in Num copyNo of PMT
        Nsipm[copyNo]++;                                                  //Plus one photon registered by Num copyNo of PMT
      }
    // collect optical photons 
      theTrack->SetTrackStatus(fStopAndKill);
    }   
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

