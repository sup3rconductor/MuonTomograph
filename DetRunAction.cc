#include "DetRunAction.hh"
#include "DetPrimaryGeneratorAction.hh"
#include "DetEventAction.hh"
#include "DetSteppingAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetRunAction::DetRunAction()
{
  ;
} 

DetRunAction::~DetRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << " start" << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetRunAction::EndOfRunAction(const G4Run* aRun)
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
