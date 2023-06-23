
#include "DetDetectorConstruction.hh"
#include "DetActionInitialization.hh"

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"

/* #ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else */
#include "G4RunManager.hh"
//#endif

#include "G4UImanager.hh"
//#include "QBBC.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include <ctime>
#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int Nevent=0;
char Title[200]={0};
char NazFile[1500], FILEName[50];
FILE *PhEffect, *CHARGED;

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //



  G4UIExecutive* ui = 0;                                                          
  if ( argc == 1 ) {                                                                
    ui = new G4UIExecutive(argc, argv);                                             
  }                                                                              

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  //
//#ifdef G4MULTITHREADED
  //G4MTRunManager* runManager = new G4MTRunManager;
//#else
  G4RunManager* runManager = new G4RunManager;
//#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new DetDetectorConstruction());
  
  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  auto opticalParams = G4OpticalParameters::Instance();

  opticalParams->SetWLSTimeProfile("delta");
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);
    
  // User action initialization
  runManager->SetUserInitialization(new DetActionInitialization());
  
  // Initialize visualization                                                      
  //                                                                               
   G4VisManager* visManager = new G4VisExecutive;                                 
   //G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.     
   //G4VisManager* visManager = new G4VisExecutive("Quiet");                        
   visManager->Initialize();                                                      

  // Get the pointer to the User Interface manager
   G4UImanager* UImanager = G4UImanager::GetUIpointer();                            

  // Process macro or start UI session
  //
  if (!ui) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }                                                                                
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  } 

  //Clear clipboard and save data to files every 15 mins
  clock_t start_time = clock();
  double ExecTime = double(clock() - start_time) / CLOCKS_PER_SEC;

  if (ExecTime >= 900)
  {
      fflush(CHARGED);
      start_time = ExecTime;
  }
  //runManager->Initialize();
  //G4int numberOfEvent = 40000;//50000
  //runManager->BeamOn(numberOfEvent);

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  if (PhEffect) fclose(PhEffect);
  if (CHARGED) fclose(CHARGED);

  delete visManager;                                                             //Switching off visualisation
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
