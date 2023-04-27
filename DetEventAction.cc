
#include "DetEventAction.hh"
#include "DetRunAction.hh"
#include "DetSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern G4int Nph[160];
extern G4int Nevent;
extern FILE *put, *pos;

extern G4double Tph[160][150000];
extern G4double Eph[160][150000], E0;
extern G4double x0, yy0, z0, teta, phi, x, y, z, ux, uy, uz;
extern G4double scx, scy, scz, scux, scuy, scuz;
extern char fpartname[7];
extern char NazFile[1500];

G4int j, p;

DetEventAction::DetEventAction(DetRunAction* runAction)
: G4UserEventAction(), fRunAction(runAction)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetEventAction::~DetEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetEventAction::BeginOfEventAction(const G4Event*)
{    
     /*for (p = 0; p < 160; p++)
     {
         for(j = 0; j < 15000; j++)
         {
             Tph[p][j] = 0;
             Eph[p][j] = 0;
         }
         Nph[p] = 0;
     } */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetEventAction::EndOfEventAction(const G4Event*)
{   
   /* G4int NumPMT;
    fprintf(pos, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Nevent, scx, scy, scz, scux, scuy, scuz);

    sprintf(NazFile, "EventM05%05d.dat", Nevent);
    put = fopen(NazFile, "w");
    if(put)    
    {
       fprintf(put, "PmtNumber\tPhotNum\tEnergy\tTime\n");     
       for(NumPMT = 0; NumPMT < 160; NumPMT++)
       {
            for (j = 0; j < Nph[NumPMT]; j++)
            {
                fprintf(put, "%d\t%d\t%lf\t%lf\n", NumPMT, j, Eph[NumPMT][j], Tph[NumPMT][j]);
            }
        }
        
        fclose(put);
        
    }
    G4cout << "Number of event: " << Nevent << G4endl;
    Nevent++; */
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
