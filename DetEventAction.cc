
#include "DetEventAction.hh"
#include "DetRunAction.hh"
#include "DetSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern long int Nph[1152];
extern G4int Nevent;
extern FILE *PhEffect, *CHARGED;

extern float Tph[1152][200000];
extern float Eph[1152][200000];
extern char NazFile[1500], FILEName[50];
//extern G4double x0, yy0, z0, teta, phi, x, y, z, ux, uy, uz;
//extern G4double strx, stry, strz, strux, struy, struz;
//extern char fpartname[7];

long int q, p;

DetEventAction::DetEventAction(DetRunAction* runAction)
: G4UserEventAction(), fRunAction(runAction)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetEventAction::~DetEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetEventAction::BeginOfEventAction(const G4Event*)
{    
     for (p = 0; p < 1152; p++)
     {
         for(q = 0; q < 200000; q++)
         {
             Tph[p][q] = 0;
             Eph[p][q] = 0;
         }
         Nph[p] = 0;
     } 

     sprintf(FILEName, "ChargedParticlesEvent%07d.dat", Nevent);
     CHARGED = fopen(FILEName, "w");
     if (CHARGED) fprintf(CHARGED, "Particle\tStripCopyNum\tEnergy, MeV\tX\tY\tZ\n");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetEventAction::EndOfEventAction(const G4Event*)
{   
    G4int NumSiPM, NumPhot;
    //fprintf(pos, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Nevent, scx, scy, scz, scux, scuy, scuz);

    sprintf(NazFile, "PhotonsEvent%05d.dat", Nevent);
    PhEffect = fopen(NazFile, "w");
    if(PhEffect)
    {
       fprintf(PhEffect, "Energy\tTime\tNcopy\n");
       for(NumSiPM = 0; NumSiPM < 1152; NumSiPM++)
       {
            for (NumPhot = 0; NumPhot < Nph[NumSiPM]; NumPhot++)
            {
                fprintf(PhEffect, "%lf\t%lf\t%d\n", Eph[NumSiPM][NumPhot], Tph[NumSiPM][NumPhot], NumSiPM);
            }
       } 
       fclose(PhEffect);   
    }
    Nevent++;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
