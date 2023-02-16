
#include "DetEventAction.hh"

#include "DetRunAction.hh"
#include "DetSteppingAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

extern G4int Nphot[8][128], Nfe[8][128], Nevent, FlMu[8], MatrMu[8][128], MatrAll[8][128];
extern FILE *st, *stMu, *stEl, *stPr, *stPi, *srStr;
extern G4double sumQ, KoordMu[8][3], MatrQ[8][128];

G4double NQ[41] = {103.1,  133.34, 167.76, 205.3,  244.36, 282.84, 318.35, 348.43, 370.79, 383.69,
				   386.06, 377.66, 359.24, 332.28, 298.84, 261.34, 222.26, 183.82, 147.88, 115.72,
				   88.12,  65.33,  47.18,  33.23,  22.87,  15.43,  10.26,  6.77,   4.5,    3.06,
				   2.18,   1.66,   1.36,   1.19,   1.1,    1.048,  1.02,   1,      0.9,    0.5,    0};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetEventAction::DetEventAction()
{ 
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetEventAction::~DetEventAction()
{ 
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetEventAction::BeginOfEventAction(const G4Event* event)
{ 
  G4int Nsl, Nst;
  
  for(Nsl = 0; Nsl<8; Nsl++)
  {
    FlMu[Nsl] = 0;
	KoordMu[Nsl][0] = 0;
	KoordMu[Nsl][1] = 0;
	KoordMu[Nsl][2] = 0;
    for(Nst = 0; Nst<128; Nst++)
    {
      Nphot[Nsl][Nst] = 0;  
      Nfe[Nsl][Nst] = 0;
	  MatrMu[Nsl][Nst] = 0; 
	  MatrAll[Nsl][Nst] = 0;
	  MatrQ[Nsl][Nst] = 0;
    }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetEventAction::EndOfEventAction(const G4Event* event)
{  
  G4double a, S, Qsl, Qsum;
  G4int k, j;
  G4int Nsl, Nst, FlSrStr;
  char FName[500];
  
  FlSrStr = 0;
  for(Nsl = 0; Nsl<8; Nsl++)
  for(Nst = 0; Nst<128; Nst++)
  {
    if(Nphot[Nsl][Nst] || MatrMu[Nsl][Nst])
	  FlSrStr = 1;

	if(Nfe[Nsl][Nst])
	{
      Qsum = 0;

      for(j=0; j<Nfe[Nsl][Nst]; j++)
      {
        Qsl = 0;
        S = sumQ*G4UniformRand();
        k = 0;
        do
        {
          S-=(NQ[k]+NQ[k+1])*0.5*0.1;
          k++;
        }
        while(S>0 && k<41);
        k--;
        S*=-1.;
  
        if(k>-1 && k<41)
        {
          a = (NQ[k+1]-NQ[k])/0.1;
          Qsl = (k+1)*0.1-(NQ[k+1]-sqrt(NQ[k+1]*NQ[k+1]-2*a*S))/a;
        }
        if(Qsl>0 && Qsl<4)
	      Qsum += Qsl;		
		//G4cout<<Nsl<<'\t'<<Nst<<'\t'<<j<<'\t'<<Qsl<<'\t'<<G4endl;
      }   

      if(Qsum>=0.25)
      {
        MatrAll[Nsl][Nst] = 1;
      } 
	  MatrQ[Nsl][Nst] = Qsum;
	}
  }
  
  if(FlSrStr)
  {
    sprintf(FName, "Ev%07dSrStr.dat", Nevent);
    srStr = fopen(FName, "w");  
    if(srStr)
    {
      fprintf(srStr, "Nevent\tNsl\tNst\tNphot\tNfe\tMatrAll\tMatrMu\tMatrQ\n"); 
	  for(Nsl = 0; Nsl<8; Nsl++)
      for(Nst = 0; Nst<128; Nst++)
	  if(Nphot[Nsl][Nst] || MatrMu[Nsl][Nst])
	  {
		fprintf(srStr,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g\n", Nevent, Nsl, Nst, Nphot[Nsl][Nst], Nfe[Nsl][Nst], MatrAll[Nsl][Nst], MatrMu[Nsl][Nst], MatrQ[Nsl][Nst]);     
	  } 
      fclose(srStr);
    }
  }

  if(stMu)   fclose(stMu);
  if(stEl)   fclose(stEl);
  if(stPr)   fclose(stPr);
  if(stPi)   fclose(stPi);
 
}  
