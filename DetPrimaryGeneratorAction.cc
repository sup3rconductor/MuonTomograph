#include "DetPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

extern G4int Nevent;
extern FILE *st, *spectr, *stMu, *stEl, *stPr, *stPi;
extern G4double mpi;

G4double zw1=3.85, xw1=9, yw1=15.5, xw2=-9, yw2=15.5, xw3=-9,
         yw3=-15.5, xw4=9, yw4=-15.5, fi1, fi2, fi3, fi4, xp, yp, zp;

G4String TypePart[9]={"mu+", "mu-", "pi+", "pi-", "e+", "e-", "proton", "neutron", "gamma"};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetPrimaryGeneratorAction::DetPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="mu-"); 
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  fParticleGun->SetParticleEnergy(3.*GeV);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetPrimaryGeneratorAction::~DetPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4int Tip;  
  G4double xv, yv, zv, xn, yn, zn, ga, yz, Len, tx, ty, tz, p; 
  G4double  xx, yy, zz, xxs, yys, zzs, zps, xps, yps, xns, yns, zns, xls,
            yls, zls, xvs, yvs, zvs;
  G4double TetaGr, FiGr, Energy;
  char FName[500], PN[50];
  G4String NameOfPart;
    
  fscanf(spectr,"%d\t%s\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
		 &Nevent, &PN, &Tip, &TetaGr, &FiGr, &Energy, &xx, &yy, &zz);
  
  stMu = NULL;
  stEl = NULL;
  stPr = NULL;
  stPi = NULL;
  
  fi1 =       atan(fabs((yw1-yy)/(xw1-xx)))*180.0/mpi;
  fi2 = 180 - atan(fabs((yw2-yy)/(xw2-xx)))*180.0/mpi;
  fi3 = 180 + atan(fabs((yw3-yy)/(xw3-xx)))*180.0/mpi;
  fi4 = 360 - atan(fabs((yw4-yy)/(xw4-xx)))*180.0/mpi;


  if(TetaGr>1)
  {    
    if(FiGr>=fi1 && FiGr<=fi2) 
    {
      zps=((yw1-yy)/(sin(TetaGr*mpi/180)*sin(FiGr*mpi/180)))*cos(TetaGr*mpi/180)+zz;
      if(zps<zw1)
      {
       	xps=((yw1-yy)/sin(FiGr*mpi/180))*cos(FiGr*mpi/180)+xx; 
        yps=yw1;
          
        xxs=xps;
        yys=yps;
        zzs=zps;
      } 
      else
      {
        xp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*cos(FiGr*mpi/180)+xx;
        yp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*sin(FiGr*mpi/180)+yy; 
        zp=zw1;
           
        xxs=xp;
        yys=yp;
        zzs=zp;
      }   
    }
     
    if(FiGr>=fi2 && FiGr<=fi3)
    {
      zns=((xw2-xx)/(sin(TetaGr*mpi/180)*cos(FiGr*mpi/180)))*cos(TetaGr*mpi/180)+zz;
     
      if(zns<zw1)
      {
        yns=((xw2-xx)/cos(FiGr*mpi/180))*sin(FiGr*mpi/180)+yy; 
        xns=xw2;
           
        xxs=xns;
        yys=yns;
        zzs=zns;
      }
      else
      {
        xp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*cos(FiGr*mpi/180)+xx;
        yp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*sin(FiGr*mpi/180)+yy; 
        zp=zw1;
           
        xxs=xp;
        yys=yp;
        zzs=zp;
      }
                 
    }    
     
    if(FiGr>=fi3 && FiGr<=fi4)
    {
      zls=((yw3-yy)/(sin(TetaGr*mpi/180)*sin(FiGr*mpi/180)))*cos(TetaGr*mpi/180)+zz;
      if(zls<zw1)
      {
        xls=((yw3-yy)/sin(FiGr*mpi/180))*cos(FiGr*mpi/180)+xx; 
        yls=yw3;
           
        xxs=xls;
        yys=yls;
        zzs=zls;
      }
      else
      {
        xp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*cos(FiGr*mpi/180)+xx;
        yp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*sin(FiGr*mpi/180)+yy; 
        zp=zw1;
           
        xxs=xp;
        yys=yp;
        zzs=zp;
      }       
    } 
    
       
    if( (FiGr>=0 && FiGr<=fi1) || (FiGr>=fi4 && FiGr<=360))//!!!!!!!!!
    {
  	  zvs=((xw4-xx)/(sin(TetaGr*mpi/180)*cos(FiGr*mpi/180)))*cos(TetaGr*mpi/180)+zz;
      if(zvs<zw1)
      {
        yvs=((xw4-xx)/cos(FiGr*mpi/180))*sin(FiGr*mpi/180)+yy; 
        xvs=xw4;
      
        xxs=xvs;
        yys=yvs;
        zzs=zvs;
      }
      else
      {
        xp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*cos(FiGr*mpi/180)+xx;
        yp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*sin(FiGr*mpi/180)+yy; 
        zp=zw1;
           
        xxs=xp;
        yys=yp;
        zzs=zp;
      }       
    }
  }
  else
  {
    xp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*cos(FiGr*mpi/180)+xx;
    yp=((zw1-zz)/cos(TetaGr*mpi/180))*sin(TetaGr*mpi/180)*sin(FiGr*mpi/180)+yy; 
    zp=zw1;
           
    xxs=xp;
    yys=yp;
    zzs=zp;
  }

  Len = sqrt( (xx-xxs)*(xx-xxs) + (yy-yys)*(yy-yys) + (zz-zzs)*(zz-zzs) );

  tx = (xx-xxs)/Len;
  ty = (yy-yys)/Len;
  tz = (zz-zzs)/Len;

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(tx, ty, tz));
  fParticleGun->SetParticlePosition(G4ThreeVector(xxs*m,yys*m,zzs*m));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
    = particleTable->FindParticle(TypePart[Tip]);
  fParticleGun->SetParticleDefinition(particle);

  fParticleGun->SetParticleEnergy(Energy*GeV);

  fprintf(st,"%d\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Nevent, TypePart[Tip].c_str(), TetaGr, FiGr, Energy, xx, yy, zz);
  fflush(st);

  fParticleGun->GeneratePrimaryVertex(anEvent); 

  G4cout<<"Nevent"<<'\t'<<Nevent<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

