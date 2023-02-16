
#include "DetDetectorConstruction.hh"
#include "DetSteppingAction.hh"

#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

extern G4int NcopyToSlStr[1024][2], SlStrToNcopy[8][128];

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetDetectorConstruction::DetDetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetDetectorConstruction::~DetDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetDetectorConstruction::Construct()
{      
  
  G4double a, z;  // atomic mass, atomic number
  G4double density, fractionmass;
  G4int ncomponents, nelements, j;
  
  G4Element *elO  = new G4Element("Oxygen"   , "O"  , z=8.,  a = 16.00*g/mole);
  G4Element *elN  = new G4Element("Nitrogen" , "N"  , z=7.,  a = 14.01*g/mole);
  G4Element *elC  = new G4Element("Carbon"   , "C"  , z=6.,  a = 12.01*g/mole);
  G4Element *elH  = new G4Element("Hydrogen" , "H"  , z=1.,  a = 1.01*g/mole); 
  G4Element *elAl = new G4Element("Aluminium", "Al" , z=13., a = 26.98*g/mole);
  G4Element *elSi = new G4Element("Silicium" , "Si" , z=14., a = 28.086*g/mole);
  G4Element *elB  = new G4Element("Borum"    , "B"  , z=5.,  a = 10.812*g/mole);
  G4Element *elNa = new G4Element("Natrium"  , "Na" , z=11., a = 22.99*g/mole);
  G4Element *elK  = new G4Element("K"        , "K"  , z=19 , a=39.1*g/mole);
  G4Element *elCa = new G4Element("Calzium"  , "Ca" , z=31 , a=69.72*g/mole);
  G4Element *elFe = new G4Element("Iron"     , "Fe" , z=26 , a=55.85*g/mole);   
  G4Element *elHg = new G4Element("Hg"       , "Hg" , z=80 , a=200.59*g/mole);

  // Бетон
  G4Material *Concrete = new G4Material("Concrete", density = 2.1*g/cm3, 10);
  Concrete->AddElement(elH  , fractionmass= 0.01);
  Concrete->AddElement(elO  , fractionmass= 0.529);
  Concrete->AddElement(elNa , fractionmass= 0.016);
  Concrete->AddElement(elHg , fractionmass= 0.002);
  Concrete->AddElement(elAl , fractionmass= 0.034);
  Concrete->AddElement(elSi , fractionmass= 0.337);
  Concrete->AddElement(elK  , fractionmass= 0.013);
  Concrete->AddElement(elCa , fractionmass= 0.044);
  Concrete->AddElement(elFe , fractionmass= 0.014);
  Concrete->AddElement(elC  , fractionmass= 0.001);

  // Керамзит
  G4Material *KeramzitSiO2 = new G4Material("keramzit", density= 0.3*g/cm3, ncomponents=2);
  KeramzitSiO2->AddElement(elSi, nelements=1);
  KeramzitSiO2->AddElement(elO , nelements=2);

  // Кирпич
  G4Material *KirpichSiO2 = new G4Material("Kirpich",density= 1.8*g/cm3, ncomponents=2);
  KirpichSiO2->AddElement(elSi, nelements=1);
  KirpichSiO2->AddElement(elO , nelements=2);

  // Алюминий
  G4Material *AlMaterial = new G4Material("MAluminium", z = 13., a = 26.98*g/mole, density = 2.8*g/cm3);

  // Воздух
  G4Material *Air = new G4Material("MAir"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  // Полиэтилен
  G4Material *Polyetilen = new G4Material("Polyetilen",density = 0.067*g/cm3,nelements=2);
  Polyetilen->AddElement(elC, nelements=2);
  Polyetilen->AddElement(elH, nelements=4);

  // Пенопласт
  G4Material *Polystyrene = new G4Material("Polystyrene",density = 0.02*g/cm3,nelements=2);
  Polystyrene->AddElement(elC, nelements=8);
  Polystyrene->AddElement(elH, nelements=8);
  	
   // Полистерол для покрытия стрипа
  G4Material *PSp = new G4Material("MpStrip", density=1.04*g/cm3, nelements=2);
  PSp->AddElement(elC, nelements=8);
  PSp->AddElement(elH, nelements=8); 

  // Полистерол для стрипа
  G4Material *PS = new G4Material("MStrip", density=1.032*g/cm3, nelements=2);
  PS->AddElement(elC, nelements=8);
  PS->AddElement(elH, nelements=8); 
  
  // Внешняя оболочка оптоволокна 
  G4Material *FP = new G4Material("MFP", density=1.43*g/cm3, nelements=3);
  FP->AddElement(elC, nelements=5);
  FP->AddElement(elH, nelements=8); 
  FP->AddElement(elO, nelements=2); 

  // Средняя часть оптоволокна 
  G4Material *PMMA = new G4Material("MPMMA", density=1.19*g/cm3, nelements=3);
  PMMA->AddElement(elC, nelements=5);
  PMMA->AddElement(elH, nelements=8); 
  PMMA->AddElement(elO, nelements=2); 

  // Ядро-полистерол для оптоволокна 
  G4Material *PSOV = new G4Material("MPSoptov", density=1.05*g/cm3, nelements=2);
  PSOV->AddElement(elC, nelements=8);
  PSOV->AddElement(elH, nelements=8);

  // фотокатод - боросиликатное стекло
  G4Material *SiO2 = new G4Material("MSiO2", density=1.8*g/cm3, nelements=2);
  SiO2->AddElement(elO,  nelements=2);
  SiO2->AddElement(elSi, nelements=1);
  G4Material *Na2O = new G4Material("MNa2O", density=2.27*g/cm3, nelements=2);
  Na2O->AddElement(elNa, nelements=2);
  Na2O->AddElement(elO,  nelements=1);
  G4Material *Al2O3 = new G4Material("MAl2O3", density=3.97*g/cm3, nelements=2);
  Al2O3->AddElement(elAl, nelements=2);
  Al2O3->AddElement(elO,  nelements=3);
  G4Material *B2O3 = new G4Material("MB2O3", density=1.844*g/cm3, nelements=2);
  B2O3->AddElement(elB, nelements=2);
  B2O3->AddElement(elO, nelements=3);
  
  G4Material *StFEU = new G4Material("MFotKat", density=2.55*g/cm3, ncomponents=4);
  StFEU->AddMaterial(SiO2,  fractionmass=80.*perCent);
  StFEU->AddMaterial(Na2O,   fractionmass=4.*perCent);
  StFEU->AddMaterial(Al2O3, fractionmass=2.*perCent);
  StFEU->AddMaterial(B2O3,  fractionmass=14.*perCent); 

  G4double EnergyStFEU[2]     = {1.9*eV, 4.0*eV};
  G4double RefractiveStFEU[2] = {1.49, 1.49};
  G4double AbsLengthStFEU[2]  = {5*m, 5*m};

  G4MaterialPropertiesTable *MPT_StFEU = new G4MaterialPropertiesTable();
  MPT_StFEU->AddProperty("RINDEX",   EnergyStFEU, RefractiveStFEU, 2);
  MPT_StFEU->AddProperty("ABSLENGTH",EnergyStFEU, AbsLengthStFEU,  2); 
  StFEU->SetMaterialPropertiesTable(MPT_StFEU);

  G4double EnergyPSp[2] = {1.9*eV, 4.0*eV};
  G4double AbsLenPSp[2] = {5.0*m,  5.0*m};
  G4double RindPSp[2]   = {1.0, 1.0};

  G4MaterialPropertiesTable *PSpPT = new G4MaterialPropertiesTable();
  PSpPT->AddProperty("RINDEX",       EnergyPSp, RindPSp,   2);
  PSpPT->AddProperty("ABSLENGTH",    EnergyPSp, AbsLenPSp, 2);
  PSp->SetMaterialPropertiesTable(PSpPT);

  G4double EnergyAir[2] = {1.9*eV, 4.0*eV};
  G4double AbsLenAir[2] = {5.0*m,  5.0*m};
  G4double RindAir[2]   = {1.0, 1.0};

  G4MaterialPropertiesTable *AirPT = new G4MaterialPropertiesTable();
  AirPT->AddProperty("RINDEX",       EnergyAir, RindAir,   2);
  AirPT->AddProperty("ABSLENGTH",    EnergyAir, AbsLenAir, 2);
  Air->SetMaterialPropertiesTable(AirPT);

  G4double EnergyStrip[22] = {1.9*eV, 2.0*eV, 2.1*eV, 2.2*eV, 2.3*eV, 2.4*eV, 2.5*eV, 2.6*eV, 2.7*eV, 2.8*eV, 2.9*eV, 3.0*eV, 3.1*eV, 3.2*eV, 3.3*eV, 3.4*eV, 3.5*eV, 3.6*eV, 3.7*eV, 3.8*eV, 3.9*eV, 4.0*eV};
  G4double SpIzlStripa[22] = {0., 0., 0.00464, 0.00882, 0.02151, 0.04781, 0.09489, 0.20273, 0.33555, 0.55049, 1.0, 1.0, 0.60731, 0.22928, 0.11889, 0.11189, 0.12626, 0.12519, 0.09297, 0.07725, 0.01796, 0.00265};
  G4double AbsLenStripa[22], RindStrip[22];
  for(G4int k=0; k<22; k++) {AbsLenStripa[k] = 17.5*cm; RindStrip[k] = 1.581;}
      
  G4MaterialPropertiesTable *StripPT = new G4MaterialPropertiesTable();
  StripPT->AddProperty("RINDEX",       EnergyStrip, RindStrip,    22);
  StripPT->AddProperty("ABSLENGTH",    EnergyStrip, AbsLenStripa, 22);
  StripPT->AddProperty("FASTCOMPONENT",EnergyStrip, SpIzlStripa,  22);
  StripPT->AddConstProperty("SCINTILLATIONYIELD", 12./keV);
  StripPT->AddConstProperty("RESOLUTIONSCALE",    1.0);
  StripPT->AddConstProperty("FASTTIMECONSTANT",   1.6*ns);
  PS->GetIonisation()->SetBirksConstant(0.126*mm/MeV); 
  PS->SetMaterialPropertiesTable(StripPT);
  
  G4double EnergyOpt[10] = {1.9*eV, 2.2*eV, 2.3*eV, 2.4*eV, 2.56*eV, 2.66*eV, 2.68*eV, 3.69*eV, 3.7*eV, 4.0*eV};
  G4double AbsLenOpt[10] = {5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m, 0.1*mm, 0.1*mm, 5.0*m, 5.0*m};
  G4double SpIzlOpt[10]  = {0.001, 0.05, 0.25, 0.7, 1., 1., 0., 0., 0., 0.};
  G4double RindexOpt[10] = {1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59};

  G4MaterialPropertiesTable *OptPT = new G4MaterialPropertiesTable();
  OptPT->AddProperty("RINDEX",       EnergyOpt, RindexOpt, 10);
  OptPT->AddProperty("WLSABSLENGTH", EnergyOpt, AbsLenOpt, 10);
  OptPT->AddProperty("WLSCOMPONENT", EnergyOpt, SpIzlOpt,  10);
  OptPT->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
  PSOV->SetMaterialPropertiesTable(OptPT);
   
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = false;
   
  // World
  G4double world_sizeX = 9.175*m;
  G4double world_sizeY = 15.65*m;
  G4double world_sizeZ = 3.9*m;
  
  G4Box *solidWorld = new G4Box("World_g", world_sizeX, world_sizeY, world_sizeZ);        
  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, Air, "World_l");                                    
  G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps); 

  // Стены    
  G4Box *solidSt = new G4Box("Steni_g", 8.975*m, 15.45*m, 3.425*m);    			     
  G4LogicalVolume *logicSt = new G4LogicalVolume(solidSt, KirpichSiO2,	"Steni_l");
  G4VPhysicalVolume *physSt = new G4PVPlacement(0,	G4ThreeVector(0*m, 0*m, -0.35*m), logicSt, "Steni", logicWorld, false, 0, checkOverlaps);
  
  // Нишы в стенах
  G4Box *solidNishaM = new G4Box("NishaM_g", 0.2*m, 14.9*m, 0.15*m);
  G4LogicalVolume *logicNishaMV= new G4LogicalVolume(solidNishaM, Air, "NishaMV_l");
  G4VPhysicalVolume *physNishaMV = new G4PVPlacement(0, G4ThreeVector( 8.625*m, 0*m, 3.275*m), logicNishaMV, "NishaMV", logicSt, false,	0, checkOverlaps);	
  G4LogicalVolume *logicNishaMN= new G4LogicalVolume(solidNishaM, Air, "NishaMN_l");
  G4VPhysicalVolume *physNishaMN = new G4PVPlacement(0, G4ThreeVector(-8.625*m, 0*m, 3.275*m), logicNishaMN, "NishaMN", logicSt, false,	0, checkOverlaps);

  G4Box *solidNisha = new G4Box("Nisha_g", 0.125*m, 14.9*m, 1.075*m);
  G4LogicalVolume *logicNishaV= new G4LogicalVolume(solidNisha, Air, "NishaV_l");
  G4VPhysicalVolume *physNishaV = new G4PVPlacement(0, G4ThreeVector( 8.55*m, 0*m, 2.05*m), logicNishaV, "NishaV", logicSt, false,	0, checkOverlaps);	
  G4LogicalVolume* logicNishaN= new G4LogicalVolume(solidNisha, Air, "NishaN_l");
  G4VPhysicalVolume* physNishaN = new G4PVPlacement(0, G4ThreeVector(-8.55*m, 0*m, 2.05*m), logicNishaN, "NishaN", logicSt, false,	0, checkOverlaps);
  
  // Воздух
  G4Box* solidExpHall = new G4Box("ExpHall_g", 8.425*m, 14.9*m, 3.425*m);                       
  G4LogicalVolume* logicExpHall= new G4LogicalVolume(solidExpHall, Air, "ExpHall_l");
  G4VPhysicalVolume* physExpHall = new G4PVPlacement(0, G4ThreeVector(0*m, 0*m, 0*m), logicExpHall, "ExpHall", logicSt, false,	0, checkOverlaps);	
   
  // Бетон - крыша
  G4Box* solidKr = new G4Box("Krisha_g", 8.975*m, 15.45*m, 0.05*m);
  G4LogicalVolume* logicKr = new G4LogicalVolume(solidKr, Concrete, "Krisha_l");      				     
  G4VPhysicalVolume* physKr = new G4PVPlacement(0, G4ThreeVector(0*m,0*m,3.125*m), logicKr, "Krisha", logicWorld, false, 0, checkOverlaps);

  // Керамзит
  G4Box*  solidKeramzit= new G4Box("Keramzit_g", 8.545*m, 15.02*m, 0.125*m);
  G4LogicalVolume* logicKeramzit = new G4LogicalVolume(solidKeramzit, KeramzitSiO2, "Keramzit_l");       			                  
  G4VPhysicalVolume* physKeramzit = new G4PVPlacement(0, G4ThreeVector(0.*m,0.*m,3.30*m), logicKeramzit, "Keramzit", logicWorld, false, 0, checkOverlaps);

  // Бортик
  G4Box* solidBorK = new G4Box("BorK_g", 8.545*m, 0.215*m, 0.3*m);    			     
  G4LogicalVolume* logicBorL = new G4LogicalVolume(solidBorK, KirpichSiO2, "BorL_l");
  G4VPhysicalVolume* physBorL = new G4PVPlacement(0, G4ThreeVector(0*m,  15.235*m, 3.475*m), logicBorL, "BorL", logicWorld, false, 0, checkOverlaps);
  G4LogicalVolume* logicBorP = new G4LogicalVolume(solidBorK, KirpichSiO2, "BorP_l");
  G4VPhysicalVolume* physBorP = new G4PVPlacement(0, G4ThreeVector(0*m, -15.235*m, 3.475*m), logicBorP, "BorP", logicWorld, false, 0, checkOverlaps);
  G4Box* solidBorD = new G4Box("BorD_g", 0.215*m, 15.45*m, 0.3*m);    			     
  G4LogicalVolume* logicBorV = new G4LogicalVolume(solidBorD, KirpichSiO2, "BorV_l");
  G4VPhysicalVolume* physBorV = new G4PVPlacement(0, G4ThreeVector( 8.76*m, 0*m, 3.475*m), logicBorV, "BorV", logicWorld, false, 0, checkOverlaps);
  G4LogicalVolume* logicBorN = new G4LogicalVolume(solidBorD, KirpichSiO2, "BorN_l");
  G4VPhysicalVolume* physBorN = new G4PVPlacement(0, G4ThreeVector(-8.76*m, 0*m, 3.475*m), logicBorN, "BorN", logicWorld, false, 0, checkOverlaps);   
  
  // СцМГ
     
  G4double DlinaStripa  = 3460*mm; // длина стрипа 
  G4double VisotaStripa = 10.6*mm; // высота стрипа
  G4double ShirinStripa = 26.3*mm; // ширина стрипа
  
  G4double RadiusKan = 0.8*mm;  // радиус нижней круглой части канавки
  G4double VisotaKan = 1.5*mm;  // высота верхней прямоуголной части канавки
  G4double ShirinKan = 1.6*mm;  // ширина верхней прямоуголной части канавки
  
  G4double TolBokPokr  = 0.05*mm; // толщина скотча и покраски сбоку  
  
  G4double RadiusOpt = 0.5*mm;  // радиус оптоволокна
  G4double DlinaOpt  = 500*mm;  // длина оптоволокна вне стрипа 15
  
  G4double ShirinStFEU = 2*mm;  // ширина  стекла ФЭУ
  G4double DlinaStFEU = 1*mm;   // толщина стекла ФЭУ
  G4double VisotaStFEU = 2*mm;  // высота  стекла ФЭУ

  G4double ShirinFKFEU = 2*mm;  // ширина  фотокатода ФЭУ
  G4double DlinaFKFEU = 0.1*mm; // толщина фотокатода ФЭУ
  G4double VisotaFKFEU = 2*mm;  // высота  фотокатода ФЭУ

  G4double TolAl = 0.8*mm;    // толщина алюминия
  G4double TolSk = 1.1*mm;    // толщина скотча
  G4double RastMStr = 0.1*mm; // расстояние между стрипами
  G4double SdvigStr = 0.4*mm; // дополнительный "сдвиг" между стрипами двух БМ
  
  
  G4double DlinaSlScMg  = 2*TolAl+2*TolSk+DlinaStripa+DlinaStFEU+DlinaFKFEU;      // длина слоя СцМГ 
  G4double VisotaSlScMg = 2*TolAl+2*TolSk+VisotaStripa;                           // высота слоя СцМГ
  G4double ShirinSlScMg = 2*TolAl+2*TolSk+128*ShirinStripa+127*RastMStr+SdvigStr; // ширина слоя СцМГ

  G4double DlinaSk  = 2*TolSk+DlinaStripa+DlinaStFEU+DlinaFKFEU;      // длина скотча СцМГ 
  G4double VisotaSk = 2*TolSk+VisotaStripa;                           // высота скотча СцМГ
  G4double ShirinSk = 2*TolSk+128*ShirinStripa+127*RastMStr+SdvigStr; // ширина скотча СцМГ

  G4double DlinaAir  = DlinaStripa+DlinaStFEU+DlinaFKFEU;      // длина воздуха СцМГ 
  G4double VisotaAir = VisotaStripa;                           // высота воздуха СцМГ
  G4double ShirinAir = 128*ShirinStripa+127*RastMStr+SdvigStr; // ширина воздуха СцМГ

  G4double MPlR = 33*cm-2*VisotaSlScMg;      // толщина пенопласта координатными плоскостями

  // положение СцМГ
  G4double ScMg_xpos = (6555-2225)*mm; 
  G4double ScMg_ypos = (22.586-13.0)*m;
  G4double ScMg_zpos = -252.5*cm;

  // Пенопласт 
  G4Box* solidPenScMg = new G4Box("PenScMg_8", ShirinSlScMg*0.5, ShirinSlScMg*0.5, MPlR*0.5);                         
  G4LogicalVolume* logicPenScMg = new G4LogicalVolume(solidPenScMg, Polystyrene, "PenScMg_l");                                   
  G4VPhysicalVolume* physPenScMgN = new G4PVPlacement(0, G4ThreeVector(ScMg_xpos, ScMg_ypos, ScMg_zpos-MPlR-2*VisotaSlScMg), logicPenScMg, "PenScMgN", logicExpHall, false, 0, checkOverlaps);
  G4VPhysicalVolume* physPenScMgS = new G4PVPlacement(0, G4ThreeVector(ScMg_xpos, ScMg_ypos, ScMg_zpos),                     logicPenScMg, "PenScMgS", logicExpHall, false, 0, checkOverlaps);
  G4VPhysicalVolume* physPenScMgV = new G4PVPlacement(0, G4ThreeVector(ScMg_xpos, ScMg_ypos, ScMg_zpos+MPlR+2*VisotaSlScMg), logicPenScMg, "PenScMgV", logicExpHall, false, 0, checkOverlaps);
  
  G4LogicalVolume *logicSlScMg[8]={NULL}, *logicSk[8]={NULL}, *logicAir[8]={NULL};
  G4VPhysicalVolume *physSlScMg[8]={NULL}, *physSk[8]={NULL}, *physAir[8]={NULL};

  G4Box *solidSlScMg = new G4Box("SlScMg_g", ShirinSlScMg*0.5, DlinaSlScMg*0.5, VisotaSlScMg*0.5);                         
  G4Box *solidSk     = new G4Box("SkScMg_g", ShirinSk*0.5, DlinaSk*0.5, VisotaSk*0.5);  
  G4Box *solidAir    = new G4Box("AirScMg_8", ShirinAir*0.5, DlinaAir*0.5, VisotaAir*0.5); 

  G4int Nsl, Nst, NCopy, NaprStr;
  G4double xsl, ysl, zsl;
  G4RotationMatrix* PovSl[8]; 

  // положение стрипа
  G4double strip_xpos, strip_ypos, strip_zpos;

  G4Box *solidPokStr[8][128]={NULL}, *solidStr[8][128]={NULL}, *solidSkotch[8][128]={NULL}, *solidStFEU[8][128]={NULL}, *solidFKFEU[8][128]={NULL};
  G4Tubs *solidOptVn[8][128]={NULL}, *solidOptPok[8][128]={NULL};
  G4LogicalVolume *logicPokStr[8][128]={NULL}, *logicStr[8][128]={NULL}, *logicSkotch[8][128]={NULL}, *logicOptVn[8][128]={NULL}, 
	              *logicOptPok[8][128]={NULL}, *logicStFEU[8][128]={NULL}, *logicFKFEU[8][128]={NULL};
  G4VPhysicalVolume *physPokStr[8][128]={NULL}, *physStr[8][128]={NULL}, *physSkotch[8][128]={NULL}, *physOptVn[8][128]={NULL}, 
	                *physOptPok[8][128]={NULL}, *physStFEU[8][128]={NULL}, *physFKFEU[8][128]={NULL};  

  NCopy = 0;
  for(Nsl = 0; Nsl<8; Nsl++)
  {
	PovSl[Nsl] = new G4RotationMatrix; 
	if(Nsl%2)
	{
	  PovSl[Nsl]->rotateZ(0*deg);
	  NaprStr = 1;
	}
	else
	{
	  PovSl[Nsl]->rotateZ(90*deg);
	  NaprStr = -1;
	}  
	xsl = ScMg_xpos; 
	ysl = ScMg_ypos; 
	zsl = ScMg_zpos-1.5*MPlR-3.5*VisotaSlScMg + Nsl*VisotaSlScMg + MPlR*int(Nsl/2);

    logicSlScMg[Nsl] = new G4LogicalVolume(solidSlScMg, AlMaterial, "SlScMg_l");                                   
    physSlScMg[Nsl] = new G4PVPlacement(PovSl[Nsl], G4ThreeVector(xsl, ysl, zsl), logicSlScMg[Nsl],"SlScMg", logicExpHall, false, 0, checkOverlaps);
                         
    logicSk[Nsl] = new G4LogicalVolume(solidSk, Polyetilen, "SkScMg_l");                                   
    physSk[Nsl] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicSk[Nsl],"SkScMg", logicSlScMg[Nsl], false, 0, checkOverlaps);
                          
    logicAir[Nsl] = new G4LogicalVolume(solidAir, Air, "AirScMg_l");                                   
    physAir[Nsl] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicAir[Nsl],"AirScMg", logicSk[Nsl], false, 0, checkOverlaps);

	G4RotationMatrix* PovOpt = new G4RotationMatrix; // поворот оптоволокна
    PovOpt->rotateX(-90*deg);

	for(Nst = 0; Nst<128; Nst++)
	{
	  strip_xpos = NaprStr*(-63.5*ShirinStripa-63.5*RastMStr-SdvigStr/2.0 + Nst*(ShirinStripa+RastMStr)+SdvigStr*int(Nst/64));
      strip_ypos = (DlinaStFEU+DlinaFKFEU)*0.5;
      strip_zpos = 0*mm;

	  // покрытие стрипа
	  solidPokStr[Nsl][Nst] = new G4Box("PokStr_g", 0.5*ShirinStripa, 0.5*DlinaStripa, 0.5*VisotaStripa);       
      logicPokStr[Nsl][Nst] = new G4LogicalVolume(solidPokStr[Nsl][Nst], PSp, "PokStr_l");               
      physPokStr[Nsl][Nst]  = new G4PVPlacement(0, G4ThreeVector(strip_xpos, strip_ypos, strip_zpos), logicPokStr[Nsl][Nst], "PokStr", logicAir[Nsl], false, 0, checkOverlaps);  

	  // стрип
      solidStr[Nsl][Nst] = new G4Box("Strip_g", 0.5*(ShirinStripa-0.3*mm), 0.5*(DlinaStripa-2*TolBokPokr), 0.5*(VisotaStripa-0.6*mm));       
      logicStr[Nsl][Nst] = new G4LogicalVolume(solidStr[Nsl][Nst], PS, "Strip_l");               
      physStr[Nsl][Nst]  = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicStr[Nsl][Nst], "Strip", logicPokStr[Nsl][Nst], false, NCopy, checkOverlaps); 

	  // Скотч в конце стрипа
      solidSkotch[Nsl][Nst] = new G4Box("Skotch_g", 0.5*ShirinStripa, 0.5*TolBokPokr, 0.5*VisotaStripa);       
      logicSkotch[Nsl][Nst] =  new G4LogicalVolume(solidSkotch[Nsl][Nst], AlMaterial, "Skotch_l");               
      physSkotch[Nsl][Nst] = new G4PVPlacement(0, G4ThreeVector(0, (DlinaStripa-TolBokPokr)/2.0, 0), logicSkotch[Nsl][Nst], "Skotch", logicPokStr[Nsl][Nst], false, NCopy, checkOverlaps);	  
	  
	  // Оптоволокно
      solidOptVn[Nsl][Nst] = new G4Tubs("OptVn_g", 0, RadiusOpt, 0.5*(DlinaStripa-2*TolBokPokr), 0.*deg, 360.*deg);         
      logicOptVn[Nsl][Nst] =  new G4LogicalVolume(solidOptVn[Nsl][Nst], PSOV,   "OptVn_l");                 
      physOptVn[Nsl][Nst] = new G4PVPlacement(PovOpt, G4ThreeVector(0, 0, 0.5*VisotaStripa-0.3*mm-VisotaKan), logicOptVn[Nsl][Nst], "OptVn", logicStr[Nsl][Nst], false, NCopy, checkOverlaps);
  
      solidOptPok[Nsl][Nst] = new G4Tubs("OptPok_g", 0, RadiusOpt, 0.5*TolBokPokr, 0.*deg, 360.*deg);       
      logicOptPok[Nsl][Nst] =  new G4LogicalVolume(solidOptPok[Nsl][Nst], PSOV, "OptPok_l");               
      physOptPok[Nsl][Nst] = new G4PVPlacement(PovOpt, G4ThreeVector(0, -(DlinaStripa-TolBokPokr)/2.0, 0.5*VisotaStripa-0.3*mm-VisotaKan), logicOptPok[Nsl][Nst], "OptPok", logicPokStr[Nsl][Nst], false, NCopy, checkOverlaps);
  
	  // Стекло ФЭУ
      solidStFEU[Nsl][Nst]  = new G4Box("StFEUg", ShirinStFEU/2.0, DlinaStFEU/2.0, VisotaStFEU/2.0);
      logicStFEU[Nsl][Nst] = new G4LogicalVolume(solidStFEU[Nsl][Nst], StFEU, "StFEUl",0,0,0);
      physStFEU[Nsl][Nst] = new G4PVPlacement(0, G4ThreeVector(strip_xpos, strip_ypos-DlinaStripa/2.0-DlinaStFEU/2.0, strip_zpos+VisotaStripa/2.0-VisotaKan-0.3*mm), logicStFEU[Nsl][Nst], "StFEU", logicAir[Nsl], false, NCopy, checkOverlaps);
	  
	  // Фотокатод
      solidFKFEU[Nsl][Nst]  = new G4Box("FKFEUg", ShirinFKFEU/2.0, DlinaFKFEU/2.0, VisotaFKFEU/2.0);
      logicFKFEU[Nsl][Nst] = new G4LogicalVolume(solidFKFEU[Nsl][Nst], AlMaterial, "FKFEUl",0,0,0);
      physFKFEU[Nsl][Nst] = new G4PVPlacement(0, G4ThreeVector(strip_xpos, strip_ypos-DlinaStripa/2.0-DlinaStFEU-DlinaFKFEU/2.0, strip_zpos+VisotaStripa/2.0-VisotaKan-0.3*mm), logicFKFEU[Nsl][Nst], "FKFEU", logicAir[Nsl], false, NCopy, checkOverlaps);
  	  
	  NcopyToSlStr[NCopy][0] = Nsl;
	  NcopyToSlStr[NCopy][1] = Nst;

	  SlStrToNcopy[Nsl][Nst]=NCopy;
	  
	  NCopy++;

	}
  }   
  
  // Поверхность стрипа: диффузное отражение
  G4double reflectivity_str[2]= {0.98, 0.98};
  G4double PhotonEnergyPov[2] = {1.9*eV, 4*eV};
  G4double reflectivity_skotch[2] = {0.83, 0.83};

  G4OpticalSurface *OpticalPovStripa = new G4OpticalSurface("PovStripaSurface");
  OpticalPovStripa->SetModel(unified);
  OpticalPovStripa->SetType(dielectric_dielectric);
  OpticalPovStripa->SetFinish(groundfrontpainted);
   
  G4MaterialPropertiesTable *PovStripaSurfacePT = new G4MaterialPropertiesTable();
  PovStripaSurfacePT->AddProperty("REFLECTIVITY", PhotonEnergyPov, reflectivity_str, 2);
  OpticalPovStripa->SetMaterialPropertiesTable(PovStripaSurfacePT);

  // Поверхность торца стрипа, покрытого скотчем: зеркальное отражение
  G4OpticalSurface* OpSurfaceSkotch = new G4OpticalSurface("SkotchSurface");
  OpSurfaceSkotch -> SetType(dielectric_metal);
  OpSurfaceSkotch -> SetFinish(ground);
  OpSurfaceSkotch -> SetModel(glisur);
    
  G4MaterialPropertiesTable *SkotchSurfacePT = new G4MaterialPropertiesTable();
  SkotchSurfacePT -> AddProperty("REFLECTIVITY", PhotonEnergyPov, reflectivity_skotch, 2);

  OpSurfaceSkotch -> SetMaterialPropertiesTable(SkotchSurfacePT);

  G4LogicalBorderSurface *ScintillatorPSTiO2[8][128]={NULL}, *ScintillatorSkotch[8][128]={NULL}, *OptovSkotch[8][128]={NULL}, *StecloFotoKat[8][128]={NULL};

  for(Nsl = 0; Nsl<8; Nsl++)
  for(Nst = 0; Nst<128; Nst++)
  {
    ScintillatorPSTiO2[Nsl][Nst] = new G4LogicalBorderSurface("StripTiO2Surface",     physStr[Nsl][Nst],   physPokStr[Nsl][Nst], OpticalPovStripa);
    ScintillatorSkotch[Nsl][Nst] = new G4LogicalBorderSurface("StripSkotchSurface",   physStr[Nsl][Nst],   physSkotch[Nsl][Nst], OpSurfaceSkotch);
    OptovSkotch[Nsl][Nst]        = new G4LogicalBorderSurface("OptovSkotchSurface",   physOptVn[Nsl][Nst], physSkotch[Nsl][Nst], OpSurfaceSkotch);
    StecloFotoKat[Nsl][Nst]      = new G4LogicalBorderSurface("StecloFotoKatSurface", physStFEU[Nsl][Nst], physFKFEU[Nsl][Nst],  OpSurfaceSkotch);
  }


  G4Colour grey     (0.5, 0.5, 0.5);
  G4Colour black    (0.0, 0.0, 0.0);
  G4Colour red      (1.0, 0.0, 0.0);
  G4Colour green    (0.0, 1.0, 0.0);
  G4Colour blue     (0.0, 0.0, 1.0);
  G4Colour cyan     (0.0, 1.0, 1.0);
  G4Colour magenta  (1.0, 0.0, 1.0);
  G4Colour yellow   (1.0, 1.0, 0.0);

  G4VisAttributes* VisAtt_green = new G4VisAttributes(green);
  VisAtt_green->SetForceSolid(true);

  G4VisAttributes* VisAtt_cyan = new G4VisAttributes(cyan);
  VisAtt_cyan->SetForceSolid(true);

  G4VisAttributes* VisAtt_yellow = new G4VisAttributes(yellow);
  VisAtt_yellow->SetForceSolid(true);

  G4VisAttributes* VisAtt_grey = new G4VisAttributes(grey);
  VisAtt_grey->SetForceSolid(true);

  G4VisAttributes* VisAtt_magenta = new G4VisAttributes(magenta);
  VisAtt_magenta->SetForceSolid(true);

  G4VisAttributes* VisAtt_blue = new G4VisAttributes(blue);
  VisAtt_blue->SetForceSolid(true);
  
  logicKr->SetVisAttributes(VisAtt_green);
  logicKeramzit->SetVisAttributes(VisAtt_magenta);

  logicNishaMV->SetVisAttributes(VisAtt_yellow);  
  logicNishaMN->SetVisAttributes(VisAtt_yellow);  
  logicNishaV ->SetVisAttributes(VisAtt_yellow);   
  logicNishaN ->SetVisAttributes(VisAtt_yellow);  
  logicExpHall->SetVisAttributes(VisAtt_yellow);  
  
  logicSt  ->SetVisAttributes(VisAtt_cyan);
  logicBorL->SetVisAttributes(VisAtt_cyan);
  logicBorP->SetVisAttributes(VisAtt_cyan);
  logicBorV->SetVisAttributes(VisAtt_cyan);
  logicBorN->SetVisAttributes(VisAtt_cyan);

  logicPenScMg->SetVisAttributes(VisAtt_green);

  for(Nsl=0; Nsl<8; Nsl++)
  {
    logicSlScMg[Nsl]->SetVisAttributes(VisAtt_grey);
    logicSk[Nsl]->SetVisAttributes(VisAtt_cyan);
    logicAir[Nsl]->SetVisAttributes(VisAtt_yellow);
	for(Nst = 0; Nst<128; Nst++)
	{
		logicOptVn[Nsl][Nst]->SetVisAttributes(VisAtt_yellow);
		logicOptPok[Nsl][Nst]->SetVisAttributes(VisAtt_yellow);
		logicFKFEU[Nsl][Nst]->SetVisAttributes(VisAtt_magenta);
		logicStFEU[Nsl][Nst]->SetVisAttributes(VisAtt_green);

		if(Nst==0 || Nst==63 || Nst==64)
		logicStr[Nsl][Nst]->SetVisAttributes(VisAtt_magenta);
        else
		logicStr[Nsl][Nst]->SetVisAttributes(VisAtt_blue);
	}
  }

  // return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
