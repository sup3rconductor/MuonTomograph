#include "DetDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Z0const, X0const, Y0const;

DetDetectorConstruction::DetDetectorConstruction()
	: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetDetectorConstruction::~DetDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetDetectorConstruction::Construct()
{
	/*	MATERIALS	*/

	G4double a, z;  //Atomic mass, atomic number
	G4double density, fractionmass;
	G4int ncomponents, nelements;


	//Chemical elements
	G4Element* elH = new G4Element("Hydrogen", "H", z = 1., a = 1.01 * g / mole);
	G4Element* elC = new G4Element("Carbon", "C", z = 6., a = 12.01 * g / mole);
	G4Element* elN = new G4Element("Nitrogen", "N", z = 7., a = 14.01 * g / mole);
	G4Element* elO = new G4Element("Oxygen", "O", z = 8., a = 16.00 * g / mole);
	G4Element* elSi = new G4Element("Silicium", "Si", z = 14., a = 28.09 * g / mole);
	G4Element* elAl = new G4Element("Aluminium", "Al", z = 13., a = 26.98 * g / mole);
	G4Element* elB = new G4Element("Boron", "B", z = 5., a = 10.812 * g / mole);
	G4Element* elFe = new G4Element("Ferrum", "Fe", z = 26., a = 55.85 * g / mole);
	G4Element* elF = new G4Element("Fluor", "F", z = 17., a = 18.99 * g / mole);


	//Air
	G4Material* Air = new G4Material("MAir", density = 1.290 * mg / cm3, ncomponents = 2);
	Air->AddElement(elN, fractionmass = 0.8);
	Air->AddElement(elO, fractionmass = 0.2);

	//Aluminium
	G4Material* AlMaterial = new G4Material("MAluminium", z = 13., a = 26.98 * g / mole, density = 2.8 * g / cm3);

	//Ferrum
	G4Material* FeMaterial = new G4Material("MFerrum", z = 26., a = 55.85 * g / mole, density = 7.9 * g / cm3);

	//Cover of fiberglass 
	G4Material* PMMA = new G4Material("MPMMA", density = 1.19 * g / cm3, ncomponents = 3);
	PMMA->AddElement(elC, nelements = 5);
	PMMA->AddElement(elH, nelements = 8);
	PMMA->AddElement(elO, nelements = 2);

	//Core of fiberglass
	G4Material* PS = new G4Material("MPS", density = 1.05 * g / cm3, ncomponents = 2);
	PS->AddElement(elC, nelements = 8);
	PS->AddElement(elH, nelements = 8);

	//Tedlar
	G4Material* PVF = new G4Material("MPVF", density = 1.39 * g / cm3, ncomponents = 3);
	PVF->AddElement(elC, nelements = 2);
	PVF->AddElement(elH, nelements = 3);
	PVF->AddElement(elF, nelements = 1);


	//Strip material
	G4Material* PFT = new G4Material("MPFT", density = 1.19 * g / cm3, ncomponents = 2);
	PFT->AddElement(elC, nelements = 18);
	PFT->AddElement(elH, nelements = 14);

	G4Material* POPOP = new G4Material("MPOPOP", density = 1.5 * g / cm3, ncomponents = 4);
	POPOP->AddElement(elC, nelements = 24);
	POPOP->AddElement(elH, nelements = 16);
	POPOP->AddElement(elN, nelements = 2);
	POPOP->AddElement(elO, nelements = 2);

	G4Material* STR = new G4Material("MSTR", density = 1.2 * g / cm3, ncomponents = 3);
	STR->AddMaterial(PS, fractionmass = 98.46 * perCent);
	STR->AddMaterial(PFT, fractionmass = 1.5 * perCent);
	STR->AddMaterial(PS, fractionmass = 0.04 * perCent);

	//Optical glue
	G4Material* Glue = new G4Material("MGlue", density = 1.02 * g / cm3, ncomponents = 4);
	Glue->AddElement(H, nelements = 108);
	Glue->AddElement(C, nelements = 65);
	Glue->AddElement(N, nelements = 20);
	Glue->AddElement(O, nelements = 7);


	/*	OPTICAL PROPERTIES	*/


		//Scintillator optical properties
	const G4int nEntries = 60;
	G4double PhotonEnergy[nEntries] = { 2.3, 2.31525, 2.33051, 2.34576, 2.36102, 2.37627, 2.39153, 2.40678, 2.42203, 2.43729, 2.45254, 2.4678, 2.48305, 2.49831, 2.51356,
		 2.52881, 2.54407, 2.55932, 2.57458, 2.58983, 2.60508, 2.62034, 2.63559, 2.65085, 2.6661, 2.68136, 2.69661, 2.71186, 2.72712, 2.74237,
		 2.75763, 2.77288, 2.78814, 2.80339, 2.81864, 2.8339, 2.84915, 2.86441, 2.87966, 2.89492, 2.91017, 2.92542, 2.94068, 2.95593, 2.97119,
		 2.98644, 3.00169, 3.01695, 3.0322, 3.04746, 3.06271, 3.07797, 3.09322, 3.10847, 3.12373, 3.13898, 3.15424, 3.16949, 3.18475, 3.2 };
	G4double RefractiveScin[nEntries];
	G4double AbsLengthScin[nEntries];
	G4double SpIzlStr[nEntries] = { 0, 0, 0.04304, 0.09311, 0.14318, 0.19325, 0.24331, 0.29338, 0.34345, 0.39352, 0.44359, 0.49365, 0.54372, 0.59379, 0.65703,
		 0.72516, 0.7829, 0.85487, 0.93619, 1.0156, 1.10002, 1.19322, 1.29936, 1.41172, 1.53233, 1.65876, 1.79893, 1.98186, 2.18771, 2.4366,
		 2.78324, 3.0698, 3.27276, 3.39218, 3.46918, 3.4941, 3.52619, 3.60856, 3.88683, 4.28688, 4.71702, 4.93565, 4.80817, 4.56821, 4.23367,
		 3.56117, 2.30136, 1.47323, 1.10353, 0.84005, 0.61903, 0.46259, 0.35545, 0.2483, 0.14115, 0.034, 0, 0, 0, 0 };

	G4int j;

	for (j = 0; j < nEntries; j++)
	{
		RefractiveScin[j] = 1.58;
		AbsLengthScin[j] = 1. * m;
		PhotonEnergy[j] = PhotonEnergy[j] * eV;
	}

	G4MaterialPropertiesTable* ScintillatorProperties = new G4MaterialPropertiesTable();
	ScintillatorProperties->AddProperty("RINDEX", PhotonEnergy, RefractiveScin, nEntries);
	ScintillatorProperties->AddProperty("ABSLENGTH", PhotonEnergy, AbsLengthScin, nEntries);
	ScintillatorProperties->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergy, SpIzlStr, nEntries);
	ScintillatorProperties->AddProperty("SCINTILLATIONCOMPONENT2", PhotonEnergy, SpIzlStr, nEntries);
	ScintillatorProperties->AddConstProperty("RESOLUTIONSCALE", 1.0);
	ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD", 1200 / MeV); // 12000
	ScintillatorProperties->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.4 * ns);
	ScintillatorProperties->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 5 * ns);
	ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
	ScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
	STR->SetMaterialPropertiesTable(ScintillatorProperties);
	STR->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);


	G4double EnergyOpt[10] = { 1.9 * eV, 2.2 * eV, 2.3 * eV, 2.4 * eV, 2.56 * eV, 2.66 * eV, 2.68 * eV, 3.69 * eV, 3.7 * eV, 4.0 * eV };
	G4double AbsLenOpt[10] = { 5.0 * m, 5.0 * m, 5.0 * m, 5.0 * m, 5.0 * m, 5.0 * m, 0.1 * mm, 0.1 * mm, 5.0 * m, 5.0 * m };
	G4double SpIzlOpt[10] = { 0.001, 0.05, 0.25, 0.7, 1., 1., 0., 0., 0., 0. };
	G4double RindexOptCore[10] = { 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59 };
	G4double RindexOptCov[10] = { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49 };

	//Core optical properties
	G4MaterialPropertiesTable* OptCore = new G4MaterialPropertiesTable();
	OptCore->AddProperty("RINDEX", EnergyOpt, RindexOptCore, 10);
	OptCore->AddProperty("WLSABSLENGTH", EnergyOpt, AbsLenOpt, 10);
	OptCore->AddProperty("WLSCOMPONENT", EnergyOpt, SpIzlOpt, 10);
	OptCore->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
	PS->SetMaterialPropertiesTable(OptCore);

	//Cover optical properties
	G4MaterialPropertiesTable* OptCov = new G4MaterialPropertiesTable();
	OptCov->AddProperty("RINDEX", EnergyOpt, RindexOptCov, 10);
	OptCov->AddProperty("WLSABSLENGTH", EnergyOpt, AbsLenOpt, 10);
	OptCov->AddProperty("WLSCOMPONENT", EnergyOpt, SpIzlOpt, 10);
	OptCov->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
	PMMA->SetMaterialPropertiesTable(OptCov);

	//Air optical properties
	G4double EnergyAir[2] = { 1.9 * eV, 4.0 * eV };
	G4double AbsLenAir[2] = { 5.0 * m,  5.0 * m };
	G4double RindAir[2] = { 1.0002926, 1.0002926 };

	G4MaterialPropertiesTable* AirPT = new G4MaterialPropertiesTable();
	AirPT->AddProperty("RINDEX", EnergyAir, RindAir, 2);
	AirPT->AddProperty("ABSLENGTH", EnergyAir, AbsLenAir, 2);
	Air->SetMaterialPropertiesTable(AirPT);



	/*	DETECTOR	*/


	G4bool checkOverlaps = true;

	//World
	G4double world_sizeX = 5 * m;
	G4double world_sizeY = 5 * m;
	G4double world_sizeZ = 5 * m;

	G4Box* solidWorld = new G4Box("World_s", 0.5 * world_sizeX, 0.5 * world_sizeY, 0.5 * world_sizeZ);
	G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Air, "World_l");
	G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

	G4double GapH = 0.1 * mm;   //Horizontal gap
	G4double GapV = 0.1 * mm;	//Vertical gap
	G4double GapFS = 0.1 * mm;	//Gap between frame and strip
	G4double GapSh = 0.1 * mm;	//Gap between shells
	G4double disloc = 5 * mm;	//Dislocation of upper layer of coordinate plate

	//Strip parameters
	G4double StrLength = 1000 * mm;
	G4double StrWidth = 10 * mm;
	G4double StrHeight = 7 * mm;

	G4double TdlrLength = StrLength;
	G4double TdlrWidth = StrWidth + 0.2 * mm;
	G4double TdlrHeight = StrHeight + 0.2 * mm;

	//Optical glue in strip
	G4double GlueLength = StrLength;
	G4double GlueWidth = 1.5 * mm;
	G4double GlueHeight = 2 * mm;

	//Optical fiber
	G4double OptRad = 0.485 * mm;
	G4double OptHeight = StrLength;
	G4double CovThickness = 0.03 * mm;
	G4RotationMatrix* OptRot = new G4RotationMatrix;
	OptRot->rotateY(90. * deg);

	//Steel shell and air hollow inside it
	G4double HollowLength = TdlrLength + 2 * GapH;
	G4double HollowWidth = NRows * (TdlrWidth + GapH) + disloc;
	G4double HollowHeight = 2 * TdlrHeight + GapV + 2 * GapFS;

	G4double ShellThickness = 1 * mm;
	G4double ShellLength = HollowLength + 2 * ShellThickness;
	G4double ShellWidth = HollowWidth + 2 * ShellThickness;
	G4double ShellHeight = HollowHeight + 2 * ShellThickness;

	//Variables for creating copies
	const G4int NRows = 96, NLvls = 2, NCoord = 3;
	G4int row, level, coord;
	G4int StrNCopy = 0, TdlrNCopy = 0, OptCoreNCopy = 0, OptCovNCopy = 0, GlueNCopy = 0, ShellNCopy = 0, HollowNCopy = 0;

	G4double XStr = 0 * mm;
	G4double YStr = -(0.5 * NRows * TdlrWidth + (0.5 * NRows - 1) * GapH + 0.5 * GapH) + 0.5 * TdlrWidth;
	G4double ZStr = 0.5 * TdlrHeight;

	G4double XGl = 0 * mm;
	G4double YGl = 0 * mm;
	G4double ZGl = 0.5 * (StrHeight - GlueHeight);

	G4double XOpt = 0 * mm;
	G4double YOpt = 0 * mm;
	G4double ZOpt = OptRad - 0.5 * GlueHeight;

	G4double Str_X, Str_Y, Str_Z, Gl_X, Gl_Y, Gl_Z, Opt_X, Opt_Y, Opt_Z;


	//Volumes
	G4Box* solidRotVolume = { NULL }, * solidShell = { NULL }, * solidHollow = { NULL }, * solidTdlr[NRows][NLvls] = {NULL}, * solidStrip[NRows][NLvls] = {NULL}, * solidGlue[NRows][NLvls] = {NULL};
	G4Tubs* solidCore[NRows][NLvls] = { NULL }, * solidCov[NRows][NLvls] = { NULL };
	G4LogicalVolume* logicRotVolume = { NULL }, * logicShell = { NULL }, * logicHollow = { NULL }, * logicTdlr[NRows][NLvls] = { NULL }, * logicStrip[NRows][NLvls] = { NULL }, * logicGlue[NRows][NLvls] = { NULL },
		* logicCov[NRows][NLvls] = { NULL }, * logicCore[NRows][NLvls] = { NULL };
	G4VPhysicalVolume* physRotVolume = { NULL }, * physShell = { NULL }, * physHollow = { NULL }, * physTdlr[NRows][NLvls] = { NULL }, * physStrip[NRows][NLvls] = { NULL }, * physGlue[NRows][NLvls] = { NULL },
		* physCov[NRows][NLvls] = { NULL }, * physCore[NRows][NLvls] = { NULL };

	//Rotating volume
	solidRotVolume = new G4Box("RotVol_s", 0.5 * m, 0.5 * m, 0.5 * m);
	logicRotVolume = new G4LogicalVolume(solidRotVolume, Air, "RotVol_l");
	physRotVolume = new G4PVPlacement(0, G4ThreeVector(), logicRotVolume, "ROTATION_VOLUME", logicWorld, false, 0, checkOverlaps);

	//Steel shell and air hollow inside it
	solidShell = new G4Box("shell_s", 0.5 * ShellLength, 0.5 * ShellWidth, 0.5 * ShellHeight);
	logicShell = new G4LogicalVolume(solidShell, FeMaterial, "shell_l");
	physShell = new G4PVPlacement(0, G4ThreeVector(0., 0., - 0.5 * m + 0.5 * ShellHeight), logicShell, "SHELL", logicRotVolume, false, 0, checkOverlaps);

	solidHollow = new G4Box("hollow_s", 0.5 * HollowLength, 0.5 * HollowWidth, 0.5 * HollowHeight);
	logicHollow = new G4LogicalVolume(solidRotVolume, Air, "hollow_l");
	physHollow = new G4PVPlacement(0, G4ThreeVector(), logicHollow, "HOLLOW", logicShell, false, 0, checkOverlaps);

	Str_X = XStr, Str_Y = YStr, Str_Z = ZStr, Gl_X = XGl, Gl_Y = YGl, Gl_Z = ZGl, Opt_X = XOpt, Opt_Y = YOpt, Opt_Z = ZOpt;
	G4double distance = TdlrWidth + GapH;

	for (row = 0; row < NRows; row++)
	{
		solidTdlr[row] = new G4Box("tdlr_s", 0.5 * TdlrLength, 0.5 * TdlrWidth, 0.5 * TdlrHeight);
		logicTdlr[row] = new G4LogicalVolume(solidTdlr[row], PVF, "tdlr_l");
		physTdlr[row] = new G4PVPlacement(0, G4ThreeVector(Str_X, Str_Y, Str_Z), logicTdlr[row], "TEDLAR", logicWorld, false, TdlrNCopy, checkOverlaps);

		logicStrip[row] = new G4LogicalVolume(solidStr, STR, "strip_l");
		physStrip[row] = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicStrip[row], "STRIP", logicTdlr[row], false, StrNCopy, checkOverlaps);

		solidCov[row] = new G4Tubs("Cov_s", 0, OptRad, 0.5 * OptHeight, 0. * deg, 360. * deg);
		logicCov[row] = new G4LogicalVolume(solidCov[row], PMMA, "Cov_l");
		physCov[row] = new G4PVPlacement(OptRot, G4ThreeVector(Opt_X, Opt_Y, Opt_Z), logicCov[row], "COVER", logicTdlr[row], false, OptCovNCopy, checkOverlaps);

		solidCore[row] = new G4Tubs("core_s", 0, OptRad - CovThickness, 0.5 * OptHeight, 0. * deg, 360. * deg);
		logicCore[row] = new G4LogicalVolume(solidCore[row], PS, "core_l");
		physCore[row] = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicCore[row], "CORE", logicCov[row], false, OptCoreNCopy, checkOverlaps);

		Str_Y += distance;

		StrNCopy++;
		TdlrNCopy++;
		OptCoreNCopy++;
		OptCovNCopy++;
	}

	//Making world invisible
	auto UniverseVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	UniverseVisAtt->SetVisibility(true);
	UniverseVisAtt->SetForceWireframe(true);
	logicWorld->SetVisAttributes(UniverseVisAtt);
	logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......