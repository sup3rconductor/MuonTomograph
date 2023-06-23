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

	//Silicium
	G4Material* SiMaterial = new G4Material("MSilicium", z = 14., a = 26.09 * g / mole, density = 2.33 * g / cm3);

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

	G4Material* STR = new G4Material("MSTR", density = 1.06 * g / cm3, ncomponents = 3);
	STR->AddMaterial(PS, fractionmass = 98.46 * perCent);
	STR->AddMaterial(PFT, fractionmass = 1.5 * perCent);
	STR->AddMaterial(POPOP, fractionmass = 0.04 * perCent);

	//Expanded polystyrene for strip cover
	G4Material* PScov = new G4Material("McovStrip", density = 1.04 * g / cm3, nelements = 2);
	PScov->AddElement(elC, nelements = 8);
	PScov->AddElement(elH, nelements = 8);

	//Optical glue
	G4Material* Glue = new G4Material("MGlue", density = 1.02 * g / cm3, ncomponents = 4);
	Glue->AddElement(elH, nelements = 108);
	Glue->AddElement(elC, nelements = 65);
	Glue->AddElement(elN, nelements = 20);
	Glue->AddElement(elO, nelements = 7);


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

	//Glue optical properties
	G4double PhotonEnergyGlue[2] = { 2.175 * eV, 3.542 * eV };
	G4double RindexGlue[2] = { 1.56, 1.56 };
	G4double AbsLengthGlue[2] = { 5 * m, 5 * m };

	G4MaterialPropertiesTable* OptGlue = new G4MaterialPropertiesTable();
	OptGlue->AddProperty("RINDEX", PhotonEnergyGlue, RindexGlue, 2);
	OptGlue->AddProperty("WLSABSLENGTH", PhotonEnergyGlue, AbsLengthGlue, 2);
	Glue->SetMaterialPropertiesTable(OptGlue);

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

	G4double GapH = 0.2 * mm;   //Horizontal gap
	G4double GapV = 0.2 * mm;	//Vertical gap
	G4double GapFS = 0.2 * mm;	//Gap between frame and strip
	G4double GapSh = 5 * mm;	//Gap between shells

	//Variables for creating copies
	const G4int NRows = 96, NLvls = 2, NPlates = 2, NCoord = 3;
	G4int row, level, plate, coord;
	G4int StrNCopy = 0, TdlrNCopy = 0, ExPSNCopy = 0, SiPMNCopy = 0, OptCoreNCopy = 0, OptCovNCopy = 0, GlueNCopy = 0, ShellNCopy = 0, HollowNCopy = 0;

	//Strip parameters
	G4double StrLength = 1000 * mm;
	G4double StrWidth = 10 * mm;
	G4double StrHeight = 7 * mm;

	G4double ExPSLength = StrLength + 0.02 * mm;
	G4double ExPSWidth = StrWidth + 0.02 * mm;
	G4double ExPSHeight = StrHeight + 0.02 * mm;

	G4double TdlrLength = StrLength + 0.2 * mm;
	G4double TdlrWidth = StrWidth + 0.2 * mm;
	G4double TdlrHeight = StrHeight + 0.2 * mm;

	G4double disloc = 0.5 * (TdlrWidth + GapH);	//Dislocation of upper layer of coordinate plate

	//Optical glue in strip
	G4double GlueLength = StrLength;
	G4double GlueWidth = 1.5 * mm;
	G4double GlueHeight = 2 * mm;

	//Optical fiber
	G4double OptRad = 0.485 * mm;
	G4double OptHeight = StrLength;
	G4double CovThickness = 0.02 * mm;
	G4RotationMatrix* OptRot = new G4RotationMatrix;
	OptRot->rotateY(90. * deg);

	//SiPM
	G4double TdlrThick = 0.1 * mm;				//Variable for placing SiPMs

	G4double SiPMLength = 0.01 * mm;
	G4double SiPMWidth = 1.3 * mm;
	G4double SiPMHeight = 1.3 * mm;

	//Rotating volume, steel shell and air hollow inside it
	G4double HollowLength = TdlrLength + 2 * GapFS;
	G4double HollowWidth = NRows * TdlrWidth + (NRows - 1) * GapH + 2 * GapFS + disloc;
	G4double HollowHeight = 2 * TdlrHeight + GapV + 2 * GapFS;

	G4double ShellThickness = 1 * mm;
	G4double ShellLength = HollowLength + 2 * ShellThickness;
	G4double ShellWidth = HollowWidth + 2 * ShellThickness;
	G4double ShellHeight = HollowHeight + 2 * ShellThickness;

	G4RotationMatrix* ShellRot[NPlates][NCoord] = { NULL };

	G4double RotVolLength = 1.05 * m;
	G4double RotVolWidth = 1.05 * m;
	G4double RotVolHeight = 1 * m;

	//Coordinates
	G4double XStr = 0 * mm;
	G4double YStr = -0.5 * HollowWidth + 0.5 * TdlrWidth + GapFS;
	G4double ZStr = -0.5 * HollowHeight + 0.5 * TdlrHeight + GapFS;

	G4double XGl = 0 * mm;
	G4double YGl = 0 * mm;
	G4double ZGl = 0.5 * (StrHeight - GlueHeight);

	G4double XOpt = 0 * mm;
	G4double YOpt = 0 * mm;
	G4double ZOpt = OptRad - 0.5 * GlueHeight;

	G4double XSiPM = 0.5 * (StrLength + SiPMLength);
	G4double YSiPM = 0 * mm;
	G4double ZSiPM = 0.5 * StrHeight - GlueHeight + OptRad;

	G4double XSh = 0 * mm;
	G4double YSh = 0 * mm;
	G4double ZSh = -0.5 * (RotVolHeight - ShellHeight);

	G4double Sh_X, Sh_Y, Sh_Z,
		Str_X, Str_Y, Str_Z,
		Gl_X, Gl_Y, Gl_Z,
		Opt_X, Opt_Y, Opt_Z,
		SiPM_X, SiPM_Y, SiPM_Z;


	//Volumes
	G4Box* solidRotVolume = { NULL }, * solidShell[NPlates][NCoord] = { NULL }, * solidHollow[NPlates][NCoord] = { NULL }, * solidTdlr[NRows][NLvls][NPlates][NCoord] = { NULL },
		* solidExPS[NRows][NLvls][NPlates][NCoord], * solidStrip[NRows][NLvls][NPlates][NCoord] = { NULL }, * solidGlue[NRows][NLvls][NPlates][NCoord] = { NULL }, * solidSiPM[NRows][NLvls][NPlates][NCoord];

	G4Tubs* solidCore[NRows][NLvls][NPlates][NCoord] = { NULL }, * solidCov[NRows][NLvls][NPlates][NCoord] = { NULL };

	G4LogicalVolume* logicRotVolume = { NULL }, * logicShell[NPlates][NCoord] = { NULL }, * logicHollow[NPlates][NCoord] = { NULL }, * logicTdlr[NRows][NLvls][NPlates][NCoord] = { NULL },
		* logicExPS[NRows][NLvls][NPlates][NCoord], * logicStrip[NRows][NLvls][NPlates][NCoord] = { NULL }, * logicGlue[NRows][NLvls][NPlates][NCoord] = { NULL }, * logicCov[NRows][NLvls][NPlates][NCoord] = { NULL },
		* logicCore[NRows][NLvls][NPlates][NCoord] = { NULL }, * logicSiPM[NRows][NLvls][NPlates][NCoord];

	G4VPhysicalVolume* physRotVolume = { NULL }, * physShell[NPlates][NCoord] = { NULL }, * physHollow[NPlates][NCoord] = { NULL }, * physTdlr[NRows][NLvls][NPlates][NCoord] = { NULL },
		* physExPS[NRows][NLvls][NPlates][NCoord], * physStrip[NRows][NLvls][NPlates][NCoord] = { NULL }, * physGlue[NRows][NLvls][NPlates][NCoord] = { NULL }, * physCov[NRows][NLvls][NPlates][NCoord] = { NULL },
		* physCore[NRows][NLvls][NPlates][NCoord] = { NULL }, * physSiPM[NRows][NLvls][NPlates][NCoord];

	//Rotating volume
	solidRotVolume = new G4Box("RotVol_s", 0.5 * RotVolLength, 0.5 * RotVolWidth, 0.5 * RotVolHeight);
	logicRotVolume = new G4LogicalVolume(solidRotVolume, Air, "RotVol_l");
	physRotVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicRotVolume, "ROTATION_VOLUME", logicWorld, false, 0, checkOverlaps);

	Sh_X = XSh, Sh_Y = YSh, Sh_Z = ZSh,
		Str_X = XStr, Str_Y = YStr, Str_Z = ZStr, Gl_X = XGl, Gl_Y = YGl, Gl_Z = ZGl, Opt_X = XOpt, Opt_Y = YOpt, Opt_Z = ZOpt, SiPM_X = XSiPM, SiPM_Y = YSiPM, SiPM_Z = ZSiPM;
	G4double distance = TdlrWidth + GapH;
	G4double interval = 0.5 * RotVolHeight - 2.5 * GapSh - 3 * ShellHeight;

	//Tomograph has 3 coordinate planes
	for (coord = 0; coord < NCoord; coord++)
	{
		for (plate = 0; plate < NPlates; plate++)
		{
			if (plate == 1)
			{
				ShellRot[plate][coord] = new G4RotationMatrix;
				ShellRot[plate][coord]->rotateZ(90. * deg);
			}

			else
			{
				ShellRot[plate][coord] = new G4RotationMatrix;
				ShellRot[plate][coord]->rotateZ(0. * deg);
			}

			solidShell[plate][coord] = new G4Box("shell_s", 0.5 * ShellLength, 0.5 * ShellWidth, 0.5 * ShellHeight);
			logicShell[plate][coord] = new G4LogicalVolume(solidShell[plate][coord], FeMaterial, "shell_l");
			physShell[plate][coord] = new G4PVPlacement(ShellRot[plate][coord], G4ThreeVector(Sh_X, Sh_Y, Sh_Z), logicShell[plate][coord], "SHELL", logicRotVolume, false, ShellNCopy, checkOverlaps);

			solidHollow[plate][coord] = new G4Box("hollow_s", 0.5 * HollowLength, 0.5 * HollowWidth, 0.5 * HollowHeight);
			logicHollow[plate][coord] = new G4LogicalVolume(solidHollow[plate][coord], Air, "hollow_l");
			physHollow[plate][coord] = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicHollow[plate][coord], "HOLLOW", logicShell[plate][coord], false, HollowNCopy, checkOverlaps);

			ShellNCopy++;
			HollowNCopy++;

			//Steel shell contains 2 scintilation strips inside
			for (level = 0; level < NLvls; level++)
			{
				//Each strip layer contains 96 scintillation strips
				for (row = 0; row < NRows; row++)
				{
					solidTdlr[row][level][plate][coord] = new G4Box("tdlr_s", 0.5 * TdlrLength, 0.5 * TdlrWidth, 0.5 * TdlrHeight);
					logicTdlr[row][level][plate][coord] = new G4LogicalVolume(solidTdlr[row][level][plate][coord], PVF, "tdlr_l");
					physTdlr[row][level][plate][coord] = new G4PVPlacement(0, G4ThreeVector(Str_X, Str_Y, Str_Z), logicTdlr[row][level][plate][coord], "TEDLAR", logicHollow[plate][coord], false, TdlrNCopy, checkOverlaps);

					solidExPS[row][level][plate][coord] = new G4Box("exPS_s", 0.5 * ExPSLength, 0.5 * ExPSWidth, 0.5 * ExPSHeight);
					logicExPS[row][level][plate][coord] = new G4LogicalVolume(solidExPS[row][level][plate][coord], PScov, "exPS_l");
					physExPS[row][level][plate][coord] = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicExPS[row][level][plate][coord], "EXPS", logicTdlr[row][level][plate][coord], false, ExPSNCopy, checkOverlaps);

					solidStrip[row][level][plate][coord] = new G4Box("strip_s", 0.5 * StrLength, 0.5 * StrWidth, 0.5 * StrHeight);
					logicStrip[row][level][plate][coord] = new G4LogicalVolume(solidStrip[row][level][plate][coord], STR, "strip_l");
					physStrip[row][level][plate][coord] = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicStrip[row][level][plate][coord], "STRIP", logicExPS[row][level][plate][coord], false, StrNCopy, checkOverlaps);

					solidGlue[row][level][plate][coord] = new G4Box("glue_s", 0.5 * GlueLength, 0.5 * GlueWidth, 0.5 * GlueHeight);
					logicGlue[row][level][plate][coord] = new G4LogicalVolume(solidGlue[row][level][plate][coord], Glue, "glue_l");
					physGlue[row][level][plate][coord] = new G4PVPlacement(0, G4ThreeVector(Gl_X, Gl_Y, Gl_Z), logicGlue[row][level][plate][coord], "GLUE", logicStrip[row][level][plate][coord], false, GlueNCopy, checkOverlaps);

					solidCov[row][level][plate][coord] = new G4Tubs("cov_s", 0, OptRad, 0.5 * OptHeight, 0. * deg, 360. * deg);
					logicCov[row][level][plate][coord] = new G4LogicalVolume(solidCov[row][level][plate][coord], PMMA, "cov_l");
					physCov[row][level][plate][coord] = new G4PVPlacement(OptRot, G4ThreeVector(Opt_X, Opt_Y, Opt_Z), logicCov[row][level][plate][coord], "COVER", logicGlue[row][level][plate][coord], false, OptCovNCopy, checkOverlaps);

					solidCore[row][level][plate][coord] = new G4Tubs("core_s", 0, OptRad - CovThickness, 0.5 * OptHeight, 0. * deg, 360. * deg);
					logicCore[row][level][plate][coord] = new G4LogicalVolume(solidCore[row][level][plate][coord], PS, "core_l");
					physCore[row][level][plate][coord] = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicCore[row][level][plate][coord], "CORE", logicCov[row][level][plate][coord], false, OptCoreNCopy, checkOverlaps);

					solidSiPM[row][level][plate][coord] = new G4Box("sipm_s", 0.5 * SiPMLength, 0.5 * SiPMWidth, 0.5 * SiPMHeight);
					logicSiPM[row][level][plate][coord] = new G4LogicalVolume(solidSiPM[row][level][plate][coord], SiMaterial, "sipm_l");
					physSiPM[row][level][plate][coord] = new G4PVPlacement(0, G4ThreeVector(SiPM_X, SiPM_Y, SiPM_Z), logicSiPM[row][level][plate][coord], "SIPM", logicExPS[row][level][plate][coord], false, SiPMNCopy, checkOverlaps);

					Str_Y += distance;

					StrNCopy++;
					TdlrNCopy++;
					ExPSNCopy++;
					GlueNCopy++;
					SiPMNCopy++;
					OptCoreNCopy++;
					OptCovNCopy++;
				}

				//Next layer of strips + its dislocation related to lower layer
				Str_Y = YStr + disloc;
				Str_Z += (TdlrHeight + GapV);
			}

			//Next shell
			Sh_Z += (ShellHeight + GapSh);

			//Returning to start position of strip and SiPM
			Str_X = XStr;
			Str_Y = YStr;
			Str_Z = ZStr;
		}

		//Moving to next coordinate plane
		Sh_Z += interval;
	}


	//Border Strip - expanded polystyrene: diffuse reflection
	G4double reflectivity_exps[2] = { 0.95, 0.95 };
	G4double PhotonEnergyExPS[2] = { 1.9 * eV, 4.0 * eV };

	G4OpticalSurface* OptPovExPS = new G4OpticalSurface("PovExpPolystyrene");
	OptPovExPS->SetType(dielectric_dielectric);
	OptPovExPS->SetFinish(groundfrontpainted);
	OptPovExPS->SetModel(unified);

	G4MaterialPropertiesTable* PovExPSPT = new G4MaterialPropertiesTable();
	PovExPSPT->AddProperty("REFLECTIVITY", PhotonEnergyExPS, reflectivity_exps, 2);
	OptPovExPS->SetMaterialPropertiesTable(PovExPSPT);

	G4int exps, lvl, plt, cord;

	G4LogicalBorderSurface* ExpPolystyreneSurface[NRows][NLvls][NPlates][NCoord] = { NULL };
	for (cord = 0; cord < NCoord; cord++)
	{
		for (plt = 0; plt < NPlates; plt++)
		{
			for (lvl = 0; lvl < NLvls; lvl++)
			{
				for (exps = 0; exps < NRows; exps++)
				{
					ExpPolystyreneSurface[exps][lvl][plt][cord] = new G4LogicalBorderSurface("ExPSStripSurface", physStrip[exps][lvl][plt][cord], physExPS[exps][lvl][plt][cord], OptPovExPS);
				}
			}
		}
	}

	//Border expanded polystyrene - Tedlar: absorbtion
	G4double reflectivity_tdlr[2] = { 0., 0. };
	G4double PhotonEnergyTedlar[2] = { 1.9 * eV, 4.0 * eV };

	G4OpticalSurface* OptPovTdlr = new G4OpticalSurface("PovTedlar");
	OptPovTdlr->SetType(dielectric_dielectric);
	OptPovTdlr->SetFinish(groundfrontpainted);
	OptPovTdlr->SetModel(unified);

	G4MaterialPropertiesTable* PovTdlrPT = new G4MaterialPropertiesTable();
	PovTdlrPT->AddProperty("REFLECTIVITY", PhotonEnergyTedlar, reflectivity_tdlr, 2);
	OptPovTdlr->SetMaterialPropertiesTable(PovTdlrPT);

	G4int tdlr;

	G4LogicalBorderSurface* TedlarSurface[NRows][NLvls][NPlates][NCoord] = { NULL };
	for (cord = 0; cord < NCoord; cord++)
	{
		for (plt = 0; plt < NPlates; plt++)
		{
			for (lvl = 0; lvl < NLvls; lvl++)
			{
				for (tdlr = 0; tdlr < NRows; tdlr++)
				{
					TedlarSurface[tdlr][lvl][plt][cord] = new G4LogicalBorderSurface("TedlarStripSurface", physExPS[tdlr][lvl][plt][cord], physTdlr[tdlr][lvl][plt][cord], OptPovTdlr);
				}
			}
		}
	}

	//Border optical fiber - SiPM: mirror reflection
	G4double reflectivity_SiPM[2] = { 0.99, 0.99 };
	G4double PhotonEnergySiPM[2] = { 1.9 * eV, 4.0 * eV };

	G4OpticalSurface* OptPovSiPM = new G4OpticalSurface("PovSiPM");
	OptPovSiPM->SetType(dielectric_metal);
	OptPovSiPM->SetFinish(ground);
	OptPovSiPM->SetModel(glisur);

	G4MaterialPropertiesTable* PovSiPM_PT = new G4MaterialPropertiesTable();
	PovSiPM_PT->AddProperty("REFLECTIVITY", PhotonEnergySiPM, reflectivity_SiPM, 2);
	OptPovSiPM->SetMaterialPropertiesTable(PovSiPM_PT);

	G4int spm;

	G4LogicalBorderSurface* CoreSiPMSurface[NRows][NLvls][NPlates][NCoord] = { NULL }, * CovSiPMSurface[NRows][NLvls][NPlates][NCoord] = { NULL };
	for (cord = 0; cord < NCoord; cord++)
	{
		for (plt = 0; plt < NPlates; plt++)
		{
			for (lvl = 0; lvl < NLvls; lvl++)
			{
				for (spm = 0; spm < NRows; spm++)
				{
					CoreSiPMSurface[spm][lvl][plt][cord] = new G4LogicalBorderSurface("CoreSiPMSurface", physCore[spm][lvl][plt][cord], physSiPM[spm][lvl][plt][cord], OptPovSiPM);
					CovSiPMSurface[spm][lvl][plt][cord] = new G4LogicalBorderSurface("CoverSiPMSurface", physCov[spm][lvl][plt][cord], physSiPM[spm][lvl][plt][cord], OptPovSiPM);
				}
			}
		}
	}



	/* VISUAL ATTRIBUTES */

	//Colours
	G4Colour grey(0.5, 0.5, 0.5);
	G4Colour black(0.0, 0.0, 0.0);
	G4Colour red(1.0, 0.0, 0.0);
	G4Colour green(0.0, 1.0, 0.0);
	G4Colour blue(0.0, 0.0, 1.0);
	G4Colour cyan(0.0, 1.0, 1.0);
	G4Colour magenta(1.0, 0.0, 1.0);
	G4Colour yellow(1.0, 1.0, 0.0);

	//Making world invisible
	auto UniverseVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	UniverseVisAtt->SetVisibility(true);
	UniverseVisAtt->SetForceWireframe(true);
	logicWorld->SetVisAttributes(UniverseVisAtt);
	logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

	//Making rotation volume invisible
	auto RotVolVisAtt = new G4VisAttributes(G4Colour(1., 1., 1.));
	RotVolVisAtt->SetVisibility(true);
	RotVolVisAtt->SetForceSolid(true);
	logicRotVolume->SetVisAttributes(RotVolVisAtt);
	logicRotVolume->SetVisAttributes(G4VisAttributes::GetInvisible());

	//SiPM visual attributes
	G4VisAttributes* SiPMVisAtt = new G4VisAttributes(magenta);
	SiPMVisAtt->SetVisibility(true);
	SiPMVisAtt->SetForceSolid(true);

	G4int i, m, k, l;
	for (i = 0; i < NCoord; i++)
	{
		for (m = 0; m < NPlates; m++)
		{
			for (k = 0; k < NLvls; k++)
			{
				for (l = 0; l < NRows; l++)
				{
					logicSiPM[l][k][m][i]->SetVisAttributes(SiPMVisAtt);
				}
			}
		}
	}

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......