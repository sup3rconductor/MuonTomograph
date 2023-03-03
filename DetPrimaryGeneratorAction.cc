#include "DetPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

G4double theta, phi, ux, uy, uz, E0;
G4double x_rand, y_rand, z, x_up, y_up, z_up;
char fpartname[7];
G4int fEvent, fpartnum;
G4double ftheta, fphi, fEkin;


//extern G4double ShellLength, ShellWidth, ShellHeight, ShellThickness, GapH, GapV, GapFP, ScrHeight;
extern G4double Z0const;
extern FILE* rdata;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetPrimaryGeneratorAction::DetPrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(),
	fParticleGun(0)
{
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);

	// default particle kinematic
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "mu+");
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
	fParticleGun->SetParticleEnergy(4. * GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetPrimaryGeneratorAction::~DetPrimaryGeneratorAction()
{
	delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//Maximal and minimal X coordinate for modelling light gathering
	//G4double x_min = -ShellLength + ShellThickness + GapFP;
	//G4double x_max = 0 * mm;
	G4double x_min = (- 200 * 2 - 100 - 0.11 * 2) * mm;
	G4double x_max = 0 * mm;

	//Maximal and minimal Y coordinate for modelling light gathering
	//G4double y_min = -ShellWidth + ShellThickness + GapFP;
	//G4double y_max = ShellWidth - ShellThickness - GapFP;
	G4double y_min = (-200 * 2 - 100 - 0.11 * 2) * mm;
	G4double y_max = (200 * 2 + 100 + 0.11 * 2) * mm;

	//Z coordinate of 5th lvl for modelling light gathering
	//G4double z_mid = ShellHeight - GapFP - 4 * (ScrHeight + GapV) - 0.5 * ScrHeight;
	G4double z_mid = - 0.5 * (0.11 + 5) * mm;

	//Random coordinares
	x_rand = x_min + (x_max - x_min) * G4UniformRand();
	y_rand = y_min + (y_max - y_min) * G4UniformRand();
	z = z_mid;

	//teta = acos(pow((1+G4UniformRand()*(pow(cos(15*twopi/360),4.2)-1)), 1/4.2));
	//phi = twopi*G4UniformRand();

	//Getting particles data from file
	fscanf(rdata, "%d\t%s\t%d\t%lf\t%lf\t%lf\n", &fEvent, &fpartname, &fpartnum, &ftheta, &fphi, &fEkin);
	G4cout << "Number of event: " << fEvent << "\tParticle type: " << fpartnum << "\tTheta angle: "<< ftheta << "\tPhi angle: " << fphi << "\tKinetic energy: " << fEkin << G4endl;

	//Angles in radians
	theta = ftheta * pi / 180.0;
	phi = fphi * pi / 180.0;

	//Launching positively charged muons
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle(fpartname);
	fParticleGun->SetParticleDefinition(particle);

	//Launching position
	x_up = x_rand + 150 * cm * sin(theta) * cos(phi);
	y_up = y_rand + 150 * cm * sin(theta) * sin(phi);
	z_up = z + 150 * cm * cos(theta); 

	fParticleGun->SetParticlePosition(G4ThreeVector(x_up, y_up, z_up));

	//Momentum direction
	ux = -sin(theta) * cos(phi);
	uy = -sin(theta) * sin(phi);
	uz = -cos(theta);

	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

	fParticleGun->SetParticleEnergy(fEkin * GeV);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

