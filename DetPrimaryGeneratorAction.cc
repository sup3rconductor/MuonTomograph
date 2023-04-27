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
G4double x, y, z, x0, yy0, z0;
char fpartname[7];
G4int fEvent, fpartnum;
G4double ftheta, fphi, fEkin;


//extern G4double ShellLength, ShellWidth, ShellHeight, ShellThickness, GapH, GapV, GapFP, ScrHeight;
extern G4double Z0const;

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
	//Particle coordinates
	x = -400 + 800 * G4UniformRand();
	y = 400 * G4UniformRand();
	z = Z0const;

	//Scanning particle data from file
	//fscanf(rdata, "%d\t%s\t%d\t%lf\t%lf\t%lf\n", &fEvent, &fpartname, &fpartnum, &ftheta, &fphi, &fEkin);

	//Converting degrees to radians
	//theta = ftheta * pi / 180.0;
	//phi = fphi * pi / 180.0;

	//Setting type of particle to a particle gun
	/* G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle(fpartname);
	fParticleGun->SetParticleDefinition(particle); */

	/* //Particle gun position
	x0 = x + 1500 * sin(theta) * cos(phi);
	yy0 = y + 1500 * sin(theta) * sin(phi);
	z0 = z + 1500 * cos(theta); */

	fParticleGun->SetParticlePosition(G4ThreeVector(0 * mm, 0 * mm, 1100 * mm));

	/* //Particle momentum direction
	ux = -sin(theta) * cos(phi);
	uy = -sin(theta) * sin(phi);
	uz = -cos(theta); */

	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, -1));

	//Particle kinetic energy
	//E0 = fEkin;

	fParticleGun->SetParticleEnergy(4. * GeV);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


