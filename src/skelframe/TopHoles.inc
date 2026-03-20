#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"

G4LogicalVolume* BuildTopHoles()
{
std::vector<G4TwoVector> outline = {
    {23.925*mm,-32.0*mm},
    {-23.925*mm,-32.0*mm},
    {-23.925*mm,32.0*mm},
    {23.925*mm,32.0*mm},
    {23.925*mm,-32.0*mm},
};

G4VSolid* solid = new G4ExtrudedSolid("plate", outline, 3.175/2*mm, G4TwoVector(),1.0, G4TwoVector(),1.0);

G4VSolid* circ0 = new G4Tubs("circ0",0,18.0*mm,3.175*mm,0,360*deg);
solid = new G4SubtractionSolid("cutc0", solid, circ0, nullptr, G4ThreeVector(1.39e-14*mm,6.9e-15*mm,0));

auto nist = G4NistManager::Instance();
auto material = nist->FindOrBuildMaterial("G4_Al");

G4LogicalVolume* plateLV =
    new G4LogicalVolume(
        solid,
        material,
        "FrontFrameLV");

return plateLV;
}
