#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"

G4LogicalVolume* BuildBottomHoles()
{
std::vector<G4TwoVector> outline = {
    {23.925*mm,-76.5*mm},
    {-23.925*mm,-76.5*mm},
    {-23.925*mm,76.5*mm},
    {23.925*mm,76.5*mm},
    {23.925*mm,-76.5*mm},
};

G4VSolid* solid = new G4ExtrudedSolid("plate", outline, 3.175/2*mm, G4TwoVector(),1.0, G4TwoVector(),1.0);

G4VSolid* circ0 = new G4Tubs("circ0",0,18.0*mm,3.175*mm,0,360*deg);
solid = new G4SubtractionSolid("cutc0", solid, circ0, nullptr, G4ThreeVector(1.39e-14*mm,-44.49999999999993*mm,0));

G4VSolid* circ1 = new G4Tubs("circ1",0,18.0*mm,3.175*mm,0,360*deg);
solid = new G4SubtractionSolid("cutc1", solid, circ1, nullptr, G4ThreeVector(1.91e-14*mm,44.50000000000006*mm,0));

G4VSolid* circ2 = new G4Tubs("circ2",0,18.0*mm,3.175*mm,0,360*deg);
solid = new G4SubtractionSolid("cutc2", solid, circ2, nullptr, G4ThreeVector(1.37e-14*mm,6.94e-14*mm,0));

auto nist = G4NistManager::Instance();
auto material = nist->FindOrBuildMaterial("G4_Al");

G4LogicalVolume* plateLV =
    new G4LogicalVolume(
        solid,
        material,
        "FrontFrameLV");

return plateLV;
}
