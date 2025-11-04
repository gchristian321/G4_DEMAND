#include "G4NistManager.hh"        //Geant4
#include "G4SystemOfUnits.hh"

#include "Materials.hh"               //local


Materials* Materials::instance = nullptr;

Materials* Materials::Instance() 
           {
            if (!instance) 
               {
                instance = new Materials();
                }
            return instance;
            }

void Materials::DefineMaterials() 
     {
      G4NistManager* nistManager = G4NistManager::Instance();
      // Material                                                                           index
      Hydrogen = nistManager->FindOrBuildMaterial("G4_H");                                   //1
      Deuterium = new G4Material("Deuterium", 1.0, 2.01*g/mole, 0.162*g/cm3);                //2
      Helium = nistManager->FindOrBuildMaterial("G4_He");                                    //3
      Lithium = nistManager->FindOrBuildMaterial("G4_Li");                                   //4
      Beryllium = nistManager->FindOrBuildMaterial("G4_Be");                                 //5
      Carbon = nistManager->FindOrBuildMaterial("G4_C");                                     //6
      Nitrogen = nistManager->FindOrBuildMaterial("G4_N");                                   //7
      Neon = nistManager->FindOrBuildMaterial("G4_Ne");                                      //8
      Aluminium = nistManager->FindOrBuildMaterial("G4_Al");                                 //9
      Iron = nistManager->FindOrBuildMaterial("G4_Fe");                                      //10 
      Copper = nistManager->FindOrBuildMaterial("G4_Cu");                                    //11
      Tungsten = nistManager->FindOrBuildMaterial("G4_W");                                   //12
      Lead = nistManager->FindOrBuildMaterial("G4_Pb");                                      //13
      Uranium = nistManager->FindOrBuildMaterial("G4_U");                                    //14
      Air = nistManager->FindOrBuildMaterial("G4_AIR");                                      //15
      Vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");    //VACUUM                  //16                               
 
      materialMap[1]  = Vacuum;
      materialMap[2]  = Vacuum;
      materialMap[3]  = Vacuum;
      materialMap[4]  = Vacuum;
      materialMap[5]  = Copper;
      materialMap[6]  = Aluminium;
      materialMap[7]  = Lead;	  
      materialMap[8]  = Air;
      materialMap[18]  = Tungsten;
 }


G4Material* Materials::GetMatByIndex(int index) const { 
    auto it = materialMap.find(index);
    if (it != materialMap.end()) {
        return it->second;
    } else {
        G4Exception("Materials::GetMaterialByIndex", "MaterialNotFound", FatalException,
                    ("No material defined for index " + std::to_string(index)).c_str());
        return nullptr;
    }
}

/*
FROM GEANT3
Command: GPMATE(0)
Output:
0===================================================     MATERIALS      ==================================================
0MATERIAL                           A         Z     DENSITY  RADIAT L  ABSORP L NMIXT

        1 HYDROGEN                 1.010     1.000     0.071 0.865E+03 0.790E+03   1
        2 DEUTERIUM                2.010     1.000     0.162 0.757E+03 0.342E+03   1
        3 HELIUM                   4.000     2.000     0.125 0.755E+03 0.478E+03   1
        4 LITHIUM                  6.940     3.000     0.534 0.155E+03 0.121E+03   1
        5 BERILLIUM                9.010     4.000     1.848 0.353E+02 0.367E+02   1
        6 CARBON                  12.010     6.000     2.265 0.188E+02 0.499E+02   1
        7 NITROGEN                14.010     7.000     0.808 0.445E+02 0.994E+02   1
        8 NEON                    20.180    10.000     1.207 0.240E+02 0.749E+02   1
        9 ALUMINIUM               26.980    13.000     2.700 0.890E+01 0.372E+02   1
       10 IRON                    55.850    26.000     7.870 0.176E+01 0.171E+02   1
       11 COPPER                  63.540    29.000     8.960 0.143E+01 0.148E+02   1
       12 TUNGSTEN               183.850    74.000    19.300 0.350E+00 0.103E+02   1
       13 LEAD                   207.190    82.000    11.350 0.560E+00 0.185E+02   1
       14 URANIUM                238.030    92.000    18.950 0.320E+00 0.120E+02   1
       15 AIR                     14.610     7.300     0.001 0.304E+05 0.675E+05   1
       16 VACUUM                   0.000     0.000     0.000 0.100E+17 0.100E+17   1
      
       17 SCINTILLATOR            11.079     5.577     1.032 0.421E+02 0.886E+02   2     A      Z     W
                                                                                        12.01   6.00  0.915
                                                                                         1.01   1.00  0.085
       18 BARIUM FLORIDE BAF2    111.656    45.812     4.890 0.203E+01 0.340E+02   2     A      Z     W
                                                                                       137.30  56.00  0.783
                                                                                        19.00   9.00  0.217
       19 CESIUM FLORIDE CSF     118.653    49.246     4.640 0.198E+01 0.377E+02   2     A      Z     W
                                                                                       132.90  55.00  0.875
                                                                                        19.00   9.00  0.125
       20 SODIUM IODIDE NAI:TL   110.958    46.556     3.670 0.259E+01 0.469E+02   2     A      Z     W
                                                                                        23.00  11.00  0.153
                                                                                       126.90  53.00  0.847
       21 CESIUM IODIDE CSI:TL   129.969    54.023     4.510 0.186E+01 0.419E+02   2     A      Z     W
                                                                                       132.90  55.00  0.512
                                                                                       126.90  53.00  0.488
       22 BGO       BI4GE3O12    155.409    62.525     7.130 0.112E+01 0.248E+02   3     A      Z     W
                                                                                       209.00  83.00  0.671
                                                                                        72.60  32.00  0.175
                                                                                        16.00   8.00  0.154
       23 LSO   LU2(SI04)O:CE    138.222    56.502     7.400 0.114E+01 0.231E+02   3     A      Z     W
                                                                                       175.00  71.00  0.764
                                                                                        28.10  14.00  0.061
                                                                                        16.00   8.00  0.175
       24 MGO (POWDER)            21.005    10.412     1.870 0.149E+02 0.611E+02   2     A      Z     W
                                                                                        24.30  12.00  0.603
                                                                                        16.00   8.00  0.397
       25 GLASS                   11.079     5.577     1.032 0.421E+02 0.886E+02   2     A      Z     W
                                                                                        12.01   6.00  0.915
                                                                                         1.01   1.00  0.085
       26 STAINLESS STEEL         55.023    25.630     7.705 0.182E+01 0.194E+02   6     A      Z     W
                                                                                        55.85  26.00  0.709
                                                                                        52.00  24.00  0.180
                                                                                        58.69  28.00  0.080
                                                                                        54.94  25.00  0.020
                                                                                        28.09  14.00  0.010
                                                                                        12.01   6.00  0.002

       WARNING: define a material with density=0 is not allowed in Geant4. 
       The material will be constructed with the default minimal density: 1e-25g/cm3
       
       27 target                   1.000     1.000     0.000 0.114E+09 0.000E+00   1
       28 centralvacuum            1.000     1.000     0.000 0.572E+11 0.000E+00   1
       29 ent1                     1.000     1.000     0.000 0.229E+10 0.000E+00   1
       30 ent2                     1.000     1.000     0.000 0.318E+10 0.000E+00   1
       31 ent3                     1.000     1.000     0.000 0.715E+10 0.000E+00   1
       32 ext1                     1.000     1.000     0.000 0.229E+10 0.000E+00   1
       33 ext2                     1.000     1.000     0.000 0.318E+10 0.000E+00   1
       34 ext3                     1.000     1.000     0.000 0.715E+10 0.000E+00   1
       35 base                     1.000     1.000     0.000 0.752E+06 0.368E-03   1
       50 SILICON                 28.080    14.000     2.330 0.270E+01 0.890E+01   1
       51 Modified SILICON        28.080    14.000     0.023 0.000E+00 0.000E+00   1
       52 MCP Carbon              12.000     6.000     0.002 0.000E+00 0.000E+00   1
       62 TARGET CARBON           12.000     6.000     0.023 0.000E+00 0.000E+00   1

FROM GEANT3
Command: GPTMED(0)
Output:
0===================================================   TRACKING MEDIA   ==================================================
0TMED                          MATERIAL ISVOL IFIELD  FIELDM  TMAXFD  STEMAX    DEEMAX   EPSIL   STMIN

      1 VACUUM ->  no field       16       0     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
      2 VACUUM -> ifield = 1      16       1     1    100.00   10.00 -1.00      -1.000   0.001  -1.000
      3 VACUUM -> ifield = 2      16       1     2    100.00   10.00 -1.00      -1.000   0.001  -1.000
      4 VACUUM -> ifield = 3      16       1     3    100.00   10.00 -1.00      -1.000   0.001  -1.000
      5 COPPER                    11       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
      6 ALUMINUM                   9       0     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
      7 LEAD                      13       0     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
      8 ATMOSPHERE (AIR)          15       0     0      0.00   10.00 -1.00      -1.000   0.100  -1.000
      9 SCINTILLATOR              17       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     10 BARIUM FLORIDE BAF2       18       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     11 CESIUM FLORIDE CSF        19       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     12 SODIUM IODIDE NAI:TL      20       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     13 CESIUM IODIDE CSI:TL      21       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     14 BGO       BI4GE3O12       22       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     15 LSO   LU2(SI04)O:CE       23       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     16 MGO (POWDER)              24       0     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     17 GLASS                     25       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     18 TUNGSTEN                  12       0     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     19 SILICON                   50       1     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     20 STAINLESS STEEL           26       0     0      0.00   10.00 -1.00      -1.000   0.001  -1.000
     21 TARGET CARBON             62       1     0      0.00   10.00 0.400E-04  -1.000   0.000  -1.000
     22 central vacuum            28       0     0      0.00   10.00 -1.00      -1.000   0.100  -1.000
     23 Modified SILICON          51       1     0      0.00   10.00 0.400E-04  -1.000   0.000  -1.000
     24 MCP Carbon                52       1     0      0.00   10.00 0.400E-04  -1.000   0.000  -1.000
     26 STAINLESS STEEL           26       0     0      0.00   10.00  10.0       0.005   0.001   0.001
     27 Gas Target                27       0     0      0.00   10.00 0.500E-01   0.100   0.000   0.000
     28 centralvacuum             28       0     0      0.00   10.00  10.0       0.100   0.000   0.000
     29 entrance1                 29       0     0      0.00   10.00  1.00       0.100   0.000   0.000
     30 entrance2                 30       0     0      0.00   10.00  1.00       0.100   0.000   0.000
     31 entrance3                 31       0     0      0.00   10.00  1.00       0.100   0.000   0.000
    
     32 exit1                     32       0     0      0.00   10.00  1.00       0.100   0.000   0.000
     33 exit2                     33       0     0      0.00   10.00  1.00       0.100   0.000   0.000
     34 exit3                     34       0     0      0.00   10.00  1.00       0.100   0.000   0.000
     35 baseline                  35       0     0      0.00   10.00  1.00       0.100   0.000   0.000

*/





















