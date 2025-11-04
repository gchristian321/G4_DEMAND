#include "DRAGONDetectorConstruction.hh"     //local
                 

namespace DRAGON
{

void DRAGONDetectorConstruction::neighborhood()
     { 
//C.
//************************************************************************
//*                                                                      *
//*                     Define neighborhood of each module               *
//*                                                                      *
//************************************************************************
//C.
	  
    n_fngr[0][0] =  1;
    n_fngr[1][0] =  2;
    n_fngr[2][0] = 11;
    n_fngr[3][0] = 12;
    n_fngr[4][0] =  3;

    n_fngr[0][1] =  2;
    n_fngr[1][1] =  1;
    n_fngr[2][1] = 13;
    n_fngr[3][1] = 14;
    n_fngr[4][1] = 11;
    n_fngr[5][1] = 12;

    n_fngr[0][2] =  3;
    n_fngr[1][2] =  4;
    n_fngr[2][2] = 15;
    n_fngr[3][2] = 16;
    n_fngr[4][2] =  5;

    n_fngr[0][3] =  4;
    n_fngr[1][3] =  5;
    n_fngr[2][3] =  3;
    n_fngr[3][3] = 15;
    n_fngr[4][3] = 16;

    n_fngr[0][4] =  5;
    n_fngr[1][4] =  4;
    n_fngr[2][4] =  6;
    n_fngr[3][4] = 15;
    n_fngr[4][4] = 16;
    n_fngr[5][4] = 21;
    n_fngr[6][4] = 22;
    n_fngr[7][4] =  3;

    n_fngr[0][5] =  6;
    n_fngr[1][5] =  5;
    n_fngr[2][5] =  7;
    n_fngr[3][5] = 21;
    n_fngr[4][5] = 22;
    n_fngr[5][5] = 27;
    n_fngr[6][5] = 28;
    n_fngr[7][5] =  9;

    n_fngr[0][6] =  7;
    n_fngr[1][6] =  6;
    n_fngr[2][6] =  9;
    n_fngr[3][6] = 27;
    n_fngr[4][6] = 28;

    n_fngr[0][7] =  8;
    n_fngr[1][7] = 10;
    n_fngr[2][7] = 25;
    n_fngr[3][7] = 26;
    n_fngr[4][7] = 29;
    n_fngr[5][7] = 30;

    n_fngr[0][8] =  9;
    n_fngr[1][8] =  7;
    n_fngr[2][8] = 27;
    n_fngr[3][8] = 28;
    n_fngr[4][8] =  6;

    n_fngr[0][9] = 10;
    n_fngr[1][9] =  8;
    n_fngr[2][9] = 29;
    n_fngr[3][9] = 30;
    n_fngr[4][9] =  9;

    n_fngr[0][10] = 11;
    n_fngr[1][10] = 17;
    n_fngr[2][10] = 13;
    n_fngr[3][10] = 15;
    n_fngr[4][10] =  1;
    n_fngr[5][10] =  2;

    n_fngr[0][11] = 12;
    n_fngr[1][11] = 18;
    n_fngr[2][11] = 14;
    n_fngr[3][11] = 16;
    n_fngr[4][11] =  1;
    n_fngr[5][11] =  2;

    n_fngr[0][12] = 13;
    n_fngr[1][12] = 17;
    n_fngr[2][12] = 19;
    n_fngr[3][12] = 11;
    n_fngr[4][12] =  2;

    n_fngr[0][13] = 14;
    n_fngr[1][13] = 18;
    n_fngr[2][13] = 20;
    n_fngr[3][13] = 12;
    n_fngr[4][13] =  2;

    n_fngr[0][14] = 15;
    n_fngr[1][14] = 21;
    n_fngr[2][14] = 11;
    n_fngr[3][14] = 17;
    n_fngr[4][14] =  3;
    n_fngr[5][14] =  4;
    n_fngr[6][14] =  5;

    n_fngr[0][15] = 16;
    n_fngr[1][15] = 22;
    n_fngr[2][15] = 12;
    n_fngr[3][15] = 18;
    n_fngr[4][15] =  3;
    n_fngr[5][15] =  4;
    n_fngr[6][15] =  5;

    n_fngr[0][16] = 17;
    n_fngr[1][16] = 21;
    n_fngr[2][16] = 23;
    n_fngr[3][16] = 19;
    n_fngr[4][16] = 11;
    n_fngr[5][16] = 13;
    n_fngr[6][16] = 15;

    n_fngr[0][17] = 18;
    n_fngr[1][17] = 22;
    n_fngr[2][17] = 24;
    n_fngr[3][17] = 20;
    n_fngr[4][17] = 12;
    n_fngr[5][17] = 14;
    n_fngr[6][17] = 16;

    n_fngr[0][18] = 19;
    n_fngr[1][18] = 23;
    n_fngr[2][18] = 17;
    n_fngr[3][18] = 13;
    n_fngr[4][18] = 25;

    n_fngr[0][19] = 20;
    n_fngr[1][19] = 24;
    n_fngr[2][19] = 18;
    n_fngr[3][19] = 14;
    n_fngr[4][19] = 26;

    n_fngr[0][20] = 21;
    n_fngr[1][20] = 23;
    n_fngr[2][20] = 15;
    n_fngr[3][20] = 17;
    n_fngr[4][20] = 27;
    n_fngr[5][20] =  5;
    n_fngr[6][20] =  6;

    n_fngr[0][21] = 22;
    n_fngr[1][21] = 24;
    n_fngr[2][21] = 16;
    n_fngr[3][21] = 18;
    n_fngr[4][21] = 28;
    n_fngr[5][21] =  5;
    n_fngr[6][21] =  6;

    n_fngr[0][22] = 23;
    n_fngr[1][22] = 21;
    n_fngr[2][22] = 25;
    n_fngr[3][22] = 17;
    n_fngr[4][22] = 27;
    n_fngr[5][22] = 29;
    n_fngr[6][22] = 19;

    n_fngr[0][23] = 24;
    n_fngr[1][23] = 22;
    n_fngr[2][23] = 26;
    n_fngr[3][23] = 18;
    n_fngr[4][23] = 28;
    n_fngr[5][23] = 30;
    n_fngr[6][23] = 20;

    n_fngr[0][24] = 25;
    n_fngr[1][24] = 23;
    n_fngr[2][24] = 19;
    n_fngr[3][24] = 29;
    n_fngr[4][24] =  8;

    n_fngr[0][25] = 26;
    n_fngr[1][25] = 24;
    n_fngr[2][25] = 20;
    n_fngr[3][25] = 30;
    n_fngr[4][25] =  8;

    n_fngr[0][26] = 27;
    n_fngr[1][26] = 21;
    n_fngr[2][26] = 23;
    n_fngr[3][26] = 29;
    n_fngr[4][26] =  9;
    n_fngr[5][26] =  6;
    n_fngr[6][26] =  7;

    n_fngr[0][27] = 28;
    n_fngr[1][27] = 22;
    n_fngr[2][27] = 24;
    n_fngr[3][27] = 30;
    n_fngr[4][27] =  9;
    n_fngr[5][27] =  6;
    n_fngr[6][27] =  7;

    n_fngr[0][28] = 29;
    n_fngr[1][28] = 23;
    n_fngr[2][28] = 25;
    n_fngr[3][28] = 27;
    n_fngr[4][28] = 10;
    n_fngr[5][28] =  8;

    n_fngr[0][29] = 30;
    n_fngr[1][29] = 24;
    n_fngr[2][29] = 26;
    n_fngr[3][29] = 28;
    n_fngr[4][29] = 10;
    n_fngr[5][29] =  8;    
   }

}