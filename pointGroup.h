#ifndef POINTGROUP_H
#define POINTGROUP_H



#include <string>
#include <vector>
#include <math.h>
#include "matrix3x3.h"
#include "vector3.h"

using namespace std;


enum MoleculeShape {
    Polygon,
    Plane,      //Ia + Ib = Ic
    Line,       //Ia = 0 Ib = Ic
    Sphere,     //Ia = Ib = Ic
    Oblate,     //Ia = Ib < Ic
    Prolate,    //Ia < Ib = Ic
    Irregular   //Ia < Ib < Ic
};

//point-group symmetry
class PointGroup
{
public:
    PointGroup(double );

public:

    uint  NAtoms; //number of atoms
    vector < vector <double> > DistsMatrix;
    vector < vector <uint> >   SubGroups;

    MoleculeShape shape;

public:
    bool    hasI,  isCoov, isDooh;
    bool    isDnd, isDnh, isDn, isS2n;
    bool    isCnv, isCn, isCnh;
    bool    isTd, isT, isTh;
    bool    isIh, isI;
    bool    isO,  isOh;
    bool    isSn, isCi, isCs;

public:
    uint    principalAxisOrder,order;
    double  angle;
    string  PGSymbol;
    string  first,last,middle;

    double tol;

public:
    matrix3x3 IMomentMatrix;
    vector3 IMoment,MassCenter;


    //coordinates and masses of atoms
    vector < vector3 >  AtomCoordinates;
    vector <double>     AtomMasses;

    //type of atom: e.g. 0=X, 1 = H, 2=He ...
    vector  < uint >    AtomTypes;
    vector  < string >  AtomSymbols;

public:
    vector3 PrincipalAxis;
    vector3 Ci;
    vector < vector3 > Cn;
    vector < vector3 > C2;
    vector < vector3 > Horizontal_C2;
    vector < vector3 > C3;
    vector < vector3 > C4;
    vector < vector3 > C5;
    vector < vector3 > C6;

    vector < vector3 > S2n;

    vector < vector3 > S12;
    vector < vector3 > S10;
    vector < vector3 > S8;
    vector < vector3 > S6;
    vector < vector3 > S4;


    vector < vector3 > SigmaD;
    vector < vector3 > SigmaV;
    vector < vector3 > SigmaH;

    std::vector<std::string> plit(const std::string& , char );
    
    void generalPerceptSymmetry();
    void detectPrincipalAxisOrder();


    void perceptMolShape();

    void loadFile(std::string );
    void writeFile(std::string );
    void refreshMol ();

    void buildDistanceMatrix();
    void buildSubGroups_perception();    //only for percept Symmetry
    void buildSubGroups_perception_cluster();    //only for percept Symmetry
    bool buildSubGroupsOnSymmetry_refine(); //only for patch
    bool buildSubGroupsOnSymmetry_patch();  //only for refine

    void CheckCenter();
    bool CheckMatrix( matrix3x3 & );

    bool CheckC2(vector3  );
    bool CheckCn(vector3, uint );
    bool CheckS2n(vector3, uint );
    bool CheckSigma(vector3);

    void SearchMirror();
    void SearchC2();
    void SearchC3();
    void SearchC4();
    void SearchC5();
    void SearchCn(uint o,vector < vector3 > & );

    bool findFirstHorizontalC2(vector3 & );
    bool perceptHorizontal_C2(bool);
    bool perceptMirrorH();

    void perceptVerticalMirror(vector < vector3 > & ); //search Mirror based on parallel to Cn
    void perceptVerticalMirror(bool); //search Mirror based on perceptHorizontal C2

    void perceptSphere();
    void perceptOblate();
    void perceptProlate();

    void perceptThMirror(bool isCheck);
    void perceptTdMirror(bool isCheck);
    void perceptOhMirror(bool isCheck);
    void perceptIhMirror(bool isCheck);

    void appendVector3 (vector < vector3 > & , vector3 & );
    void refineSymmetryElements();
    void report(string pgsymbol="");
    void summary();

public:
    //build point group operation
    vector<matrix3x3 > PGOperation;
    vector<matrix3x3 > InvPGOperation;

    void buildInvOperation();
    void buildOperation();

    void buildOperation_Cs();
    void buildOperation_Ci();
    void buildOperation_Cn(uint);
    void buildOperation_S2n(uint);
    void buildOperation_Cnv(uint);
    void buildOperation_Cnh(uint);

    void buildOperation_Dn(uint);
    void buildOperation_Dnd(uint);
    void buildOperation_Dnh(uint);

    void buildOperation_T();
    void buildOperation_Th();
    void buildOperation_Td();

    void buildOperation_O();
    void buildOperation_Oh();

    void buildOperation_I();
    void buildOperation_Ih();

    void buildOperation_D00h();
    void buildOperation_C00v();


    // Patch molecule
    uint addingAtoms;
    bool patchMolecule(uint & NAtomsofPatched);
    void patchC00v();
    void patchD00h();


    //regroup the atoms according to the point-group operation;
    //vector < vector < uint > > AtomGroupList;

    //refine the cartesian data
    bool refine();  //based on symmetry elements
    void ClearAll();


public:
    void PerceptSymmetry(bool);
    void StandardOrientation();
    bool isOrientated;

    matrix3x3  rotMatrix,rotMatrix1,rotMatrix2;

    //rotate one vector to Z, and the second to X axis;
    void Orientation2ZX(vector3, vector3);

    //rotate one vector to Z, and the second to Y axis;
    void Orientation2ZY(vector3, vector3);

    void Orientation_Ci();
    void Orientation_Cs();
    void Orientation_Cn();
    void Orientation_S2n();
    void Orientation_Cnv();
    void Orientation_Cnh();
    void Orientation_Dn();
    void Orientation_Dnd();
    void Orientation_Dnh();
    void Orientation_Line();

    void Orientation_T();
    void Orientation_O();
    void Orientation_I();

    matrix3x3 OrientationFromTo(vector3,vector3);


//  tools
    bool isEquivalent (vector3, vector3 );
    bool isOrthogonal (vector3, vector3 );
    bool isRegularPolygon( vector < vector3 >); //using atoms' coordinates
    bool isRegularPolygon( vector < uint > ); // using atoms' index

    void Translate(vector3 & );

    void Centralize();
    void CalcInertialMoment();

    bool isInside (vector < vector3> & , vector3 & );
    void CalcGCD (vector <uint > & , uint );

    void refinePrincipalAxis ();
    void setTolerance (double );
    void setSymmetry (string );
    uint detectAxisOrder(vector3 );

    void generate_O_from_2C3(vector3, vector3);
    void generate_O_from_2C4(vector3, vector3);
    void generate_O_from_C4C2(vector3, vector3);

    void generate_I_from_2C2_orth(vector3, vector3);
    void generate_I_from_2C5(vector3, vector3); //neighbored C5

    void printInertialMoment();
    void ParseSymmetrySybmol (string );

    bool CombinationCn(vector < uint > & arr, uint data[], uint start, uint end,
                     uint index, uint r,
                     vector < vector3 > & container,
                     uint order0);


    double elementMasses[119] = {0.0, 1.00794, 4.002602, 6.941,
        9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
        22.98976928, 24.3050, 26.9815386, 28.0855, 30.973762, 32.065, 35.453, 39.948,
        39.0983, 40.078, 44.955912, 47.867, 50.9415, 51.9961, 54.938045, 55.845,
        58.933195, 58.6934, 63.546, 65.38, 69.723, 72.64, 74.92160, 78.96, 79.904,
        83.798, 85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.96, 98, 101.07,
        102.90550, 106.42, 107.8682, 112.411, 114.818, 118.710, 121.760, 127.60,
        126.90447, 131.293, 132.9054519, 137.327, 138.90547, 140.116, 140.90765,
        144.242, 145, 150.36, 151.964, 157.25, 158.92535, 162.500, 164.93032,
        167.259, 168.93421, 173.054, 174.9668, 178.49, 180.94788, 183.84, 186.207,
        190.23, 192.217, 195.084, 196.966569, 200.59, 204.3833, 207.2, 208.98040,
        209, 210, 222, 223, 226, 227, 232.03806, 231.03588, 238.02891, 237, 244, 243,
        247, 247, 251, 252, 257, 258, 259, 262, 265, 268, 271, 272, 270, 276, 281,
        280, 285, 284, 289, 288, 293, 292, 294
    };

    string ElementNames[119] = {"X", "H", "He", "Li", "Be", "B", "C", "N",
        "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
        "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
        "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
        "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
        "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg",
        "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"
    };

    uint sybmol2AtomicNum (string  s)
    {
        uint _atomicNum = 0;
        for (int i=0;i<119;i++)
        {
            if(ElementNames[i]==s)
            {
                _atomicNum=i;
                break;
            }
        }
        return _atomicNum;
    }
};


#endif // POINTGROUP_H
