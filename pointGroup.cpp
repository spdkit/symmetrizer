#include "pointGroup.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cctype>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <string>
#include <cstdlib>

#include <filesystem>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <filesystem>

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <cstdio>

using namespace std;


bool PointGroup::isRegularPolygon( vector < uint > idx)
{
    uint size0=idx.size();

    vector < vector3 > coords;
    coords.resize(size0);

    for (uint i=0;i< size0;i++)
        coords[i]=AtomCoordinates[idx[i]];
    return isRegularPolygon(coords);
}


bool PointGroup::isRegularPolygon( vector < vector3 > coords)
{
    double eps=0.1;

    //the same plane or not?
    uint size0=coords.size();
    if ( size0< 3) return false;

    vector3 center0=VZero;
    for (uint i=0; i < size0; i++) center0+=coords[i];
    center0/=size0;

    //centralize
    for (uint i=0;i < size0; i++)  coords[i]-=center0;

    //cheap check based on average distances
    double dist=0.0;
    for (uint i=0;i < size0; i++)
    {
        dist+=coords[i].length();
    }
    dist/=size0;

    double err=0.0;
    for (uint i=0;i < size0; i++)
    {
        err+=abs(dist-coords[i].length());
    }
    err/=size0;

    //scale
    //cout << "How to scale the err here inside isRegularPolygon?"<<endl;

    if ( err > eps )
        return false;

    //triangle is enough
    if(size0==3) return true;

    //size0 >= 4
    vector3 normal_vec=cross(coords[0],coords[1]);
    for (uint i = 0; i < size0; ++ i)
    {
        if (i<=1) continue;
        if  (!isOrthogonal(coords[0],normal_vec))
            return false;
    }

    return true;
}



bool PointGroup::isEquivalent (vector3  vec1, vector3 vec2)
{
   if (isEqual  (vec1.normalize(), vec2.normalize()) ) return true;
   if (isNegativeEqual(vec1.normalize(), vec2.normalize()) ) return true;
   return false;
}


bool PointGroup::isOrthogonal (vector3  vec1, vector3  vec2)
{
    double dotvalue=dot(vec1,vec2);
    if (abs(dotvalue) < 0.1) return true;

    /*double angle=abs(vectorAngle(vec1,vec2)-90.0);
    if (angle < 2.0) return true;
    */
    return false;
}




void PointGroup::setTolerance(double tole)
{
    cout << "Set tolerance " << tole << endl;
    tol=tole/0.529177;
}


bool PointGroup::isInside (vector < vector3> & container, vector3 & v1)
{
    bool find=false;
    uint size=container.size();

    for (uint i=0; i < size; i++ )
    {
        if (isEquivalent (container[i],v1)) {
            find = true;
            break;
        }
    }
    return find;
}


uint findMax(vector3 & vec)
{
    double a = abs(vec[0]), b = abs(vec[1]), c = abs(vec[2]);

    double maxone=a;
    uint id=0;
    if (b > maxone) {maxone=b;id=1;}
    if (c > maxone) {maxone=c;id=2;}
    return id;
}

void printVector3(vector3 & vec)
{
    printf ("   %10.4f %10.4f %10.4f\n",vec[0],vec[1],vec[2]);
}

void refineVector3 (vector3 & vec){
    if (vec[findMax(vec)] < 0.0)  vec=-vec;
    if ( abs(vec[0]) < 0.00001) vec.SetX(0.0);
    if ( abs(vec[1]) < 0.00001) vec.SetY(0.0);
    if ( abs(vec[2]) < 0.00001) vec.SetZ(0.0);
}


PointGroup::PointGroup(double tolerance)
{
    ClearAll();
    tol = tolerance/0.529177;
}





void PointGroup::loadFile(std::string fileName)
{
    fstream file(fileName);
    if (!file.good())  {
        cout << fileName << " is not existed "<<endl;
        return;
    }

    char abs_path_buff[100];
    if(realpath(fileName.c_str(), abs_path_buff))
        cout << abs_path_buff<<endl;

    std::string fileOut;
    fileOut=fileName+"_new";

    std::ofstream out(fileOut);


    string s,s1,sybmol;
    std::vector<std::string> single;
    uint nAtoms=0,_num;
    string title;
    double x,y,z;
    string::size_type idx;
    char buffer [128];

    vector <string> symmetryList;

    uint nMolNum=0;
    while(getline(file,s))
    { 
        ClearAll();
        string lines="";

        std::istringstream iss (s);
        iss >> nAtoms;
        if (iss.fail()) {
            cout << s<<endl;
            break;
        }

       NAtoms=nAtoms;
       nMolNum++;

       //cout << nAtoms <<endl;

       //skip one line
       getline(file,title);
       //Energy =  -5430.59490500 a.u.  #Cluster:    1   Ni3CO8-neutral-0-singlet-B3LYP-def2SVP-opt-freq-130   C1   NO   IMG_FREQ   1-A.


       while(nAtoms--)
       {
           getline(file,s);
           std::istringstream str(s);

           str>>sybmol>> x >> y >>z;
           vector3 v(x,y,z);

           AtomSymbols.push_back(sybmol);
           AtomCoordinates.push_back(v);
           _num=sybmol2AtomicNum(sybmol);
           AtomMasses.push_back(elementMasses[_num]);
           AtomTypes.push_back(_num);

           cout << _num << "  " << elementMasses[_num]<<endl;
       }

       cout << " the " <<  nMolNum << " molecule "<<endl;
       PerceptSymmetry(true);
       report();


       if(refine())
           cout << "refine OK!"<<endl;
       else
           cout << "cannot refine OK!"<<endl;

       lines=to_string(NAtoms)+"\n";
       lines+= (PGSymbol+title+"\n");

       symmetryList.push_back(PGSymbol);

       for(uint i=0;i<NAtoms;i++)
       {
           sprintf(buffer,"%s  %12.8f  %12.8f  %12.8f\n",AtomSymbols[i].c_str(),
                   AtomCoordinates[i][0],
                   AtomCoordinates[i][1],
                   AtomCoordinates[i][2]);
           lines+=buffer;

           //cout << buffer;
       }

       out.write(lines.c_str(),lines.size());

       //goto label;
    }

    file.close();
    out.close();

    for(uint i=0;i<symmetryList.size();i++)
        cout << i+1 << " : " << symmetryList[i]<<endl;
}

void PointGroup::writeFile(std::string fileName)
{

      //  ofstream out("b.txt");
}

/*
void PointGroup:: GetData()
{  
    NAtoms=mol-> numAtoms();

    if (NAtoms<1) return;
    for (uint i=0;i<mol->numAtoms();i++) {
        AtomCoordinates.push_back(mol->atomPos(i));
        AtomMasses.push_back(mol->atomicMass(i));
        AtomSymbols.push_back(mol->atomSymbol(i).toStdString());
        AtomTypes.push_back(mol->atomicNum(i));
        //cout << i+1<<": "<< mol->atomSymbol(i).toStdString().c_str()<<" "<< mol->atomPos(i) <<" "<< mol->atomicMass(i)<< endl;
    }
}
*/


void PointGroup::ClearAll()
{
    addingAtoms=0;
    NAtoms=0;

    MassCenter=VZero;

    PGOperation.clear();
    InvPGOperation.clear();
    AtomCoordinates.clear();
    AtomMasses.clear();
    AtomTypes.clear();
    DistsMatrix.clear();
    SubGroups.clear();


    Cn.clear();
    C2.clear();
    C3.clear();
    C4.clear();
    C5.clear();
    C6.clear();

    SigmaD.clear();
    SigmaV.clear();
    SigmaH.clear();

    hasI=false;
    isCoov=isDooh=false;
    isDnd=isDnh=isDn=false;
    isS2n=false;
    isCnv=isCn=isCnh=false;
    isTd=isT=isTh=false;
    isIh=isI=false;
    isO=isOh=false;
    isS2n=isCi=isCs=false;
}



void PointGroup::Translate(vector3 & vec)
{
    for (uint i = 0; i < NAtoms; ++ i) {
        AtomCoordinates[i] =  AtomCoordinates[i] + vec;
    }
}


void PointGroup::Centralize()
{
        MassCenter=VZero;
        double  TotalMass = 0.0;
        for (uint i = 0; i < NAtoms; i++){
            MassCenter += AtomCoordinates[i]*AtomMasses[i];
            TotalMass += AtomMasses[i];
        }

        MassCenter /= TotalMass;
        cout << "\nMass Center:"<< MassCenter <<endl;

        for (uint i = 0; i < NAtoms; i++){
            AtomCoordinates[i] -= MassCenter;
        }
}

void  PointGroup::printInertialMoment()
{
    std::cout.setf(ios::fixed);
    std::cout << setprecision(6);

    cout << endl;
    cout <<"       Molecular shape based on the eigenvalues"<<endl;
    cout <<"                 of its inertia Moment:"<<endl;
    cout <<"-----------------------------------------------------"<<endl;
    cout <<"             Ia + Ib = Ic Ia = Ib :  Polygon"<<endl;
    cout <<"             Ia + Ib = Ic Ia != Ib:  Plane"<<endl;
    cout <<"             Ia = 0 Ib = Ic :    Line"<<endl;
    cout <<"             Ia = Ib = Ic   :    Sphere"<<endl;
    cout <<"             Ia = Ib < Ic   :    Oblate"<<endl;
    cout <<"             Ia < Ib = Ic   :    Prolate"<<endl;
    cout <<"             Ia < Ib < Ic   :    Irregular shape"<<endl;
    cout <<"-----------------------------------------------------"<<endl;

    std::cout << endl <<"        InertialMoment Matrix:" << endl;
    std::cout << "     "<<IMomentMatrix.GetColumn(0) << endl;
    std::cout << "     "<<IMomentMatrix.GetColumn(1) << endl;
    std::cout << "     "<<IMomentMatrix.GetColumn(2) << endl<<endl;;

    cout <<"     "<<IMoment<<endl<<endl;
}

//mass-weighted
void  PointGroup::CalcInertialMoment()
{

    double Ixx=0.0, Ixy=0.0, Ixz=0.0,
                    Iyy=0.0, Iyz=0.0,
                             Izz=0.0;
    double x, y, z;
    //                  | Ixx     Ixy     Ixz |
    //          I   =   | Iyx     Iyy     Iyz |
    //                  | Izx     Izy     Izz |
    for (uint i = 0; i < NAtoms; ++ i)
    {
        x=AtomCoordinates[i][0];
        y=AtomCoordinates[i][1];
        z=AtomCoordinates[i][2];

        //cout << x << " "<< y << " "<< z << " "<<AtomMasses[i]<<endl;

        Ixx += (y*y + z*z)*AtomMasses[i];
        Iyy += (x*x + z*z)*AtomMasses[i];
        Izz += (x*x + y*y)*AtomMasses[i];
        Ixy -= (x*y)*AtomMasses[i];
        Ixz -= (x*z)*AtomMasses[i];
        Iyz -= (y*z)*AtomMasses[i];
    }

    double I[9] = {Ixx,Ixy,Ixz,Ixy,Iyy,Iyz,Ixz,Iyz,Izz};

    matrix3x3 InertialMomentMatrix(I);

    IMomentMatrix = InertialMomentMatrix.findEigenvectorsIfSymmetric(IMoment);

    printInertialMoment();

    // Diagonalize the matrix
    //matrix3x3 diagonalMatrix = eigenmatrix.inverse() * InertialMomentMatrix * eigenmatrix;
}


void PointGroup::perceptMolShape()
{
    /*
     Polygon:      Ia = Ib and 2*Ia = Ic
     Line:         Ia = 0  and Ib = Ic
     Plane:        Ia + Ib = Ic
     Sphere:       Ia = Ib = Ic
     Oblate:       Ia = Ib < Ic
     Prolate:      Ia < Ib = Ic
     Irregular:    Ia < Ib < Ic
    */

    CalcInertialMoment();
    double a=IMoment[0],
           b=IMoment[1],
           c=IMoment[2];

    if (tol <0.0001) tol = 0.001;
    double eps=2.0*c/1000.0; //0.3% error
    cout << "The tolerance of InertialMoment:  " << eps << endl;

    //Ia < Ib < Ic
    shape=Irregular;

    //Sphere: Ia = Ib = Ic
    if ( (abs(a-b) < eps) && (abs(b-c) < eps) && (abs(a-c) < eps) )
        shape=Sphere;

    //Line   : Ia =0.0 and Ia < Ib = Ic
    else if ( abs(a)< eps)
        shape=Line;

    //Prolate: Ia < Ib = Ic
    else if ( abs(b-c)< eps && (a < c) && abs(a-c) > eps )
         shape=Prolate;


   //Polygon an Plane: Ia=Ib and Ia+Ib=Ic
    else if ( abs(a+b-c) < eps)
    {
        if ( abs(a-b)< eps ) // Ia=Ib
            shape=Polygon;
        else
            shape=Plane;  //Ia!=Ib
    }


    //Oblate:  Ia = Ib < Ic
    else if ( abs(a-b)< eps && (b<c) )
               shape=Oblate;


    switch (shape) {
        case Sphere:
            cout << "Sphere Shape"<<endl;
            break;
        case Line:
            cout << "Line Shape"<<endl;
            break;
        case Plane:
            cout << "Plane Shape"<<endl;
            break;
        case Polygon:
            cout << "Polygon Shape"<<endl;
            break;
        case Prolate:
            cout << "Prolate Shape"<<endl;
            break;
        case Oblate:
            cout << "Oblate Shape"<<endl;
            break;
        default:
            cout << "Irregular Shape"<<endl;
            break;
    }

}

void PointGroup:: buildDistanceMatrix()
{
    DistsMatrix.clear();
    DistsMatrix.resize(NAtoms);
    for (uint i=0; i < NAtoms ; i++)
        DistsMatrix[i].resize(NAtoms);

    cout << "\nBuilding DistanceMatrix " <<endl;
    std::cout.setf(ios::fixed);
    std::cout << setprecision(2);

    for (uint i = 0; i < NAtoms; ++ i) {
        for (uint j = 0; j < NAtoms; ++ j){
            DistsMatrix[i][j]=(AtomCoordinates[i]-AtomCoordinates[j]).length();
            cout << "   " << DistsMatrix[i][j] ;
        }
        cout << endl;
    }
    cout << "Done\n" <<endl;
}



void PointGroup::buildSubGroups_perception()
{
    buildDistanceMatrix();

    cout << "Regroup atoms based on distance matrix:\n";
    cout << "Atom Type:  Atom index \n";

    //initialize
    SubGroups.clear();

    uint idx=0;

    //label the atoms checked
    vector <uint> Labels (NAtoms,1);
    vector <uint> single;

    for (uint i=0; i < NAtoms; i++)
    {
        //check it first
        if (Labels[i]==1) Labels[i]=0;
        else continue;

        vector <double> dists_i (DistsMatrix [i]);
        vector <uint> each;

        each.clear();
        each.push_back(i);

        for (uint j=i+1; j < NAtoms; j++)
        {
            //check it first
            if (Labels[j]<1) continue;

            //not the same type
            if (AtomTypes [i] != AtomTypes [j]) continue;

            vector <double> dists_j (DistsMatrix [j]);

            //check average distance first, cheap
            double dist1=0.0,dist2=0.0;
            for (uint k = 0;k < dists_j.size(); k++)
            {
                dist1+=dists_i[k];
                dist2+=dists_j[k];
            }
            dist1/=dists_j.size();
            dist2/=dists_j.size();
            if (abs(dist1-dist2) > tol) continue;

            std::sort (dists_i.begin(),dists_i.end());
            std::sort (dists_j.begin(),dists_j.end());

            //for (uint k = 0;k< dists_i.size(); k++)
            //    cout << " " <<dists_i[k];
            //cout << endl;

            bool ok=true;
            for (uint k = 0;k < dists_j.size(); k++)
                if (abs(dists_i[k]-dists_j[k]) > tol )
                { ok=false; break;}

            if (!ok) continue;

            each.push_back(j);
            Labels[j]=0;
        }

        if (each.size() < 1) continue;

        //store
        SubGroups.push_back(each);
        //for (uint k = 0;k< each.size(); k++)     cout <<  1+each[k] << "   ";      cout <<endl;

        /*
        SubGroups.resize(idx+1);
        for (uint k = 0;k< each.size(); k++)
             SubGroups[idx].push_back(each[k]);
        idx=idx+1;
        */

    }


    //cout << "\nThere are "<<SubGroups.size()<< " groups" <<endl;
    uint subSize=0;

    for (uint ii=0; ii < SubGroups.size(); ii++)
    {
        cout <<  setw(5) << "   " <<AtomTypes[SubGroups[ii][0]]<<":";
        subSize=SubGroups[ii].size();

        for (uint jj=0; jj < subSize; jj++)
        {
            cout << setw(5)<<SubGroups[ii][jj]+1;

            if(jj+1<15)
                continue;
            if ((jj+1)%15==0)
                cout << "\n        ";
        }

        //vector<uint>::iterator iter;
       // for( iter= SubGroups[i].begin();iter != SubGroups[i].end(); ++iter)
         //   cout << setw(5)<< (*iter) +1;
        cout << endl;
    }
    cout << endl;

    cout << "buildSubGroups done!"<<endl;
}




//Cluster algorithms:
//Marcus Johansson and Valera Veryazov, J Cheminform 2017,9,8 (doi 10.1186/s13321-017-0193-3)
struct data0 {
    uint idx;
    vector3 pos;
};

bool cmpX (data0 a,data0 b){
    return (a.pos.x() < b.pos.x());
}

void PointGroup::buildSubGroups_perception_cluster()
{
    //unit sphere, plane
    vector < vector3 > S,P;
    S.resize(NAtoms); //Si=Sigmaj(Uij*Dij/|Dij|);
    P.resize(NAtoms); //Pi=Sigmaj{Uij*(aj-(ai.aj)*ai/|ai|^2)}

    //reduced mass
    vector < vector <double> > U,D;
    U.clear();
    U.resize(NAtoms);
    for (uint i=0; i < NAtoms ; i++) U[i].resize(NAtoms);

    D.clear();
    D.resize(NAtoms);
    for (uint i=0; i < NAtoms ; i++)  D[i].resize(NAtoms);

    //S.resize(NAtoms); //Si=Sigmaj(Uij*Dij/|Dij|);
    //P.resize(NAtoms); //Pi=Sigmaj{Uij*(aj-(ai.aj)*ai/|ai|^2)}
    for (uint i = 0; i < NAtoms; ++ i)
    {
        vector3 s=VZero,
                p=VZero;

        for (uint j = 0; j < NAtoms; ++ j)
        {
            if (j==i) continue;

            double miu=1.0;   //ghost atom
            if (AtomMasses[i]>0.0 && AtomMasses[j]>0.0)
                miu=AtomMasses[i]*AtomMasses[j]/(AtomMasses[i]+AtomMasses[j]);

            U[i][j]=miu;
            vector3 vec = AtomCoordinates[i]-AtomCoordinates[j];

            D[i][j]=miu*vec.length();
            if (vec.length()>0.0)
                s+=miu*vec.normalize();

            double length2=AtomCoordinates[i].length_2();
            double dot0=dot(AtomCoordinates[i],AtomCoordinates[j]);

            //Pi=Uij*(aj-(ai.aj)*ai/|ai|^2)
            p+=miu*(AtomCoordinates[j]-dot0*AtomCoordinates[i]/length2);
        }
        P[i]=p;
        S[i]=s;
        //cout <<s<<"  " <<p<<endl;

    }

    //mean Di and its length
    //vector <vector3 > Di;
    //vector <double > LDi;
    //Di.resize(NAtoms);  LDi.resize(NAtoms);


    vector <data0> par;
    par.resize(NAtoms);

    double x_bar=0.0,y_bar=0.0,z_bar=0.0;
    double x_min=10000.0,y_min=10000.0,z_min=10000.0;
    double x,y,z;
    for (uint i = 0; i < NAtoms; ++ i) {
        double d=0.0,dl=0.0;
        for (uint j = 0; j < NAtoms; ++ j) {
            d+=D[i][j];
            dl+=D[i][j]*D[i][j];
        }

        x=P[i].length();
        y=d;
        z=sqrt(dl);

        //Di[i]=d/NAtoms; LDi[i]=sqrt(dl);
        par[i].idx=i;
        par[i].pos.Set(x,y,z);

        if (x_min > x) {
            if (x>0.00001) x_min=x;
        }
        if (y_min >y) y_min=y;
        if (z_min >z) z_min=z;

        x_bar+=x;
        y_bar+=y;
        z_bar+=z;

        //cout << i+1 <<":  "<<P[i].length() <<" " <<d<< " " <<sqrt(dl)<<endl;
    }

    //sort
    //std::sort (Di.begin(),Di.end());
    //for (uint i = 0; i < NAtoms; ++ i)        cout << Di[i] << " ";

    std::sort (par.begin(),par.end(),cmpX);

    std::cout.setf(ios::fixed);
    std::cout << setprecision(6);

    cout << "\n\nRegroup atoms based on modified cluster algorithms:\n";
    cout << "  Marcus Johansson and Valera Veryazov, J. Cheminform. 2017,9,8. (doi 10.1186/s13321-017-0193-3)"<<endl;
    cout << "Atom Type:  Atom index \n";

    //initialize
    uint idx=0;
    SubGroups.clear();
    SubGroups.resize(idx);

    x_bar/=NAtoms;
    y_bar/=NAtoms;
    z_bar/=NAtoms;

    double eps_x=tol*x_min/5.0;
    if (eps_x<0.001) eps_x=0.001;

    double eps_y=tol*y_min/50.0;
    double eps_z=tol*z_min/50.0;

    cout << "Tolerances : "<< tol<<endl;
    cout << "   "<<eps_x <<"   "<<eps_y << "   "<<eps_z<<endl;
    for (uint i = 0; i < NAtoms; ++ i)  cout <<par[i].idx +1 << " " <<par[i].pos<<endl;

    //label the atoms checked
    vector <uint> Labels (NAtoms,1);

    for (uint i=0; i < NAtoms; i++)
    {
        //check it first
        if (Labels[i]==1) Labels[i]=0;
        else continue;

        vector <uint> each;
        each.clear();
        each.push_back(par[i].idx);

        double x1,y1,z1;
        double x=par[i].pos[0],
               y=par[i].pos[1],
               z=par[i].pos[2];

        //group based on X
        for (uint j=i+1; j < NAtoms; j++)
        {
            if (AtomTypes [par[i].idx] != AtomTypes [par[j].idx]) continue;

            //check it first
            if (Labels[j]<1) continue;

            x1=par[j].pos[0];
            if (fabs(x-x1) > eps_x) continue;

            //check Y
            y1=par[j].pos[1];
            if (fabs(y-y1) > eps_y) continue;

            //check Z
            z1=par[j].pos[2];
            if (fabs(z-z1) > eps_z) continue;
            each.push_back(par[j].idx);
            Labels[j]=0;
        }


        if (each.size() < 1) continue;

        //store
        SubGroups.resize(idx+1);
        for (uint k = 0;k< each.size(); k++)
             SubGroups[idx].push_back(each[k]);
        idx=idx+1;
    }



    //print
    for (uint i=0; i < SubGroups.size(); i++) {
        cout <<  setw(5) << "   " <<AtomTypes[SubGroups[i][0]]<<":";
        for (uint j=0; j < SubGroups[i].size(); j++){
            cout << setw(5)<<SubGroups[i][j]+1;
            if ((j+1)%15==0) cout << "\n        ";
        }
        cout << endl;
    }
    cout << "Done\n\n";
}


bool PointGroup::CheckMatrix( matrix3x3 & m)
{
/*
    cout << "\nCheck Matrix:\n";
    cout << "tolerance is: "<<tol<<endl;
    std::cout << setprecision(4);
    std::cout << m.GetRow(0) << endl;
    std::cout << m.GetRow(1) << endl;
    std::cout << m.GetRow(2) << endl<<endl;

    for (uint i=0 ; i < SubGroups.size() ;i++)
    {
        uint size=SubGroups[i].size();

        for (uint j=0 ; j < size ; j++)
        {
            cout << AtomTypes[SubGroups[i][j]] <<" ";
            vector3 vec0 = AtomCoordinates[SubGroups[i][j]]*0.529177;
            vector3 vec1 = m * vec0;
            cout << vec0 << " " << vec1 << endl;
        }
    }
*/

    vector3  vec1, vec2, vec3;
    bool find=false;
    for (uint i=0 ; i < SubGroups.size() ;i++)
    {
        uint size=SubGroups[i].size();

        for (uint j=0 ; j < size ; j++)
        {
            vec1 =  m * AtomCoordinates[SubGroups[i][j]]*0.529177;

            find=false;
            for (uint k=0 ; k < size ; k++)
            {
                vec2 =  AtomCoordinates[SubGroups[i][k]]*0.529177;
                //cout << vec2 << endl;
                vec3=vec1-vec2;

                if ( vec3.length() < tol ) {
                    find = true;
                    break;
                }
            }

            if (!find) return false;
        }
    }
    return true;
}



void PointGroup::CheckCenter()
{
    matrix3x3 m;
    m.Set(0,0,-1.0);
    m.Set(1,1,-1.0);
    m.Set(2,2,-1.0);
    hasI=CheckMatrix(m);
}

bool PointGroup::CheckSigma(vector3  vec)
{
    matrix3x3 OperateMatrix;
    vector3 vec0=vec.normalize();
    OperateMatrix.SetupReflectiontMatrix(vec0);
/*
    std::cout << OperateMatrix.GetRow(0) << endl;
    std::cout << OperateMatrix.GetRow(1) << endl;
    std::cout << OperateMatrix.GetRow(2) << endl<<endl;
*/
    return CheckMatrix(OperateMatrix);
}

bool PointGroup::CheckS2n(vector3 vec, uint order)
{
    if (vec.length() < 0.1) return false;
    vector3 vec0=vec.normalize();

    matrix3x3 OperateMatrix;
    OperateMatrix=SetupRotationReflectiontMatrix(vec0, 2*order);

    if  (CheckMatrix(OperateMatrix)) {
        return true;
    }

    return false;
}

bool PointGroup::CheckCn(vector3 vec, uint  order)
{
    if (vec.length() < 0.1) return false;

    vector3 vec0=vec.normalize();

    //cout << "Check order of axis: " <<order << " " << angle*180.0/3.141592654<<endl;
    //cout << "Check order of axis: " <<order << " " << angle<<endl;
    //cout << "axis : " << vec0 << endl;

    matrix3x3 OperateMatrix;
    OperateMatrix.SetupRotateMatrix(vec0, order);
    //std::cout << OperateMatrix.GetRow(0) << endl;
    //std::cout << OperateMatrix.GetRow(1) << endl;
    //std::cout << OperateMatrix.GetRow(2) << endl<<endl;

    if  (CheckMatrix(OperateMatrix)) {
        //cout << "succeed! " <<endl <<endl;
        return true;
    }
    //cout << "failed " <<endl <<endl;
    return false;
}



uint PointGroup::detectAxisOrder(vector3 vec)
{
    uint order0=1;
    uint maxOrder=1;

    for (uint i=0 ;i < SubGroups.size() ; i++){

        if (SubGroups[i].size()==1) continue;
        if (maxOrder < SubGroups[i].size() )
            maxOrder = SubGroups[i].size() ;
    }

    if(maxOrder<2) return 1;

    order0=maxOrder;
    while (order0)
    {
         //cout << order0 <<endl;
         if (CheckCn ( vec, order0 ) )
             break;

          order0=order0-1;
         if(order0==1) break;
    }

    return order0;
}

void PointGroup::detectPrincipalAxisOrder()
{
    order=detectAxisOrder(PrincipalAxis);
    refinePrincipalAxis();
    //cout << order <<endl;
}


void PointGroup::SearchMirror()
{
    double eps = tol;

    vector3 One;

    for (uint i=0 ; i < SubGroups.size() ;i++)
    {
        uint kh = SubGroups[i].size();

        if (kh==1)
            continue;

        else if (kh==2)
        {
            uint idx1=SubGroups[i][0],
                 idx2=SubGroups[i][1];

            One=AtomCoordinates[idx1]+AtomCoordinates[idx2];

            //pass through the center
            if (!isEquivalent(One, MassCenter ) ) {
                if (CheckSigma(One))
                    appendVector3 (SigmaV,One);
            }

            One=AtomCoordinates[idx1]-AtomCoordinates[idx2];
            // check
            if (!CheckSigma(One)) continue;
            appendVector3 (SigmaV,One);
        }

        else if (kh > 2)
        {

            for (uint j = 0; j < SubGroups[i].size(); j++)
            {
                for (uint k = j+1; k < SubGroups[i].size(); k++)
                {
                     uint idx1 = SubGroups[i][j],
                          idx2 = SubGroups[i][k];

                    One=AtomCoordinates[idx1]+AtomCoordinates[idx2];

                    //pass through the center
                    if (!isEquivalent(One, MassCenter ) ) {
                        // check
                        if (CheckSigma(One))
                            appendVector3 (SigmaV,One);
                    }

                    One=AtomCoordinates[idx1]-AtomCoordinates[idx2];

                    // check
                    if (!CheckSigma(One)) continue;
                    appendVector3 (SigmaV,One);
                }
            }
        }
    }
}


bool PointGroup::CombinationCn(vector < uint > & arr, uint data[],
                 uint start, uint end, uint index, uint r,
                 vector < vector3 > & container,
                 uint order0)
{
    // Current combination is ready to be printed

    if (index == r)
    {
        std::vector <uint > idx;
        idx.clear();

        for (uint j=0; j<r; j++)   idx.push_back(data[j]);

        //it is a regular polygon or not?
        if(!isRegularPolygon(idx))  return false;

        vector3 One=VZero;
        for (uint k=0;k<order0;k++)
            One +=AtomCoordinates[idx[k]];
         One /= double (order0);

        if (CheckCn(One,order0))
            appendVector3(container, One);

        //10 C3
        if (order0==3)
            if (container.size()==10) return true;

        //3 C4
        if (order0==4)
            if (container.size()==4)  return true;

        //6 C5
        if (order0==5)
            if (container.size()==6)  return true;

        return false;
    }

    bool isOK=false;
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        if ( CombinationCn(arr, data, i+1, end, index+1, r, container,order0) )
        {
            isOK=true;
            break;
        }
    }
    return isOK;
}


void PointGroup::SearchCn(uint order0, vector < vector3 > & container)
{
    double eps = tol;

    vector < uint > idx;
    idx.resize(order0);
    uint k;
    vector3 One;

    container.clear();

    //cout << "SubGroups size: " << SubGroups.size()<<endl;

    for (uint i=0 ; i < SubGroups.size() ;i++)
    {
        uint kh = SubGroups[i].size();

        if (kh==1)  {
            One=SubGroups[i][0];
            if (!CheckCn(One,order0)) continue;
            appendVector3(container, One);
            continue;
        }

        if (kh==order0)
        {
            //10 C3
            if (order0==3)
                if (container.size()==10) break;

            //3 C4
            if (order0==4)
                if (container.size()==4) break;

            //6 C5
            if (order0==5)
                if (container.size()==6) break;


            idx.clear();

            //it is a regular polygon or not?
            if(!isRegularPolygon(SubGroups[i])) continue;

            for (k=0;k<order0;k++)
                idx[k]=SubGroups[i][k];

            One=VZero;
            for ( k=0; k<order0; k++)
                One +=AtomCoordinates[idx[k]];

            One /= order0;

            if (!CheckCn(One,order0)) continue;

            appendVector3(container, One);
            continue;
        }

        if ( kh > order0 )
        {
            uint m=order0, n=kh-1;
            uint data[m];

            //selected m from n each time, and then check/store
            CombinationCn(SubGroups[i], data, 0, n, 0, m, container,order0);
        }
    }
}

//----------------------------------------------------------
void PointGroup::SearchC5()
{
    SearchCn(5,C5);
}

void PointGroup::SearchC3()
{
    SearchCn(3,C3);
}

void PointGroup::SearchC4()
{
    SearchCn(4,C4);
}


void PointGroup::SearchC2()
{
    double eps = tol;

    uint idx1,idx2;
    vector3 One;

    for (uint i=0 ; i < SubGroups.size() ;i++)
    {
        uint kh = SubGroups[i].size();

        if (kh==1)
        {
            One=SubGroups[i][0];

            //pass through the center
            if (isEquivalent(One, MassCenter ) )
                continue;

            if (!CheckCn(One,2)) continue;
            appendVector3(C2, One);

            continue;
        }

        if (kh==2)
        {
           idx1=SubGroups[i][0];
           idx2=SubGroups[i][1];

           //pass through the center or not?
           if (isEquivalent(AtomCoordinates[idx1] + AtomCoordinates[idx2], MassCenter) )
               One = (AtomCoordinates[idx1] - AtomCoordinates[idx2]);
           else
               One = (AtomCoordinates[idx1] + AtomCoordinates[idx2]);

           if (!CheckCn(One,2)) continue;
           appendVector3(C2, One);
           continue;
        }

        if (kh > 2)
        {
            for (uint j = 0; j < SubGroups[i].size(); j++)
            {
                for (uint k = j+1; k < SubGroups[i].size(); k++)
                {
                    idx1 = SubGroups[i][j];
                    idx2 = SubGroups[i][k];

                    //pass through the center or not?
                    if (isEquivalent(AtomCoordinates[idx1] + AtomCoordinates[idx2], MassCenter ) )
                        One = AtomCoordinates[idx1] - AtomCoordinates[idx2];
                    else
                        One = AtomCoordinates[idx1] + AtomCoordinates[idx2];

                    if (!CheckCn(One,2)) continue;
                    appendVector3(C2, One);
                }
            }
        }
    }

}


bool PointGroup::findFirstHorizontalC2(vector3 & firstOne)
{
    bool find=false;

    //cout << "MassCenter is: " << MassCenter <<endl;

    for (uint i=0 ; i < SubGroups.size() ;i++)
    {
        uint kh = SubGroups[i].size();

        //cout << "Kh: " <<kh<<endl;

        if (kh==1)
            continue;

        else if (kh==2)
        {
            uint idx1=SubGroups[i][0],
                 idx2=SubGroups[i][1];

            //pass through the center or not?
            if (isEquivalent(AtomCoordinates[idx1] + AtomCoordinates[idx2], MassCenter ) )
                firstOne = (AtomCoordinates[idx1] - AtomCoordinates[idx2]);
            else
                firstOne = (AtomCoordinates[idx1] + AtomCoordinates[idx2]);

            //orthnogal to PrincipalAxis or not?
            if ( !isOrthogonal(firstOne, PrincipalAxis)  ) continue;

            // check
            if (!CheckCn(firstOne,2)) continue;
            find=true;
            break;

            //cout << "inside kh==2"<<endl;
            //cout << idx1+1 << " " <<idx2+2<< " "<< firstOne << " SSSSSSSS" <<endl;
        }

        else if (kh > 2)
        {

            //selected two for n each time
            for (uint j = 0; j < SubGroups[i].size(); j++)
            {
                for (uint k = j+1; k < SubGroups[i].size(); k++)
                {
                    uint idx1=SubGroups[i][j],
                         idx2=SubGroups[i][k];


                    //pass through the center or not?
                    if (isEquivalent(AtomCoordinates[idx1] + AtomCoordinates[idx2], MassCenter ) )
                        firstOne = (AtomCoordinates[idx1] - AtomCoordinates[idx2]);
                    else
                        firstOne = (AtomCoordinates[idx1] + AtomCoordinates[idx2]);

                    //orthnogal to PrincipalAxis or not?
                    if ( !isOrthogonal(firstOne, PrincipalAxis)  ) continue;

                    // check
                    if (!CheckCn(firstOne,2)) continue;

                    //cout << idx1+1 << " " <<idx2+2<< " "<< firstOne << " SSSSSSSS" <<endl;

                    find=true;
                    //break;
                    return true;

                }
             //if (find) break;
            }
        }
        //if (find) break;
    }
    return false;
}



bool PointGroup::perceptHorizontal_C2(bool isCheck=true)
{
    if (order < 2) return false;

    Horizontal_C2.clear();

    //find the first one and then check the orther
    vector3 firstOne;

    //Now, we obtained the first C2
    if (findFirstHorizontalC2(firstOne) )
          firstOne.normalize();
    else{
          cout << "cannot find the first C2!\n\n";
          return false;
    }

    //check the others
    matrix3x3 rotatation;
    rotatation.SetupRotateMatrix(PrincipalAxis,order*2);

/*
    std::cout << rotatation.GetColumn(0) << endl;
    std::cout << rotatation.GetColumn(1) << endl;
    std::cout << rotatation.GetColumn(2) << endl<<endl;;
    cout << "Search horizontal C2 based on order " << order<<endl;
    cout << firstOne <<endl;
*/

    vector3 nextOne=firstOne;

    //cout << "order is "<< order<<endl;
    //if (isCheck)  cout << "check them one by one "<<endl;
    //else cout << "do not check"<<endl;

    for (int i=0 ; i < order ;i++)
    {
        if (i<1) {
            appendVector3(Horizontal_C2,firstOne);
            continue;
        }

        nextOne=rotatation*nextOne;

        //do not check
        if (!isCheck)
        {
            appendVector3(Horizontal_C2,nextOne);
            cout << " ok !" << Horizontal_C2.size () <<endl;
            cout << nextOne <<endl;
            continue;
        }

        if  (CheckCn(nextOne,2)) {
            appendVector3(Horizontal_C2,nextOne);
            //cout << " ok !" << Horizontal_C2.size () <<endl;
            //cout << nextOne <<endl;
            continue;
        }

        else {//wrong!
            //cout << "Something is wrong inside perceptHorizontal-C2!"<<endl;
            //cout << "Please set larger tolerance!"<<endl;
        }
    }

    return true;
}



void PointGroup::perceptVerticalMirror(vector < vector3 > & container)
{
    double eps = tol;

    container.clear();

    //find the first one and then check the orther
    vector3 firstOne;
    bool find=false;


    for (uint i=0 ; i < SubGroups.size() ;i++)
    {
        uint kh = SubGroups[i].size();

        if (kh==1)
            continue;

        else if (kh==2)
        {
            uint idx1=SubGroups[i][0],
                 idx2=SubGroups[i][1];

            firstOne = AtomCoordinates[idx1] - AtomCoordinates[idx2];

            //orthnogal to PrincipalAxis or not?
            if (dot(firstOne, PrincipalAxis) > eps) continue;


            // check
            if (!CheckSigma(firstOne)) continue;

            find=true;
        }

        else if (kh > 2)
        {

            //selected two for n each time
            for (uint j = 0; j < SubGroups[i].size(); j++)
            {
                for (uint k = j+1; k < SubGroups[i].size(); k++)
                {

                    uint idx1 = SubGroups[i][j],
                         idx2 = SubGroups[i][k];

                    firstOne=AtomCoordinates[idx1]-AtomCoordinates[idx2];
                    if (firstOne.length() < eps) continue;

                    //orthnogal to PrincipalAxis or not?
                    if (dot(firstOne, PrincipalAxis) > eps) continue;


                    // check
                    if (!CheckSigma(firstOne)) continue;

                    find=true;
                    break;
                }
                if (find) break;
            }
        }
        if (find) break;
    }

    if (!find) {
        cout << "Cannot find perceptVerticalMirror\n";
        return;
    }


    //only one mirror existed, Cs symmetry
    if (find && order <2)
    {
        appendVector3(container,firstOne);
        return;
    }

    //Now, we obtained the first sigma
    firstOne.normalize();

    //check the others
    matrix3x3 rotatation;
    rotatation.SetupRotateMatrix(PrincipalAxis,order*2);

    vector3 nextOne=firstOne;
    for (uint i=0 ; i < order ;i++)
    {
        if (i<1) {
            appendVector3(container,firstOne);
            continue;
        }

        nextOne=rotatation*nextOne;
        if  (CheckSigma(nextOne)) {
            appendVector3(container,nextOne);
            continue;
        }

        //wrong!
        cout << "\nSomething is wrong inside perceptMirrorV!\n\n";
    }
}


bool PointGroup::perceptMirrorH()
{
    return (CheckSigma(PrincipalAxis)) ;
}



void PointGroup::perceptVerticalMirror(bool isCheck=true)
{
    //general cases
    if (Horizontal_C2.size() < 1 ) {
        SigmaV.clear();
        perceptVerticalMirror(SigmaV);
        return;
    }

    printf("\nSearch vertical mirrors based on %d Horizontal C2 ...\n",Horizontal_C2.size());

    isDnh=false;
    isDnd=false;
    isDn=false;
    isCn=false;
    isCnv=false;
    isCi=false;
    isCs=false;

    //based on Horizontal C2
    vector3 firstOne;

    if (CheckSigma(PrincipalAxis)) {
       isDnh=true;
       firstOne = cross ( Horizontal_C2[0] , PrincipalAxis);
    }


    if (!isDnh){
        firstOne = cross ( (Horizontal_C2[0] + Horizontal_C2[1]), PrincipalAxis);

        if  (CheckSigma(firstOne))
            isDnd=true;
    }

    if ( (!isDnh) &&  (!isDnd) )
        isDn=true;


    if (isDn)  {
        cout << "Dn  PointGoup" <<endl;
        return;
    }

    if (isDnh) cout << "Dnh PointGoup" <<endl;
    if (isDnd) cout << "Dnd PointGoup" <<endl;

    matrix3x3 rotatation;
    rotatation.SetupRotateMatrix(PrincipalAxis,order*2);

    vector3 nextOne=firstOne;
    for (uint i=0 ; i < order ;i++)
    {
        if (i<1) {
            if (isDnh)  appendVector3(SigmaV,firstOne);
            if (isDnd)  appendVector3(SigmaD,firstOne);
            continue;
        }

        nextOne=rotatation*nextOne;

        if (!isCheck) {
            if (isDnh)  appendVector3(SigmaV,nextOne);
            if (isDnd)  appendVector3(SigmaD,nextOne);
            continue;
        }

        if  (CheckSigma(nextOne)) {
            if (isDnh)  appendVector3(SigmaV,nextOne);
            if (isDnd)  appendVector3(SigmaD,nextOne);
            continue;
        }
    }
    cout << "Done "<<endl;
}





void PointGroup::appendVector3(vector < vector3 > & container, vector3 & vec)
{
    vector3 vec0=vec.normalize();

    if (vec0.length() < 0.01) return;

    //dangerous
    //if (vec0[findMax(vec0)] < 0.0)   vec0=-vec0;

    if (container.size() < 1) {
        container.push_back(vec0);
        return;
    }

    if (isInside(container,vec))
        return;
    container.push_back(vec0);
}

void PointGroup::refinePrincipalAxis ()
{
    if (isEquivalent(PrincipalAxis,VX) ) PrincipalAxis=VX;
    else if (isEquivalent(PrincipalAxis,VY) ) PrincipalAxis=VY;
    else if (isEquivalent(PrincipalAxis,VZ) ) PrincipalAxis=VZ;

    if (PrincipalAxis[findMax(PrincipalAxis)] < 0.0)
        PrincipalAxis=-PrincipalAxis;
}




void PointGroup::PerceptSymmetry (bool isStandardOrientation)
{
    order=1;
    cout << "Percept symmetry based on tolerance of " << tol << " Bohr" <<endl;
    isOrientated=isStandardOrientation;

    Centralize();

    perceptMolShape();

    //buildSubGroups_perception_cluster();


    buildSubGroups_perception();


    vector <vector3> V;
    vector3 vec;


    CheckCenter();



    switch (shape)
    {
    case Line:
        PrincipalAxis = IMomentMatrix.GetColumn(0); //get the first vector
        PrincipalAxis.normalize();

        if (hasI) {
            PGSymbol="D00h";
            isDooh=true;
            SigmaH.clear();
            SigmaH.push_back(PrincipalAxis);
            cout << "D00h PointGoup !" <<endl;
            break;
        }
        else{
            PGSymbol="C00v";
            isCoov=true;
            cout << "C00v PointGoup !" <<endl;
            break;
        }
        break;

    case Plane:
        cout << "Plane shape!\n\n";
        if (hasI)
        { //Cnh,Dnh
            isCnh=true;
            PrincipalAxis = IMomentMatrix.GetColumn(2); //get the last vector
            PrincipalAxis.normalize();

            detectPrincipalAxisOrder();

            perceptHorizontal_C2(false);

            if (perceptMirrorH()){
                SigmaH.clear();
                SigmaH.push_back(PrincipalAxis);
            }
            perceptVerticalMirror(false);

            if (Horizontal_C2.size()>0){
                isDnh=true;
                cout << "Dnh Symmetry!\n\n";
            }
            else
            {
                isCnh=true;
                cout << "Cnh Symmetry!\n\n";
            }
            break;
        }
        else { //Cnv,Cs
            SigmaV.clear();
            for (uint i=0;i<3;i++)
                if (CheckSigma(IMomentMatrix.GetColumn(i)))
                    SigmaV.push_back(IMomentMatrix.GetColumn(i));


            if (SigmaV.size() == 1) {  //only one mirror
                isCs=true;
                cout <<"Cs symmetry\n\n";
            }



            if (SigmaV.size() ==2) {//C2v
                order=detectAxisOrder( cross (SigmaV[0],SigmaV[1]));
                isCnv=true;
                cout << "\nC2v PointGoup" <<endl;
            }

            cout << order << endl;
              cout << SigmaV.size()  <<endl;

              if(order<1)
               return;

            if (order >1) {
                PrincipalAxis = cross (SigmaV[0],SigmaV[1]);
                PrincipalAxis.normalize();
            }

        }
        break;

    case Polygon:

        PrincipalAxis = IMomentMatrix.GetColumn(2); //get the last vector
        PrincipalAxis.normalize();

        detectPrincipalAxisOrder();
        cout << " Polygon shape   :    Ia + Ib = Ic and Ia=Ib "<<endl;
        cout << " Principal axis : " << PrincipalAxis << endl;
        cout << " Order of principal axis : " <<order << endl<<endl;

        perceptHorizontal_C2(false);


        if (perceptMirrorH()){
            SigmaH.clear();
            SigmaH.push_back(PrincipalAxis);
        }
        perceptVerticalMirror(false);

        if (Horizontal_C2.size()>1) return;

        //Cnh
        if (SigmaH.size() ==1)
        {
            isCnh=true;
            cout << "Cnh Symmetry!\n\n";
            break;
        }

        isCn=true;

        break;

    case Oblate:                 //Ia = Ib < Ic
        PrincipalAxis = IMomentMatrix.GetColumn(2); //get the last vector
        PrincipalAxis.normalize();
        detectPrincipalAxisOrder();

        cout << " Oblate shape   :    Ia = Ib < Ic"<<endl;
        cout << " Principal axis : " << PrincipalAxis << endl;
        cout << " Order of principal axis : " <<order << endl<<endl;



        perceptOblate();




        if (isCnv){
            if (SigmaV.size()==order) break;
            cout << "Wrong result! The given tolerance is too large!"<<endl;
            order=SigmaV.size();
            if (order==1)
            {
                isCnv=false;
                isCs=true;
            }
        }
        break;

      case Prolate:               //Ia < Ib = Ic
        PrincipalAxis = IMomentMatrix.GetColumn(0); //get the first vector
        PrincipalAxis.normalize();
        detectPrincipalAxisOrder();

        cout << "Prolate shape  :   Ia < Ib = Ic"<<endl;
        cout << "Principal axis : " <<PrincipalAxis << endl;
        cout << "Order of principal axis : " <<order << endl<<endl;
        perceptProlate();

        if (isCnv){
            if (SigmaV.size()==order) break;
            cout << "Wrong result! The given tolerance is too large!"<<endl;
            order=SigmaV.size();
            if (order==1)
            {
                isCnv=false;
                isCs=true;
            }
        }
        break;

    case Sphere:

        /*
        //  test ok, please try!
        cout << "Perfrom brute-force algorithms to search symmetry elements ..." <<endl;
        SearchC2();cout << "there is "<<C2.size() << " C2"<<endl;
        SearchC3();cout << "there is "<<C3.size() << " C3"<<endl;
        SearchC4();cout << "there is "<<C4.size() << " C4"<<endl;
        SearchC5();cout << "there is "<<C5.size() << " C5"<<endl;
        SearchMirror();cout << "there is "<<SigmaV.size() << " SigmaV"<<endl;
        return;
        */


        perceptSphere();

        if (isT && hasI) {
            isTh=true;
            isT=false;
            perceptThMirror(false);
        }

        if (isT)
        {
           perceptTdMirror(true);
           if (SigmaD.size()>0) {
               isT=false;
               isTd=true;
           }
        }

        if (isO && hasI) {
            isOh=true;
            isO=false;
            perceptOhMirror(false);
        }

        if (isI && hasI) {
            isIh=true; isI=false;
            perceptIhMirror(false);
        }

        break;

    case Irregular:
        order = 1;
        for (uint i =0;i<3;i++)
        {
            PrincipalAxis = IMomentMatrix.GetColumn(i);
            PrincipalAxis.normalize();
            detectPrincipalAxisOrder();
            if (order == 1) continue;
            break;
        }

        //cout << "order inside irregular "<<order <<endl;

        if (order < 2) { //Cs,Ci,C1
            if (hasI) {
                isCi=true;
                Ci=VZero;
                cout << "\nCi point group\n\n";
                break;
            }

            SigmaV.clear();
            perceptVerticalMirror(SigmaV);
            if (SigmaV.size() >0)
            {
                isCs=true;
                cout << "There is only " << SigmaV.size()  << " Sigma V"<<endl;
                cout << "\nCs point group\n\n";
            }
            break;
        }


        cout << "Order of irregular shape: "<< order<<endl;

        isCnv=false;
        if (order > 1) {// find one vertical SigmaV
            perceptVerticalMirror(SigmaV);
            if (SigmaV.size() == order) {
                cout << "Cnv Point Group"<<endl;
                isCnv=true;
                break;
            }
            else { //special cases, the molecule is distorted
                generalPerceptSymmetry();
            }
        }


        SigmaH.clear();
        if (order >1){ // check it is Horizantal SigmaH
            if ( CheckSigma(PrincipalAxis) ) {
                 appendVector3(SigmaH,PrincipalAxis);
                 isCnh=true;
                 cout << "Cnh Point Group"<<endl;
                 break;
            }
        }

        //Cn, Sn, Cs, Ci are left
        isCs=false;
        SigmaH.clear();
        for (uint i = 0;i < 3; i++)
        {
            if ( CheckSigma(PrincipalAxis) ) {
                isCs = true;
                appendVector3(SigmaH,PrincipalAxis);
                break;
            }
        }

        if (C2.size() == 1){
            isCn=true;
            break;
        }


        if (hasI) {
            isCi=true;
            break;
        }

        isCn=true;

        //cout << " Irregular shape  :   Ia < Ib < Ic"<<endl;
        //cout << " Principal axis : " <<PrincipalAxis << endl;
        //cout << " Order of principal axis : " <<order << endl<<endl;

        break;
    }

    if (isCn)
        if (CheckS2n(PrincipalAxis,order))
        {
            isCn=false;
            isS2n=true;
            order=2*order;
        }

    //append S2n
    if (isDnd) {
        if (order==2) S4.push_back(PrincipalAxis);
        if (order==3) S6.push_back(PrincipalAxis);
        if (order==4) S8.push_back(PrincipalAxis);
        if (order==5) S10.push_back(PrincipalAxis);
        if (order==6) S12.push_back(PrincipalAxis);
        if (order>6)  S2n.push_back(PrincipalAxis);
    }



}

void PointGroup::perceptOblate()
{
    isCi=false;
    isCs=false;
    isCn=false;

    if (order < 2) {
        //Cs,Ci
        if (hasI) {
                isCi=true;
                Ci=VZero;
                cout << "\nCi point group\n\n";
                return;
        }

        SigmaV.clear();
        perceptVerticalMirror(SigmaV);

        if (SigmaV.size() > 0)
        {
             isCs=true;
             cout << "There is only " << SigmaV.size()  << " Sigma V"<<endl;
             cout << "\nCs point group\n\n";
             return;
        }

        return;
    }


    perceptHorizontal_C2(false);


    if (Horizontal_C2.size() > 0 )
    {
        isCnv=false;
        isCnh=false;
        isCn=false;

        if (perceptMirrorH()){
            SigmaH.clear();
            SigmaH.push_back(PrincipalAxis);
        }
        perceptVerticalMirror(false);
        return;
    }


    //Cnv and Cnh
    isCnv=false;
    isCnh=false;

    if (order > 1) { // find vertical SigmaV
        perceptVerticalMirror(SigmaV);
        if (SigmaV.size()>0)
        {
            isCnv=true;
            isCn=false;
            isCnh=false;
            return;
        }
    }


   //Cnh is left
   SigmaH.clear();
   if (order >1){ // check it is Horizantal SigmaH
            if ( CheckSigma(PrincipalAxis) ) {
                 appendVector3(SigmaH,PrincipalAxis);
                 isCnh=true;
                 isCnv=false;
                 isCn=false;
                 return;
                 cout << "Cnh Point Group"<<endl;
            }
    }

    //Cn is left
    isCn=true;
}


void PointGroup::perceptProlate()
{
    cout<<"PointGroup::perceptProlate()"<<endl;
    perceptOblate();
}


void PointGroup::perceptThMirror(bool isCheck=false)
{
    //has one center
    if (!hasI) {
      cout << "There one center in Th symmetry!\n";
      return;
    }

    uint size = C2.size();
    if ( size<3) {
        cout << "Something is wrong when check the three mirror of Th!\n";
        cout << "Only "<<size<<" is found!\n";
    }

    SigmaH.clear();
    for (uint i =0 ; i < size ; i ++)
    {
        if (isCheck)
            if ( CheckSigma(C2[i])) appendVector3(SigmaH,C2[i]);
        else
           appendVector3(SigmaH,C2[i]);
    }
}

void PointGroup::generalPerceptSymmetry()
{
    SearchC2();

    if (C2.size()==1) { //maybe C2
        PrincipalAxis=C2[0].normalize();
        order=2;
        perceptVerticalMirror(SigmaV);
        if (SigmaV.size()==2)
        {
            isCnv=true;
        }

        if (CheckSigma(C2[0])){
            isCnh=true;
        }

        if ( (!isCnv) && (!isCnh))
            isCn=true;

        return;
    }


    //SearchMirror();
    if (C2.size()==15) {
        //find two orth C2 out of 15
        vector3 c2_1=C2[0],c2_2;
        double ang=0.0;
        bool isfind=false;
        for (uint k=1;k<C2.size();k++)
        {
            c2_2=C2[k];
            ang=vectorAngle(c2_1,c2_2);
            if ( abs (ang-90.0) < 3.0) {
                 isfind = true;break;
            }
        }

        if (isfind) {
             generate_I_from_2C2_orth(c2_1, c2_2);
        }

        isI=true;
        if (SigmaV.size()>10) {
            isIh=true;
            isI=false;
            return;
        }
    }
}

void PointGroup::perceptTdMirror(bool isCheck=false)
{
    uint size = C3.size();

    if (size <4) {
        cout << "Something is wrong when check the three mirror of Th!\n";
        cout << "Only "<<size<<" is found!\n";
    }

    SigmaD.clear();
    vector3 vec;
    for (uint i =0 ; i < size ; i ++)
    {
        //all the combination of four C3
        for (uint j = i+1 ; j < size ; j ++)
        {
            vec = cross(C3[i],C3[j]);
            if (isCheck){
                if ( CheckSigma(vec)) appendVector3(SigmaD,vec);
            }
            else
               appendVector3(SigmaD,vec);
        }
    }

    if (SigmaD.size() <1) return;


    //append S4
    if (C2.size()>3) return; //maybe Oh symmetry

    S4.clear();
    for (uint i =0 ; i < C2.size(); i ++)
    {
        S4.push_back(C2[i]);
    }


}


void PointGroup::perceptOhMirror(bool isCheck=false)
{
    perceptTdMirror(isCheck);

    uint size=C4.size();

    if ( size<3) {
        cout << "Something is wrong when check the three mirror of Oh!\n";
        cout << "Only "<<C4.size()<<" is found!\n";
        return;
    }

    SigmaH.clear();
    for (uint i =0 ; i <size ; i ++)
    {
        if (isCheck) {
            if ( CheckSigma(C4[i]))
                appendVector3(SigmaH,C4[i]);
        }
        else
           appendVector3(SigmaH,C4[i]);
    }


    //append S4
    S4.clear();
    for (uint i =0 ; i < C4.size(); i ++)
    {
        S4.push_back(C4[i]);
    }

    //append S6
    S6.clear();
    for (uint i =0 ; i < C3.size(); i ++)
    {
        S6.push_back(C3[i]);
    }

}


void PointGroup::perceptIhMirror(bool isCheck=false)
{
    uint size = C5.size();
    if ( size < 6) {
        cout << "Something is wrong when check the three mirror of Ih!\n";
        cout << "Only "<<size<<" is found!\n";
    }


    SigmaD.clear();
    vector3 vec;
    for (uint i =0 ; i < size ; i ++)
    {
        //all the combination of four C3
        for (uint j = i+1 ; j < size ; j ++)
        {
            vec = cross(C5[i],C5[j]);
            if (isCheck) {
                if ( CheckSigma(vec))
                    appendVector3(SigmaD,vec);
            }
            else
                appendVector3(SigmaD,vec);
        }
    }

    //append S10
    S10.clear();
    for (uint i =0 ; i < C5.size(); i ++)
    {
        S10.push_back(C5[i]);
    }
}


void PointGroup::perceptSphere()
{
  //    I/Ih=[36,60,72,90,108,120,144,180]
  //    O/Oh=[45,60,90]
  //    T/Td/Th=[90]

  isTd = isT  = isTh=false;
  isIh = isI  = false;
  isO  = isOh = false;

  //CheckCenter();

  double eps=tol;

  //find two c2
  vector3 One;
  bool find=false;

  C2.clear();

  for (uint i=0 ; i < SubGroups.size() ;i++)
  {
      uint kh = SubGroups[i].size();

      if (kh<3)  continue;

      for (uint j = 0; j < SubGroups[i].size(); j++)
      {
          for (uint k = j+1; k < SubGroups[i].size(); k++)
          {
              uint idx1 = SubGroups[i][j],
                   idx2 = SubGroups[i][k];

              One=(AtomCoordinates[idx1]+AtomCoordinates[idx2])/2.0;

              if (!CheckCn(One,2)) continue;
              appendVector3(C2,One);

              if (C2.size() ==2) break;
          }
          if (C2.size() ==2) break;
      }
      if (C2.size() ==2) break;
  }


  if (C2.size() <2) {
      cout << " Something is wrong when search of C2 inside the sphere molecule!"<<endl;
      return;
  }

  //dectect based on the two c2
  vector3 c2_1=C2[0],
          c2_2=C2[1];


  double angle = vectorAngle(c2_1,c2_2);
  if (angle > 91.0)
  {
      angle=180.0-angle;
      c2_2=-C2[1];
  }

  //cout << "Angle of two C2 "<< angle <<endl;
 // I/Ih=[36,60,72,90,108,120,144,180]
 // O/Oh=[45,60,90]
 // T/Td/Th=[90]

  //Now, we obtained the two C2

  matrix3x3 rotatation;
  double err=3.0;

  //36 degree:  I/Ih=[36,60,72,90]  I/Ih=[36,60,72,90,108,120,144,180]
  if (abs (angle - 72.0) < err)  angle=36.0;
  if (abs (angle - 36.0) < err)
  {
        angle=36.0;

        isI=true;

        //C5 ^( two C2)
        vector3 vec=cross(c2_1, c2_2);
        vector3 c5=vec.normalize();

        //two neighbored C5
        vector3 mid = (c2_1 + c2_2)/2.0;
        vector3 axis=cross(c5, mid);
        rotatation.SetupRotateMatrix(axis,-63.435);
        vector3 c5_1= rotatation*c5;

        if (!CheckCn(c5_1,5)){
            rotatation.SetupRotateMatrix(axis,63.435);
            c5_1= rotatation*c5;
        }
        generate_I_from_2C5 (c5,c5_1);
        return;
  }


  //45 degree:  O/Oh=[45,60,90]
  if (abs (angle - 135.0) < err) { angle=45.0;
      c2_2=-c2_2;
  }

  if (abs (angle - 45.0) < err)
  {
        isO=true;
        angle=45.0;

        //generate all the symmetry elements (one of them is C4)
        if (CheckCn(c2_1,4))
            generate_O_from_C4C2(c2_1,c2_2);

        if (CheckCn(c2_2,4))
            generate_O_from_C4C2(c2_2,c2_1);

        return;
  }




  //60 degree
  // I/Ih=[36,60,72,90]  and O/Oh=[45,60,90]
  if (abs (angle - 120.0) < err) {
      angle = 60.0;
      c2_2=-c2_2;
  }

  if (abs (angle - 60.0) < err)
   {
        angle = 60.0;

        C3.clear();
        //generate the first C3
        vector3 c3_1=cross(c2_1, c2_2);
        C3.push_back(c3_1.normalize());


        //find the second C3 of Oh (C3^C3=70.529)
        vector3 mid = (c2_1 + c2_2)/2.0;
        vector3 axis = cross(mid , c3_1);


        rotatation.SetupRotateMatrix(axis, 70.529);
        vector3 c3_2 = rotatation*c3_1;

        if (!CheckCn(c3_2,3))
        {
            rotatation.SetupRotateMatrix(axis, -70.529);
            c3_2 = rotatation*c3_1;
        }

        isO=false;
        if (CheckCn(c3_2,3)) {
            C3.push_back(c3_2);
            isO=true;
        }

        //generall all the C2, C3 and C4 of O symmetry
        if (isO) {
            generate_O_from_2C3(C3[0],C3[1]);
            return;
        }

        //only I/Ih is left
        //the nearest C2 around C5, C2^C2=60
        isI=true; //at least

        mid=(c2_1 + c2_2)/2.0;
        vector3 c3=C3[0];
        axis = cross(mid , c3);

        //generate two C2
        rotatation.SetupRotateMatrix(axis,20.905);
        c2_1 = rotatation*c3;

        rotatation.SetupRotateMatrix(axis,-20.905);
        c2_2=rotatation*c3;

        if (CheckCn(c2_1,2))
            rotatation.SetupRotateMatrix(axis,-37.377368);

        if (CheckCn(c2_2,2))
            rotatation.SetupRotateMatrix(axis,37.377368);

        vector3 c5_1=rotatation*c3;

        rotatation.SetupRotateMatrix(c3,120.0);
        vector3 c5_2=rotatation*c5_1;

        c5_1.normalize();c5_2.normalize();
        generate_I_from_2C5(c5_1,c5_2); //neighbored C5

   }// end 60 degree


  //90 degree:  I/Ih=[36,60,72,90]  O/Oh=[45,60,90]  T/Td/Th=[90]
  if (abs (angle - 90.0) < err)
  {
        angle=90.0;

        //cout << "H " << c2_1 <<endl;
        //cout << "H " << c2_2 <<endl;

        //two C4
        if ( CheckCn(c2_1,4) && CheckCn(c2_2,4) ){
             generate_O_from_2C4(c2_1,c2_2);

             isO=true;
             return;
        }

        //one C4
        if ( CheckCn(c2_1,4)){
             generate_O_from_C4C2(c2_1,c2_2);
             isO=true;
             return;
        }

        //one C4
        if ( CheckCn(c2_2,4)){
             generate_O_from_C4C2(c2_2,c2_1);
             isO=true;
             return;
        }


        //two C2 left
        //Generate the third C2
        vector3 vec = cross(c2_1, c2_2);
        vector3 c2_3=vec.normalize();


        //it is C4
        if ( CheckCn(c2_3,4)){
            isO=true;
            generate_O_from_C4C2(c2_3,c2_1);
            return;
        }


        // T or I ?
        // test one c5^c2=58.283, 90-58.283
        rotatation.SetupRotateMatrix(c2_1,58.283);
        vec=rotatation*c2_3;
        if ( CheckCn(c2_3,5)){
            isI=true;
            generate_I_from_2C2_orth(c2_1,c2_2);
        }

        //only only T/Td/Th is left
        isT=true;
        appendVector3(C2,c2_1);
        appendVector3(C2,c2_2);
        appendVector3(C2,c2_3);

        //generate all C3
        vector3 c3_1 = ( c2_1 + c2_2 + c2_3).normalize();
        vector3 c3_2 = ( c2_1 + c2_2 - c2_3).normalize();
        vector3 c3_3 = ( c2_1 - c2_2 + c2_3).normalize();
        vector3 c3_4 = (-c2_1 + c2_2 + c2_3).normalize();
        appendVector3(C3,c3_1);
        appendVector3(C3,c3_2);
        appendVector3(C3,c3_3);
        appendVector3(C3,c3_4);
  }

}



void PointGroup::generate_I_from_2C2_orth(vector3 c2_1, vector3 c2_2)
{
    vector3 vec=cross(c2_1, c2_2);
    vector3 c2_3=vec.normalize();

    vector < vector3 > c5;
    c5.clear();

    matrix3x3 rotatation;

    rotatation.SetupRotateMatrix(c2_1,31.717);
    vec=rotatation*c2_3;
    if (CheckCn(vec,5))
        c5.push_back(vec);

    rotatation.SetupRotateMatrix(c2_1,-31.717);
    vec=rotatation*c2_3;
    if (CheckCn(vec,5))
        c5.push_back(vec);


    if (c5.size()<2)
    {
        rotatation.SetupRotateMatrix(c2_1,58.283);
        vec=rotatation*c2_3;
        if (CheckCn(vec,5))
            c5.push_back(vec);
    }

    if (c5.size()<2)
    {
        rotatation.SetupRotateMatrix(c2_1,-58.283);
        vec=rotatation*c2_3;
        if (CheckCn(vec,5))
            c5.push_back(vec);
    }


    if (c5.size()==2) {
        generate_I_from_2C5(c5[0],c5[1]);
    }
    else
    {
        cout << "Something is wrong inside generate_I_from_2C2_orth"<<endl;
    }


}


void PointGroup::generate_I_from_2C5(vector3 c5_1, vector3 c5_2)
{
    C2.clear();
    C3.clear();
    C4.clear();
    C5.clear();

    vector3 c2=(c5_1 + c5_2).normalize();

    matrix3x3 rotatation;
    rotatation.SetupRotateMatrix(c5_1,72.0);
    vector3 vec=c2;

    for (uint i=0;i<5;i++){
        vec=rotatation*vec;
        appendVector3(C2,vec);
    }

    //generate c5 around c5 , c5^c2=63.4349488
    appendVector3(C5,c5_1);
    appendVector3(C5,c5_2);

    //generate all the six c5 based on two c5
    rotatation.SetupRotateMatrix(c5_1,72.0);
    vec=C5[1];

    for (uint i=0;i<5;i++)
    {
        if (C5.size()==6)  break;
        vec=rotatation*vec;

        bool isfind=false;

        for (uint j=0;j<C5.size();j++)
        {
            if (isEquivalent(C5[j],vec))
            {
                isfind=true;break;
            }
        }

        if (isfind) continue;
        if (CheckCn(vec,5))  appendVector3(C5,vec);
    }

    if (C5.size() <6){
        cout << "Something is wrong inside generation of two neighbored C5 : generate_I_from_2C5"<<endl;
        return;
    }


    //generate fifteen c2 based on six c5, each mid-point of C5-C5 edge
    vector < vector3 > B12;
    B12.clear();
    for (uint i=0;i<6;i++){
        B12.push_back(C5[i]);
        B12.push_back(-C5[i]);
    }

    vector3 v1,v2,v;
    C2.clear();
    for (uint i=0; i < 12; i++)
    {
        if (C2.size()==15 )  break;
        v1=B12[i];

        for (uint j=0;j<12;j++)
        {
            if (C2.size()==15 )  break;

            v2=B12[j];

            //find the neighbored
            if ( (v1-v2).length() > 1.2)   continue;

            //mid-point of C5-C5 edge
            v = (v1+v2)/2.0;

            if (isInside(C2,v)) continue;

            if (CheckCn(v,2))
                appendVector3(C2,v);
         }
    }

    //generate c3 based on six c5
    C3.clear();
    vector3 c3_1=(c5_1+c5_2+C5[2]).normalize();
    rotatation.SetupRotateMatrix(c5_1,72.0);

    vec=c3_1;
    for (uint k=0;k<5;k++){
        vec=rotatation*vec; appendVector3(C3,vec);
    }


    rotatation.SetupRotateMatrix(c5_2,72.0);
    vec=c3_1;
    for (uint k=0;k<5;k++){
        vec=rotatation*vec; appendVector3(C3,vec);
    }

    rotatation.SetupRotateMatrix(c5_1,72.0);
    vec=C3[C3.size()-1];

    for (uint k=0;k<5;k++){
        if (k<1) continue;
        vec=rotatation*vec;
        appendVector3(C3,vec);
    }

    vec=C3[C3.size()-2];
    for (uint k=0 ; k<5; k++){
        if (k<1) continue;
        vec=rotatation*vec;
        appendVector3(C3,vec);
    }

    /*
    cout <<"C2:"<<endl;
    for (uint k=0;k<C2.size();k++)
        cout << "H " << 4.5* C2[k] <<endl;

      cout <<"\n\nC3:"<<endl;
    for (uint k=0;k<C3.size();k++)
        cout << "H " << 4.5* C3[k] <<endl;

        cout <<"\n\nC5:"<<endl;
    for (uint k=0;k<C5.size();k++)
        cout << "H " << 4.5* C5[k] <<endl;
        */
}

void PointGroup::generate_O_from_2C3(vector3 c3_1, vector3 c3_2)
{
    C2.clear();
    C3.clear();
    C4.clear();

    c3_1.normalize();
    c3_2.normalize();

    appendVector3(C3,c3_1);
    appendVector3(C3,c3_2);


    // generate all C3
    matrix3x3 rotatation;
    rotatation.SetupRotateMatrix(c3_1, 120.0);

    vector3 c3_3=rotatation*c3_2;
    vector3 c3_4=rotatation*c3_3;
    appendVector3(C3,c3_3);
    appendVector3(C3,c3_4);

    // generate all C4
    C4.clear();
    vector3 c4_1=(c3_2 + c3_3);
    vector3 c4_2=(c3_3 + c3_4);
    vector3 c4_3=(c3_2 + c3_4);
    appendVector3(C4,c4_1);
    appendVector3(C4,c4_2);
    appendVector3(C4,c4_3);

    // generate three C2 based on neighbored C3
    vector3 c2_1=(c3_1 + c3_2);
    vector3 c2_2=(c3_1 + c3_3);
    vector3 c2_3=(c3_1 + c3_4);

    appendVector3(C2,c2_1);
    appendVector3(C2,c2_2);
    appendVector3(C2,c2_3);

    rotatation.SetupRotateMatrix(c3_2,120.0);
    vector3 vec=c2_1;
    for (uint i=0 ; i<2;i++) {
        vec = rotatation*vec;
        appendVector3(C2,vec);
    }

    rotatation.SetupRotateMatrix(c3_3,120.0);
    vec=c2_2;
    for (uint i=0 ; i<2;i++) {
        vec = rotatation*vec;
        appendVector3(C2,vec);
    }
}

void PointGroup::generate_O_from_2C4(vector3 c4_1, vector3 c4_2)
{
    C2.clear();
    C3.clear();
    C4.clear();

    c4_1.normalize();
    c4_2.normalize();

    vector3 c4_3=cross(c4_1,c4_2);
    c4_3.normalize();

    appendVector3(C4,c4_1);
    appendVector3(C4,c4_2);
    appendVector3(C4,c4_3);

    // generate all C3
    vector3  c3_1=(c4_1 + c4_2 + c4_3).normalize();
    vector3  c3_2=(c4_1 + c4_2 - c4_3).normalize();
    vector3  c3_3=(c4_1 - c4_2 + c4_3).normalize();
    vector3  c3_4=(-c4_1 + c4_2 + c4_3).normalize();

    appendVector3(C3,c3_1);
    appendVector3(C3,c3_2);
    appendVector3(C3,c3_3);
    appendVector3(C3,c3_4);

    //generate all C2/45
    matrix3x3 rotatation;
    rotatation.SetupRotateMatrix(c4_1,45.0);
    vector3 vec = rotatation*c3_1;
    appendVector3(C2,vec);

    /*
    cout <<"C2:"<<endl;
    for (uint k=0;k<C2.size();k++)
        cout << "H " << 4.5* C2[k] <<endl;

      cout <<"\n\nC3:"<<endl;
    for (uint k=0;k<C3.size();k++)
        cout << "H " << 4.5* C3[k] <<endl;

    cout <<"C4:"<<endl;
    for (uint k=0;k<C4.size();k++)
        cout << "H " << 4.5* C4[k] <<endl;
        */
}

void PointGroup::generate_O_from_C4C2(vector3 c4, vector3 c2)
{
    C2.clear();
    C3.clear();
    C4.clear();

    c4.normalize();
    c2.normalize();

    appendVector3(C2,c2);

    // produce another C2 around C4
    matrix3x3 rotatation;
    rotatation.SetupRotateMatrix(c4,90.0);
    vector3 c2_2 = rotatation*c2;
    appendVector3(C2,c2_2);

    vector3 c4_2,c4_3;

    appendVector3(C4,c4);
    c4_2=(c2+c2_2).normalize();
    if (CheckCn(c4_2,4)) {
        c4_3=(c2-c2_2).normalize();
        appendVector3(C4,c4_2);

        appendVector3(C4,c4_3);
    }

    else{// produce the other two C4 around C2
        rotatation.SetupRotateMatrix(c2,180.0);
        vector3 c4_2=rotatation*c4;
        appendVector3(C4,c4_2);

        rotatation.SetupRotateMatrix(c2_2,180.0);
        vector3 c4_3=rotatation*c4;
        appendVector3(C4,c4_3);
    }

/*
    cout << "C2\n";
    cout << "H " << c2<<endl;
    cout << "H " << c2_2<<endl;
    cout << C2.size()<<endl;

    cout << "C4\n";
    cout << "H " << c4<<endl;
    cout << "H " << c4_2<<endl;
    cout << "H " << c4_3<<endl;

    cout << C4.size()<<endl;

*/


    // generate all C3 based on C4
    vector3 c3_1=(c4 + c4_2 + c4_3);
    vector3 c3_2=(c4 + c4_2 - c4_3);
    vector3 c3_3=(c4 - c4_2 + c4_3);
    vector3 c3_4=(-c4 + c4_2 + c4_3);

    appendVector3(C3,c3_1);
    appendVector3(C3,c3_2);
    appendVector3(C3,c3_3);
    appendVector3(C3,c3_4);

    // generate another two c2
    rotatation.SetupRotateMatrix(c3_1,-120.0);
    vector3 vec=rotatation*C2[1];
    appendVector3(C2,vec);

    rotatation.SetupRotateMatrix(c3_2,-120.0);
    vec=rotatation*C2[0];
    appendVector3(C2,vec);

    rotatation.SetupRotateMatrix(c4,90.0);
    for (uint k=0;k<4;k++){
        vec=rotatation*vec;
        if (!CheckCn(vec,2)) continue;
        appendVector3(C2,vec);
    }

/*
    cout <<"C2:"<<endl;
    for (uint k=0;k<C2.size();k++)
        cout << "H " << 2* C2[k] <<endl;

    cout <<"\n\nC3:"<<endl;
    for (uint k=0;k<C3.size();k++)
        cout << "H " << 2* C3[k] <<endl;

    cout <<"C4:"<<endl;
    for (uint k=0;k<C4.size();k++)
        cout << "H " << 2* C4[k] <<endl;
*/
}





void PointGroup::summary()
{
    stringstream ss;
    string sorder;

    if (isTd){
        PGSymbol = "Td";return;
    }
    else if (isTh){
        PGSymbol = "Th";return;
    }
    else if (isT){
        PGSymbol = "T";return;
    }

    else if (isOh){
        PGSymbol = "Oh";return;
    }

    else if (isO){
        PGSymbol = "O";return;
    }

    else if (isIh){
        PGSymbol = "Ih";
    }

    else if (isI){
        PGSymbol = "I";return;
    }


    else if (isDnh){
        ss<<"D"<<order<<"h"; PGSymbol=ss.str();
        return;
    }
    else if (isDnd){
        ss<<"D"<<order<<"d";  PGSymbol=ss.str();return;
    }

    else if (isDn){
        ss<<"D"<<order;  PGSymbol=ss.str();return;
    }


    else if (isCnv){
        ss<<"C"<<order<<"v";  PGSymbol=ss.str();return;
    }

    else if (isCnh){
        ss<<"C"<<order<<"h";  PGSymbol=ss.str();return;
    }

    else if (isCn){
        ss<<"C"<<order;  PGSymbol=ss.str();return;
    }


    else if (isCoov)
    {
        PGSymbol = "Coov";
        return;
    }
    else if (isDooh) {
        PGSymbol = "Dooh";return;
    }


    else if (isCi){
        PGSymbol = "Ci";return;
    }

    else if (isCs){
        PGSymbol = "Cs";return;
    }

    else if (isS2n){
       ss<<"S"<<order;  PGSymbol=ss.str();
    }

    else
        PGSymbol = "C1";

}



void PointGroup::report(string  pgsymbol)
{
    if(pgsymbol.length()<1)
        pgsymbol=PGSymbol;

    refineSymmetryElements();
    if (isOrientated)
        StandardOrientation();

    summary();
    pgsymbol=PGSymbol;

    std::cout.setf(ios::fixed);
    std::cout << setprecision(4);

    cout << "\n\n                   "<<PGSymbol<<endl;
    cout << "     Summary of Symmetric Elements\n"<<endl;

    if (PrincipalAxis.length()>0.1){
        cout << "-----------  Principal Axis  -----------"<<endl;
        printVector3 (PrincipalAxis);
        cout<< endl;
    }
    if (Horizontal_C2.size() >0) {

        cout << "-----------  Horizontal  C2  -----------"<<endl;
        for (uint i=0; i < Horizontal_C2.size(); i++){
            printVector3 (Horizontal_C2[i]);
        }
        cout << endl;
    }

    if (C2.size() >0) {
        cout << "---------------    C2    ---------------"<<endl;
        for (uint i=0; i < C2.size(); i++){
            printVector3 (C2[i]);
        }
        cout << endl;
    }

    if (C3.size() >0) {
        cout << "---------------    C3    ---------------"<<endl;
        for (uint i=0; i < C3.size(); i++){
            printVector3 (C3[i]);
        }
        cout << endl;
    }

    if (C4.size() >0) {
        cout << "---------------    C4    ---------------"<<endl;
        for (uint i=0; i < C4.size(); i++){
            printVector3 (C4[i]);
        }
        cout << endl;
    }

    if (C5.size() >0) {
        cout << "---------------    C5    ---------------"<<endl;
        for (uint i=0; i < C5.size(); i++){
            printVector3 (C5[i]);
        }
        cout << endl;
    }

    if (S4.size() >0) {
        cout << "---------------    S4    ---------------"<<endl;
        for (uint i=0; i < S4.size(); i++){
            printVector3 (S4[i]);
        }
        cout << endl;
    }

    if (S6.size() >0) {
        cout << "---------------    S6    ---------------"<<endl;
        for (uint i=0; i < S6.size(); i++) {
            printVector3 (S6[i]);
        }
        cout << endl;
    }

    if (S8.size() >0) {
        cout << "---------------    S8    ---------------"<<endl;
        for (uint i=0; i < S8.size(); i++){
            printVector3 (S8[i]);
        }
        cout << endl;
    }

    if (S10.size() >0) {
        cout << "---------------   S10    ---------------"<<endl;
        for (uint i=0; i < S10.size(); i++) {
            printVector3 (S10[i]);
        }
        cout << endl;
    }

    if (S12.size() >0) {
        cout << "---------------   S12    ---------------"<<endl;
        for (uint i=0; i < S12.size(); i++) {
            printVector3 (S12[i]);
        }
        cout << endl;
    }

    if (S2n.size() >0) {
        cout << "---------------   S2n    ---------------"<<endl;
        for (uint i=0; i < S2n.size(); i++){
            printVector3 (S2n[i]);
        }
        cout << endl<< endl;
    }

    if (SigmaH.size() >0) {
        cout << "--------------- Mirror-H ---------------"<<endl;
        for (uint i=0; i < SigmaH.size(); i++)
            printVector3 (SigmaH[i]);

        cout << endl;
    }

    if (SigmaV.size() >0) {
        cout << "--------------- Mirror-V ---------------"<<endl;
        for (uint i=0; i < SigmaV.size(); i++) {
            printVector3 (SigmaV[i]);
        }
        cout << endl;
    }

    if (SigmaD.size() >0) {
        cout << "--------------- Mirror-D ---------------"<<endl;
        for (uint i=0; i < SigmaD.size(); i++){
            printVector3 (SigmaD[i]);
        }
        cout << endl;
    }

    cout << "More infomation about symmetry, please go to http://symmetry.otterbein.edu/gallery/index.html"<<endl;
}

void PointGroup::refineSymmetryElements()
{
    refineVector3(PrincipalAxis);

    uint i=0;
    if (C2.size() >0)
        for ( i=0; i < C2.size(); i++)
            refineVector3(C2[i]);

    if (Horizontal_C2.size() >0)
        for ( i=0; i < Horizontal_C2.size(); i++)
            refineVector3(Horizontal_C2[i]);

    if (C2.size() >0)
        for ( i=0; i < C2.size(); i++)
            refineVector3(C2[i]);

    if (C3.size() >0)
        for ( i=0; i < C3.size(); i++)
            refineVector3(C3[i]);

    if (C4.size() >0)
        for (uint i=0; i < C4.size(); i++)
            refineVector3(C4[i]);

    if (C5.size() >0)
        for ( i=0; i < C5.size(); i++)
            refineVector3(C5[i]);

    if (S4.size() >0)
        for ( i=0; i < S4.size(); i++)
            refineVector3(S4[i]);

    if (S6.size() >0)
        for ( i=0; i < S6.size(); i++)
            refineVector3(S6[i]);

    if (S8.size() >0)
        for ( i=0; i < S8.size(); i++)
            refineVector3(S8[i]);

    if (S10.size() >0)
        for ( i=0; i < S10.size(); i++)
            refineVector3(S10[i]);

    if (S12.size() >0)
        for ( i=0; i < S12.size(); i++)
            refineVector3(S12[i]);

    if (S2n.size() >0)
        for ( i=0; i < S2n.size(); i++)
            refineVector3(S2n[i]);

    if (SigmaH.size() >0)
        for ( i=0; i < SigmaH.size(); i++)
            refineVector3(SigmaH[i]);

    if (SigmaV.size() >0)
        for ( i=0; i < SigmaV.size(); i++)
            refineVector3(SigmaV[i]);

    if (SigmaD.size() >0)
        for ( i=0; i < SigmaD.size(); i++)
            refineVector3(SigmaD[i]);

}







void PointGroup::StandardOrientation()
{
    rotMatrix.UnitMatrix();
    rotMatrix1.UnitMatrix();
    rotMatrix2.UnitMatrix();

    if (isT || isTd || isTh)
        Orientation_T();

    else if (isI || isIh )
        Orientation_I();

    else if (isO ||isOh )
        Orientation_O();

    else  if (isDn || isDnd || isDnh )
        Orientation_Dn();

    else  if (isS2n)
        Orientation_S2n();

    else  if (isCn)
        Orientation_Cn();

    else  if (isCnh)
        Orientation_Cnh();

    else  if (isCnv)
        Orientation_Cnv();

    else  if (isCoov )
        Orientation_Line();

    else  if (isDooh )
        Orientation_Line();

    else  if (isCs)
        Orientation_Cs();


    //rotMatrix
    for (uint i = 0; i < NAtoms; ++ i) {
        vector3 vec = rotMatrix*AtomCoordinates[i];
        AtomCoordinates[i]=vec;
        //cout << AtomTypes[i]<<" " << vec*0.529177 <<endl;
     }
}

void PointGroup::Orientation_Line()
{
    Orientation_Cn();
}

void PointGroup::Orientation_Ci()
{
    return;
}

void PointGroup::Orientation_Cs()
{
    vector3 vec0;
    if (SigmaH.size() > 0 )
        vec0=SigmaH[0];
    else if (SigmaV.size() > 0 )
        vec0=SigmaV[0];

    double angle=vectorAngle(VY,vec0);

    if( (fabs(angle)>179.0) || (fabs(angle) < 1.0)) return;
    vector3 vec=cross(vec0,VY);
    rotMatrix.SetupRotateMatrix(vec,angle);
}

matrix3x3 PointGroup::OrientationFromTo(vector3 Vec1, vector3 Vec2)
{
    matrix3x3 mat;
    mat.UnitMatrix();

    double angle=vectorAngle(Vec1,Vec2);

    if( angle > 178.0 || angle < 2.0 ) return mat;

     vector3 vec=cross(Vec1,Vec2);
     mat.SetupRotateMatrix(vec,angle);
     return mat;
}



void PointGroup::Orientation_Cn()
{
    rotMatrix = OrientationFromTo(PrincipalAxis,VZ);
}


void PointGroup::Orientation_S2n()
{
    Orientation_Cn();
}


void PointGroup::Orientation_Cnv()
{
    cout << "Orientation_Cnv" <<endl;

    vector3 sigma;
    bool isfind=false;
    if (isEquivalent (PrincipalAxis,VZ)) { //check it has done or not?
        for (uint i=0;i<SigmaV.size();i++)
        {
            sigma=SigmaV[i];
            if (isEquivalent (sigma,VY)) {
                isfind = true;
                break;
            }
        }
    }

    if (isfind) return;

    Orientation2ZY(PrincipalAxis,SigmaV[0]);
}


void PointGroup::Orientation_Cnh()
{
    //cout << "Orientation_Cnh"<<endl;
    Orientation_Cn();
}


void PointGroup::Orientation_Dn()
{
    vector3 c2;
    bool isfind=false;
    if (isEquivalent (PrincipalAxis,VZ)) { //check it has done or not?
        for (uint i=0;i<Horizontal_C2.size();i++)
        {
            c2=Horizontal_C2[i];
            if (isEquivalent (c2,VY)) {
                isfind = true;
                break;
            }
        }
    }

    if (isfind) return;
    Orientation2ZX(PrincipalAxis,Horizontal_C2[0]);
}

void PointGroup::Orientation_Dnd()
{
    Orientation_Dn();
}

void PointGroup::Orientation_Dnh()
{
    Orientation_Dn();
}



void PointGroup::PointGroup::Orientation_T()
{
    vector3 vec1=C2[0],vec2=C2[1];
    Orientation2ZX(vec1,vec2);
}


void PointGroup::Orientation_O()
{
    vector3 vec1=C4[0],vec2=C4[1];
    Orientation2ZX(vec1,vec2);
}

void PointGroup::Orientation_I()
{
    rotMatrix.UnitMatrix();
    rotMatrix2.UnitMatrix();
    rotMatrix1.UnitMatrix();

    cout << "Orientation_I"<<endl;
    vector3 vec1=C5[0],vec2=C5[1];

    rotMatrix1 = OrientationFromTo(vec1,VZ);
    vec2 = rotMatrix1 * vec2;


    angle=vectorAngle(vec2,VX);

    vector3 vec=VX;
    if(vec2.z()<0) vec=-VX;

    vec2.SetZ(0.0);
    rotMatrix2 = OrientationFromTo(vec2,vec);

    rotMatrix=rotMatrix2*rotMatrix1;
}

void PointGroup::Orientation2ZY(vector3 vec1, vector3 vec2)
{
    double angle=vectorAngle(vec1,VZ);

    //cout << angle<<"-------------"<<endl;

    //std::cout << rotMatrix1.GetRow(0) << endl;
    //std::cout << rotMatrix1.GetRow(1) << endl;
    //std::cout << rotMatrix1.GetRow(2) << endl<<endl;

    vector3 vec;
    if(angle > 1.0 && angle < 179.0 )
    {
        vec = cross (vec1,VZ);
        rotMatrix1.SetupRotateMatrix(vec,angle);
        vec2 = rotMatrix1 * vec2;
    }
/*
    std::cout << rotMatrix1.GetRow(0) << endl;
    std::cout << rotMatrix1.GetRow(1) << endl;
    std::cout << rotMatrix1.GetRow(2) << endl<<endl;
*/

    angle=vectorAngle(vec2,VY);
    if(angle>1.0 && angle<179.0 )
    {
        vec=cross (vec2,VY);
        rotMatrix2.SetupRotateMatrix(vec,angle);
    }

    /*
    cout << angle<<"-------------"<<endl;
    std::cout << rotMatrix2.GetRow(0) << endl;
    std::cout << rotMatrix2.GetRow(1) << endl;
    std::cout << rotMatrix2.GetRow(2) << endl<<endl;
*/
    rotMatrix=rotMatrix2*rotMatrix1;

/*
    for (uint i = 0; i < NAtoms; ++ i) {
        vector3 vec=rotMatrix*AtomCoordinates[i];
        AtomCoordinates[i]=vec;
        cout << AtomTypes[i]<<" " << vec*0.529177 <<endl;
     }
*/
}

void PointGroup::Orientation2ZX(vector3 vec1, vector3 vec2)
{
    double angle=vectorAngle(vec1,VZ);
    vector3 vec;
    if(angle > 1.0 && angle < 179.0 )
    {
        vec = cross (vec1,VZ);
        rotMatrix1.SetupRotateMatrix(vec,angle);
        vec2 = rotMatrix1 * vec2;
    }

    angle=vectorAngle(VX,vec2);
    if(angle>1.0 && angle<179.0 )
    {
        vec=cross (vec2,VX);
        rotMatrix2.SetupRotateMatrix(vec,angle);
    }
    rotMatrix=rotMatrix2*rotMatrix1;

    /*
    for (uint i = 0; i < NAtoms; ++ i) {
        vector3 vec=rotMatrix*AtomCoordinates[i];
        AtomCoordinates[i]=vec;
        cout << AtomTypes[i]<<" " << vec*0.529177 <<endl;
     }
     */
}



















//--------------------------------------------
void PointGroup::buildOperation()
{
    PGOperation.clear();

    if (isT) {
        buildOperation_T();
        return;
    }

    else if (isTd) {
        buildOperation_Td();
        return;
    }

    else if (isTh) {
        buildOperation_Th();
        return;
    }

    else if (isO) {
        buildOperation_O();
        return;
    }

    else if (isOh){
        buildOperation_Oh();
        return;
    }

    else if (isI) {
        buildOperation_I();
        return;
    }

    else if (isIh) {
        buildOperation_Ih();
        //buildOperation_I();
        return;
    }

    else  if (isDn) {
        buildOperation_Dn(order);
        return;
    }

    else  if (isDnd){
        buildOperation_Dnd(order);
        return;
    }

    else  if (isDnh) {
        buildOperation_Dnh(order);
        return;
    }


    else  if (isS2n) {
        buildOperation_S2n(order);
        return;
    }

    else  if (isCn) {
        buildOperation_Cn(order);
        return;
    }

    else  if (isCnv) {
        buildOperation_Cnv(order);
        return;
    }

    else  if (isCnh) {
        buildOperation_Cnh(order);
        return;
    }


    else  if (isCoov) {
        buildOperation_C00v();
        return;
    }

    else  if (isDooh ) {
        buildOperation_D00h();
        return;
    }

    else  if (isCi) {
       buildOperation_Ci();
       return;
    }

    else if (isCs)
        buildOperation_Cs();
}

void PointGroup::buildOperation_C00v()
{

}


void PointGroup::buildOperation_Cs()
{
    //E
    matrix3x3  matrix(1.0);
    PGOperation.push_back(matrix);

    //XOY
    matrix.Set(0,0,1.0);
    matrix.Set(1,1,-1.0);
    matrix.Set(2,2,01.0);
    PGOperation.push_back(matrix);
}

//----------------------------------------------
void PointGroup::buildOperation_Ci()
{
    //E
    matrix3x3 matrix(1.0);
    PGOperation.push_back(matrix);

    //I
    matrix.Set(0,0,-1.0);
    matrix.Set(1,1,-1.0);
    matrix.Set(2,2,-1.0);
    PGOperation.push_back(matrix);
}

//---------------------------------------
void PointGroup::buildOperation_Cn(uint n)
{
    matrix3x3 matrix;
    uint i;
    angle=360.0/double(n);
    double ang=0.0;

    for(i=0;i<n;i++)
    {
       matrix.rotateZAngle(ang);
       PGOperation.push_back(matrix);
       ang+=angle;
    }
}


//---------------------------------------
void PointGroup::buildOperation_S2n(uint n)
{
    matrix3x3 matrix;
    int i;
    angle=360.0/double(2.0*n);
    double ang=0.0;

    int j=1;
    for(i=0;i<n;i++)
    {
       matrix.rotateZAngle(ang);
       matrix.Set(2,2,j);

       PGOperation.push_back(matrix);
       ang+=angle;
       j=-j;
    }
}

//---------------------------------------
void PointGroup::buildOperation_Cnv(uint n)
{
    matrix3x3 matrix;
    matrix3x3 mirror; mirror.reflectionXOZ();

    angle=360.0/double(n);
    double ang=0.0;

    //!{E,Cn(1),Cn(2),Cn(3)...  Cn(n-1)}
    uint i;
    for(i=0;i<n;i++)
    {
       matrix.rotateZAngle(ang);
       PGOperation.push_back(matrix);
       ang+=angle;
    }

    //!{Sigma,sigma*Cn(1),sigma*Cn(2)...  sigma*Cn(n-1)}
    for(i=0;i<n;i++)
    {
       matrix=PGOperation[i]*mirror;
       PGOperation.push_back(matrix);
    }

}

//---------------------------------------
void PointGroup::buildOperation_Cnh(uint n)
{
    matrix3x3 matrix;
    angle=360.0/double(n);
    double ang=0.0;

    //!{E,Cn(1),Cn(2),Cn(3)...  Cn(n-1)}
    uint i;
    for(i=0;i<n;i++)
    {
       matrix.rotateZAngle(ang);

       PGOperation.push_back(matrix);
       ang+=angle;
    }

    //!{Sigma,sigma*Cn(1),sigma*Cn(2)...  sigma*Cn(n-1)}
    for(i=0;i<n;i++)
    {
       matrix=PGOperation[i];
       matrix.Set(2,2,-1.0);

       PGOperation.push_back(matrix);
    }

}

//---------------------------------------
void PointGroup::buildOperation_Dn(uint n)
{
    matrix3x3 matrix;
    uint i;
    angle=360.0/double(n);
    double ang=0.0;

    //{E,Cn(1),Cn(2),Cn(3)...  Cn(n-1)}
    for(i=0;i<n;i++)
    {
       matrix.rotateZAngle(ang);
       PGOperation.push_back(matrix);
       ang+=angle;
    }

    matrix3x3 c2;
    c2.rotateXAngle(180.0);

    //!{C2,C2*Cn(1),C2*Cn(2)...  C2*Cn(n-1)}
    for(i=0;i<n;i++)
    {
       matrix=PGOperation[i]*c2;
       PGOperation.push_back(matrix);
    }
}

//--------------------------------------------------
void PointGroup::buildOperation_Dnd(uint n)
{
    matrix3x3 matrix;
    int i;
    angle=360.0/double(n);
    double	ang=180.0/double(n);

    buildOperation_Dn(n);

    matrix3x3 mirror;
    mirror.reflectionXOZwithAngle(ang/2.0);
    //!{Sigma,sigma*Cn(1),sigma*Cn(2)...  sigma*Cn(n-1)}
    for(i=0;i<n;i++)
    {
       matrix=PGOperation[i]*mirror;
       PGOperation.push_back(matrix);
    }


     //!{E,S2n(1),S2n(3),S2n(5)...  S2n(2n-1)}
    for(i=0;i<n;i++)
    {
       matrix.rotateZReflection(ang);
       PGOperation.push_back(matrix);
       ang+=angle;
    }
}

//--------------------------------------------
void PointGroup::buildOperation_Dnh(uint n)
{
    matrix3x3 matrix;
    int i;
    angle=360.0/double(n);
    double ang=0.0;

    buildOperation_Dn(n);

    matrix3x3 mirror;mirror.reflectionXOZ();
    //!{SigmaV,SigmaV*Cn(1),SigmaV*Cn(2)...  SigmaV*Cn(n-1)}
    for(i=0;i<n;i++)
    {
       matrix=PGOperation[i]*mirror;
       PGOperation.push_back(matrix);
    }

    //!{SigmaH,SigmaH*Cn(1),SigmaH*Cn(2)...  SigmaH*Cn(n-1)}
    for(i=0;i<n;i++)
    {
       matrix=PGOperation[i];
       matrix.Set(2,2,-1.0);

       PGOperation.push_back(matrix);
    }
}

//-----------------------------------------------------------
void PointGroup::buildOperation_T()
{
    matrix3x3 matrix,C2;
    int i;
    double ang=0.0;
    vector3 C3(1.0,1.0,1.0);

    //!E
    matrix.Set(0,0,1.0);matrix.Set(1,1,1.0);matrix.Set(2,2,1.0);
        PGOperation.push_back(matrix);

    //!three C2
    C2.rotateZAngle(180.0);PGOperation.push_back(C2);
    C2.rotateXAngle(180.0);PGOperation.push_back(C2);
    C2.rotateYAngle(180.0);PGOperation.push_back(C2);

        //!{E,C3(1),C3(2)}=(1 1 1)
        C3.Set(1.0,1.0,1.0);
    for(int j=0;j<2;j++)
    {
            ang+=120.0;
            matrix.SetupRotateMatrix(C3,ang);
        for(i=0;i<4;i++)
        {
           PGOperation.push_back(matrix*PGOperation[i]);
        }
    }
}
//-------------------------------------------------------
void PointGroup::buildOperation_Th()
{
    //!Th=T*Ci;
    int i=0;
    buildOperation_T();

    matrix3x3 mirror,matrix,Ci;
    Ci.Set(0,0,-1.0);Ci.Set(1,1,-1.0);Ci.Set(2,2,-1.0);
    for(i=0;i<12;i++)
    {
       matrix=PGOperation[i];
       PGOperation.push_back(Ci*matrix);
    }
}

//-------------------------------------------------------
void PointGroup::buildOperation_Td()
{
   buildOperation_T();
    vector3 norm(1.0,1.0,0.0);
    matrix3x3 matrix, mirror;
    mirror.SetupReflectiontMatrix(norm);

    for(int i=0;i<12;i++)
    {
       matrix=PGOperation[i];
       PGOperation.push_back(mirror*matrix);
    }
}

//--------------------------------------------
void PointGroup::buildOperation_O()
{
    matrix3x3 matrix;
    int i,j;
    double ang=0.0;

    //c3=(1 1 1)
    vector3 C3(1.0,1.0,1.0);
    for(i=0;i<3;i++)
    {
      matrix.SetupRotateMatrix(C3,ang);
      PGOperation.push_back(matrix);
      ang+=120.0;
    }

    //C4=(0 0 1)
    ang=0.0;
    for(i=0;i<3;i++)
    {
       ang+=90.0;
       matrix.rotateZAngle(ang);
       for(j=0;j<3;j++)   PGOperation.push_back(matrix*PGOperation[j]);
    }

    //C2=(1 0 0);
    matrix3x3 C2;
    C2.rotateXAngle(180.0);
    for(i=0;i<12;i++)
    {
       matrix=PGOperation[i];
       PGOperation.push_back(C2*matrix);
    }
}
//--------------------------------------------
void PointGroup::buildOperation_Oh()
{
    //Oh=O*Ci;
    buildOperation_O();

    matrix3x3 Ci,matrix;
    Ci.Set(0,0,-1.0);
    Ci.Set(1,1,-1.0);
    Ci.Set(2,2,-1.0);

    for(int i=0;i<24;i++)
    {
       matrix=PGOperation[i];
       PGOperation.push_back(Ci*matrix);
    }
}
//--------------------------------------------
void PointGroup::buildOperation_I()
{
    matrix3x3 matrix;
    unsigned int i,j;
    double ang,ang1;

    //C5=(0  0  1)
    ang=0.0;
    for(i=0;i<5;i++)
    {
      matrix.rotateZAngle(ang);
      PGOperation.push_back(matrix);
      ang+=72.0;
    }

    //C2={sin31.717474413001757,0.0,cos31.717474413001757}
    ang1=31.717474413001757*DEG_TO_RAD;
    vector3 C2(sin(ang1),0.0, cos(ang1));
    matrix.SetupRotateMatrix(C2,180.0);
    for(i=0;i<5;i++)
        PGOperation.push_back(matrix*PGOperation[i]);

    //C5=(0  0  1)
    for(i=0;i<4;i++)
    {
      ang+=72.0;
          matrix.rotateZAngle(ang);
      for(j=0;j<5;j++)
            PGOperation.push_back(matrix*PGOperation[j+5]);
    }

    //C2
    matrix.rotateYAngle(180.0);
    for(i=0;i<30;i++)
        PGOperation.push_back(matrix*PGOperation[i]);
}
//--------------------------------------------
void PointGroup::buildOperation_Ih()
{
    //!Ih=I*Ci;
    buildOperation_I();

    matrix3x3 Ci,matrix;
    Ci.Set(0,0,-1.0);
    Ci.Set(1,1,-1.0);
    Ci.Set(2,2,-1.0);
    for(unsigned int i=0;i<60;i++)
    {
       matrix=PGOperation[i];
       PGOperation.push_back(Ci*matrix);

       //one is enough
       //break;
    }
}

void PointGroup::buildInvOperation()
{
    InvPGOperation.clear();

    unsigned int size=PGOperation.size();
    unsigned int i;

    matrix3x3 matrix;

    for(i=0;i<size;i++){
       matrix=PGOperation[i];
       InvPGOperation.push_back( matrix.transpose() );
    }
}


void PointGroup::buildOperation_D00h()
{
    //E
    matrix3x3  matrix;
    matrix.Set(0,0,1.0);
    matrix.Set(1,1,1.0);
    matrix.Set(2,2,1.0);
    PGOperation.push_back(matrix);

    //xoy
    matrix.Set(0,0,1.0);
    matrix.Set(1,1,1.0);
    matrix.Set(2,2,-1.0);
    PGOperation.push_back(matrix);
}



bool PointGroup::refine()
{
    cout << "PointGroup::refine()"<<endl;

    if ( PGSymbol=="C1" || PGSymbol=="c1"){
        cout << "Cannot refine based on C1 symmetry!\n\n";
        return false;
    }

       cout << "There are " << NAtoms <<" atoms !!!"<<endl;

    //Be careful! This funcations contains "clearAll" commands
    ParseSymmetrySybmol(PGSymbol);
    cout << "There are " << NAtoms <<" atoms !!!"<<endl;

    if (NAtoms < 1)      {
        cout << "No atoms!\n\n";
        return false;
    }




    vector3 vec;
    if (isCoov || isDooh )
    {
        for (uint i = 0; i < NAtoms; ++ i) {
            vec.Set(0.0,0.0,AtomCoordinates[i][2]);
            AtomCoordinates[i]=vec;
        }

        if (isCoov) return true;
    }


    if (!buildSubGroupsOnSymmetry_refine())
    {
        cout << "\n\nCannot regroup based on symmetry of "<<PGSymbol<<endl;
        return false;
    }

    //check the SubGroups
    for (uint i=0;i < SubGroups.size();i++)
    {
        if (SubGroups.at(i).size()<PGOperation.size() ){
            cout << "\n\nCannot regroup based on symmetry of "<<PGSymbol<<endl;
            return false;
        }
    }


    buildInvOperation();

    uint operationCount= PGOperation.size();

    /*
    cout <<"========================================"<<endl;
    for(uint i=0; i < operationCount ; i++)
    {
        cout << PGOperation[i];
    }

    for(uint i=0; i < operationCount ; i++)
    {
        //cout << InvPGOperation[i]<<endl;
    }
    cout <<"========================================"<<endl;
*/
    double eps=1e-15,eps0;
    vector3 sum0,sum1;
    uint idx,times=0;


   std::cout << setprecision(16);
   cout << "Begin refine with presion of "<<eps<<"  ......"<<endl;
   times=0;
   for (uint i=0;i < SubGroups.size();i++)
   {
       cout << "\nGroup : "<<i+1<<endl;
       idx=SubGroups[i][0];
       sum0=AtomCoordinates[idx];
       times=0;

       while(1)
       {
           if (times>999)
           {
               cout << "refine 1000 times!"<<endl;
               cout <<"the error of " << eps0<<" "<< " cannot reach "<<eps<<endl;
               break;
           }

           //the first one
           sum1=VZero;

           //collect
           for(uint j=0; j < operationCount ; j++)
           {
               idx=SubGroups[i][j];
               if (j==0){
                   sum1+=AtomCoordinates[idx];
                   continue;
               }
               sum1+=InvPGOperation[j]*AtomCoordinates[idx];
           }
           sum1/=operationCount;

           //reach aim
           eps0=(sum0-sum1).length();
           if( eps0 < eps)  break;

           sum0=sum1;
           times=times+1;
           cout << "times  " << times<< " error:  "<<eps0<<endl;


           for(uint j=0; j < operationCount ; j++) //distribute
           {
               idx=SubGroups[i][j];

               //E
               if(j==0) {
                   AtomCoordinates[idx]=sum0;
                   continue;
               }

               AtomCoordinates[idx]=PGOperation[j]*sum0;
           }

       }
   }

   cout << "Refine done!"<<endl;
   cout << "\nRefined coordinates in angstrom:\n";
   std::cout.setf(ios::fixed);
   std::cout << setprecision(10);
   double x,y,z;
   for (uint k = 0; k < NAtoms; ++ k) {
       x=AtomCoordinates[k].x()*0.5291772;
       y=AtomCoordinates[k].y()*0.5291772;
       z=AtomCoordinates[k].z()*0.5291772;
       printf ("%3s  %12.8f %12.8f %12.8f\n",AtomSymbols[k].c_str(),x,y,z);
       //cout << AtomSymbols[k]<<" "<< AtomCoordinates[k]*0.5291772<<endl;
   }
   cout << endl;
   return true;
}


void PointGroup::setSymmetry(string pgsybol)
{
    PGSymbol = pgsybol;
    transform(PGSymbol.begin(),PGSymbol.end(),PGSymbol.begin(),::tolower);
    //str2Lower(PGSymbol);

    PGSymbol[0]=toupper(PGSymbol[0]);
    cout << "Set symmetry: "<< PGSymbol<<endl;
}




void PointGroup::ParseSymmetrySybmol (string schoefiles)
{
    //ClearAll();
    string first,last,middle;

    transform(schoefiles.begin(),schoefiles.end(),schoefiles.begin(),::tolower);
    //str2Lower(schoefiles);

    int length=schoefiles.length();
    first=schoefiles.at(0);
    last=schoefiles.at(length-1);

    string symbol=schoefiles;

    //cout << symbol << " inside ParseSymmetrySybmol " <<endl;

    if (schoefiles=="coov") {
        order = 1;
        isCoov=true;
        return;
    }

    else if (schoefiles=="dooh") {
        order = 1;
        //PGSymbol="Dooh";
        isDooh=true;
        return;
    }

    else if (schoefiles=="cs") {
        order = 1;
        isCs=true;
        return;
    }

    else if (schoefiles=="ci") {
        order = 1;
        isCi=true;
        return;
    }

    else if (schoefiles=="td") {
        order = 2;
        isTd=true;
        return;
    }

    else if (schoefiles=="th") {
        order = 2;
        isTh=true;
        return;
    }

    else if (schoefiles=="t") {
        order = 2;
        isT=true;
        return;
    }


    else if (schoefiles=="o") {
        order = 4;
        isO=true;
        return;
    }

    else if (schoefiles=="oh") {
        order = 2;
        isOh=true;
        return;
    }


    else if (schoefiles=="i") {
        order = 5;
        isI=true;
        return;
    }

    else if (schoefiles=="ih") {
        order = 5;
        isIh=true;
        return;
    }

    else if(first.find("c")==0 ) //Cnv  Cnh  Cn
    {
       if(last.find("v")==0 ){
           middle=symbol.substr(1,length-2);
           order=atoi(middle.c_str());
           isCnv=true;
           return;

       }
       else if(last.find("h")==0){
           middle=symbol.substr(1,length-2);
           order=atoi(middle.c_str());
           isCnh=true;
           return;
       }
       else{
           middle=symbol.substr(1,length-1);
           order=atoi(middle.c_str());
           isCn=true;
           return;
       }
    }

    else if(first.find("d")==0 ) //Dn  Dnd  Dnh
    {
        if(last.find("d")==0 ){
           middle=symbol.substr(1,length-2);
           order=atoi(middle.c_str());
           isDnd=true;
           return;
         }
         else if(last.find("h")==0){
           middle=symbol.substr(1,length-2);
           order=atoi(middle.c_str());
           isDnh=true;
           return;
         }
         else{
           middle=symbol.substr(1,length-1);
           order=atoi(middle.c_str());
           isDn=true;
           return;
         }
    }
    else if(first.find("s")==0 ){ //S2n
           middle=symbol.substr(1,length);
           order=atoi(middle.c_str())/2;
           isS2n=true;
    }
}




bool PointGroup::patchMolecule(uint & NAtomsofPatched)
{
    ParseSymmetrySybmol(PGSymbol);

    addingAtoms=0;
    //cout << PGSymbol<<endl;

    vector3 vec;

    if (NAtoms <1)  {
        cout << "No atoms!\n\n";
        return false;
    }

    cout << NAtoms <<endl;

    //                            true means for patch
    if (!buildSubGroupsOnSymmetry_patch()) return false;


   uint operationCount= PGOperation.size();

   //distance tolerance is 0.5 Bohr
   double tol0=0.5;

   uint idx=0;


   std::cout.setf(ios::fixed);
   std::cout << setprecision(12);
   cout << "Begin patch ... ..."<<endl;
   vector3 vec1,vec2,vec3;

   for (uint i=0;i < SubGroups.size();i++)
   {
       cout << "\nGroup : "<<i+1<<endl;

       //check the first one is enough
       idx=SubGroups[i][0];

       for(uint j=1; j < operationCount ; j++)
       {
            vec1=PGOperation[j]*AtomCoordinates[idx];

            bool find=false;
            for (uint k=0 ; k < SubGroups[i].size (); k++)
            {
                vec2 =  AtomCoordinates[SubGroups[i][k]];
                if ( (vec1-vec2).length() < tol0 ) {
                    find = true;
                    break;
                }
            }
            if (find) continue;

            //cannot find, store
            cout << AtomTypes[idx] <<"   " <<AtomSymbols[idx] <<" "<< vec1*0.529177<<endl;

            AtomCoordinates.push_back(vec1);
            AtomTypes.push_back(AtomTypes[idx]);
            AtomMasses.push_back(AtomMasses[idx]);
            AtomSymbols.push_back(AtomSymbols[idx]);

            uint size0=SubGroups[i].size();
            SubGroups[i].resize(size0+1);
            SubGroups[i][size0]=NAtoms+addingAtoms;
            addingAtoms=addingAtoms+1;
       }
   }

   cout << "Patch done!"<<endl;
   cout << "\nPatched coordinates in angstrom:\n";

   double x,y,z;
   for (uint k = 0; k < NAtoms+addingAtoms; ++ k) {
           x=AtomCoordinates[k].x()*0.5291772;
           y=AtomCoordinates[k].y()*0.5291772;
           z=AtomCoordinates[k].z()*0.5291772;
           printf ("%3s  %12.8f %12.8f %12.8f\n",AtomSymbols[k].c_str(),x,y,z);
           //cout << AtomSymbols[k]<<" "<< AtomCoordinates[k]*0.5291772<<endl;
   }
   cout << endl;

   NAtomsofPatched=addingAtoms;
   if (addingAtoms<1) {
       cout << "\nNo atom needs patch!\n";
       NAtomsofPatched=0;
   }
   cout << endl<<endl;

   return true;
}


bool PointGroup::buildSubGroupsOnSymmetry_patch()
{
    buildOperation ();

    uint operationCount= PGOperation.size();

    cout << "\nRegroup atoms based on Symmetry "<<PGSymbol<<":"<<endl;
    cout << "with size of "<<operationCount<<endl;


    SubGroups.clear();
    SubGroups.resize(0);

    uint idx=0;
    vector3 vec1,vec2;

    double tol0=tol;
    //label the atoms checked


    vector <uint> Labels (NAtoms,1);


    vector <uint> each;

    for (uint i=0; i < NAtoms; i++)
    {

        //check it first
        if (Labels[i]==0)   continue;

        //cout << i+1 << endl;
        Labels[i]=0;

        each.clear();
        each.push_back(i);


        for(uint j=1; j < operationCount ; j++)
        {
             vec1=PGOperation[j]*AtomCoordinates[i];

             bool find=false;

             //if refine, we need find all coordinates of each element of poting-group (k=0)
             //if patch, only the same group (k=i+1)
             uint k=i+1;
             for (; k < NAtoms; k++)
             {
                //not the same atom' type
                if (AtomTypes [i] != AtomTypes [k]) continue;

                vec2=AtomCoordinates[k];

                if ( (vec1-vec2).length() < tol0)
                {
                    find=true;
                    break;
                }
             }

              if (!find) continue;

              each.push_back(k);
              Labels[k]=0;
        }

        if (each.size() < 1) continue;


        //store
        SubGroups.resize(idx+1);
        for (uint kk = 0;kk< each.size(); kk++)
             SubGroups[idx].push_back(each[kk]);
        idx=idx+1;
        continue;



        /*//print
        for (uint kk=0; kk < Labels.size(); kk++) {
            cout <<Labels[kk]<<" ";
        }
        cout << endl;*/

    }



    //print
    cout << "Atom Type:  Atom index \n";
    for (uint i=0; i < SubGroups.size(); i++) {
        cout << "      " <<AtomTypes[SubGroups[i][0]]<<":   ";
        for (uint j=0; j < SubGroups[i].size(); j++)
            cout << "  " << SubGroups[i][j]+1;
        cout << endl;
    }

    cout << "Done\n\n";
    return true;
}





bool PointGroup::buildSubGroupsOnSymmetry_refine()
{
    buildOperation ();

    uint operationCount= PGOperation.size();

    cout << "\nRegroup atoms based on Symmetry "<<PGSymbol<<":"<<endl;
    cout << "with size of "<<operationCount<<endl;


    SubGroups.clear();
    SubGroups.resize(0);

    uint idx=0;
    vector3 vec1,vec2;

    double tol0=tol;
    //label the atoms checked


    vector <uint> Labels (NAtoms,1);
    vector <uint> each;

    bool find=false;
    uint k,i;

    for (i=0; i < NAtoms; i++)
    {
        //check it first
        if (Labels[i]==0)   continue;

        //cout << i+1 << endl;
        Labels[i]=0;

        each.clear();
        each.push_back(i);


        for(uint j=1; j < operationCount ; j++)
        {
             vec1=PGOperation[j]*AtomCoordinates[i];
Label:
             find=false;
             //if refine, we need find all coordinates of each element of poting-group (k=0)
             for ( k=0; k < NAtoms; k++)
             {
                //not the same atom' type
                if (AtomTypes [i] != AtomTypes [k]) continue;

                vec2=AtomCoordinates[k];

                if ( (vec1-vec2).length() < tol0)
                {
                    find=true;
                    break;
                }
             }

             //cannot find, maybe the tolerance is too small
             if (!find)
             {
                 if(tol0 > 0.8)  return false;
                 cout << "\nSomething is wrong when regoup atoms based on PGSymbol:\n";
                 cout << "The number of atoms inside this group is less than the order of the Symmtry Group\n";
                 cout << AtomTypes[i] <<": ";
                 for (uint l=0; l < each.size(); l++)
                     cout << each[l]+1<<" ";

                 cout << "\nEnhance the tolerance "<< tol0 << " to "<< tol0+0.1 << " to try!\n";
                 tol0=tol0+0.1;
                 goto Label;
              }

              tol0=tol; //recovery
              each.push_back(k);
              Labels[k]=0;
        }

        if (each.size() < 1) continue;

        //store
        SubGroups.resize(idx+1);
        for (uint kk = 0;kk< each.size(); kk++)
             SubGroups[idx].push_back(each[kk]);
        idx=idx+1;
        continue;

        //print
        //for (uint kk=0; kk < Labels.size(); kk++)  cout <<Labels[kk]<<" ";
        //cout << endl;
    }

    //print
    cout << "\nAtom Type:  Atom index \n";
    for (uint i=0; i < SubGroups.size(); i++) {
        cout << "      " <<AtomTypes[SubGroups[i][0]]<<":   ";
        for (uint j=0; j < SubGroups[i].size(); j++){
            cout << "  " << SubGroups[i][j]+1;
            if ((j+1)%15==0) cout << endl << "          ";
        }
        cout << endl;
    }

    cout << "Done\n\n";
    return true;

}


void  PointGroup::patchC00v()
{
    cout<< "Cannot patch based on C00v symmetry!\n";
    return;
}

void  PointGroup::patchD00h()
{
    addingAtoms=0;
    vector3 vec1, vec2;

    for (uint i=0; i < NAtoms; i++)
    {
        uint idx=i;
        vec1 =AtomCoordinates[i];

        bool find=false;
        for (uint j=i+1 ; j < NAtoms; j++)
        {
            vec2 =  AtomCoordinates[j];
            if ( (vec1-vec2).length() < tol ) {
                find = true;
                break;
            }
        }

        if (find) continue;

        //cannot find, store
        cout << AtomTypes[idx] <<"   " <<AtomSymbols[idx] <<" "<< vec1*0.529177<<endl;

        AtomCoordinates.push_back(vec1);
        AtomTypes.push_back(AtomTypes[idx]);
        AtomMasses.push_back(AtomMasses[idx]);
        AtomSymbols.push_back(AtomSymbols[idx]);

        addingAtoms=addingAtoms+1;
    }
}



