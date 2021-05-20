/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CLTool.h"
#include "CLToolRegister.h"
#include "core/PlumedMain.h"
#include "tools/Vector.h"
#include "tools/Random.h"
#include "tools/OpenMP.h"
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include "tools/Angle.h"
#include "tools/Torsion.h"

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS simplemd
/*
simplemd allows one to do molecular dynamics on systems of Lennard-Jones atoms.

The input to simplemd is specified in an input file. Configurations are input and
output in xyz format. The input file should contain one directive per line.
The directives available are as follows:

\par Examples

You run an MD simulation using simplemd with the following command:
\verbatim
plumed simplemd < in
\endverbatim

The following is an example of an input file for a simplemd calculation. This file
instructs simplemd to do 50 steps of MD at a temperature of 0.722
\verbatim
nputfile input.xyz
outputfile output.xyz
temperature 0.722
tstep 0.005
friction 1
forcecutoff 2.5
listcutoff  3.0
nstep 50
nconfig 10 trajectory.xyz
nstat   10 energies.dat
\endverbatim

If you run the following a description of all the directives that can be used in the
input file will be output.
\verbatim
plumed simplemd --help
\endverbatim

*/
//+ENDPLUMEDOC

class SimpleMD:
  public PLMD::CLTool
{
  string description()const override {
    return "run lj code";
  }

  bool write_positions_first;
  bool write_statistics_first;
  int write_statistics_last_time_reopened;
  FILE* write_statistics_fp;

  //To include bonds effects
  std::vector<double> bonddo;
  std::vector<int> Bpair; 

  //To include angle effects
  std::vector<double> thetao;
  std::vector<int> Atrip; 

  //To include torsion effects
  std::vector<double> phio;
  std::vector<int> Tquad;

public:
  static void registerKeywords( Keywords& keys ) {
    keys.add("compulsory","nstep","The number of steps of dynamics you want to run");
    keys.add("compulsory","temperature","NVE","the temperature at which you wish to run the simulation in LJ units");
    keys.add("compulsory","friction","off","The friction (in LJ units) for the Langevin thermostat that is used to keep the temperature constant");
    keys.add("compulsory","tstep","0.005","the integration timestep in LJ units");
    keys.add("compulsory","epsilon","1.0","LJ parameter");
    keys.add("compulsory","sigma","1.0","LJ parameter");
    keys.add("compulsory","inputfile","An xyz file containing the initial configuration of the system");
    keys.add("compulsory","forcecutoff","2.5","");
    keys.add("compulsory","listcutoff","3.0","");
    keys.add("compulsory","outputfile","An output xyz file containing the final configuration of the system");
    keys.add("compulsory","nconfig","10","The frequency with which to write configurations to the trajectory file followed by the name of the trajectory file");
    keys.add("compulsory","nstat","1","The frequency with which to write the statistics to the statistics file followed by the name of the statistics file");
    keys.add("compulsory","maxneighbours","10000","The maximum number of neighbors an atom can have");
    keys.add("compulsory","idum","0","The random number seed");
    keys.add("compulsory","ndim","3","The dimensionality of the system (some interesting LJ clusters are two dimensional)");
    keys.add("compulsory","wrapatoms","false","If true, atomic coordinates are written wrapped in minimal cell");

//    keys.add("option","bondsfile","Bond effects added to LJ");
  }

  explicit SimpleMD( const CLToolOptions& co ) :
    CLTool(co),
    write_positions_first(true),
    write_statistics_first(true),
    write_statistics_last_time_reopened(0),
    write_statistics_fp(NULL)
  {
    inputdata=ifile;
  }

private:

  void
  read_input(double& temperature,
             double& tstep,
             double& friction,
             double& forcecutoff,
             double& listcutoff,
             int&    nstep,
             int&    nconfig,
             int&    nstat,
             bool&   wrapatoms,
             string& inputfile,
             string& outputfile,
             string& trajfile,
             string& statfile,
             int&    maxneighbours,
             int&    ndim,
             int&    idum,
             double& epsilon,
             double& sigma)
  {

    // Read everything from input file
    char buffer1[256];
    std::string tempstr; parse("temperature",tempstr);
    if( tempstr!="NVE" ) Tools::convert(tempstr,temperature);
    parse("tstep",tstep);
    std::string frictionstr; parse("friction",frictionstr);
    if( tempstr!="NVE" ) {
      if(frictionstr=="off") { fprintf(stderr,"Specify friction for thermostat\n"); exit(1); }
      Tools::convert(frictionstr,friction);
    }
    parse("forcecutoff",forcecutoff);
    parse("listcutoff",listcutoff);
    parse("nstep",nstep);
    parse("maxneighbours",maxneighbours);
    parse("idum",idum);
    parse("epsilon",epsilon);
    parse("sigma",sigma);

    // Read in stuff with sanity checks
    parse("inputfile",inputfile);
    if(inputfile.length()==0) {
      fprintf(stderr,"Specify input file\n");
      exit(1);
    }
    parse("outputfile",outputfile);
    if(outputfile.length()==0) {
      fprintf(stderr,"Specify output file\n");
      exit(1);
    }
    std::string nconfstr; parse("nconfig",nconfstr);
    sscanf(nconfstr.c_str(),"%100d %255s",&nconfig,buffer1);
    trajfile=buffer1;
    if(trajfile.length()==0) {
      fprintf(stderr,"Specify traj file\n");
      exit(1);
    }
    std::string nstatstr; parse("nstat",nstatstr);
    sscanf(nstatstr.c_str(),"%100d %255s",&nstat,buffer1);
    statfile=buffer1;
    if(statfile.length()==0) {
      fprintf(stderr,"Specify stat file\n");
      exit(1);
    }
    parse("ndim",ndim);
    if(ndim<1 || ndim>3) {
      fprintf(stderr,"ndim should be 1,2 or 3\n");
      exit(1);
    }
    std::string w;
    parse("wrapatoms",w);
    wrapatoms=false;
    if(w.length()>0 && (w[0]=='T' || w[0]=='t')) wrapatoms=true;

//    parse("bondsfile", bondsfile);
  }

  void read_natoms(const string & inputfile,int & natoms) {
// read the number of atoms in file "input.xyz"
    FILE* fp=fopen(inputfile.c_str(),"r");
    if(!fp) {
      fprintf(stderr,"ERROR: file %s not found\n",inputfile.c_str());
      exit(1);
    }

// call fclose when fp goes out of scope
    auto deleter=[](FILE* f) { fclose(f); };
    std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);

    int ret=fscanf(fp,"%1000d",&natoms);
    if(ret==0) plumed_error() <<"Error reading number of atoms from file "<<inputfile;
  }

  void read_positions(const string& inputfile,int natoms,vector<Vector>& positions,double cell[3]) {
// read positions and cell from a file called inputfile
// natoms (input variable) and number of atoms in the file should be consistent
    FILE* fp=fopen(inputfile.c_str(),"r");
    if(!fp) {
      fprintf(stderr,"ERROR: file %s not found\n",inputfile.c_str());
      exit(1);
    }
// call fclose when fp goes out of scope
    auto deleter=[](FILE* f) { fclose(f); };
    std::unique_ptr<FILE,decltype(deleter)> fp_deleter(fp,deleter);

    char buffer[256];
    char atomname[256];
    char* cret=fgets(buffer,256,fp);
    if(cret==nullptr) plumed_error() <<"Error reading buffer from file "<<inputfile;
    int ret=fscanf(fp,"%1000lf %1000lf %1000lf",&cell[0],&cell[1],&cell[2]);
    if(ret==0) plumed_error() <<"Error reading cell line from file "<<inputfile;
    for(int i=0; i<natoms; i++) {
      ret=fscanf(fp,"%255s %1000lf %1000lf %1000lf",atomname,&positions[i][0],&positions[i][1],&positions[i][2]);
// note: atomname is read but not used
      if(ret==0) plumed_error() <<"Error reading atom line from file "<<inputfile;
    }
  }

  void randomize_velocities(const int natoms,const int ndim,const double temperature,const vector<double>&masses,vector<Vector>& velocities,Random&random) {
// randomize the velocities according to the temperature
    for(int iatom=0; iatom<natoms; iatom++) for(int i=0; i<ndim; i++)
        velocities[iatom][i]=sqrt(temperature/masses[iatom])*random.Gaussian();
  }

  void pbc(const double cell[3],const Vector & vin,Vector & vout) {
// apply periodic boundary condition to a vector
    for(int i=0; i<3; i++) {
      vout[i]=vin[i]-floor(vin[i]/cell[i]+0.5)*cell[i];
    }
  }

  void check_list(const int natoms,const vector<Vector>& positions,const vector<Vector>&positions0,const double listcutoff,
                  const double forcecutoff,bool & recompute)
  {
// check if the neighbour list have to be recomputed
    recompute=false;
    auto delta2=(0.5*(listcutoff-forcecutoff))*(0.5*(listcutoff-forcecutoff));
// if ANY atom moved more than half of the skin thickness, recompute is set to .true.
    for(int iatom=0; iatom<natoms; iatom++) {
      if(modulo2(positions[iatom]-positions0[iatom])>delta2) recompute=true;
    }
  }


  void compute_list(const int natoms,const vector<Vector>& positions,const double cell[3],const double listcutoff,
                    vector<vector<int> >& list) {
    double listcutoff2=listcutoff*listcutoff; // squared list cutoff
    list.assign(natoms,vector<int>());
#   pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(static,1)
    for(int iatom=0; iatom<natoms-1; iatom++) {
      for(int jatom=iatom+1; jatom<natoms; jatom++) {
        auto distance=positions[iatom]-positions[jatom];
        Vector distance_pbc; // minimum-image distance of the two atoms
        pbc(cell,distance,distance_pbc);
// if the interparticle distance is larger than the cutoff, skip
        if(modulo2(distance_pbc)>listcutoff2)continue;
        list[iatom].push_back(jatom);
      }
    }
  }

  void compute_forces(const int natoms,double epsilon, double sigma,const vector<Vector>& positions,const double cell[3],
                      double forcecutoff,const vector<vector<int> >& list,vector<Vector>& forces,double & engconf)
  {
    double forcecutoff2=forcecutoff*forcecutoff; // squared force cutoff
    engconf=0.0;
    for(int i=0; i<natoms; i++)for(int k=0; k<3; k++) forces[i][k]=0.0;
    double engcorrection=4.0*epsilon*(1.0/pow(forcecutoff2/(sigma*sigma),6.0)-1.0/pow(forcecutoff2/(sigma*sigma),3)); // energy necessary shift the potential avoiding discontinuities
#   pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      std::vector<Vector> omp_forces(forces.size());
      #pragma omp for reduction(+:engconf) schedule(static,1) nowait
      for(int iatom=0; iatom<natoms-1; iatom++) {
        for(int jlist=0;jlist<list[iatom].size();jlist++){
          const int jatom=list[iatom][jlist];
          auto distance=positions[iatom]-positions[jatom];
          Vector distance_pbc;    // minimum-image distance of the two atoms
          pbc(cell,distance,distance_pbc);
          auto distance_pbc2=modulo2(distance_pbc);   // squared minimum-image distance
// if the interparticle distance is larger than the cutoff, skip
          if(distance_pbc2>forcecutoff2) continue;
          auto distance_pbcm2=sigma*sigma/distance_pbc2;
          auto distance_pbcm6=distance_pbcm2*distance_pbcm2*distance_pbcm2;
          auto distance_pbcm8=distance_pbcm6*distance_pbcm2;
          auto distance_pbcm12=distance_pbcm6*distance_pbcm6;
          auto distance_pbcm14=distance_pbcm12*distance_pbcm2;
          engconf+=4.0*epsilon*(distance_pbcm12 - distance_pbcm6) - engcorrection;
          auto f=24.0*distance_pbc*(2.0*distance_pbcm14-distance_pbcm8)*epsilon/sigma;
          omp_forces[iatom]+=f;
          omp_forces[jatom]-=f;
        }
      }
//Adding bonds effects
#     pragma omp for reduction(+: engconf) schedule(static, 1) nowait
      for (int i=0; i < bonddo.size(); ++i) {
      	   int index1 = 2*i;
	   int index2 = index1 + 1;
	   int atom1  = Bpair[index1] - 1;
	   int atom2  = Bpair[index2] - 1;

	   auto r21 = positions[atom2] - positions[atom1];
	   Vector r21_pbc; //OJO PLUMED Vector NO std::vector
	   pbc(cell, r21, r21_pbc);

	   const double r = r21_pbc.modulo();
	   const double rdo = r - bonddo[i];

	   const double KAPPAB = 1000.0; //KAPPA bonds effect

	   engconf += 0.5 * KAPPAB * rdo * rdo;
	   auto f   = KAPPAB * (rdo/r) * r21_pbc;

	   omp_forces[atom1] += f;
	   omp_forces[atom2] -= f;	   
      }

//Adding angles effects
#     pragma omp for reduction(+: engconf) schedule(static, 1) nowait
      for (int i=0; i < thetao.size(); ++i) {
      	   int index1 = 3*i;
	   int index2 = index1 + 1;
	   int index3 = index1 + 2;
	   int atom1a  = Atrip[index1] - 1; //atom1a = atom1 of angle effect
	   int atom2a  = Atrip[index2] - 1; //atom2a = atom2 of angle effect
	   int atom3a  = Atrip[index3] - 1; //atom3a = atom3 of angle effect

	   auto r21a = positions[atom2a] - positions[atom1a];
	   auto r23a = positions[atom2a] - positions[atom3a];

	   Vector r21a_pbc, r23a_pbc; //OJO PLUMED Vector NO std::vector
	   pbc(cell, r21a, r21a_pbc);
	   pbc(cell, r23a, r23a_pbc);

	   Angle theta_plmd;
	   double theta = theta_plmd.compute(r21a, r23a, r21a_pbc, r23a_pbc);
	   const double deltatheta = theta - thetao[i];

	   const double KAPPAA = 1000.0; //KAPPA angle effect

	   engconf += 0.5 * KAPPAA * deltatheta * deltatheta;
	   auto f21   = KAPPAA * deltatheta * r21a_pbc;
	   auto f23   = KAPPAA * deltatheta * r23a_pbc;

	   omp_forces[atom1a] += f21;
	   omp_forces[atom2a] -= f21 + f23;
	   omp_forces[atom3a] += f23;	   
      }
        
//Adding torsions effects
#     pragma omp for reduction(+: engconf) schedule(static, 1) nowait
      for (int i=0; i < phio.size(); ++i) {
      	   int index1t = 4*i;
	   int index2t = index1t + 1;
	   int index3t = index1t + 2;
	   int index4t = index1t + 3;
	   int atom1t  = Tquad[index1t] - 1; //atom1t = atom1 of torsion effect
	   int atom2t  = Tquad[index2t] - 1; //atom2t = atom2 of torsion effect
	   int atom3t  = Tquad[index3t] - 1; //atom3t = atom3 of torsion effect
	   int atom4t  = Tquad[index4t] - 1; //atom4t = atom4 of torsion effect

	   auto r12t = positions[atom1t] - positions[atom2t];
	   auto r23t = positions[atom2t] - positions[atom3t];
	   auto r34t = positions[atom3t] - positions[atom4t];

	   Vector r12t_pbc, r23t_pbc, r34t_pbc; //OJO PLUMED Vector NO std::vector
	   pbc(cell, r12t, r12t_pbc);
	   pbc(cell, r23t, r23t_pbc);
	   pbc(cell, r34t, r34t_pbc);

	   Torsion phi_plmd;
	   double phi = phi_plmd.compute(r12t, r23t, r34t, r12t_pbc, r23t_pbc, r34t_pbc);
	   const double deltaphi = phi - phio[i];

	   const double KAPPAT = 100.0; //KAPPA torsion effect

	   engconf += 0.5 * KAPPAT * ( 1.0 + cos(deltaphi) );
	   auto f12   =  0.5 * KAPPAT * sin(deltaphi) * r12t_pbc;
	   auto f23   =  0.5 * KAPPAT * sin(deltaphi) * r23t_pbc;
	   auto f34   =  0.5 * KAPPAT * sin(deltaphi) * r34t_pbc;

	   omp_forces[atom1t] -= f12;
	   omp_forces[atom2t] +=  f12 - f23;
	   omp_forces[atom3t] +=  f23 - f34;
	   omp_forces[atom4t] += f34;
      }

#     pragma omp critical
      for(unsigned i=0; i<omp_forces.size(); i++) forces[i]+=omp_forces[i];
    }
    }

//  }

  void compute_engkin(const int natoms,const vector<double>& masses,const vector<Vector>& velocities,double & engkin)
  {
// calculate the kinetic energy from the velocities
    engkin=0.0;
    for(int iatom=0; iatom<natoms; iatom++) {
        engkin+=0.5*masses[iatom]*modulo2(velocities[iatom]);
      }
  }


  void thermostat(const int natoms,const int ndim,const vector<double>& masses,const double dt,const double friction,
                  const double temperature,vector<Vector>& velocities,double & engint,Random & random) {
// Langevin thermostat, implemented as described in Bussi and Parrinello, Phys. Rev. E (2007)
// it is a linear combination of old velocities and new, randomly chosen, velocity,
// with proper coefficients
    double c1=exp(-friction*dt);
    for(int iatom=0; iatom<natoms; iatom++) {
      double c2=sqrt((1.0-c1*c1)*temperature/masses[iatom]);
      for(int i=0; i<ndim; i++) {
        engint+=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
        velocities[iatom][i]=c1*velocities[iatom][i]+c2*random.Gaussian();
        engint-=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
      }
    }
  }

  void write_positions(const string& trajfile,int natoms,const vector<Vector>& positions,const double cell[3],const bool wrapatoms)
  {
// write positions on file trajfile
// positions are appended at the end of the file
    Vector pos;
    FILE*fp;
    if(write_positions_first) {
      fp=fopen(trajfile.c_str(),"w");
      write_positions_first=false;
    } else {
      fp=fopen(trajfile.c_str(),"a");
    }
    fprintf(fp,"%d\n",natoms);
    fprintf(fp,"%f %f %f\n",cell[0],cell[1],cell[2]);
    for(int iatom=0; iatom<natoms; iatom++) {
// usually, it is better not to apply pbc here, so that diffusion
// is more easily calculated from a trajectory file:
      if(wrapatoms) pbc(cell,positions[iatom],pos);
      else pos=positions[iatom];
      fprintf(fp,"Ar %10.7f %10.7f %10.7f\n",pos[0],pos[1],pos[2]);
    }
    fclose(fp);
  }

  void write_final_positions(const string& outputfile,int natoms,const vector<Vector>& positions,const double cell[3],const bool wrapatoms)
  {
// write positions on file outputfile
    Vector pos;
    FILE*fp;
    fp=fopen(outputfile.c_str(),"w");
    fprintf(fp,"%d\n",natoms);
    fprintf(fp,"%f %f %f\n",cell[0],cell[1],cell[2]);
    for(int iatom=0; iatom<natoms; iatom++) {
// usually, it is better not to apply pbc here, so that diffusion
// is more easily calculated from a trajectory file:
      if(wrapatoms) pbc(cell,positions[iatom],pos);
      else pos=positions[iatom];
      fprintf(fp,"Ar %10.7f %10.7f %10.7f\n",pos[0],pos[1],pos[2]);
    }
    fclose(fp);
  }


  void write_statistics(const string & statfile,const int istep,const double tstep,
                        const int natoms,const int ndim,const double engkin,const double engconf,const double engint) {
// write statistics on file statfile
    if(write_statistics_first) {
// first time this routine is called, open the file
      write_statistics_fp=fopen(statfile.c_str(),"w");
      write_statistics_first=false;
    }
    if(istep-write_statistics_last_time_reopened>100) {
// every 100 steps, reopen the file to flush the buffer
      fclose(write_statistics_fp);
      write_statistics_fp=fopen(statfile.c_str(),"a");
      write_statistics_last_time_reopened=istep;
    }
    fprintf(write_statistics_fp,"%d %f %f %f %f %f\n",istep,istep*tstep,2.0*engkin/double(ndim*natoms),engconf,engkin+engconf,engkin+engconf+engint);
  }



  int main(FILE* in,FILE*out,PLMD::Communicator& pc) override {
    int            natoms;       // number of atoms
    vector<Vector> positions;    // atomic positions
    vector<Vector> velocities;   // velocities
    vector<double> masses;       // masses
    vector<Vector> forces;       // forces
    double         cell[3];      // cell size
    double         cell9[3][3];  // cell size

// neighbour list variables
// see Allen and Tildesey book for details
    vector< vector<int> >  list; // neighbour list
    vector<Vector> positions0;   // reference atomic positions, i.e. positions when the neighbour list

// input parameters
// all of them have a reasonable default value, set in read_input()
    double      tstep;             // simulation timestep
    double      temperature;       // temperature
    double      friction;          // friction for Langevin dynamics (for NVE, use 0)
    double      listcutoff;        // cutoff for neighbour list
    double      forcecutoff;       // cutoff for forces
    int         nstep;             // number of steps
    int         nconfig;           // stride for output of configurations
    int         nstat;             // stride for output of statistics
    int         maxneighbour;      // maximum average number of neighbours per atom
    int         ndim;              // dimensionality of the system (1, 2, or 3)
    int         idum;              // seed
    int         plumedWantsToStop; // stop flag
    bool        wrapatoms;         // if true, atomic coordinates are written wrapped in minimal cell
    string      inputfile;         // name of file with starting configuration (xyz)
    string      outputfile;        // name of file with final configuration (xyz)
    string      trajfile;          // name of the trajectory file (xyz)
    string      statfile;          // name of the file with statistics

    double      epsilon, sigma;    // LJ parameters

    double engkin;                 // kinetic energy
    double engconf;                // configurational energy
    double engint;                 // integral for conserved energy in Langevin dynamics

    bool recompute_list;           // control if the neighbour list have to be recomputed

    Random random;                 // random numbers stream

    std::unique_ptr<PlumedMain> plumed;


// Commenting the next line it is possible to switch-off plumed
    plumed.reset(new PLMD::PlumedMain);

    if(plumed) {
      int s=sizeof(double);
      plumed->cmd("setRealPrecision",&s);
    }

    read_input(temperature,tstep,friction,forcecutoff,
               listcutoff,nstep,nconfig,nstat,
               wrapatoms,inputfile,outputfile,trajfile,statfile,
               maxneighbour,ndim,idum,epsilon,sigma);

// number of atoms is read from file inputfile
    read_natoms(inputfile,natoms);

// write the parameters in output so they can be checked
    fprintf(out,"%s %s\n","Starting configuration           :",inputfile.c_str());
    fprintf(out,"%s %s\n","Final configuration              :",outputfile.c_str());
    fprintf(out,"%s %d\n","Number of atoms                  :",natoms);
    fprintf(out,"%s %f\n","Temperature                      :",temperature);
    fprintf(out,"%s %f\n","Time step                        :",tstep);
    fprintf(out,"%s %f\n","Friction                         :",friction);
    fprintf(out,"%s %f\n","Cutoff for forces                :",forcecutoff);
    fprintf(out,"%s %f\n","Cutoff for neighbour list        :",listcutoff);
    fprintf(out,"%s %d\n","Number of steps                  :",nstep);
    fprintf(out,"%s %d\n","Stride for trajectory            :",nconfig);
    fprintf(out,"%s %s\n","Trajectory file                  :",trajfile.c_str());
    fprintf(out,"%s %d\n","Stride for statistics            :",nstat);
    fprintf(out,"%s %s\n","Statistics file                  :",statfile.c_str());
    fprintf(out,"%s %d\n","Max average number of neighbours :",maxneighbour);
    fprintf(out,"%s %d\n","Dimensionality                   :",ndim);
    fprintf(out,"%s %d\n","Seed                             :",idum);
    fprintf(out,"%s %s\n","Are atoms wrapped on output?     :",(wrapatoms?"T":"F"));
    fprintf(out,"%s %f\n","Epsilon                          :",epsilon);
    fprintf(out,"%s %f\n","Sigma                            :",sigma);

// Setting the seed
    random.setSeed(idum);

// allocation of dynamical arrays
    positions.resize(natoms);
    positions0.resize(natoms);
    velocities.resize(natoms);
    forces.resize(natoms);
    masses.resize(natoms);
    list.resize(natoms);

//Reading the files Bpair.dat and bonddo.dat file and store the parameters in order to include bonds effects

    //Bpair.dat
    ifstream inbp("preprocdata/Bpair.dat");
    string line;

    if (!inbp) {
		cerr << "Cannot open the file: Bpair.dat" << endl;
		return false;
    }

    while (getline(inbp, line)) {
	if (line.size() > 0)
		Bpair.push_back(std::stoi(line));
    }
    inbp.close();

    //bonddo.dat
    ifstream inbd("preprocdata/bonddo.dat");

    if (!inbd) {
		cerr << "Cannot open the file: bonddo.dat" << endl;
		return false;
    }

    while (getline(inbd, line)) {
	if (line.size() > 0)
		bonddo.push_back(std::stod(line));
    }
    inbd.close();

//Reading the files Atrip.dat and thetao.dat file and store the parameters in order to include angles effects

    //Atrip.dat
    ifstream inat("preprocdata/Atrip.dat");
    string line2;

    if (!inat) {
		cerr << "Cannot open the file: Atrip.dat" << endl;
		return false;
    }

    while (getline(inat, line2)) {
	if (line2.size() > 0)
		Atrip.push_back(std::stoi(line2));
    }
    inat.close();

    //thetao.dat
    ifstream inatheta("preprocdata/thetao.dat");

    if (!inatheta) {
		cerr << "Cannot open the file: theta.dat" << endl;
		return false;
    }

    while (getline(inatheta, line2)) {
	if (line2.size() > 0)
		thetao.push_back(std::stod(line2));
    }
    inatheta.close();

//Reading the files Tquad.dat and phio.dat file and store the parameters in order to include torsions effects

    //Tquad.dat
    ifstream intq("preprocdata/Tquad.dat");
    string line3;

    if (!intq) {
		cerr << "Cannot open the file: Tquad.dat" << endl;
		return false;
    }

    while (getline(intq, line3)) {
	if (line3.size() > 0)
		Tquad.push_back(std::stoi(line3));
    }
    intq.close();

    //phio.dat
    ifstream intphi("preprocdata/phio.dat");

    if (!intphi) {
		cerr << "Cannot open the file: phio.dat" << endl;
		return false;
    }

    while (getline(intphi, line3)) {
	if (line3.size() > 0)
		phio.push_back(std::stod(line3));
    }
    intphi.close();

/*    
    for (int i=0; i < bonddo.size(); ++i){
	cout << " " << bonddo[i] << endl;   //Uncoment this if you want to be sure you are reading the file correctly
    }
    
    for (int i=0; i < Bpair.size(); ++i){
	cout << " " << Bpair[i] << endl;   //Uncoment this if you want to be sure you are reading the file correctly
    }

    for (int i=0; i < Atrip.size(); ++i){
	cout << " " << Atrip[i] << endl;   //Uncoment this if you want to be sure you are reading the file correctly
    }

    for (int i=0; i < thetao.size(); ++i){
	cout << " " << thetao[i] << endl;   //Uncoment this if you want to be sure you are reading the file correctly
    }

    for (int i=0; i < Tquad.size(); ++i){
	cout << " " << Tquad[i] << endl;   //Uncoment this if you want to be sure you are reading the file correctly
    }

    for (int i=0; i < phio.size(); ++i){
	cout << " " << phio[i] << endl;   //Uncoment this if you want to be sure you are reading the file correctly
    }

    string line;
    int k=1;
    ifstream file("bonds.dat");
    while (getline (file, line) ) {
	    if (k%2 != 0) {
		    size_t pos = line.find("=");
		    string ATOMS=line.substr(pos+1);
		    stringstream ss(ATOMS);
		    while (ss.good()) {
			    string substr;
			    getline(ss, substr, ',');
			    Bpair.push_back(std::stoi(substr));
		    }
	    } else {
		    size_t pos = line.find("AT=");
		    size_t p=15;
		    string at=line.substr(pos+3,p);
		    do_b.push_back(std::stod(at,&p));

	    }
	    ++k;
    }
    file.close();

    cout << do_b.size() << " " << endl;
*/


// masses are hard-coded to 1
    for(int i=0; i<natoms; ++i) masses[i]=1.0;

// energy integral initialized to 0
    engint=0.0;

// positions are read from file inputfile
    read_positions(inputfile,natoms,positions,cell);

// velocities are randomized according to temperature
    randomize_velocities(natoms,ndim,temperature,masses,velocities,random);

    if(plumed) {
      plumed->cmd("setNoVirial");
      plumed->cmd("setNatoms",&natoms);
      plumed->cmd("setMDEngine","simpleMD");
      plumed->cmd("setTimestep",&tstep);
      plumed->cmd("setPlumedDat","plumed.dat");
      int pversion=0;
      plumed->cmd("getApiVersion",&pversion);
// setting kbT is only implemented with api>1
// even if not necessary in principle in SimpleMD (which is part of plumed)
// we leave the check here as a reference
      if(pversion>1) {
        plumed->cmd("setKbT",&temperature);
      }
      plumed->cmd("init");
    }

// neighbour list are computed, and reference positions are saved
    compute_list(natoms,positions,cell,listcutoff,list);

    int list_size=0;
    for(int i=0;i<list.size();i++) list_size+=list[i].size();
    fprintf(out,"List size: %d\n",list_size);
    for(int iatom=0; iatom<natoms; ++iatom) positions0[iatom]=positions[iatom];

// forces are computed before starting md
    compute_forces(natoms,epsilon,sigma,positions,cell,forcecutoff,list,forces,engconf);

// remove forces if ndim<3
    if(ndim<3)
      for(int iatom=0; iatom<natoms; ++iatom) for(int k=ndim; k<3; ++k) forces[iatom][k]=0.0;

// here is the main md loop
// Langevin thermostat is applied before and after a velocity-Verlet integrator
// the overall structure is:
//   thermostat
//   update velocities
//   update positions
//   (eventually recompute neighbour list)
//   compute forces
//   update velocities
//   thermostat
//   (eventually dump output informations)
    for(int istep=0; istep<nstep; istep++) {
      thermostat(natoms,ndim,masses,0.5*tstep,friction,temperature,velocities,engint,random);

      for(int iatom=0; iatom<natoms; iatom++)
          velocities[iatom]+=forces[iatom]*0.5*tstep/masses[iatom];

      for(int iatom=0; iatom<natoms; iatom++)
          positions[iatom]+=velocities[iatom]*tstep;

// a check is performed to decide whether to recalculate the neighbour list
      check_list(natoms,positions,positions0,listcutoff,forcecutoff,recompute_list);
      if(recompute_list) {
        compute_list(natoms,positions,cell,listcutoff,list);
        for(int iatom=0; iatom<natoms; ++iatom) positions0[iatom]=positions[iatom];
        int list_size=0;
        for(int i=0;i<list.size();i++) list_size+=list[i].size();
        fprintf(out,"List size: %d\n",list_size);
      }

      compute_forces(natoms,epsilon,sigma,positions,cell,forcecutoff,list,forces,engconf);

      if(plumed) {
        int istepplusone=istep+1;
        plumedWantsToStop=0;
        for(int i=0; i<3; i++)for(int k=0; k<3; k++) cell9[i][k]=0.0;
        for(int i=0; i<3; i++) cell9[i][i]=cell[i];
        plumed->cmd("setStep",&istepplusone);
        plumed->cmd("setMasses",&masses[0]);
        plumed->cmd("setForces",&forces[0]);
        plumed->cmd("setEnergy",&engconf);
        plumed->cmd("setPositions",&positions[0]);
        plumed->cmd("setBox",cell9);
        plumed->cmd("setStopFlag",&plumedWantsToStop);
        plumed->cmd("calc");
        if(plumedWantsToStop) nstep=istep;
      }
// remove forces if ndim<3
      if(ndim<3)
        for(int iatom=0; iatom<natoms; ++iatom) for(int k=ndim; k<3; ++k) forces[iatom][k]=0.0;

      for(int iatom=0; iatom<natoms; iatom++)
          velocities[iatom]+=forces[iatom]*0.5*tstep/masses[iatom];

      thermostat(natoms,ndim,masses,0.5*tstep,friction,temperature,velocities,engint,random);

// kinetic energy is calculated
      compute_engkin(natoms,masses,velocities,engkin);

// eventually, write positions and statistics
      if((istep+1)%nconfig==0) write_positions(trajfile,natoms,positions,cell,wrapatoms);
      if((istep+1)%nstat==0)   write_statistics(statfile,istep+1,tstep,natoms,ndim,engkin,engconf,engint);

    }

// call final plumed jobs
    plumed->cmd("runFinalJobs");

// write final positions
    write_final_positions(outputfile,natoms,positions,cell,wrapatoms);

// close the statistic file if it was open:
    if(write_statistics_fp) fclose(write_statistics_fp);

    return 0;
  }


};

PLUMED_REGISTER_CLTOOL(SimpleMD,"simplemd")

}
}




