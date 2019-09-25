/*
 * Copyright (c) 2018, Aleksi Kurkela, Aleksas Mazeliauskas, Jean-Francois
 * Paquet, Soeren Schlichting and Derek Teaney
 * All rights reserved.
 *
 * KoMPoST is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/KMPST/KoMPoST/
 */
#include <cstring>
#include <omp.h>

//input file parser
#include "INIReader.h"

// ENERGY-MOMENTUM TENSOR DATA HOLDER //
#include "EnergyMomentumTensor.h"

// Code-for evolution
#include "KineticEvolution.h"

#include "GreensFunctions.h"

// ENERGY-MOMENTUM TENSOR IO //
#include "EnergyMomentumTensorIO_music.inc"

int main(int argc, char **argv) {
  // Try to open input file
  if (argc < 2) {
    std::cerr << "Error: ** " << argv[0] << " ** no input filename " << std::endl;
    std::cerr << "USAGE:" << std::endl;
    std::cerr << "    " << argv[0] << " setup.ini" << std::endl;
    exit(EXIT_FAILURE);
  }

  INIReader reader(argv[1]) ;
  if (reader.ParseError()){
    std::cerr << "Error: ** " << argv[0] << " ** failed to open " << argv[1] << std::endl;
    exit(EXIT_FAILURE);
  }

  //  Get KoMPoST input/output parameters from input file
  double default_tIn = 0.2005;
  double default_tOut = 1.2005;
  std::string default_InputFile="examples/ipglasma/EnergyMomentumTensorT0.2005fmID1483491970.txt_sigma02_big_smeared_coarse";
  std::string default_OutputFileTag="EnergyMomentumTensorT1.2005fmID1483491970" ;

  double tIn = reader.GetReal("KoMPoSTInputs","tIn",default_tIn) ;
  double tOut = reader.GetReal("KoMPoSTInputs","tOut",default_tOut);
  std::string InputFile = reader.GetString("KoMPoSTInputs", "InputFile", default_InputFile) ;
  std::string OutputFileTag =reader.GetString("KoMPoSTInputs", "OutputFileTag", default_OutputFileTag) ;

  //  Get grid parameters from input file
  EventInput::Setup(reader);
  //  Get KoMPoST parameters from input file
  // Smearing width of response functions
  // Hard-coded

  // SETUP OpenMP
  int NumberOfOpenMPThreads = omp_get_max_threads();
  std::cerr << "#CALCULATING WITH " << NumberOfOpenMPThreads
            << " OPEN-MP THREADS" << std::endl;

  // ALLOCATE INITIAL ENERGY-MOMENTUM TENSOR //
  EnergyMomentumTensorMap *Tmunu_In = new EnergyMomentumTensorMap(tIn);
  
  // LOAD INITIAL ENERGY MOMENTUM TENSOR //
  EnergyMomentumTensorMapLoad(Tmunu_In, InputFile);

  // OUTPUT INITIAL ENERGY DENSITY AND FLOW PROFILE //
  //std::string tmp_file_name;
  bool use_pimunu_NS=false;
  std::stringstream tmp_file_name;
  tmp_file_name << OutputFileTag << ".tOut" << tIn << ".music_init_flowNonLinear_pimunuTransverse.txt";
  write_initial_conditions_MUSIC(tmp_file_name.str(), use_pimunu_NS, EventInput::afm, EventInput::afm, EventInput::Ns, EventInput::Ns, Tmunu_In, Tmunu_In);

  const double delta_tau_aim=0.1;
  const double tau_steps=int((tOut-tIn)/delta_tau_aim)+1;
  const double delta_tau=(tOut-tIn)/tau_steps;
  std::cout << "Outputting the energy density & flow every " << delta_tau << " fm\n";
  // LOOP OVER TIME //
  for(int itau=1;itau<=tau_steps;itau++) {

          const double tau=tIn+itau*delta_tau;

          KoMPoSTParameters::Sigma=0.1/(tau-tIn);
          KoMPoSTParameters::Setup(reader);

          // ALLOCATE FINAL AND BACKGROUND ENERGY-MOMENTUM TENSOR //
          EnergyMomentumTensorMap *Tmunu_OutFull = new EnergyMomentumTensorMap(tau);
          EnergyMomentumTensorMap *Tmunu_OutBG   = new EnergyMomentumTensorMap(tau);

          // LOAD RESPONSE FUNCTIONS //
          KoMPoST::Setup() ;
          // COMPUTE EVOLVE ENERGY-MOMENTUM TENSOR //
          KoMPoST::Run(Tmunu_In, Tmunu_OutBG, Tmunu_OutFull) ;
          // END OF KOMPOST EVOLUTION //

          // WRITE OUT EVOLVED ENERGY-MOMENTUM TENSOR //
          tmp_file_name.str("");
          //std::iostream 
          tmp_file_name << OutputFileTag << ".tOut" << tau << ".music_init_flowNonLinear_pimunuTransverse.txt";
          void write_initial_conditions_MUSIC(std::string outfile_name, bool use_sigmamunu_NavierStokes, double dx, double dy, int Nx, int Ny, EnergyMomentumTensorMap *Tmunu_Out_Full, EnergyMomentumTensorMap *Tmunu_Out_BG);
          //tmp_file_name = std::string(OutputFileTag) + "tau" + str(tau) + ".music_init_flowNonLinear_pimunuTransverse.txt";
          //use_pimunu_NS=false;
          write_initial_conditions_MUSIC(tmp_file_name.str(), use_pimunu_NS, EventInput::afm, EventInput::afm, EventInput::Ns, EventInput::Ns, Tmunu_OutFull, Tmunu_OutBG);
          //tmp = std::string(OutputFileTag) + ".music_init_flowNonLinear_pimunuTransverse_pimunuNS.txt";
          //use_pimunu_NS=true;
          //write_initial_conditions_MUSIC(tmp, use_pimunu_NS, EventInput::afm, EventInput::afm, EventInput::Ns, EventInput::Ns, Tmunu_OutFull, Tmunu_OutBG);


          //tmp = std::string(OutputFileTag) + ".input.txt";
          //EnergyMomentumTensorMapSave(Tmunu_In, tmp);
          //tmp = std::string(OutputFileTag) + ".txt";
          //EnergyMomentumTensorMapSave(Tmunu_OutFull, tmp);
          //tmp = std::string(OutputFileTag) + ".background.txt";
          //EnergyMomentumTensorMapSave(Tmunu_OutBG, tmp);
          }

}
