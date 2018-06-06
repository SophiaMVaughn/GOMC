#include "InputFileReader.h"
#include <string>
#include <stdlib.h>
#include <iostream>
#include <args.hxx>

std::string EnsembleSearch(std::string filename);
void CPUversion(std::string ensType, std::string fileName) //GOMC::CPUversion? ***************************************
void GPUfunct(std::string ensType, std::string fileName) //GOMC::GPUfunct?
//Do we define the parsing flags here or in the main? Is it define the same way as the test case?

int main(int argc, char* argv[])
{
  // Check if the user provided input file as an argument
  if(argc < 2) {
    std::cerr << "Error: Input file name is required!" << std::endl;
    std::cerr << "Example: GOMC in.conf" << std::endl;
    exit(0);
  }

  // Read the input file and find the Ensemble keyword
  // and assign the ensembleType to that value
  std::string filename(argv[1]);
	std::string ensembleType = EnsembleSearch(filename);

  // If the returned value was NAN, it means it couldn't find the ensemble keyword
  // So we will exit here
  if(ensembleType == "NAN") {
    std::cerr << "Error: Ensemble type is required in your input file!" << std::endl;
    std::cerr << "Example: Ensemble GCMC" << std::endl;
    exit(0);
  }
	
//set the parser here with the input line **************************************************************************
//if --gpu flag is true, run GPUfunct
	GPUfunct(ensembleType, filename);
//if --gpu flag is false, run CPUversion
	CPUversion(ensembleType, filename);
//return 0
	
	return 0;
}



// This function will return the ensemble type by reading the input file
// @filename: Input file name that passed by user
// @return value: The ensemble
// NOTE: It will return "NAN" if it couldn't find the ensemble
std::string EnsembleSearch(std::string filename)
{
  // line will contain each line, broken into tokens
  // So line[0] will hold the first word (token/keyword) and 
  // line[1] will be the second word in the line
  std::vector<std::string> line;

  // The InputFileReader is designed to automatically read each line and return
  // the tokens to you by using readNextLine function
  // reader.Open will open the file for you
  InputFileReader reader;
  reader.Open(filename);

  // We will loop here and ignore the other lines until we find "Ensemble" keyword
	while (reader.readNextLine(line))
	{
		if (line.size() == 0)
			continue;
		if (line[0] == "Ensemble" || line[0] == "ENSEMBLE" || line[0] == "ensemble")
		{
			return line[1]; //return the word after ensemble
		}
	}

  // return "NAN" if we come out of the loop. AKA couldn't find the keyword
  return "NAN";
}


//this will run the CPU compiled version of the ensemble specified
void CPUversion(std::string ensType, std::string fileName) //GOMC::CPUversion? ******************************
{
#ifdef _WIN32
  // Generate the command string based on the ensemble
  // It should look somewhat close to :
  // GOMC_CPU_GCMC.exe in.conf
  std::string Executable_To_Run = "GOMC_CPU_";
  Executable_To_Run += ensembleType;
  Exectuable_To_Run += ".exe";
  Executable_To_Run += " ";
  Executable_To_Run += filename;

  // Call the system function to actually run the simulation
  system(Executable_To_Run.c_str());
#endif

#if defined(__linux__) || defined(__APPLE__)
  // Generate the command string based on the ensemble
  // It should look somewhat close to :
  // GOMC_CPU_GCMC in.conf
  std::string Executable_To_Run = "GOMC_CPU_";
  Executable_To_Run += ensembleType;
  Executable_To_Run += " ";
  Executable_To_Run += filename;

  // Call the system function to actually run the simulation
  system(Executable_To_Run.c_str());
#endif
}


//this will run the GPU specified version of the executable (with # of threads)
void GPUfunct(std::string ensType, std::string fileName) //GOMC::GPUfunct? ***********************************
{
	
//if --thread true, record the # in std::string threadInt (or cast the int to string so it can be concatenated) *******************
//std::string pNum = "+p";
//pNum += threadInt;
	
#ifdef _WIN32
  // Generate the command string based on the ensemble
  // It should look somewhat close to :
  // GOMC_GPU_GCMC.exe +p# in.conf
  std::string Executable_To_Run = "GOMC_GPU_";
  Executable_To_Run += ensType;
  Exectuable_To_Run += ".exe";
  Executable_To_Run += " ";
	
  //Executable_To_Run += pNum; ****************************************
  //Executable_To_Run += " ";
	
  Executable_To_Run += fileName;

  // Call the system function to actually run the simulation
  system(Executable_To_Run.c_str());
#endif

#if defined(__linux__) || defined(__APPLE__)
  // Generate the command string based on the ensemble
  // It should look somewhat close to :
  // GOMC_GPU_GCMC +p# in.conf
  std::string Executable_To_Run = "GOMC_GPU_";
  Executable_To_Run += ensembleType;
  Executable_To_Run += " ";
	
  //Executable_To_Run += pNum; **************************************
  //Executable_To_Run += " ";
	
  Executable_To_Run += filename;

  // Call the system function to actually run the simulation
  system(Executable_To_Run.c_str());
#endif
}
