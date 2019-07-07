//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: FSORunAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file FSORunAction.cc
/// \brief Implementation of the FSORunAction class

#include "FSORunAction.hh"
#include "FSOPrimaryGeneratorAction.hh"
#include "FSODetectorConstruction.hh"
// #include "FSORun.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4Run.hh"

#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"

#include <fstream>
#include <iostream>
#include <string>
#ifdef _WIN32
#include <Windows.h>
#include <ShObjIdl.h>
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


FSORunAction::FSORunAction()
: G4UserRunAction()
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 
  accumulablesRegistered = false;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FSORunAction::~FSORunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void FSORunAction::BeginOfRunAction(const G4Run*)
{ 

  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  if(!accumulablesRegistered)
  {
	  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
	  SDCount = SDMan->GetHCtable()->entries();

	  for (G4int i = 0; i < SDCount && i > -1; i++)
	  {
		  //sensitiveDetectors.push_back(SDMan->GetHCtable()->GetSDname(i)); //part of working out how to decide at the time of a run action  which sensitive detectors need event actions
		  sdEdep.push_back(0);
		  sdEdep2.push_back(0);
	  }
	  for (G4int i = 0; i < SDCount && i > -1; i++)
	  {
		  accumulableManager->RegisterAccumulable(sdEdep[i]); //Do this in a second separate loop as vector pointers may become invalid after push_back alters the size of the vector?
		  accumulableManager->RegisterAccumulable(sdEdep2[i]);
	  }
	  accumulablesRegistered = true;
  }
	  for (G4int i = 0; i < SDCount && i > -1; i++)
	  {
		  sdEdep[i].Reset();
		  sdEdep2[i].Reset();
	  }
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FSORunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEventToBeProcessed();
  if (nofEvents == 0) return;


  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge(); //Thankfully Merge() still works with vectorized accumulables

  
  const FSODetectorConstruction* detectorConstruction
   = static_cast<const FSODetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  vector <G4double> masses = detectorConstruction->GetMasses();
  string outputFilename = detectorConstruction->GetOutputFilename();
  G4double particleEnergy;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  
  const FSOPrimaryGeneratorAction* generatorAction
   = static_cast<const FSOPrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition = "";
  if (generatorAction)
  {
    const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += ",";
	particleEnergy = particleGun->GetCurrentSource()->GetEneDist()->GetMonoEnergy();
	std::ostringstream strtemp;
	strtemp << particleEnergy / keV;
	runCondition += strtemp.str();
	runCondition += ",keV,";
	if (particleGun->GetCurrentSource()->GetPosDist()->GetConfined())
	{
		runCondition += particleGun->GetCurrentSource()->GetPosDist()->GetConfineVolume(); //Append which volume the source is confined to
	}
	else runCondition += "false"; //No volume confinement for this run
  } 
#ifdef G4MULTITHREADED
  if (G4Threading::G4GetThreadId() == 0) 
	  //In multithreaded mode, only the worker threads will have access to run conditions
	  //This makes sure that just one of the worker threads outputs the necessary table information
  {
	  stringstream runResults;
	  runResults << nofEvents << "," << runCondition;
	  ofstream myfile;
	  myfile.open(outputFilename, ios::app);
	  myfile << runResults.rdbuf();
	  myfile.close();
  }
  if (IsMaster()) //Master thread time!
  {
	  ifstream myFile;
	  myFile.open(outputFilename, ios::in);
	  myFile.seekg(-1, ios_base::end); //Seek to the character just before EOF
	  bool keepLookingBack = true;
	  while (keepLookingBack)
	  {
		  char tempchar;
		  myFile.get(tempchar);
		  if((int)myFile.tellg() <= 1)
		  {
			  myFile.seekg(0); //Very first line is the very last line, so go to the start
			  keepLookingBack = false;
		  }
		  else if (tempchar == '\n') keepLookingBack = false; //Found a newline!
		  else myFile.seekg(-2, ios_base::cur); //Keep looking backwards
	  }
	  string finalLine;
	  std::getline(myFile, finalLine);
	  myFile.close();
	  char delimiter = ',';
		  vector <string> tokenized_line;
		  stringstream check1(finalLine);
		  string intermediate;
		  while (getline(check1, intermediate, delimiter))
		  {
			  if (intermediate != "") //Merge delimiters style, a bunch of sequential delimiters shouldn't add blank entries into the token list for our purposes
			  {
				  tokenized_line.push_back(intermediate);
			  }
		  }
		  std::istringstream i(tokenized_line[2]);
		  i >> particleEnergy;
		  particleEnergy = particleEnergy * keV;
  }
#endif
  if (IsMaster()) {
  // Print
  //  
  
    G4cout
     << G4endl
	 << "--------------------End of Global Run-----------------------" << G4endl;
	stringstream runResults;
#ifndef G4MULTITHREADED
	runResults
		<< nofEvents << "," << runCondition;
#endif
	for (G4int i = 0; i < SDCount; i++)
	{
		G4double dose = sdEdep[i].GetValue() / nofEvents;
		G4double absFraction = dose * masses[i] / particleEnergy;
		G4String doseReport = G4BestUnit(dose, "Dose");
		doseReport = doseReport.substr(0, doseReport.length() - 1); //Unfortunately G4BestUnit appends a space onto the end, we'll remove that
		std::replace(doseReport.begin(), doseReport.end(), ' ', ','); //Replace the space in the middle with a comma to keep our output in CSV format
		runResults
			<< ","
			<< doseReport << "," << absFraction << "," << std::sqrt((sdEdep2[i].GetValue() / nofEvents - dose * dose) / (nofEvents - 1)) / dose;

	}


	runResults << G4endl;

	G4cout << runResults.str();



	ofstream myfile;
	myfile.open(outputFilename, ios::app);
	myfile << runResults.rdbuf();
	myfile.close();
  }
  else {
    G4cout
     << G4endl
	 << "--------------------End of Local Run: " << runCondition << "------------------------" << G4endl;
  }
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FSORunAction::AddEdeps(std::vector <G4double> energyDepositions)
{
	for (G4int i = 0; i < SDCount; i++)
	{
		sdEdep[i] += energyDepositions[i];
		sdEdep2[i] += energyDepositions[i] * energyDepositions[i];
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

