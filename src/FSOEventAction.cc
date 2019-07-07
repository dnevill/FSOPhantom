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
// $Id: FSOEventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file FSOEventAction.cc
/// \brief Implementation of the FSOEventAction class

#include "FSOEventAction.hh"
#include "FSORunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4THitsMap.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FSOEventAction::FSOEventAction(FSORunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{
	
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FSOEventAction::~FSOEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FSOEventAction::BeginOfEventAction(const G4Event*)
{   
  SDCount = fRunAction->getSDCount();
  sdEdeps.clear();
  for (G4int i = 0; i < SDCount; i++)
  {
	  sdEdeps.push_back(0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FSOEventAction::EndOfEventAction(const G4Event* evt)
{   
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if (!HCE) return;

  for (G4int i = 0; i < SDCount; i++)
  {
	  G4THitsMap<G4double>* eventMap = (G4THitsMap<G4double>*) (HCE->GetHC(i));

	  std::map<G4int, G4double*>::iterator itr = eventMap->GetMap()->begin();
	  for (; itr != eventMap->GetMap()->end(); itr++) {
		  sdEdeps[i] += *(itr->second);
	  }
	  //G4double testResult = *(itr->second);
	  //sdEdeps.push_back(0);
  }
  fRunAction->AddEdeps(sdEdeps);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
