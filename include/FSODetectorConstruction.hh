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
// $Id: FSODetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file FSODetectorConstruction.hh
/// \brief Definition of the FSODetectorConstruction class

#ifndef FSODetectorConstruction_h
#define FSODetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <vector>

using namespace std;

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class FSODetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    FSODetectorConstruction();
    virtual ~FSODetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();
	
	std::vector <G4double> GetMasses() const { return masses; }
	string GetOutputFilename() const { return outputFilename; }
  protected:
	std::vector <G4double> masses;
	string outputFilename;
	std::vector <G4LogicalVolume*> stlLogicalVolumes;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

