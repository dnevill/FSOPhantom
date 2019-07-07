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
// $Id: FSODetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file FSODetectorConstruction.cc
/// \brief Implementation of the FSODetectorConstruction class

#include "FSODetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4TriangularFacet.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4TessellatedSolid.hh"
#include "G4TessellatedGeometryAlgorithms.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <codecvt>
#include <locale>
#ifdef _WIN32
#include <Windows.h>
#include <ShObjIdl.h>
#endif
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//These includes are for the steady_clock used for testing loading times
#include <ctime>
#include <ratio>
#include <chrono>


FSODetectorConstruction::FSODetectorConstruction()
	: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FSODetectorConstruction::~FSODetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Takes in an std::string, returns a string vector of each token in the line, parsed using the char delimiter
vector <string> tokenizeLine(string line, char delimiter)
{
	vector <string> tokenized_line;
	stringstream check1(line);
	string intermediate;
	while (getline(check1, intermediate, delimiter))
	{
		if (intermediate != "") //Merge delimiters style, a bunch of sequential delimiters shouldn't add blank entries into the token list for our purposes
		{
			tokenized_line.push_back(intermediate);
		}
	}
	return tokenized_line;
}

vector <wstring> tokenizePWSTR(wstring line, wchar_t delimiter)
{
	vector <wstring> tokenized_line;
	wstringstream check1(line);
	wstring intermediate;
	while (getline(check1, intermediate, delimiter))
	{
		if (intermediate != L"") //Merge delimiters style, a bunch of sequential delimiters shouldn't add blank entries into the token list for our purposes
		{
			tokenized_line.push_back(intermediate);
		}
	}
	return tokenized_line;
}

//Takes a string representation of a double and returns it as a double
inline double convertToDouble(const std::string& s)
{
	std::istringstream i(s);
	double x;
	i >> x;
	return x;
}

//Takes a string representation of a double and returns it as a double
inline G4int convertToInt(const std::string& s)
{
	std::istringstream i(s);
	G4int x;
	i >> x;
	return x;
}

inline bool isNormalInverted(double vertex1[3], double vertex2[3], double vertex3[3], double normal[3])
{
	//Calculate the cross product of the differences between vertex 1 & 2, and 2 & 3 to get the normal for order 1->2->3
	//This will then be compared to the normal given in the STL to see if the facet points this way or if it should be 3->2->1
	double vx1[3] = { vertex2[0] - vertex1[0], vertex2[1] - vertex1[1], vertex2[2] - vertex1[3] };
	double vx2[3] = { vertex3[0] - vertex2[0], vertex3[1] - vertex2[1], vertex3[2] - vertex2[3] };
	double normal_calc[3] = { vx1[1] * vx2[2] - vx1[2] * vx2[1],
		-1 * vx1[0] * vx2[2] + vx1[2] * vx2[0],
		vx1[0] * vx2[1] - vx1[1] * vx2[0] };

	//Sum of the products of X, Y, Z of the normal given in the STL file and the normal calculated using our assumed orientation
	double testDirection = (double)normal_calc[0] * normal[0];
	testDirection += (double)normal_calc[1] * normal[1];
	testDirection += (double)normal_calc[2] * normal[2];
	return testDirection > 0.0;
}

#ifdef _WIN32
inline G4Material* loadMaterialFromFile(G4NistManager* nist, LPCWSTR windowTitle)
{
	G4Material* testTargetMaterial;
	PWSTR pszFilePath2;
	HRESULT hr2 = CoInitializeEx(NULL, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);
	if (SUCCEEDED(hr2))
	{
		IFileOpenDialog *pFileOpen;

		hr2 = CoCreateInstance(CLSID_FileOpenDialog, NULL, CLSCTX_ALL, IID_IFileOpenDialog, reinterpret_cast<void**>(&pFileOpen));
		if (SUCCEEDED(hr2))
		{
			pFileOpen->SetTitle(windowTitle);
			COMDLG_FILTERSPEC rgSpec[] =
			{
				{ L"Material Definition File", L"*.mdef" },
				{ L"All files", L"*.*" },
			};
			pFileOpen->SetFileTypes(2, rgSpec);
			hr2 = pFileOpen->Show(NULL);
			if (SUCCEEDED(hr2))
			{
				IShellItem *pItem;
				hr2 = pFileOpen->GetResult(&pItem);
				if (SUCCEEDED(hr2))
				{

					hr2 = pItem->GetDisplayName(SIGDN_FILESYSPATH, &pszFilePath2);
					if (SUCCEEDED(hr2))
					{

						//Alright, we have a path to a valid file
						
						ifstream myfile;
						myfile.open(pszFilePath2, ios::in);
						if (myfile.is_open())
						{
							string line;

							//int materialCount = count(istreambuf_iterator<char>(myfile), istreambuf_iterator<char>(), '\m'); //This will return the number of lines - 1
							if (getline(myfile, line))
							{
								vector <string> thisToken = tokenizeLine(line, ' ');
								testTargetMaterial = new G4Material(thisToken[0], convertToDouble(thisToken[1])*g / cm3, convertToInt(thisToken[2]));
								//First line of the materialdef file should give the name and density in g/cm^3
							}
							int actual_material_count = 0;
							while (getline(myfile, line))
							{
								vector <string> thisToken = tokenizeLine(line, ' ');
								//Each following line should give the element and the weight fraction
								testTargetMaterial->AddElement(nist->FindOrBuildElement(thisToken[0]), convertToDouble(thisToken[1]));
								actual_material_count++;

							}
							if (actual_material_count != testTargetMaterial->GetNumberOfElements()) //Well, they messed that up, at least let them know that they did
							{
								cout << "Error: MDEF file declared " << testTargetMaterial->GetNumberOfElements() << " elements but " << actual_material_count << " were found." << endl;
							}
							myfile.close();
						}
						else std::cout << "Oops, couldn't open that file";

						//MessageBox(NULL, pszFilePath, L"File Path", MB_OK);
						CoTaskMemFree(pszFilePath2);
					}
					pItem->Release();
				}
			}
			pFileOpen->Release();
		}
		CoUninitialize();
	}
	return testTargetMaterial;
}
#endif


inline G4Material* loadMaterialFromFileGeomCfg(G4NistManager* nist, string pathToMDEF)
{
	G4Material* testTargetMaterial;

	ifstream myfile;
	myfile.open(pathToMDEF, ios::in);
	if (myfile.is_open())
	{
		string line;
		int num_components = 0;
		//int materialCount = count(istreambuf_iterator<char>(myfile), istreambuf_iterator<char>(), '\m'); //This will return the number of lines - 1
		if (getline(myfile, line))
		{
			vector <string> thisToken = tokenizeLine(line, ' ');
			testTargetMaterial = new G4Material(thisToken[0], convertToDouble(thisToken[1])*g / cm3, convertToInt(thisToken[2]));
			num_components = convertToInt(thisToken[2]);
			//First line of the materialdef file should give the name and density in g/cm^3 followed by the number of elements
		}
		int actual_material_count = 0;
		while (getline(myfile, line))
		{
			vector <string> thisToken = tokenizeLine(line, ' ');
			//Each following line should give the element and the weight fraction
			testTargetMaterial->AddElement(nist->FindOrBuildElement(thisToken[0]), convertToDouble(thisToken[1]));
			actual_material_count++;

		}
		if (actual_material_count != num_components) //Well, they messed that up, at least let them know that they did
		{
			cout << "Error: MDEF file declared " << num_components << " elements but " << actual_material_count << " were found." << endl;
		}
		myfile.close();
	}
	else
	{
		std::cout << "Error: Could not open " << pathToMDEF << endl
			<< "Using air at STP instead." << endl;
		testTargetMaterial = nist->FindOrBuildMaterial("G4_AIR");
	}
	return testTargetMaterial;
}

inline bool checkIfFileExists(const std::string& name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	}
	else {
		return false;
	}
}


G4VPhysicalVolume* FSODetectorConstruction::Construct()
{
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();


	// Option to switch on/off checking of volumes overlaps
	//
	G4bool checkOverlaps = true;
	G4bool readFromGeomCfg = false;
	vector <G4double>	x_min, x_max, y_min, y_max, z_min, z_max; //Used while iterating over STLs to ensure we can redefine the world & envelope to fully contain the models and give boundings for source confinement

	stringstream setupResults; //Used to propagate configuration information to the results file

	//vector<G4LogicalVolume*> stlLogicalVolumes; //Used to store all create logical volumes while iterating over STL files, but this is now in the header as a private variable for re-use in constructsdandfield
	G4double xmult = 1., ymult = 1., zmult = 1.;
	string pathToEnvelopeMDEF;
	if (checkIfFileExists("geometry.cfg"))
	{
		readFromGeomCfg = true;
		ifstream geomFile;
		geomFile.open("geometry.cfg", ios::in);
		std::getline(geomFile, outputFilename); //Use readline here so they aren't forced to use quote marks around a long filename
cout << outputFilename << " is output filename" << endl;
		//G4String convertBuffer;
		//std::getline(geomFile, convertBuffer);

		int num_of_bodies;
		geomFile >> num_of_bodies;
		std::cout << num_of_bodies << " STLs to load." << endl;


		double x_offset = 0.0;
		geomFile >> x_offset;

		double y_offset = 0.0;
		geomFile >> y_offset;

		double z_offset = 0.0;
		geomFile >> z_offset;

		std::cout << "XYZ offset in mm: " << x_offset << " " << y_offset << " " << z_offset << endl;
		string temp;
		std::getline(geomFile, temp); //It appears that mixing >> with getline will result in an extra 0 length line the first getline used after >> assignments.
		vector<string> stlPaths, mdefPaths;
		for (int i = 0; i < num_of_bodies; i++)
		{
			string tempBuffer;
			std::getline(geomFile, tempBuffer);
			stlPaths.push_back(tempBuffer);
			std::getline(geomFile, tempBuffer);
			mdefPaths.push_back(tempBuffer);

		}
		geomFile >> xmult;
		if (xmult < 1.0) { xmult = 1.0; std::cout << "Warning: Envelope must be at least large enough to fully contain all loaded STLs." << endl; }
		geomFile >> ymult;
		if (ymult < 1.0) { ymult = 1.0; std::cout << "Warning: Envelope must be at least large enough to fully contain all loaded STLs." << endl; }
		geomFile >> zmult;
		if (zmult < 1.0) { zmult = 1.0; std::cout << "Warning: Envelope must be at least large enough to fully contain all loaded STLs." << endl; }
		std::getline(geomFile, temp); //It appears that mixing >> with getline will result in an extra 0 length line the first getline used after >> assignments.
		std::getline(geomFile, pathToEnvelopeMDEF);


		G4TessellatedSolid* testTarget = new G4TessellatedSolid("Test");
		for (int body_creation_iteration = 0; body_creation_iteration < num_of_bodies; body_creation_iteration++)
		{

			bool good_file = true;
			//Alright, we have a path to a valid file, lets make the geometry

			ifstream myfile;
			myfile.open(stlPaths[body_creation_iteration], ios::in | ios::binary);
			std::cout << "Now loading " << stlPaths[body_creation_iteration] << endl;
#ifdef WIN32
			vector <string> filepath = tokenizeLine(stlPaths[body_creation_iteration], '\\'); //This probably needs an #ifdef by platform for linux machines 
#elif __linux__
			vector <string> filepath = tokenizeLine(stlPaths[body_creation_iteration], '\/');
#endif
			string filename = filepath[filepath.size() - 1];
			string meshname = filename; //Left this way to keep the following code comparable to the UI version
			std::replace(meshname.begin(), meshname.end(), ',', '_'); //The output file's table is in comma separated format, and some filesystems might allow them
			std::replace(meshname.begin(), meshname.end(), ' ', '_'); //Macro commands can have trouble with spaces in volume names, not all accept "" to bracket a long name

			if (myfile.is_open())
			{
				char file_header[80]; //Presently I only actually use the first few characters of the header to make sure this isn't an ASCII STL

cout << "File successfully opened." << endl;

				if (myfile.read(file_header, sizeof(file_header)))
				{
					file_header[6] = 0;
					string checkIfAscii(file_header);
					if (checkIfAscii == "solid ")
					{

						std::cout << "Error: This is an ASCII STL, not a Binary STL.";
						good_file = false;

					}
					else
					{
						uint32_t triangleCount;

						if (myfile.read(reinterpret_cast<char *>(&triangleCount), sizeof(triangleCount)))
						{

							testTarget = new G4TessellatedSolid(meshname);
							for (uint32_t j = 0; j < triangleCount; j++)
							{

								G4float normalRaw[3]; //Normal is given by 3 little-endian 32-bit floating point numbers
								myfile.read(reinterpret_cast<char *>(&normalRaw), sizeof(normalRaw));
								G4float vertex1Raw[3]; //Next comes each of the three vertices as sets of 3 32-bit floating point numbers
								myfile.read(reinterpret_cast<char *>(&vertex1Raw), sizeof(vertex1Raw));
								G4float vertex2Raw[3];
								myfile.read(reinterpret_cast<char *>(&vertex2Raw), sizeof(vertex2Raw));
								G4float vertex3Raw[3];
								myfile.read(reinterpret_cast<char *>(&vertex3Raw), sizeof(vertex3Raw));
								G4ThreeVector normal(normalRaw[0], normalRaw[1], normalRaw[2]); //Convert from a raw array to the Geant4 ThreeVector expected most places
								G4ThreeVector vertex1(vertex1Raw[0] + x_offset, vertex1Raw[1] + y_offset, vertex1Raw[2] + z_offset); //Here we'll also add any offset requested by the user

								if (j == 0) //If this is the first set of vertices loaded, set the minimum and maximum seen so far to the values for this first vertex
								{
									x_max.push_back(vertex1[0]);
									x_min.push_back(vertex1[0]);
									y_max.push_back(vertex1[1]);
									y_min.push_back(vertex1[1]);
									z_max.push_back(vertex1[2]);
									z_min.push_back(vertex1[2]);
								}
								else //Otherwise, check whether any of the maximum or minimum x/y/z seen so far has been exceeded
								{
									if (vertex1[0] > x_max[body_creation_iteration]) x_max[body_creation_iteration] = vertex1[0]; //Track the minimum and maximum boundaries on each axis for world volume and source confinement
									else if (vertex1[0] < x_min[body_creation_iteration]) x_min[body_creation_iteration] = vertex1[0];

									if (vertex1[1] > y_max[body_creation_iteration]) y_max[body_creation_iteration] = vertex1[1];
									else if (vertex1[1] < y_min[body_creation_iteration]) y_min[body_creation_iteration] = vertex1[1];

									if (vertex1[2] > z_max[body_creation_iteration]) z_max[body_creation_iteration] = vertex1[2];
									else if (vertex1[2] < z_min[body_creation_iteration]) z_min[body_creation_iteration] = vertex1[2];
								}


								G4ThreeVector vertex2(vertex2Raw[0] + x_offset, vertex2Raw[1] + y_offset, vertex2Raw[2] + z_offset);

								if (vertex2[0] > x_max[body_creation_iteration]) x_max[body_creation_iteration] = vertex2[0]; //Track the minimum and maximum boundaries on each axis for world volume and source confinement
								else if (vertex2[0] < x_min[body_creation_iteration]) x_min[body_creation_iteration] = vertex2[0];

								if (vertex2[1] > y_max[body_creation_iteration]) y_max[body_creation_iteration] = vertex2[1];
								else if (vertex2[1] < y_min[body_creation_iteration]) y_min[body_creation_iteration] = vertex2[1];

								if (vertex2[2] > z_max[body_creation_iteration]) z_max[body_creation_iteration] = vertex2[2];
								else if (vertex2[2] < z_min[body_creation_iteration]) z_min[body_creation_iteration] = vertex2[2];



								G4ThreeVector vertex3(vertex3Raw[0] + x_offset, vertex3Raw[1] + y_offset, vertex3Raw[2] + z_offset);
								if (vertex3[0] > x_max[body_creation_iteration]) x_max[body_creation_iteration] = vertex3[0]; //Track the minimum and maximum boundaries on each axis for world volume and source confinement
								else if (vertex3[0] < x_min[body_creation_iteration]) x_min[body_creation_iteration] = vertex3[0];

								if (vertex3[1] > y_max[body_creation_iteration]) y_max[body_creation_iteration] = vertex3[1];
								else if (vertex3[1] < y_min[body_creation_iteration]) y_min[body_creation_iteration] = vertex3[1];

								if (vertex3[2] > z_max[body_creation_iteration]) z_max[body_creation_iteration] = vertex3[2];
								else if (vertex3[2] < z_min[body_creation_iteration]) z_min[body_creation_iteration] = vertex3[2];

								uint16_t attributes; //I don't presently use this, some software will store facet colors here, not sure I can color by facet in Geant4
								myfile.read(reinterpret_cast<char *>(&attributes), sizeof(attributes));


								//Calculate the cross product of the differences between vertex 1 & 2, and 2 & 3 to get the normal for order 1->2->3
								//This will then be compared to the normal given in the STL to see if the facet points this way or if it should be 3->2->1
								G4ThreeVector vx1(vertex2.getX() - vertex1.getX(), vertex2.getY() - vertex1.getY(), vertex2.getZ() - vertex1.getZ());
								G4ThreeVector vx2(vertex3.getX() - vertex2.getX(), vertex3.getY() - vertex2.getY(), vertex3.getZ() - vertex2.getZ());
								G4ThreeVector normal_calc(vx1.getY()*vx2.getZ() - vx1.getZ()*vx2.getY(),
									-1 * vx1.getX()*vx2.getZ() + vx1.getZ()*vx2.getX(),
									vx1.getX()*vx2.getY() - vx1.getY()*vx2.getX());

								//Sum of the products of X, Y, Z of the normal given in the STL file and the normal calculated using our assumed orientation
								G4double testDirection = (G4double)normal_calc.getX()*normal.getX();
								testDirection += (G4double)normal_calc.getY()*normal.getY();
								testDirection += (G4double)normal_calc.getZ()*normal.getZ();



								G4double vtx1_d[3] = { vertex1.getX(), vertex1.getY(), vertex1.getZ() };
								G4double vtx2_d[3] = { vertex2.getX(), vertex2.getY(), vertex2.getZ() };
								G4double vtx3_d[3] = { vertex3.getX(), vertex3.getY(), vertex3.getZ() };
								G4double nrml_d[3] = { normal.getX(), normal.getY(), normal.getZ() };
								bool checkNormalDirection = isNormalInverted(vtx1_d, vtx2_d, vtx3_d, nrml_d); //I could re-write this to natively handle G4ThreeVectors
								if (testDirection > 0)
								{
									G4TriangularFacet *thisFacet = new G4TriangularFacet(vertex1, vertex2, vertex3, (G4FacetVertexType)0); //0 = Absolute
									testTarget->AddFacet((G4VFacet*)thisFacet);
								}
								else
								{
									G4TriangularFacet *thisFacetBackwards = new G4TriangularFacet(vertex3, vertex2, vertex1, (G4FacetVertexType)0); //0 = Absolute
									testTarget->AddFacet((G4VFacet*)thisFacetBackwards);
								}

								string derp = ""; //just here to insert a breakpoint after the above have finished.
							}
							//Okay, now this line has been vectorized when I todo put the new vectorize function here


						}
					}
				}
				else
				{
					std::cout << "Error: File does not start with a solid declaration. Are you sure this is an ASCII STL file?" << endl;
				}


			}
			myfile.close();
			if (good_file)
			{
				std::cout << "Setting solid closed...this may take a while!" << endl;
				std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
				testTarget->SetSolidClosed(true);
				std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
				std::chrono::duration<double> time_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
				std::cout << time_elapsed.count() << " seconds elapsed while closing this solid." << endl;
			}



			if (good_file)
			{
				G4Material* testTargetMaterial = loadMaterialFromFileGeomCfg(nist, mdefPaths[body_creation_iteration]);

				//G4cout << *(G4Material::GetMaterialTable()) << endl;

				G4LogicalVolume* testTarget_logical = new G4LogicalVolume(testTarget, testTargetMaterial, testTarget->GetName(), 0, 0, 0);
				masses.push_back(testTarget_logical->GetMass());
				stlLogicalVolumes.push_back(testTarget_logical);
			}
			else body_creation_iteration--; //Bad STL file was loaded, so this iteration didn't actually get to add a new STL.  Decrement the counter and retry.
		}
	}
#ifdef _WIN32
	else
	{
		std::cout << "Name for run-results output file: " << endl;
		std::getline(std::cin, outputFilename);

		

		std::cout << "How many STLs will be loaded?" << endl;
		int num_of_bodies;
		std::cin >> num_of_bodies;

		std::cout << "Transform needed?" << endl << "X offset (mm): ";
		double x_offset = 0.0;
		std::cin >> x_offset;
		std::cout << endl << "Y offset (mm): ";
		double y_offset = 0.0;
		std::cin >> y_offset;

		std::cout << endl << "Z offset (mm): ";
		double z_offset = 0.0;
		std::cin >> z_offset;

		G4TessellatedSolid* testTarget = new G4TessellatedSolid("Test");
		for (int body_creation_iteration = 0; body_creation_iteration < num_of_bodies; body_creation_iteration++)
		{

			bool good_file = true;
			PWSTR pszFilePath;
			HRESULT hr = CoInitializeEx(NULL, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);
			if (SUCCEEDED(hr))
			{
				IFileOpenDialog *pFileOpen;

				hr = CoCreateInstance(CLSID_FileOpenDialog, NULL, CLSCTX_ALL, IID_IFileOpenDialog, reinterpret_cast<void**>(&pFileOpen));
				if (SUCCEEDED(hr))
				{
					pFileOpen->SetTitle(L"Select Mesh (STL) Definition File");
					COMDLG_FILTERSPEC rgSpec[] =
					{
						{ L"STL Mesh", L"*.stl" },
						{ L"All files", L"*.*" },
					};
					pFileOpen->SetFileTypes(2, rgSpec);
					hr = pFileOpen->Show(NULL);
					if (SUCCEEDED(hr))
					{
						IShellItem *pItem;
						hr = pFileOpen->GetResult(&pItem);
						if (SUCCEEDED(hr))
						{

							hr = pItem->GetDisplayName(SIGDN_FILESYSPATH, &pszFilePath);
							if (SUCCEEDED(hr))
							{

								//Alright, we have a path to a valid file, lets make the geometry

								ifstream myfile;
								myfile.open(pszFilePath, ios::in | ios::binary);
								std::wcout << "Now loading " << pszFilePath << endl;
								vector <wstring> filepath = tokenizePWSTR(pszFilePath, L'\\');
								wstring filename = filepath[filepath.size() - 1];
								using convert_type = std::codecvt_utf8<wchar_t>;
								std::wstring_convert<convert_type, wchar_t> converter;
								string meshname = converter.to_bytes(filename);
								std::replace(meshname.begin(), meshname.end(), ',', '_'); //The output file's table is in comma separated format, and some filesystems might allow them
								std::replace(meshname.begin(), meshname.end(), ' ', '_'); //Macro commands can have trouble with spaces in volume names, not all accept "" to bracket a long name

								if (myfile.is_open())
								{
									char file_header[80]; //Presently I only actually use the first few characters of the header to make sure this isn't an ASCII STL



									if (myfile.read(file_header, sizeof(file_header)))
									{
										file_header[6] = 0;
										string checkIfAscii(file_header);
										if (checkIfAscii == "solid ")
										{

											std::cout << "Error: This is an ASCII STL, not a Binary STL.";
											good_file = false;

										}
										else
										{
											uint32_t triangleCount;

											if (myfile.read(reinterpret_cast<char *>(&triangleCount), sizeof(triangleCount)))
											{

												testTarget = new G4TessellatedSolid(meshname);
												for (uint32_t j = 0; j < triangleCount; j++)
												{

													G4float normalRaw[3]; //Normal is given by 3 little-endian 32-bit floating point numbers
													myfile.read(reinterpret_cast<char *>(&normalRaw), sizeof(normalRaw));
													G4float vertex1Raw[3]; //Next comes each of the three vertices as sets of 3 32-bit floating point numbers
													myfile.read(reinterpret_cast<char *>(&vertex1Raw), sizeof(vertex1Raw));
													G4float vertex2Raw[3];
													myfile.read(reinterpret_cast<char *>(&vertex2Raw), sizeof(vertex2Raw));
													G4float vertex3Raw[3];
													myfile.read(reinterpret_cast<char *>(&vertex3Raw), sizeof(vertex3Raw));
													G4ThreeVector normal(normalRaw[0], normalRaw[1], normalRaw[2]); //Convert from a raw array to the Geant4 ThreeVector expected most places
													G4ThreeVector vertex1(vertex1Raw[0] + x_offset, vertex1Raw[1] + y_offset, vertex1Raw[2] + z_offset); //Here we'll also add any offset requested by the user

													if (j == 0) //If this is the first set of vertices loaded, set the minimum and maximum seen so far to the values for this first vertex
													{
														x_max.push_back(vertex1[0]);
														x_min.push_back(vertex1[0]);
														y_max.push_back(vertex1[1]);
														y_min.push_back(vertex1[1]);
														z_max.push_back(vertex1[2]);
														z_min.push_back(vertex1[2]);
													}
													else //Otherwise, check whether any of the maximum or minimum x/y/z seen so far has been exceeded
													{
														if (vertex1[0] > x_max[body_creation_iteration]) x_max[body_creation_iteration] = vertex1[0]; //Track the minimum and maximum boundaries on each axis for world volume and source confinement
														else if (vertex1[0] < x_min[body_creation_iteration]) x_min[body_creation_iteration] = vertex1[0];

														if (vertex1[1] > y_max[body_creation_iteration]) y_max[body_creation_iteration] = vertex1[1];
														else if (vertex1[1] < y_min[body_creation_iteration]) y_min[body_creation_iteration] = vertex1[1];

														if (vertex1[2] > z_max[body_creation_iteration]) z_max[body_creation_iteration] = vertex1[2];
														else if (vertex1[2] < z_min[body_creation_iteration]) z_min[body_creation_iteration] = vertex1[2];
													}


													G4ThreeVector vertex2(vertex2Raw[0] + x_offset, vertex2Raw[1] + y_offset, vertex2Raw[2] + z_offset);

													if (vertex2[0] > x_max[body_creation_iteration]) x_max[body_creation_iteration] = vertex2[0]; //Track the minimum and maximum boundaries on each axis for world volume and source confinement
													else if (vertex2[0] < x_min[body_creation_iteration]) x_min[body_creation_iteration] = vertex2[0];

													if (vertex2[1] > y_max[body_creation_iteration]) y_max[body_creation_iteration] = vertex2[1];
													else if (vertex2[1] < y_min[body_creation_iteration]) y_min[body_creation_iteration] = vertex2[1];

													if (vertex2[2] > z_max[body_creation_iteration]) z_max[body_creation_iteration] = vertex2[2];
													else if (vertex2[2] < z_min[body_creation_iteration]) z_min[body_creation_iteration] = vertex2[2];



													G4ThreeVector vertex3(vertex3Raw[0] + x_offset, vertex3Raw[1] + y_offset, vertex3Raw[2] + z_offset);
													if (vertex3[0] > x_max[body_creation_iteration]) x_max[body_creation_iteration] = vertex3[0]; //Track the minimum and maximum boundaries on each axis for world volume and source confinement
													else if (vertex3[0] < x_min[body_creation_iteration]) x_min[body_creation_iteration] = vertex3[0];

													if (vertex3[1] > y_max[body_creation_iteration]) y_max[body_creation_iteration] = vertex3[1];
													else if (vertex3[1] < y_min[body_creation_iteration]) y_min[body_creation_iteration] = vertex3[1];

													if (vertex3[2] > z_max[body_creation_iteration]) z_max[body_creation_iteration] = vertex3[2];
													else if (vertex3[2] < z_min[body_creation_iteration]) z_min[body_creation_iteration] = vertex3[2];

													uint16_t attributes; //I don't presently use this, some software will store facet colors here, not sure I can color by facet in Geant4
													myfile.read(reinterpret_cast<char *>(&attributes), sizeof(attributes));


													//Calculate the cross product of the differences between vertex 1 & 2, and 2 & 3 to get the normal for order 1->2->3
													//This will then be compared to the normal given in the STL to see if the facet points this way or if it should be 3->2->1
													G4ThreeVector vx1(vertex2.getX() - vertex1.getX(), vertex2.getY() - vertex1.getY(), vertex2.getZ() - vertex1.getZ());
													G4ThreeVector vx2(vertex3.getX() - vertex2.getX(), vertex3.getY() - vertex2.getY(), vertex3.getZ() - vertex2.getZ());
													G4ThreeVector normal_calc(vx1.getY()*vx2.getZ() - vx1.getZ()*vx2.getY(),
														-1 * vx1.getX()*vx2.getZ() + vx1.getZ()*vx2.getX(),
														vx1.getX()*vx2.getY() - vx1.getY()*vx2.getX());

													//Sum of the products of X, Y, Z of the normal given in the STL file and the normal calculated using our assumed orientation
													G4double testDirection = (G4double)normal_calc.getX()*normal.getX();
													testDirection += (G4double)normal_calc.getY()*normal.getY();
													testDirection += (G4double)normal_calc.getZ()*normal.getZ();



													G4double vtx1_d[3] = { vertex1.getX(), vertex1.getY(), vertex1.getZ() };
													G4double vtx2_d[3] = { vertex2.getX(), vertex2.getY(), vertex2.getZ() };
													G4double vtx3_d[3] = { vertex3.getX(), vertex3.getY(), vertex3.getZ() };
													G4double nrml_d[3] = { normal.getX(), normal.getY(), normal.getZ() };
													bool checkNormalDirection = isNormalInverted(vtx1_d, vtx2_d, vtx3_d, nrml_d); //I could re-write this to natively handle G4ThreeVectors
													if (testDirection > 0)
													{
														G4TriangularFacet *thisFacet = new G4TriangularFacet(vertex1, vertex2, vertex3, (G4FacetVertexType)0); //0 = Absolute
														testTarget->AddFacet((G4VFacet*)thisFacet);
													}
													else
													{
														G4TriangularFacet *thisFacetBackwards = new G4TriangularFacet(vertex3, vertex2, vertex1, (G4FacetVertexType)0); //0 = Absolute
														testTarget->AddFacet((G4VFacet*)thisFacetBackwards);
													}

													string derp = ""; //just here to insert a breakpoint after the above have finished.
												}
												//Okay, now this line has been vectorized when I todo put the new vectorize function here


											}
										}
									}
									else
									{
										std::cout << "Error: File does not start with a solid declaration. Are you sure this is an ASCII STL file?" << endl;
									}


								}
								myfile.close();
								if (good_file)
								{
									std::cout << "Setting solid closed...this may take a while!" << endl;
									std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
									testTarget->SetSolidClosed(true);
									std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
									std::chrono::duration<double> time_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
									std::cout << time_elapsed.count() << " seconds elapsed while closing this solid." << endl;
								}
							}
							else G4cout << "Oops, couldn't open that file";

							//MessageBox(NULL, pszFilePath, L"File Path", MB_OK);
							CoTaskMemFree(pszFilePath);
						}
						pItem->Release();
					}
				}
				pFileOpen->Release();
			}
			CoUninitialize();

			if (good_file)
			{
				G4Material* testTargetMaterial = loadMaterialFromFile(nist, L"Select Material Definition File");

				//G4cout << *(G4Material::GetMaterialTable()) << endl;

				G4LogicalVolume* testTarget_logical = new G4LogicalVolume(testTarget, testTargetMaterial, testTarget->GetName(), 0, 0, 0);
				masses.push_back(testTarget_logical->GetMass());
				stlLogicalVolumes.push_back(testTarget_logical);
			}
			else body_creation_iteration--; //Bad STL file was loaded, so this iteration didn't actually get to add a new STL.  Decrement the counter and retry.
		}
	}
#elif __linux__
else {
cout << "Error: FSOPhantom currently does not support interactive file loading. Please use geometry.cfg." << endl;
}
#endif
	//Get the maximum extent from all loaded STLs
	G4double envx_max = *max_element(x_max.begin(), x_max.end());
	G4double envx_min = *min_element(x_min.begin(), x_min.end());
	G4double envy_max = *max_element(y_max.begin(), y_max.end());
	G4double envy_min = *min_element(y_min.begin(), y_min.end());
	G4double envz_max = *max_element(z_max.begin(), z_max.end());
	G4double envz_min = *min_element(z_min.begin(), z_min.end());


	//Define a box that would perfectly contain the maximum extent of all the loaded STLs
	G4double envhalfx = (envx_max - envx_min) / 2.0;
	G4double envhalfy = (envy_max - envy_min) / 2.0;
	G4double envhalfz = (envz_max - envz_min) / 2.0;
	G4double envcenterx = (envx_max + envx_min) / 2.0;
	G4double envcentery = (envy_max + envy_min) / 2.0;
	G4double envcenterz = (envz_max + envz_min) / 2.0;

	// Envelope parameters
	//
	if (!readFromGeomCfg) 
	{
	std::cout << "An envelope that is " << envhalfx * 2 << "mm (X) by " << envhalfy * 2 << "mm (Y) by " << envhalfz * 2 << "mm (Z) is the minimum to contain the models" << endl
		<< "How many times larger should the envelope be per axis (where 1.0 = exactly the same, 2.0 = double, etc., and larger multipliers increase accuracy" << endl
		<< "X multiplier: ";

	std::cin >> xmult;
	if (xmult < 1.0) xmult = 1.0;
	std::cout << "Y multiplier: ";
	std::cin >> ymult;
	if (ymult < 1.0) ymult = 1.0;
	std::cout << "Z multiplier: ";
	std::cin >> zmult;
	if (zmult < 1.0) zmult = 1.0;
	}



	G4double env_sizeX = envhalfx * xmult, env_sizeY = envhalfy * ymult, env_sizeZ = envhalfz * zmult; //Envelope will be mult times larger in each dimension than the maximum extent of all STLs
	G4Material* env_mat;
	if (readFromGeomCfg) 
	{
		env_mat = loadMaterialFromFileGeomCfg(nist, pathToEnvelopeMDEF);
	}
	else 
	{
#ifdef _WIN32
		env_mat = loadMaterialFromFile(nist, L"Select Material Definition File for surrounding environment");
#elif __linux__
                //env_mat = loadMaterialFromFile(nist); Add a function to ask from the terminal interactively for the file?
#endif
	}
	//     
	// World
	//
	//Everything in GEANT4 must be within the world volume, so we make it slightly larger than the envelope
	//(We don't just use the envelope as the world volume to support, in the future, a separate ground/air or sediment/water geometry option)
	//Note that these sizes are half-lengths per box specifications
	G4double world_sizeX = 1.01*env_sizeX + abs(envcenterx);

	/*
	If the user did not center their volume before segmenting, the resulting model might be very far away from the origin.
	In that case, extend the size of the world volume on that axis to include it.  Since the World volume *must* be centered
	on the origin, this means the World will also extend an equal distance on the other side of the origin from the loaded STLs as well*/

	G4double world_sizeY = 1.01*env_sizeY + abs(envcentery);

	G4double world_sizeZ = 1.01*env_sizeZ + abs(envcenterz);


	G4Box* solidWorld =
		new G4Box("World",                       //its name
		world_sizeX, world_sizeY, world_sizeZ);     //its size

	G4LogicalVolume* logicWorld =
		new G4LogicalVolume(solidWorld,          //its solid
		env_mat,           //its material
		"World");            //its name

	G4VPhysicalVolume* physWorld =
		new G4PVPlacement(0,                     //no rotation
		G4ThreeVector(),       //World volume MUST be centered on 0,0,0
		logicWorld,            //its logical volume
		"World",               //its name
		0,                     //its mother  volume
		false,                 //no boolean operation
		0,                     //copy number
		checkOverlaps);        //overlaps checking

	//     
	// Envelope
	//  
	G4Box* solidEnv =
		new G4Box("Envelope",                    //its name
		env_sizeX, env_sizeY, env_sizeZ); //its size

	G4LogicalVolume* logicEnv =
		new G4LogicalVolume(solidEnv,            //its solid
		env_mat,             //its material
		"Envelope");         //its name
	std::cout << "Envelope volume Mass: " << logicEnv->GetMass() << G4endl;

	setupResults << "Envelope: " << G4endl
		<< "Boundaries (mm): " << G4endl
		<< "x_max " << envcenterx + env_sizeX << " x_min: " << envcenterx - env_sizeX
		<< " y_max " << envcentery + env_sizeY << " y_min: " << envcentery - env_sizeY
		<< " z_max " << envcenterz + env_sizeZ << " z_min: " << envcenterz - env_sizeZ << G4endl
		<< "Mass: " << G4BestUnit(logicEnv->GetMass(), "Mass") << G4endl
		<< logicEnv->GetMaterial() << G4endl;

	new G4PVPlacement(0,                       //no rotation
		G4ThreeVector(envcenterx, envcentery, envcenterz),         //centered on the max extent of the loaded STLs
		logicEnv,                //its logical volume
		"Envelope",              //its name
		logicWorld,              //its mother  volume
		false,                   //no boolean operation
		0,                       //copy number
		checkOverlaps);          //overlaps checking


	std::cout << "To confine source particles to be isometrically emitted from the entire envelope of surrounding material execute the following commands: " << G4endl
		<< "/gps/pos/type Volume" << G4endl
		<< "/gps/pos/shape Para" << G4endl
		<< "/gps/pos/halfx " << env_sizeX / 10. << G4endl
		<< "/gps/pos/halfy " << env_sizeY / 10. << G4endl
		<< "/gps/pos/halfz " << env_sizeZ / 10. << G4endl
		<< "/gps/pos/centre " << envcenterx / 10. << " " << envcentery / 10. << " " << envcenterz / 10. << G4endl
		<< "/gps/pos/confine Envelope" << G4endl
		<< "/gps/ang/type iso" << G4endl;


	stringstream tableHeader;
	tableHeader << "NoOfEvents,Particle,Energy,unit";

	for (uint32_t i = 0; i < stlLogicalVolumes.size(); i++)
	{
		G4String volume_name = stlLogicalVolumes[i]->GetName();
		tableHeader << "," << volume_name << " Mean Dose"
			<< "," << "Unit"
			<< "," << "AbsorbedFraction"
			<< "," << "Mean Dose Fractional Standard Deviation";
		G4VPhysicalVolume* testTarget_physical = new G4PVPlacement(0, G4ThreeVector(-1 * envcenterx, -1 * envcentery, -1 * envcenterz), stlLogicalVolumes[i], volume_name, logicEnv, false, 0, checkOverlaps);
		cout << volume_name << " boundaries are: " << G4endl
			<< "x_max " << x_max[i] << " x_min: " << x_min[i] << " y_max: " << y_max[i]
			<< " y_min: " << y_min[i] << " z_max: " << z_max[i] << " z_min: " << z_min[i] << G4endl;
		setupResults << volume_name << G4endl
			<< "Boundaries: " << G4endl
			<< "x_max " << x_max[i] << " x_min: " << x_min[i] << " y_max: " << y_max[i]
			<< " y_min: " << y_min[i] << " z_max: " << z_max[i] << " z_min: " << z_min[i] << G4endl;
		cout << "Mass: " << G4BestUnit(stlLogicalVolumes[i]->GetMass(), "Mass") << G4endl;
		setupResults << "Mass: " << G4BestUnit(stlLogicalVolumes[i]->GetMass(), "Mass") << G4endl;
		setupResults << stlLogicalVolumes[i]->GetMaterial() << G4endl;

		G4double halfx = (x_max[i] - x_min[i]) / 20.0; //Dividing by 2 to get the half-length, then by another 10 since these values are in mm and the user commands default to cm
		G4double halfy = (y_max[i] - y_min[i]) / 20.0;
		G4double halfz = (z_max[i] - z_min[i]) / 20.0;

		//Dividing by 2 to get the average, then by another 10 since these values are in mm and the user commands default to cm
		//Then add the offset from world coordinates that the envelope itself will be at (it can be offset as part of handling STLs that are not centered with respect to their own coordinate system)
		G4double centerx = (x_max[i] + x_min[i]) / 20.0;
		G4double centery = (y_max[i] + y_min[i]) / 20.0;
		G4double centerz = (z_max[i] + z_min[i]) / 20.0;



		cout << "To confine source particles to be isometrically emitted from this volume, execute the following commands: " << G4endl
			<< "/gps/pos/type Volume" << G4endl
			<< "/gps/pos/shape Para" << G4endl
			<< "/gps/pos/halfx " << halfx << G4endl
			<< "/gps/pos/halfy " << halfy << G4endl
			<< "/gps/pos/halfz " << halfz << G4endl
			<< "/gps/pos/centre " << centerx << " " << centery << " " << centerz << G4endl
			<< "/gps/pos/confine " << volume_name << G4endl
			<< "/gps/ang/type iso" << G4endl;

	}


	ofstream myfile;
	myfile.open(outputFilename, ios::app);
	myfile << setupResults.rdbuf();
	myfile << tableHeader.rdbuf() << G4endl;
	myfile.close();
	//return the physical World
	return physWorld;
}

void FSODetectorConstruction::ConstructSDandField()
{
	for (int i = 0; i < stlLogicalVolumes.size(); i++)
	{
		G4MultiFunctionalDetector* testTarget_MFD = new G4MultiFunctionalDetector("Test Target MFD " + to_string(i)); //Create a new MFD
		G4SDManager::GetSDMpointer()->AddNewDetector(testTarget_MFD); //Alert the SD tracker of the new detector
		stlLogicalVolumes[i]->SetSensitiveDetector(testTarget_MFD); //Alert the logical volume it is a sensitive volume using this detector

		G4VPrimitiveScorer* testTarget_PSDose = new G4PSDoseDeposit("Test Target PSDose " + to_string(i)); //Create a total dose primitive scorer
		testTarget_MFD->RegisterPrimitive(testTarget_PSDose); //Alert the MFD of the new primitive scorer
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
