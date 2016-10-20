/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3DAll.h>
#include "Mesh.h"
#include "Tree.h"
#include "BranchDimensions.h"
#include "FruitDimensions.h"

/** \brief Application framework. */
class App : public GApp {
public:
    class Options {//stores parameter options to provide flexible specification of parameters
	public: 
		
    //L-system options
    int maxRecursionDepthL = 4;
	float initialHeightL = 1.0f;
	int circlePtsL = 3;
	int branchSectionsL = 1;
	const Array<String> phenotypesL = Array<String>("Normal", "Random", "Bush", "Pine");
	int phenotypesIndexL = 0;
	const Array<String> branchCallbackL = Array<String>("Straight", "Curvy", "Corkscrew");
	int branchCallbackIndexL = 0;

    //Space colonization options
    int anchorCountSC = 1000;
    float heightSC = 20.0f;
    float radiusSC = 10.0f;
    int circlePointsSC = 10;
    float treeDistanceSC = 0.5f;
    float killDistanceSC = 2.0f;
    float branchRadiusSC = 0.01f;
    float radiusGrowthSC = 2.0f;
    float attractionRadiusSC = 100.0f;
    float leafinessSC = 10.0f;
    float discountRateSC = 0.9f;

    //Orchard options
	int numRows = 2;
	int numTrees = 3;
	const Array<String> types = Array<String>("L-System", "Space Colonization");
	int typesIndex = 0;
	const Array<String> fruits = Array<String>("Apple", "Lemon", "Pear", "Banana", "Money", "Teapot");
	int fruitsIndex = 0;
	Array<FruitDimensions> fruitDims;   
    };

protected:
    Options m_options = Options();

    /** Called from onInit */
    void makeGUI();
    void initializeFruitDims();

    void makeLTree(String filename, Array<Point3>& fruitLocations);
    void generateOrchard();
	void customOrchard();
  
public:
    App(const GApp::Settings& settings = GApp::Settings());
    virtual void onInit() override;
};