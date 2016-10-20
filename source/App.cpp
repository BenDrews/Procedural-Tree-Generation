/** \file App.cpp */
#include "App.h"
#include "Mesh.h"
#include "BranchDimensions.h"
#include "Tree.h"
#include "FruitDimensions.h"
#include "SCGenerator.h"
#include "LGenerator.h"
#include <cmath>
#include <map>
#include <tuple>
#include <stdlib.h> 
#include <iostream>
#include <sstream>

// Tells C++ to invoke command-line main() function even on OS X and Win32.
G3D_START_AT_MAIN();

int main(int argc, const char* argv[]) {
    {
        G3DSpecification g3dSpec;
        g3dSpec.audio = false;
        initGLG3D(g3dSpec);
    }

    GApp::Settings settings(argc, argv);

    settings.window.caption             = argv[0];
    settings.window.width               = 1280; settings.window.height       = 720; settings.window.fullScreen          = false;
    settings.window.resizable           = ! settings.window.fullScreen;
    settings.window.framed              = ! settings.window.fullScreen;
    settings.window.asynchronous        = false;

    settings.hdrFramebuffer.depthGuardBandThickness = Vector2int16(64, 64);
    settings.hdrFramebuffer.colorGuardBandThickness = Vector2int16(0, 0);
    settings.dataDir                    = FileSystem::currentDirectory();
    settings.screenshotDirectory        = "../journal/";

    settings.renderer.deferredShading = true;
    settings.renderer.orderIndependentTransparency = false;

    return App(settings).run();
}


App::App(const GApp::Settings& settings) : GApp(settings) {
}

// Called before the application loop begins.  Load data here and
// not in the constructor so that common exceptions will be
// automatically caught.
void App::onInit() {
    GApp::onInit();
    setFrameDuration(1.0f / 120.0f);

    showRenderingStats      = false;

    makeGUI();
    developerWindow->cameraControlWindow->moveTo(Point2(developerWindow->cameraControlWindow->rect().x0(), 0));

    loadScene("Tree Testing");

	initializeFruitDims();
}

void App::initializeFruitDims() {
    m_options.fruitDims = Array<FruitDimensions>();
	
	FruitDimensions appleDims = FruitDimensions("apple textured obj.obj", m_options.initialHeightL / 1000.0f, 0.08f, 0.0f);
	FruitDimensions lemonDims = FruitDimensions("lemon whole.obj", m_options.initialHeightL / 1000.0f, 0.04f, 0.0f);
	FruitDimensions pearDims = FruitDimensions("pear_export.obj", m_options.initialHeightL / 10.0f, 0.04f, 0.0f);
	FruitDimensions bananaDims = FruitDimensions("banana.obj", m_options.initialHeightL / 2000.0f, 0.04f, -90.0f);
	FruitDimensions moneyDims = FruitDimensions("Dollar stack wild.obj", m_options.initialHeightL / 100.0f, 0.08f, 90.0f);
	FruitDimensions teapotDims = FruitDimensions("glassTeapot.obj", m_options.initialHeightL / 1000.0f, 0.08f, 0.0f);

	m_options.fruitDims.push(appleDims);
	m_options.fruitDims.push(lemonDims);
	m_options.fruitDims.push(pearDims);
	m_options.fruitDims.push(bananaDims);
	m_options.fruitDims.push(moneyDims);
	m_options.fruitDims.push(teapotDims);
}


void App::makeGUI() {
    // Initialize the developer HUD
    createDeveloperHUD();

    debugWindow->setVisible(true);
    developerWindow->videoRecordDialog->setEnabled(true);

    debugPane->beginRow(); {
        GuiTabPane* containerPane = debugPane->addTabPane();
	    
        // L-System tree generation GUI
        GuiPane* treeLPane = containerPane->addTab("L-System Tree");
        treeLPane->setNewChildSize(500, -1, 300);
	    treeLPane->addNumberBox("Recursion depth:", &m_options.maxRecursionDepthL, "", GuiTheme::LINEAR_SLIDER, 2, 10);
	    treeLPane->addNumberBox("Initial height:", &m_options.initialHeightL, "", GuiTheme::LOG_SLIDER, 0.5f, 10.0f);
	    treeLPane->addNumberBox("Circle points:", &m_options.circlePtsL, "", GuiTheme::LOG_SLIDER, 3, 100);
	    treeLPane->addNumberBox("Branch sections:", &m_options.branchSectionsL, "", GuiTheme::LOG_SLIDER, 1, 100);
        treeLPane->addDropDownList("Phenotype", m_options.phenotypesL, &m_options.phenotypesIndexL);
        treeLPane->addDropDownList("Branch callback function", m_options.branchCallbackL, &m_options.branchCallbackIndexL);
        treeLPane->addCheckBox("Autumn", &m_options.fall);


        treeLPane->addButton("Generate tree", [this](){
	    	drawMessage("Generating tree...");
	    	Array<Point3> fruitLocations = Array<Point3>();

            Stopwatch sw;
    	    sw.tick(); //start the timer
            makeLTree("tree", fruitLocations);
            //generateForest();
         
            sw.tock();

            debugPrintf("Elapsed Time: %f\n", sw.elapsedTime());
	    	ArticulatedModel::clearCache();
	    	GApp::loadScene("Tree Testing");
	    });
		treeLPane->pack();

        // Space tree generation GUI
        GuiPane* treeSCPane = containerPane->addTab("Space Col Tree");
        treeSCPane->setNewChildSize(500, -1, 300);
	    treeSCPane->addNumberBox("Anchor points:", &m_options.anchorCountSC, "", GuiTheme::LOG_SLIDER, 1, 10000);
	    treeSCPane->addNumberBox("Height:", &m_options.heightSC, "", GuiTheme::LINEAR_SLIDER, 10.0f, 100.0f);
	    treeSCPane->addNumberBox("Radius:", &m_options.radiusSC, "", GuiTheme::LINEAR_SLIDER, 10.0f, 100.0f);
	    treeSCPane->addNumberBox("Circle points:", &m_options.circlePointsSC, "", GuiTheme::LOG_SLIDER, 1, 100);
	    treeSCPane->addNumberBox("Tree node distance:", &m_options.treeDistanceSC, "", GuiTheme::LOG_SLIDER, 0.01f, 1.0f);
	    treeSCPane->addNumberBox("Kill distance:", &m_options.killDistanceSC, "", GuiTheme::LINEAR_SLIDER, 1.0f, 10.0f);
	    treeSCPane->addNumberBox("Branch initial radius:", &m_options.branchRadiusSC, "", GuiTheme::LOG_SLIDER, 0.01f, 10.0f);
	    treeSCPane->addNumberBox("Radius growth:", &m_options.radiusGrowthSC, "", GuiTheme::LINEAR_SLIDER, 2.0f, 3.0f);
	    treeSCPane->addNumberBox("Attraction radius:", &m_options.attractionRadiusSC, "", GuiTheme::LOG_SLIDER, 1.0f, 100.0f);
	    treeSCPane->addNumberBox("Leafiness:", &m_options.leafinessSC, "", GuiTheme::LOG_SLIDER, 1.0f, 100.0f);
	    treeSCPane->addNumberBox("Discount Rate:", &m_options.discountRateSC, "", GuiTheme::LINEAR_SLIDER, 0.0f, 1.0f);
        treeSCPane->addButton("Generate tree", [this](){
	    	drawMessage("Generating tree...");
            
            Array<Point3> fruitLocations;
            SCGenerator genSC;
            shared_ptr<Tree> skeleton = genSC.makeSCTreeSkeleton(m_options.anchorCountSC, [this](float y) {return SCGenerator::bulbEnvelope(y);}, m_options.heightSC, m_options.radiusSC, m_options.killDistanceSC, m_options.treeDistanceSC, m_options.attractionRadiusSC, m_options.discountRateSC, Point3(0,0,0));
            
	        Stopwatch sw;
    	    sw.tick(); //start the timer
            genSC.skeletonToMeshSC(m_options.circlePointsSC, m_options.branchRadiusSC, m_options.radiusGrowthSC, m_options.leafinessSC, "tree", skeleton, fruitLocations);
            sw.tock();

            debugPrintf("Elapsed Time: %f\n", sw.elapsedTime());
	    	ArticulatedModel::clearCache();
	    	GApp::loadScene("Tree Testing");
	    });
	    treeSCPane->pack();

		GuiPane* orchardPane = containerPane->addTab("Orchard");
		orchardPane->setNewChildSize(500, -1, 300);
		orchardPane->addNumberBox("Number of rows:", &m_options.numRows, "", GuiTheme::LOG_SLIDER, 1, 20);
	    orchardPane->addNumberBox("Trees per row:", &m_options.numTrees, "", GuiTheme::LOG_SLIDER, 1, 20);
		orchardPane->addDropDownList("Generation type", m_options.types, &m_options.typesIndex);
        orchardPane->addDropDownList("Fruit type", m_options.fruits, &m_options.fruitsIndex);
		orchardPane->addButton("Generate orchard", [this]() {
	    	drawMessage("Generating orchard...");
	    	generateOrchard();
	    	ArticulatedModel::clearCache();
	    	GApp::loadScene("Orchard");
	    });
		orchardPane->addButton("Generate custom orchard", [this]() {
	    	drawMessage("Generating custom orchard...");
	    	customOrchard();
	    	ArticulatedModel::clearCache();
	    	GApp::loadScene("Custom Orchard");
	    });
	    orchardPane->pack();
        }
        debugPane->endRow();
        debugWindow->pack();
        debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}

void App::makeLTree(String filename, Array<Point3>& fruitLocations) {
    LGenerator genL;
    std::function<void(Array<BranchDimensions>&, float, const CoordinateFrame&, Point3&, int, int)> phenotype;
    Mesh treeMesh = Mesh(filename);
    Mesh leafMesh = Mesh("leaf");

	if (m_options.phenotypesIndexL == 0) {
		phenotype = [this](Array<BranchDimensions>& nextBranches, float initialLength, const CoordinateFrame& initial, const Point3& branchEnd, int maxRecursionDepth, int currentRecursionDepth)
			{return LGenerator::normalTree(nextBranches, initialLength, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);};
	}else if (m_options.phenotypesIndexL == 1) {
        phenotype = [this](Array<BranchDimensions>& nextBranches, float initialLength, const CoordinateFrame& initial, const Point3& branchEnd, int maxRecursionDepth, int currentRecursionDepth)
			{return LGenerator::randomTree(nextBranches, initialLength, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);};
	}else if (m_options.phenotypesIndexL == 2) {
        phenotype = [this](Array<BranchDimensions>& nextBranches, float initialLength, const CoordinateFrame& initial, const Point3& branchEnd, int maxRecursionDepth, int currentRecursionDepth)
			{return LGenerator::bushTree(nextBranches, initialLength, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);};
	}else if (m_options.phenotypesIndexL == 3) {
        phenotype = [this](Array<BranchDimensions>& nextBranches, float initialLength, const CoordinateFrame& initial, const Point3& branchEnd, int maxRecursionDepth, int currentRecursionDepth)
			{return LGenerator::pineTree(nextBranches, initialLength, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);};
	}


    std::function<Vector3(float t)> branchCallback;
	if (m_options.branchCallbackIndexL == 0) {
		branchCallback = [this](float t)
		{ return LGenerator::straight(t); };
	} else if (m_options.branchCallbackIndexL == 1) {
		branchCallback = [this](float t)
		{ return LGenerator::curvy(t); };
	} else if (m_options.branchCallbackIndexL == 2) {
		branchCallback = [this](float t)
		{ return LGenerator::corkscrew(t); };
	}

	Random& rand = Random::threadCommon();
    int barkNum = rand.integer(0,3);
    String bark = "bark" + (String)(std::to_string(barkNum));

    shared_ptr<Tree> tree = genL.makeLTreeSkeleton(CoordinateFrame(), phenotype, [this](float t) {return LGenerator::straight(t);},m_options.initialHeightL, m_options.maxRecursionDepthL, m_options.maxRecursionDepthL);
    genL.skeletonToMeshL(treeMesh, leafMesh, tree, [this](float t) {return LGenerator::straight(t);}, [this](float t, shared_ptr<Tree> tree) {return LGenerator::branchRadius(t, tree);}, fruitLocations,  m_options.circlePtsL, m_options.branchSectionsL, m_options.initialHeightL, m_options.fall, bark);
    
    treeMesh.addMesh(leafMesh);
    treeMesh.toOBJ();
}

/**
	Generates a Scene.Any file that contains an orchard of trees
*/
void App::generateOrchard() {
	Array<Point3> fruitLocations = Array<Point3>();
    SCGenerator genSC;
	float height;

	if (m_options.typesIndex == 0) {
		makeLTree("firstTree", fruitLocations);
		height = m_options.initialHeightL;
	}
	else {
		shared_ptr<Tree> skeleton = genSC.makeSCTreeSkeleton(m_options.anchorCountSC, [this](float y) {return SCGenerator::bulbEnvelope(y);}, m_options.heightSC, m_options.radiusSC, m_options.killDistanceSC, m_options.treeDistanceSC, m_options.attractionRadiusSC, m_options.discountRateSC, Point3(0,0,0));
		genSC.skeletonToMeshSC(m_options.circlePointsSC, m_options.branchRadiusSC, m_options.radiusGrowthSC, m_options.leafinessSC, "firstTree", skeleton, fruitLocations);
		height = m_options.heightSC;
	}

	FruitDimensions fDims = m_options.fruitDims[m_options.fruitsIndex];

    TextOutput writer = TextOutput("scene/orchard.Scene.Any");
	Random& rand = Random::threadCommon();

    writer.printf("{");
    writer.writeNewline();
    writer.printf("name = \"Orchard\";");
    writer.writeNewlines(2);

    // models section
    writer.printf("models = {");
    writer.writeNewline();

    writer.printf("treeModel = ArticulatedModel::Specification {");
	writer.writeNewline();
	writer.printf("filename = \"tree.OBJ\"; };");
	writer.writeNewlines(2);

	writer.printf("fruitModel = ArticulatedModel::Specification {");
	writer.writeNewline();
	writer.printf("filename = \"models/" + m_options.fruits[m_options.fruitsIndex] + "/" + fDims.filename + "\";");
	writer.writeNewline();
	writer.printf("preprocess = { transformGeometry(all(), Matrix4::scale(%f) ); } };", fDims.scale );
	writer.writeNewline();

    writer.printf("};");
    writer.writeNewlines(2);

    // entities section
    writer.printf("entities = {");
    writer.writeNewline();

    writer.printf("skybox = Skybox { texture = \"cubemap/plainsky/null_plainsky512_*.jpg\"; };");
    writer.writeNewlines(2);
    writer.printf("light = Light {");
	writer.writeNewline();
	writer.printf("attenuation = (0, 0, 1); bulbPower = Power3(4e+06); castsShadows = true; frame = CFrame::fromXYZYPRDegrees(-15, 500, -41, -164, -77, 77); shadowMapSize = Vector2int16(2048, 2048); spotHalfAngleDegrees = 5; type = \"SPOT\"; };");
    writer.writeNewlines(2);
    writer.printf("camera = Camera {");
	writer.writeNewline();
	writer.printf("depthOfFieldSettings = DepthOfFieldSettings { enabled = true; farBlurRadiusFraction = 0.005; farBlurryPlaneZ = -100; farSharpPlaneZ = -40; focusPlaneZ = -10; lensRadius = 0.01; model = \"NONE\"; nearBlurRadiusFraction = 0.015; nearBlurryPlaneZ = -0.25; nearSharpPlaneZ = -1; };");
    writer.writeNewline();
	writer.printf("frame = CFrame::fromXYZYPRDegrees(0, 1, 4); };");

    for (int i = 0; i < m_options.numRows; ++i) {
		float xOffset = height * 3.0f * i;
		
		for (int j = 0; j < m_options.numTrees; ++j) {
			int zVar = rand.integer(0, 10);
			float zOffset = height * 3.0f * j;

			writer.writeNewlines(2);
			writer.printf("tree%d%d = VisibleEntity { model = \"treeModel\";", i, j);
			writer.writeNewline();
			writer.printf("frame = CFrame::fromXYZYPRDegrees(%f, 0, %f, %d, 0, 0); };", xOffset, zOffset, rand.integer(0,360));
			writer.writeNewlines(2);

			CoordinateFrame frame = CoordinateFrame(Point3(xOffset, 0.0f, zOffset));

			for (int k = 0; k < fruitLocations.length(); ++k) {
				int placeFruit = rand.integer(-2, 1);
				if (placeFruit == 1) {
					Point3 location = frame.pointToWorldSpace(fruitLocations[k]);
					writer.printf("fruit%d%d%d = VisibleEntity { model = \"fruitModel\";", i, j, k);
					writer.writeNewline();
					writer.printf("frame = CFrame::fromXYZYPRDegrees(%f, %f, %f, %d, %d, %f); };", location.x, location.y, location.z, rand.integer(0,360), 0, fDims.roll);
					writer.writeNewline();
				}
			}
		}
    }

    writer.printf("};");
    writer.writeNewlines(2);

    writer.printf("};");

    writer.commit();

}

/**
	Generates custom orchard scene
*/
void App::customOrchard() {	
	// options for L-system tree generation
	m_options.maxRecursionDepthL = 5;
	//m_options.initialHeightL = 1.0f;
	m_options.circlePtsL = 12;
	m_options.branchSectionsL = 10;
	m_options.phenotypesIndexL = 2;
	m_options.branchCallbackIndexL = 0;

	// options for size of orchard
	m_options.numRows = 4;
	m_options.numTrees = 3;

	// number of different tree models generated, and type of fruit
	int numPhenotypes = m_options.numRows * m_options.numTrees;
	int fruitIndex = 4;
	FruitDimensions fDims = m_options.fruitDims[fruitIndex];

    TextOutput writer = TextOutput("scene/customOrchard.Scene.Any");
	Random& rand = Random::threadCommon();

    writer.printf("{");
    writer.writeNewline();
    writer.printf("name = \"Custom Orchard\";");
    writer.writeNewlines(2);
    
	// models section
    writer.printf("models = {");
    writer.writeNewline();

	writer.printf("fruitModel = ArticulatedModel::Specification {");
	writer.writeNewline();
	writer.printf("filename = \"models/" + m_options.fruits[fruitIndex] + "/" + fDims.filename + "\";");
	writer.writeNewline();
	writer.printf("preprocess = { transformGeometry(all(), Matrix4::scale(%f) ); } };", fDims.scale);
	writer.writeNewline();
	

	// generate tree models
	//Array<Point3> fruitLocations = Array<Point3>();
    //SCGenerator genSC;
    //shared_ptr<Tree> skeleton = genSC.makeSCTreeSkeleton(m_options.anchorCountSC, [this](float y) {return SCGenerator::bulbEnvelope(y);}, m_options.heightSC, m_options.radiusSC, m_options.killDistanceSC, m_options.treeDistanceSC, m_options.attractionRadiusSC, m_options.discountRateSC, Point3(0,0,0));
    //genSC.skeletonToMeshSC(m_options.circlePointsSC, m_options.branchRadiusSC, m_options.radiusGrowthSC, m_options.leafinessSC, "tree0", skeleton, fruitLocations);

	Array<Array<Point3>> fruitLocations = Array<Array<Point3>>();
	for (int i = 0; i < numPhenotypes; ++i) {
		m_options.initialHeightL = 1.0f + rand.uniform(0.0f, 0.5f);
		fruitLocations.push(Array<Point3>());
		std::ostringstream strs;
		strs << "tree" << i;
		String filename = (String) strs.str();
		makeLTree(filename, fruitLocations[i]);
	
		writer.printf("treeModel%d = ArticulatedModel::Specification {", i);
		writer.writeNewline();
		writer.printf("filename = \"tree%d.OBJ\"; };", i);
		writer.writeNewlines(2);
	}

    writer.printf("};");
    writer.writeNewlines(2);

    // entities section
    writer.printf("entities = {");
    writer.writeNewline();

    writer.printf("skybox = Skybox { texture = \"cubemap/plainsky/null_plainsky512_*.jpg\"; };");
    writer.writeNewlines(2);
    writer.printf("light = Light {");
	writer.writeNewline();
	writer.printf("attenuation = (0, 0, 1); bulbPower = Power3(4e+06); castsShadows = true; frame = CFrame::fromXYZYPRDegrees(-15, 500, -41, -164, -77, 77); shadowMapSize = Vector2int16(2048, 2048); spotHalfAngleDegrees = 5; type = \"SPOT\"; };");
    writer.writeNewlines(2);
    writer.printf("camera = Camera {");
	writer.writeNewline();
	writer.printf("depthOfFieldSettings = DepthOfFieldSettings { enabled = true; farBlurRadiusFraction = 0.005; farBlurryPlaneZ = -100; farSharpPlaneZ = -40; focusPlaneZ = -10; lensRadius = 0.01; model = \"NONE\"; nearBlurRadiusFraction = 0.015; nearBlurryPlaneZ = -0.25; nearSharpPlaneZ = -1; };");
    writer.writeNewline();
	writer.printf("frame = CFrame::fromXYZYPRDegrees(0, 1, 4); };");

    for (int i = 0; i < m_options.numRows; ++i) {
		float xOffset = m_options.initialHeightL * 4.5f * i;
		
		for (int j = 0; j < m_options.numTrees; ++j) {
			int zVar = rand.integer(0, 10);
			float zOffset = m_options.initialHeightL * 3.5f * j;

			int phenotype = rand.integer(0, numPhenotypes - 1);
			writer.writeNewlines(2);
			writer.printf("tree%d%d = VisibleEntity { model = \"treeModel%d\";", i, j, phenotype);
			writer.writeNewline();
			writer.printf("frame = CFrame::fromXYZYPRDegrees(%f, 0, %f, 0, 0, 0); };", xOffset, zOffset);
			writer.writeNewlines(2);

			CoordinateFrame frame = CoordinateFrame(Point3(xOffset, 0.0f, zOffset));
			Array<Point3> currentFruits = fruitLocations[phenotype];

			for (int k = 0; k < currentFruits.length(); ++k) {
				int placeFruit = rand.integer(-2, 1);
				if (placeFruit == 1) {
					Point3 location = frame.pointToWorldSpace(currentFruits[k]);
					writer.printf("fruit%d%d%d = VisibleEntity { model = \"fruitModel\";", i, j, k);
					writer.writeNewline();
					writer.printf("frame = CFrame::fromXYZYPRDegrees(%f, %f, %f, %d, 0, %f); };", location.x, location.y, location.z, rand.integer(0,360), fDims.roll);
					writer.writeNewline();
				}
			}
		}
    }

    writer.printf("};");
    writer.writeNewlines(2);

    writer.printf("};");

    writer.commit();
};

void App::generateForest() {
    Array<Point3> fruitLocations = Array<Point3>();
    float density = 50;

	makeLTree("firstTree", fruitLocations);


    TextOutput writer = TextOutput("scene/forest.Scene.Any");
	Random& rand = Random::threadCommon();

    writer.printf("{");
    writer.writeNewline();
    writer.printf("name = \"Forest\";");
    writer.writeNewlines(2);

    // models section
    writer.printf("models = {");
    writer.writeNewline();

    for(int i = 0; i < 30; ++i){
        String filename = "tree" + (String)(std::to_string(i));
        m_options.maxRecursionDepthL = rand.integer(5, 7);

        if(m_options.maxRecursionDepthL == 5){
             m_options.initialHeightL = rand.uniform(0.5f, 1.0f);
        }else if(m_options.maxRecursionDepthL == 6){
             m_options.initialHeightL = rand.uniform(1.5f, 2.5f);
        }else{
            m_options.initialHeightL = rand.uniform(2.0f, 3.0f);
        }

       // m_options.phenotypesIndexL = rand.integer(1,2);
        makeLTree(filename, fruitLocations);
        writer.printf("treeModel%d = ArticulatedModel::Specification {", i);
	    writer.writeNewline();
	    writer.printf("filename = \"tree%d.OBJ\"; };", i);
	    writer.writeNewlines(2);
    }
    writer.printf("};");
    writer.writeNewlines(2);
    // entities section
    writer.printf("entities = {");
    writer.writeNewline();

    writer.printf("skybox = Skybox { texture = \"cubemap/plainsky/null_plainsky512_*.jpg\"; };");
    writer.writeNewlines(2);
    writer.printf("light = Light {");
	writer.writeNewline();
	writer.printf("attenuation = (0, 0, 1); bulbPower = Power3(4e+06); castsShadows = true; frame = CFrame::fromXYZYPRDegrees(-15, 500, -41, -164, -77, 77); shadowMapSize = Vector2int16(2048, 2048); spotHalfAngleDegrees = 5; type = \"SPOT\"; };");
    writer.writeNewlines(2);
    writer.printf("camera = Camera {");
	writer.writeNewline();
	writer.printf("depthOfFieldSettings = DepthOfFieldSettings { enabled = true; farBlurRadiusFraction = 0.005; farBlurryPlaneZ = -100; farSharpPlaneZ = -40; focusPlaneZ = -10; lensRadius = 0.01; model = \"NONE\"; nearBlurRadiusFraction = 0.015; nearBlurryPlaneZ = -0.25; nearSharpPlaneZ = -1; };");
    writer.writeNewline();
	writer.printf("frame = CFrame::fromXYZYPRDegrees(0, 1, 4); };");

    for (int i = 0; i < 100; ++i) {
		float xPos = rand.uniform(-(density/2.0f),(density/2.0f));				
		float zPos = rand.uniform(-(density/2.0f),(density/2.0f));

		writer.writeNewlines(2);
        int modelNumber = i%30;
		writer.printf("tree%d = VisibleEntity { model = \"treeModel%d\";", i, modelNumber);
		writer.writeNewline();
		writer.printf("frame = CFrame::fromXYZYPRDegrees(%f, %f, %f, 0, 0, 0); };", xPos, 0.0f, zPos);
		writer.writeNewline();
    }

    writer.printf("};");
    writer.writeNewlines(2);

    writer.printf("};");

    writer.commit();
}
