/** \file App.cpp */
#include "App.h"
#include "Mesh.h"
#include "BranchDimensions.h"
#include "Tree.h"
#include "FruitDimensions.h"
#include <cmath>
#include <map>
#include <tuple>
#include <stdlib.h> 
#include "SCGenerator.h"
#include "LGenerator.h"

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
	FruitDimensions pearDims = FruitDimensions("pear_export.obj", m_options.initialHeightL / 50.0f, 0.04f, 0.0f);
	FruitDimensions bananaDims = FruitDimensions("banana.obj", m_options.initialHeightL / 2000.0f, 0.04f, -90.0f);
	FruitDimensions moneyDims = FruitDimensions("Dollar stack wild.obj", m_options.initialHeightL / 100.0f, 0.05f, 90.0f);
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
        treeLPane->addButton("Generate tree", [this](){
	    	drawMessage("Generating tree...");
	    	Array<Point3> fruitLocations = Array<Point3>();
            
	        Stopwatch sw;
    	    sw.tick(); //start the timer

            makeLTree("tree", fruitLocations);
            
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
            
	        Stopwatch sw;
    	    sw.tick(); //start the timer

            shared_ptr<Tree> skeleton = genSC.makeSCTreeSkeleton(m_options.anchorCountSC, [this](float y) {return SCGenerator::bulbEnvelope(y);}, m_options.heightSC, m_options.radiusSC, m_options.killDistanceSC, m_options.treeDistanceSC, m_options.attractionRadiusSC, m_options.discountRateSC, Point3(0,0,0));
            genSC.skeletonToMeshSC(m_options.circlePointsSC, m_options.branchRadiusSC, m_options.radiusGrowthSC, m_options.leafinessSC, "tree", skeleton, fruitLocations);

            sw.tock();
            debugPrintf("Elapsed Time: %f\n", sw.elapsedTime());
	    	ArticulatedModel::clearCache();
	    	GApp::loadScene("Tree Testing");
	    });
	    treeSCPane->pack();

		GuiPane* orchardPane = containerPane->addTab("Orchard");
		orchardPane->addNumberBox("Number of rows:", &m_options.numRows, "", GuiTheme::LOG_SLIDER, 3, 20);
	    orchardPane->addNumberBox("Trees per row:", &m_options.numTrees, "", GuiTheme::LOG_SLIDER, 1, 20);
		orchardPane->addDropDownList("Generation type", m_options.types, &m_options.typesIndex);
        orchardPane->addDropDownList("Fruit type", m_options.fruits, &m_options.fruitsIndex);
		orchardPane->addButton("Generate orchard", [this]() {
	    	drawMessage("Generating orchard...");
	    	generateOrchard();
	    	ArticulatedModel::clearCache();
	    	GApp::loadScene("Orchard");
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

    shared_ptr<Tree> tree = genL.makeLTreeSkeleton(CoordinateFrame(), phenotype, [this](float t) {return LGenerator::curvy(t);},m_options.initialHeightL, m_options.maxRecursionDepthL, m_options.maxRecursionDepthL);
    genL.skeletonToMeshL(treeMesh, leafMesh, tree, [this](float t) {return LGenerator::curvy(t);}, [this](float t, shared_ptr<Tree> tree) {return LGenerator::branchRadius(t, tree);},
			fruitLocations,  m_options.circlePtsL, m_options.branchSectionsL, m_options.initialHeightL);
    
    treeMesh.addMesh(leafMesh);
    treeMesh.toOBJ();
}

/**
	Generates a Scene.Any file that contains an orchard of trees
*/
void App::generateOrchard() {
	Array<Point3> fruitLocations = Array<Point3>();
    SCGenerator genSC;

	if (m_options.typesIndex == 0) {
		makeLTree("firstTree", fruitLocations);
	}
	else {
		shared_ptr<Tree> skeleton = genSC.makeSCTreeSkeleton(m_options.anchorCountSC, [this](float y) {return SCGenerator::bulbEnvelope(y);}, m_options.heightSC, m_options.radiusSC, m_options.killDistanceSC, m_options.treeDistanceSC, m_options.attractionRadiusSC, m_options.discountRateSC, Point3(0,0,0));
		genSC.skeletonToMeshSC(m_options.circlePointsSC, m_options.branchRadiusSC, m_options.radiusGrowthSC, m_options.leafinessSC, "firstTree", skeleton, fruitLocations);
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
		float xOffset = m_options.maxRecursionDepthL * i;
		
		for (int j = 0; j < m_options.numTrees; ++j) {
			int zVar = rand.integer(0, 10);
			float zOffset = m_options.maxRecursionDepthL * j;

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