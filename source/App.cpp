/** \file App.cpp */
#include "App.h"
#include "Mesh.h"
#include "BranchDimensions.h"
//#include "Tree.h"
#include <cmath>
#include <map>
#include <tuple>
#include <stdlib.h> 

// Tells C++ to invoke command-line main() function even on OS X and Win32.
G3D_START_AT_MAIN();

int main(int argc, const char* argv[]) {
    {
        G3DSpecification g3dSpec;
        g3dSpec.audio = false;
        initGLG3D(g3dSpec);
    }

    GApp::Settings settings(argc, argv);

    // Change the window and other startup parameters by modifying the
    // settings class.  For example:
    settings.window.caption             = argv[0];

    // Set enable to catch more OpenGL errors
    // settings.window.debugContext     = true;

    // Some common resolutions:
    // settings.window.width            =  854; settings.window.height       = 480;
    // settings.window.width            = 1024; settings.window.height       = 768;
    settings.window.width               = 1280; settings.window.height       = 720;
    //settings.window.width             = 1920; settings.window.height       = 1080;
    // settings.window.width            = OSWindow::primaryDisplayWindowSize().x; settings.window.height = OSWindow::primaryDisplayWindowSize().y;
    settings.window.fullScreen          = false;
    settings.window.resizable           = ! settings.window.fullScreen;
    settings.window.framed              = ! settings.window.fullScreen;

    // Set to true for a significant performance boost if your app can't render at 60fps, or if
    // you *want* to render faster than the display.
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

    // Call setScene(shared_ptr<Scene>()) or setScene(MyScene::create()) to replace
    // the default scene here.
    
    showRenderingStats      = false;

    makeGUI();
    // For higher-quality screenshots:
    // developerWindow->videoRecordDialog->setScreenShotFormat("PNG");
    // developerWindow->videoRecordDialog->setCaptureGui(false);
    developerWindow->cameraControlWindow->moveTo(Point2(developerWindow->cameraControlWindow->rect().x0(), 0));
    makeTree();

    //makeTreeSkeleton(5, [this](float y) {return App::envelopePerimeter(y);}, 10.0f, 10.0f, 2.0f, 0.5f, Point3(0,0,0));

    loadScene(
        //"G3D Sponza"
        "Tree Testing" // Load something simple
        //developerWindow->sceneEditorWindow->selectedSceneName()  // Load the first scene encountered 
        );
}


void App::makeGUI() {
    // Initialize the developer HUD
    createDeveloperHUD();

    debugWindow->setVisible(true);
    developerWindow->videoRecordDialog->setEnabled(true);

	// Tree generation GUI
    GuiPane* treePane = debugPane->addPane("Tree");
    treePane->setNewChildSize(500, -1, 300);
	treePane->addNumberBox("Recursion depth:", &m_maxRecursionDepth, "", GuiTheme::LINEAR_SLIDER, 2, 10);
	treePane->addNumberBox("Initial height:", &m_initialHeight, "", GuiTheme::LOG_SLIDER, 0.5f, 10.0f);
	treePane->addNumberBox("Circle points:", &m_circlePts, "", GuiTheme::LOG_SLIDER, 3, 100);
	treePane->addNumberBox("Branch sections:", &m_branchSections, "", GuiTheme::LOG_SLIDER, 1, 100);
    treePane->addDropDownList("Phenotype", m_phenotypes, &m_phenotypesIndex);
    treePane->addButton("Generate tree", [this](){
		drawMessage("Generating tree...");
        makeTree();
		ArticulatedModel::clearCache();
		GApp::loadScene("Tree Testing");
	});
	treePane->addButton("Generate orchard", [this]() {
		drawMessage("Generating orchard...");
		generateOrchard();
		ArticulatedModel::clearCache();
		GApp::loadScene("Orchard");
	});
	treePane->pack();

    debugWindow->pack();
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}


void App::makeTree() {
    Mesh tree = Mesh("tree");
    Mesh leafMesh = Mesh("leaf");

	if (m_phenotypesIndex == 0) {
		makeBranch(tree, leafMesh, CoordinateFrame() * CoordinateFrame::fromXYZYPRDegrees(0,0,0,0,0,0), m_initialHeight, [this](float t) {return App::spineCurve(t);}, [this](float t, int branchFactor, int depth) {return App::branchRadius(t, branchFactor, depth);},
			[this](Array<BranchDimensions>& nextBranches, float initialLength, const CoordinateFrame& initial, const Point3& branchEnd, int maxRecursionDepth, int currentRecursionDepth)
			{return App::normalTree(nextBranches, initialLength, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);},
			m_maxRecursionDepth, m_maxRecursionDepth, m_circlePts, m_branchSections);
	}
	if (m_phenotypesIndex == 1) {
		makeBranch(tree, leafMesh, CoordinateFrame() * CoordinateFrame::fromXYZYPRDegrees(0,0,0,0,0,0), m_initialHeight, [this](float t) {return App::spineCurve(t);}, [this](float t, int branchFactor, int depth) {return App::branchRadius(t, branchFactor, depth);},
			[this](Array<BranchDimensions>& nextBranches, float initialLength, const CoordinateFrame& initial, const Point3& branchEnd, int maxRecursionDepth, int currentRecursionDepth)
			{return App::randomTree(nextBranches, initialLength, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);},
			m_maxRecursionDepth, m_maxRecursionDepth, m_circlePts, m_branchSections);
	}
    
    tree.addMesh(leafMesh);
    tree.toOBJ();
}


void App::makeBranch(Mesh& mesh, Mesh& leafMesh, const CoordinateFrame& initial, float& length, std::function<Vector3(float)> spineCurve, std::function<float(float, int, int)> branchRadius,
    std::function<void(Array<BranchDimensions>&, float, const CoordinateFrame&, Point3&, int, int)> phenotype, int maxRecursionDepth, int currentRecursionDepth, int circlePoints, int branchSections) const {
    if (currentRecursionDepth != 0) {
        int index = mesh.numVertices();
	    //float sectionHeight;
	    float sectionRadius;
        float distanceAlongBranch = 1.0f; //Ranges from 0.0f to 1.0f
        Point3 branchEnd = initial.pointToWorldSpace(Point3(0,length*(19.0f/20.0f),0));
        Point3 branchMid = initial.pointToWorldSpace(Point3(0,length*(6.0f/10.0f),0));


        // callback function to decide how to recurse, populates an array of BranchDimensions which contain coordinate frames and lengths for the next branches
        Array<BranchDimensions> nextBranches;
        phenotype(nextBranches, length, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);
        debugAssert(nextBranches.size() > 0);


        // Add vertices of intial circle at the bottom of the branch we are currently making to mesh
	    for(int i = 0; i < circlePoints; ++i) {
	    	float angle = (i * 2.0f * pif()) / circlePoints;
            sectionRadius = branchRadius(distanceAlongBranch, nextBranches.size(), currentRecursionDepth);
            Vector3 vec = Vector3(cos(angle) * sectionRadius, 0.0f, sin(angle) * sectionRadius);
            vec = initial.pointToWorldSpace(vec);
	    	mesh.addVertex(vec.x, vec.y, vec.z);  
	    }
	    
        // Add vertices of circles on top of the initial circle to mesh
	    for(int i = 1; i <= branchSections; ++i) {
            distanceAlongBranch =  (float)(i) / float(branchSections);
	    	sectionRadius = branchRadius(distanceAlongBranch, nextBranches.size(), currentRecursionDepth);
	    	// TODO:: pass a coordinate frame that is returned by space curve function (instead of initial)
            CoordinateFrame section = initial;
            section.translation = initial.pointToWorldSpace(length * spineCurve(distanceAlongBranch));
            addCylindricSection(mesh, circlePoints, section, sectionRadius);
	    }

        //Add the leaves
        if (currentRecursionDepth == 1){
            //Add one leaf on the end of the branch
            CoordinateFrame leaf = initial;
            leaf.translation = branchEnd;
            float leafSize = 0.15f;
            addLeaves(leafMesh, leafSize, leaf);
            
            //Add random leaves along branch
            int leafNumber = 5;
            float rPitch;
            float rYaw;
            float rDisplacement;
            for(int i = 0; i < 5; ++i){
                leaf = initial;
                rDisplacement = ( (float)(rand() % 100) / 100.0f ) * length;
                rYaw = rand() % 360;
                rPitch = rand()%40 + 30.0f;  
                leaf.translation = initial.pointToWorldSpace(Point3(0.0f, rDisplacement, 0.0f));
                leaf = leaf * CoordinateFrame::fromXYZYPRDegrees(0.0f,0.0f,0.0f,rYaw);
                leaf = leaf * CoordinateFrame::fromXYZYPRDegrees(0.0f,0.0f,0.0f,0.0f,rPitch);
                addLeaves(leafMesh, leafSize, leaf);
            }
        }
        else {
            for (int i = 0; i < nextBranches.length(); ++i) {
                BranchDimensions nextBranch = nextBranches[i];
                CoordinateFrame branch = nextBranch.frame;
                float newLength = nextBranch.length;
                
                makeBranch(mesh, leafMesh, branch, newLength, spineCurve, branchRadius, phenotype, maxRecursionDepth, currentRecursionDepth - 1, circlePoints, branchSections);   
            }
        }
    }
}


void App::addLeaves(Mesh& leafMesh, float& leafSize, const CoordinateFrame& leafFrame) const{
    int index = leafMesh.numVertices();
    //float leafSize = length*1.5f;
    Vector3 vec1 = Vector3(leafSize / 2.0f, leafSize, 0.0f);
    vec1 = leafFrame.pointToWorldSpace(vec1);
    Vector3 vec2 = Vector3(-leafSize / 2.0f, leafSize, 0.0f);
    vec2 = leafFrame.pointToWorldSpace(vec2);
    Vector3 vec3 = Vector3(0.0f, 0.0f, 0.0f);
    vec3 = leafFrame.pointToWorldSpace(vec3);
    leafMesh.addVertex(vec1);
    leafMesh.addVertex(vec2);
    leafMesh.addVertex(vec3);
    leafMesh.addFace(index, index+1, index+2, 4, 3, 5, "leaf");
    leafMesh.addFace(index+2, index+1, index, 5, 3, 4, "leaf");
}


void App::addCylindricSection(Mesh& mesh, const int& pts, const CoordinateFrame& origin, const float& radius) const {
	int index = mesh.numVertices();

    for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
        Vector3 vec = Vector3(cos(angle) * radius, 0.0f, sin(angle) * radius);
        vec = origin.pointToWorldSpace(vec);

		mesh.addVertex(vec);
	}
	int offset = index - pts;
	for(int i = 0; i < pts; ++i) {
		mesh.addFace( offset + i, offset + i + (pts),offset + ((i + 1) % pts), 2, 3, 1, "bark");
		mesh.addFace(offset + i + (pts), offset + ((i + 1) % pts) + (pts), offset + ((i + 1) % pts), 2, 4, 3, "bark");
	}
}


/**
	Callback functions for the curve of the tree
*/
Vector3 App::spineCurve(float t) {
    return Vector3(0.0f, t, 0.0f);
}


/**
	Callback functions for the radii of the branches
*/
float App::branchRadius(float t, int branchingNumber, int recursionDepth) {
    if (recursionDepth == 0){
        return 0.005f;
    }
    float end = branchRadius(0.0f, branchingNumber, recursionDepth - 1);
    float base = branchingNumber * (end*end);
    base = sqrt(base);
    float diff = base - end;
    
    end += diff/4.0f;
    //base -= (diff/10.0f);

    debugAssert(t >= 0.0f && t <= 1.0f);
    //a, b, and c from the quadratic eq. ax^2 + bx + c = y
    float a = (base - end);
    float b = 2.0f*(end - base);
    float c = base;

    return (t*t*a) + (b*t) + c;
}


/**
	Callback functions for the phenotype of the tree
*/
void App::randomTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth) {  

    float newBranchLength = 3.0f * initialLength / 5.0f;
    float rPitch = rand() % 40 + 15; //This will give a random pitch within the range 25 - 45
    float rYaw;
    CoordinateFrame branch1 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, rPitch, 0.0f);
    branch1.translation = branchEnd;
    
    rPitch = rand() % 40 + 15; 
    rYaw = rand() % 40 + 70; //This will give a random yaw fom 70 to 110
    CoordinateFrame branch2 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, rYaw, 0.0f, 0.0f);
    branch2 = branch2 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, rPitch, 0.0f);
    branch2.translation = branchEnd;
    
    rPitch = rand() % 40 + 15; 
    rYaw = rand() % 40 + 70;
    CoordinateFrame branch3 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, rYaw, 0.0f, 0.0f);
    branch3 = branch3 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -rPitch, 0.0f);
    branch3.translation = branchEnd;
    
    rPitch = rand() % 40 + 15; 
    CoordinateFrame branch4 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -rPitch, 0.0f);
    branch4.translation = branchEnd;

    
    //These random values will determine how many child branches will be added (Fix realistic Radii with random child branches)
    bool makeBranch1 = true;//(rand() % 100) > 30;
    bool makeBranch2 = true;//(rand() % 100) > 30;
    bool makeBranch3 = true;//(rand() % 100) > 30;
    bool makeBranch4 = true;//(rand() % 100) > 30;

    if(makeBranch1){
        nextBranches.push(BranchDimensions(branch1, newBranchLength));
    }
    if(makeBranch2){
        nextBranches.push(BranchDimensions(branch2, newBranchLength));
    }
    if(makeBranch3){
        nextBranches.push(BranchDimensions(branch3, newBranchLength));
    }
    if(makeBranch4){
        nextBranches.push(BranchDimensions(branch4, newBranchLength));
    }
    debugAssert(nextBranches.size() > 0);
}

void App::normalTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth) {  

    float newBranchLength = 3.0f * initialLength / 5.0f;

    CoordinateFrame branch1 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, 35.0f, 0.0f);
    branch1.translation = branchEnd;

    CoordinateFrame branch2 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -35.0f, 0.0f);
    branch2.translation = branchEnd;

    CoordinateFrame branch3 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 90.0f, 0.0f, 0.0f);
    branch3 = branch3 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, 35.0f, 0.0f);
    branch3.translation = branchEnd;
    
    CoordinateFrame branch4 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 90.0f, 0.0f, 0.0f);
    branch4 = branch4 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -35.0f, 0.0f);
    branch4.translation = branchEnd;
    
    

    nextBranches.push(BranchDimensions(branch1, newBranchLength));
    nextBranches.push(BranchDimensions(branch2, newBranchLength));
    nextBranches.push(BranchDimensions(branch3, newBranchLength));
    nextBranches.push(BranchDimensions(branch4, newBranchLength));
 
    debugAssert(nextBranches.size() > 0);
}


/**
	Generates a Scene.Any file that contains an orchard of trees
*/
void App::generateOrchard() {
	makeTree();

    TextOutput writer = TextOutput("scene/orchard.Scene.Any");

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
	writer.writeNewline();
	writer.writeNewlines(2);

    // continuing the entities section, use a for loop to write in 3 shelves on the bookcase
    for (int i = 0; i < 2; ++i) {
		float xOffset = 3.0 * i;
		for (int j = 0; j < 10; ++j) {
			float zOffset = 2.0 * j;
			writer.printf("tree%d%d = VisibleEntity { model = \"treeModel\";", i, j);
			writer.writeNewline();
			writer.printf("frame = CFrame::fromXYZYPRDegrees(%f, 0, %f); };", xOffset, zOffset);
			writer.writeNewlines(2);
		}
    };

    writer.printf("};");
    writer.writeNewlines(2);

    writer.printf("};");

    writer.commit();
}


float App::envelopePerimeter(float y) {
    if(y < 0.5f) {
        return 0.0f;
    } else {
        return 2.0f - (2.0f * y);
    }
}


//Tree App::makeTreeSkeleton(int anchorPoints, std::function<float(float)> envelopePerimeter, float height, float radius, float killDistance, float nodeDistance, Point3 initTreeNode) {
//
//    Tree result = Tree(initTreeNode);
//
//    Tree testNode = Tree(initTreeNode);
//
//    Array<Point3> anchorPointsList = generateAnchorPoints(anchorPoints, height, radius, envelopePerimeter);
//    
//    while(anchorPointsList.size() > 0) {
//        Array<Tree*> selectedNodes = Array<Tree*>();
//        
//        //For each anchor, select the closest tree node.
//        for(int i = 0; i < anchorPointsList.size(); ++i) {
//            Point3 currentAnchor = anchorPointsList[i];
//            Array<Tree*> stack = Array<Tree*>(&result);
//            Tree* closestNode = &result;
//
//            while(stack.size() > 0) {
//                Tree* challenger = stack.pop();
//                if((challenger->getContents() - currentAnchor).magnitude() < (closestNode->getContents() - currentAnchor).magnitude()) {
//                    closestNode = challenger;
//                }
//                Array<Tree*>* tempArray = challenger->getChildren();
//                stack.append(*challenger->getChildren());
//            }
//            selectedNodes.push(closestNode);
//        }
//
//        //Assign the directions for new tree nodes.
//        std::map<Tree*, Vector3> nodeDirections;
//        for(int i = 0; i < anchorPointsList.size(); ++i) {
//            Tree* tempNode = selectedNodes[i];
//            if(nodeDirections.find(tempNode) == nodeDirections.end()) {
//                nodeDirections[tempNode] = (anchorPointsList[i] - tempNode->getContents()).direction();
//            } else {
//                nodeDirections[tempNode] = (nodeDirections[tempNode] + (anchorPointsList[i] - tempNode->getContents()).direction()).direction();
//            }
//        }
//
//        //Iterate over the selected nodes and spawn new tree nodes.
//        for(auto it = nodeDirections.begin(); it != nodeDirections.end(); ++it) {
//            Tree* parent = it->first;
//            Point3 pos1 = parent->getContents();
//            Point3 pos2 = nodeDistance * nodeDirections[parent];
//            Tree child = Tree(pos1 + pos2, parent);
//            Array<Tree*>* children = parent->getChildren();
//            children->push(&child);
//        }
//
//        //Kill anchors too close to the tree.
//        for(int i = 0; i < anchorPointsList.size(); ++i) {
//            Point3 currentAnchor = anchorPointsList[i];
//
//            Array<Tree*> stack = Array<Tree*>(&result);
//
//            while(stack.size() > 0) {
//                Tree* challenger = stack.pop();
//                if((challenger->getContents() - currentAnchor).magnitude() < killDistance * nodeDistance) {
//                    anchorPointsList.remove(i);
//                    break;
//                }
//                stack.append(*challenger->getChildren());
//            }
//        }
//    }
//    return result;
//}


//Array<Point3> App::generateAnchorPoints(int count, float height, float radius, std::function<float(float)> radiusCurve) {
//    Random& rng = Random::threadCommon();
//    Array<Point3> result = Array<Point3>();
//    for(int i = 0; i < count; ++i) {
//        Point3 newPoint;
//        do {
//            newPoint = Point3(rng.uniform(), rng.uniform(), rng.uniform());
//        } while(newPoint.xz().length() >= radiusCurve(newPoint.y));
//        result.push(Point3(newPoint.x * radius, newPoint.y * height, newPoint.z * radius));
//    }
//    return result;
//}
//
//void App::skeletonToMesh(float initRadius, float radiusGrowth, String filename, Tree skeleton) {
//    
//}