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

    loadScene(
        //"G3D Sponza"
        "Tree Testing" // Load something simple
        //developerWindow->sceneEditorWindow->selectedSceneName()  // Load the first scene encountered 
        );


	// initialize fruitDims
	fruitDims = Array<FruitDimensions>();
	//FruitDimensions appleDims = FruitDimensions("apple textured obj.obj", 0.001f, 0.08f);
	//FruitDimensions moneyDims = FruitDimensions("Dollar stack wild.obj", 0.01f, 0.05f);
	//FruitDimensions teapotDims = FruitDimensions("glassTeapot.obj", 0.001f, 0.08f);
	
	FruitDimensions appleDims = FruitDimensions("apple textured obj.obj", m_initialHeight / 1000.0f, 0.08f);
	FruitDimensions moneyDims = FruitDimensions("Dollar stack wild.obj", m_initialHeight / 100.0f, 0.05f);
	FruitDimensions teapotDims = FruitDimensions("glassTeapot.obj", m_initialHeight / 1000.0f, 0.08f);


	fruitDims.push(appleDims);
	fruitDims.push(moneyDims);
	fruitDims.push(teapotDims);
}


void App::makeGUI() {
    // Initialize the developer HUD
    createDeveloperHUD();

    debugWindow->setVisible(true);
    developerWindow->videoRecordDialog->setEnabled(true);

    debugPane->beginRow(); {


        GuiTabPane* containerPane = debugPane->addTabPane();
	    
        // L-System tree generation GUI
        GuiPane* treePane = containerPane->addTab("L-System Tree");
        treePane->setNewChildSize(500, -1, 300);
	    treePane->addNumberBox("Recursion depth:", &m_maxRecursionDepth, "", GuiTheme::LINEAR_SLIDER, 2, 10);
	    treePane->addNumberBox("Initial height:", &m_initialHeight, "", GuiTheme::LOG_SLIDER, 0.5f, 10.0f);
	    treePane->addNumberBox("Circle points:", &m_circlePts, "", GuiTheme::LOG_SLIDER, 3, 100);
	    treePane->addNumberBox("Branch sections:", &m_branchSections, "", GuiTheme::LOG_SLIDER, 1, 100);
        treePane->addDropDownList("Phenotype", m_phenotypes, &m_phenotypesIndex);
        treePane->addButton("Generate tree", [this](){
	    	drawMessage("Generating tree...");
	    	Array<Point3> fruitLocations = Array<Point3>();
            makeLTree("tree", fruitLocations);
	    	ArticulatedModel::clearCache();
	    	GApp::loadScene("Tree Testing");
	    });
		treePane->pack();

        // Space tree generation GUI
        GuiPane* spaceTreePane = containerPane->addTab("Space Col Tree");
        spaceTreePane->setNewChildSize(500, -1, 300);
	    spaceTreePane->addNumberBox("Anchor Points:", &m_spaceAnchorCount, "", GuiTheme::LOG_SLIDER, 1, 10000);
	    spaceTreePane->addNumberBox("Height:", &m_spaceHeight, "", GuiTheme::LINEAR_SLIDER, 10.0f, 100.0f);
	    spaceTreePane->addNumberBox("Radius:", &m_spaceRadius, "", GuiTheme::LINEAR_SLIDER, 10.0f, 100.0f);
	    spaceTreePane->addNumberBox("Circle Points:", &m_spaceCirclePoints, "", GuiTheme::LOG_SLIDER, 1, 100);
	    spaceTreePane->addNumberBox("Tree Node Distance:", &m_spaceTreeDistance, "", GuiTheme::LOG_SLIDER, 0.01f, 1.0f);
	    spaceTreePane->addNumberBox("Kill Distance:", &m_spaceKillDistance, "", GuiTheme::LINEAR_SLIDER, 1.0f, 10.0f);
	    spaceTreePane->addNumberBox("Branch Initial Radius:", &m_spaceBranchRadius, "", GuiTheme::LOG_SLIDER, 0.01f, 10.0f);
	    spaceTreePane->addNumberBox("Radius Growth:", &m_spaceRadiusGrowth, "", GuiTheme::LINEAR_SLIDER, 2.0f, 3.0f);
	    spaceTreePane->addNumberBox("Attraction Radius:", &m_spaceAttractionRadius, "", GuiTheme::LOG_SLIDER, 1.0f, 100.0f);
        spaceTreePane->addButton("Generate tree", [this](){
	    	drawMessage("Generating tree...");

			shared_ptr<Tree> skeleton = makeSCTreeSkeleton(m_spaceAnchorCount, [this](float y) {return App::envelopePerimeter(y);}, m_spaceHeight, m_spaceRadius, m_spaceKillDistance, m_spaceTreeDistance, m_spaceAttractionRadius, Point3(0,0,0));
			Array<Point3> fruitLocations;
			skeletonToMesh(m_spaceCirclePoints, m_spaceBranchRadius, m_spaceRadiusGrowth, "tree", skeleton, fruitLocations);

	    	ArticulatedModel::clearCache();
	    	GApp::loadScene("Tree Testing");
	    });
	    spaceTreePane->pack();

		GuiPane* orchardPane = containerPane->addTab("Orchard");
		orchardPane->addNumberBox("Number of rows:", &m_numRows, "", GuiTheme::LOG_SLIDER, 3, 20);
	    orchardPane->addNumberBox("Trees per row:", &m_numTrees, "", GuiTheme::LOG_SLIDER, 1, 20);
		orchardPane->addDropDownList("Generation type", m_types, &m_typesIndex);
        orchardPane->addDropDownList("Fruit type", m_fruits, &m_fruitsIndex);
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

shared_ptr<Tree> App::makeLTreeSkeleton(const CoordinateFrame& initial, std::function<void(Array<BranchDimensions>&, float, const CoordinateFrame&, Point3&, int, int)> phenotype, std::function<Vector3(float)> spineCurve, float length, int maxRecursionDepth, int currentRecursionDepth, shared_ptr<Tree> parent){
    shared_ptr<Tree> tree = Tree::create(initial, parent);
    Point3 branchEnd = initial.pointToWorldSpace(length * spineCurve(19.0f / 20.0f));

    // callback function to decide how to recurse, populates an array of BranchDimensions which contain coordinate frames and lengths for the next branches
    Array<BranchDimensions> nextBranches;
    phenotype(nextBranches, length, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);

    if (currentRecursionDepth > 0) {
        for (int i = 0; i < nextBranches.length(); ++i) {
            BranchDimensions nextBranch = nextBranches[i];
            CoordinateFrame branch = nextBranch.frame;
            float newLength = nextBranch.length;
            
            shared_ptr<Tree> childTree = makeLTreeSkeleton(branch, phenotype, spineCurve, newLength, maxRecursionDepth, currentRecursionDepth - 1, tree);
            if(notNull(childTree)){
                tree->addChild(childTree);
            }
        }
        return tree;
    } else {
        return nullptr;
    }
}

void App::makeLTree(String filename, Array<Point3>& fruitLocations) {
    Mesh mesh = Mesh(filename);
    Mesh leafMesh = Mesh("leaf");
     std::function<void(Array<BranchDimensions>&, float, const CoordinateFrame&, Point3&, int, int)> phenotype;
	if (m_phenotypesIndex == 0) {
		phenotype = [this](Array<BranchDimensions>& nextBranches, float initialLength, const CoordinateFrame& initial, const Point3& branchEnd, int maxRecursionDepth, int currentRecursionDepth)
			{return App::normalTree(nextBranches, initialLength, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);};
	}
	if (m_phenotypesIndex == 1) {
        phenotype = [this](Array<BranchDimensions>& nextBranches, float initialLength, const CoordinateFrame& initial, const Point3& branchEnd, int maxRecursionDepth, int currentRecursionDepth)
			{return App::randomTree(nextBranches, initialLength, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);};
	}

    shared_ptr<Tree> tree = makeLTreeSkeleton(CoordinateFrame(), phenotype, [this](float t) {return App::straight(t);},m_initialHeight, m_maxRecursionDepth, m_maxRecursionDepth);
    buildTree(mesh, leafMesh, tree, [this](float t) {return App::straight(t);}, [this](float t, shared_ptr<Tree> tree) {return App::branchRadius(t, tree);},
			fruitLocations,  m_circlePts, m_branchSections, m_initialHeight);


    mesh.addMesh(leafMesh);
    mesh.toOBJ();
}

void App::buildTree(Mesh& mesh, Mesh& leafMesh, const shared_ptr<Tree> tree, std::function<Vector3(float)> spineCurve, std::function<float(float, shared_ptr<Tree>)> branchRadius, Array<Point3>& fruitLocations, int circlePoints, int branchSections, float initialLength){
    float distanceAlongBranch = 0.0f;
    shared_ptr<Array<shared_ptr<Tree>>> children = tree->getChildren();
    float sectionRadius;
    CoordinateFrame initial = tree->getContents();
    
    // Add vertices of intial circle at the bottom of the branch we are currently making to mesh
	    for(int i = 0; i < circlePoints; ++i) {
	    	float angle = (i * 2.0f * pif()) / circlePoints;
            sectionRadius = branchRadius(distanceAlongBranch, tree);
            Vector3 vec = Vector3(cos(angle) * sectionRadius, 0.0f, sin(angle) * sectionRadius);
            vec = initial.pointToWorldSpace(vec);
	    	mesh.addVertex(vec.x, vec.y, vec.z);  
	    }

         // Add vertices of circles on top of the initial circle to mesh
	    for(int i = 1; i <= branchSections; ++i) {
            distanceAlongBranch =  (float)(i) / float(branchSections);
	    	sectionRadius = branchRadius(distanceAlongBranch, tree);
	    	// TODO:: pass a coordinate frame that is returned by space curve function (instead of initial)
            CoordinateFrame section = initial;
            section.translation = initial.pointToWorldSpace(initialLength * spineCurve(distanceAlongBranch));
            addCylindricSection(mesh, circlePoints, section, sectionRadius);
	    }

        if(children->size() == 0){
            addLeaves(initial, initialLength, leafMesh, fruitLocations);
            
        }else{
            float newLength = ((3.0f/5.0f)*initialLength);
                
            for(int i = 0; i < children->size(); ++i){
                buildTree(mesh, leafMesh, children->operator[](i), spineCurve, branchRadius, fruitLocations, circlePoints, branchSections, newLength);
            }
        
        }

}

void App::addLeaves(CoordinateFrame& initial, float length, Mesh& leafMesh, Array<Point3>& fruitLocations){
    Random& rand = Random::threadCommon();
    //Add one leaf on the end of the branch
    CoordinateFrame leafFrame = initial;
    leafFrame.translation = initial.pointToWorldSpace(Point3(0.0f, (19.0f/20.0f)*length, 0.0f));
    float leafSize = 0.15f;
    addLeaf(leafMesh, leafSize, leafFrame);
	addFruits(fruitLocations, leafFrame);

    
    //Add random leaves along branch
    int leafNumber = 5;
    float rPitch;
    float rYaw;
    float rDisplacement;
    for(int i = 0; i < 5; ++i){
        leafFrame = initial;
        rDisplacement = ( (float)rand.integer(0,100) / 100.0f ) * length;
        rYaw = rand.integer(0, 360);
        rPitch = rand.integer(0, 40) + 30.0f;  
        leafFrame.translation = initial.pointToWorldSpace(Point3(0.0f, rDisplacement, 0.0f));
        leafFrame = leafFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, rYaw);
        leafFrame = leafFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, rPitch);
        addLeaf(leafMesh, leafSize, leafFrame);
    }

}
void App::addLeaf(Mesh& leafMesh, float& leafSize, const CoordinateFrame& leafFrame) const {
    int index = leafMesh.numVertices();
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


void App::addFruits(Array<Point3>& fruitLocations, const CoordinateFrame& fruitFrame) {
	fruitLocations.push( fruitFrame.pointToWorldSpace(Vector3(0.0f, 0.0f, 0.0f)) );
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

//Version that takes a hard coded index
void App::addCylindricSection(Mesh& mesh, const int& parentIndex, const int& currentIndex, const int& pts, const CoordinateFrame& origin, const float& radius) const {

    for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
        Vector3 vec = Vector3(cos(angle) * radius, 0.0f, sin(angle) * radius);
        vec = origin.pointToWorldSpace(vec);

		mesh.addVertex(vec);
	}

	for(int i = 0; i < pts; ++i) {
		mesh.addFace(parentIndex + ((i + 1) % pts), parentIndex + i, currentIndex + i, 2, 3, 1, "bark");
		mesh.addFace(currentIndex + ((i + 1) % pts), parentIndex + ((i + 1) % pts), currentIndex + i, 2, 4, 3, "bark");
	}
}

/**
	Callback functions for the curve of the tree
*/
Vector3 App::straight(float t) {
    // x:cost z:sint
    return Vector3(0.0f, t, 0.0f);
}

Vector3 App::curvy(float t) { 
	return Vector3(-sqrt(t), t, 0.0f);
}


/**
	Callback functions for the radii of the branches
*/
float App::branchRadius(float t, shared_ptr<Tree> tree) {
    shared_ptr<Array<shared_ptr<Tree>>> children = tree->getChildren();
    //float end = branchRadius(0.0f, branchingNumber, recursionDepth - 1);
    if(children->size() == 0){
        return 0.005f;
    }
    float squareSum = 0.0f;
    float sum = 0.0f;

    for(int i = 0; i < children->size(); ++i){
        float childRadius = branchRadius(0.0f, children->operator[](i));
        squareSum += pow(childRadius, 2);
        sum += childRadius;
        
    }
    float end = (sum/(float)(children->size()));
    float base = sqrt(squareSum);


    debugAssert(t >= 0.0f && t <= 1.0f);
    //a, b, and c from the quadratic eq. ax^2 + bx + c = y
    float a = (base - end);
    float b = 2.0f*(end - base);
    float c = base;

    return (t*t*a) + (b*t) + c;
}

float App::envelopePerimeter(float y) {
    if(y < 0.20f) {
        return 0.0f;
    } else {
        return sin((5*pif() * y / 4) - pif()/4.0f);
    }
}

/**
	Callback functions for the phenotype of the tree
*/
void App::randomTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth) {  

    float newBranchLength = 3.0f * initialLength / 5.0f;
    float rPitch = rand() % 40 + 15.0f; //This will give a random pitch within the range 25 - 45
    float rYaw;
    CoordinateFrame branch1 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, rPitch, 0.0f);
    branch1.translation = branchEnd;
    
    rPitch = rand() % 40 + 15.0f; 
    rYaw = rand() % 40 + 70.0f; //This will give a random yaw fom 70 to 110
    CoordinateFrame branch2 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, rYaw, 0.0f, 0.0f);
    branch2 = branch2 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, rPitch, 0.0f);
    branch2.translation = branchEnd;
    
    rPitch = rand() % 40 + 15.0f; 
    rYaw = rand() % 40 + 70.0f;
    CoordinateFrame branch3 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, rYaw, 0.0f, 0.0f);
    branch3 = branch3 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -rPitch, 0.0f);
    branch3.translation = branchEnd;
    
    rPitch = rand() % 40 + 15.0f; 
    CoordinateFrame branch4 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -rPitch, 0.0f);
    branch4.translation = branchEnd;

    
    //These random values will determine how many child branches will be added (Fix realistic Radii with random child branches)
    bool makeBranch1 = (rand() % 100) > 30;
    bool makeBranch2 = (rand() % 100) > 30;
    bool makeBranch3 = (rand() % 100) > 30;
    bool makeBranch4 = (rand() % 100) > 30;

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
    

shared_ptr<Tree> App::makeSCTreeSkeleton(int anchorPointsCount, std::function<float(float)> envelopePerimeter, float height, float radius, float killDistance, float nodeDistance, float attractionRadius, Point3 initTreeNode) {

    shared_ptr<Tree> result = Tree::create(initTreeNode, nullptr);
    
	//Points in space
    Array<Point3> anchorPoints;
    generateAnchorPoints(anchorPoints, anchorPointsCount, height, radius, envelopePerimeter);
    Array<Point3> newAnchors;
    
    anchorPointsCount *= 0.01f;

    //Tree nodes about to grow
    Array<shared_ptr<Tree>> growingNodes;
    
    //Number of consecutive iterations that don't modify the tree.
    int staticIterations = 0;

    //Keep generating more tree nodes until all anchor points are killed.
    while(anchorPoints.size() > 0) {
        debugPrintf("Anchor nodes remaining: %d\n", anchorPoints.size());
        growingNodes.fastClear();   
        
        generateAnchorPoints(newAnchors, anchorPointsCount, height, radius, envelopePerimeter);
        anchorPoints.append(newAnchors);
        anchorPointsCount *= 0.9f;

        //For each anchor, select the closest tree node. Make this a method.
        for(const Point3& currentAnchor : anchorPoints) {
            
            Array<shared_ptr<Tree>> stack;
            stack.push(result);  
            shared_ptr<Tree> closestNode = nullptr;
            shared_ptr<Array<shared_ptr<Tree>>> children;

            while(stack.size() > 0) {
                shared_ptr<Tree> challenger = stack.pop();
                float distanceToNode = (challenger->getContents().translation - currentAnchor).magnitude();
                if(distanceToNode < attractionRadius && (isNull(closestNode) || distanceToNode < (closestNode->getContents().translation - currentAnchor).magnitude())) {
                    closestNode = challenger;
                }
                children = challenger->getChildren();

                stack.append(*children);
            }
            if(!isNull(closestNode)) {
                growingNodes.push(closestNode);
            }
        }  
        
        //Assign the directions for new tree nodes. Make a method
        std::map<shared_ptr<Tree>, Vector3> growthDirections;
        growthDirections.clear();

        for(int i = 0; i < growingNodes.size(); ++i) {
            shared_ptr<Tree> tempNode = growingNodes[i];
            if(growthDirections.find(tempNode) == growthDirections.end()) {
                growthDirections[tempNode] = (anchorPoints[i] - tempNode->getContents().translation).direction();
            } else {
                growthDirections[tempNode] = (growthDirections[tempNode] + (anchorPoints[i] - tempNode->getContents().translation).direction()).direction();
            }
        }

        //Iterate over the selected nodes and spawn new tree nodes. Make a method
        bool noNewNodes = true;
        for(auto it = growthDirections.begin(); it != growthDirections.end(); ++it) {
            shared_ptr<Tree> parent = it->first;
        
            Point3 newPos = parent->getContents().translation + (nodeDistance * growthDirections[parent]);
            debugPrintf("New tree node at: %f, %f, %f\n", newPos.x, newPos.y, newPos.z);
            shared_ptr<Tree> child = Tree::create(newPos, parent);
            shared_ptr<Array<shared_ptr<Tree>>> children = parent->getChildren();
            bool isNewNode = true;
            for(int i = 0; i < children->length(); ++i) {
                if(children->operator[](i)->getContents().fuzzyEq(newPos)) {
                    isNewNode = false;
                }
            }
            if(isNewNode) {
               children->push(child);
               noNewNodes = false;
            }
        }
        if(noNewNodes) {
            anchorPoints.clear();
        }

        //Kill anchors too close to the tree.
        for(Point3 currentAnchor : anchorPoints) {
            Array<shared_ptr<Tree>> stack;
            stack.push(result); while(stack.size() > 0) {
                shared_ptr<Tree> challenger = stack.pop();
                if((challenger->getContents().translation - currentAnchor).magnitude() < killDistance * nodeDistance) {
                    if(anchorPoints.findIndex(currentAnchor) > -1) {
                        anchorPoints.remove(anchorPoints.findIndex(currentAnchor));
                    }
                    staticIterations = 0;
                    break;
                }
                stack.append(*challenger->getChildren());
            }
        }
        staticIterations++;
        if(staticIterations > 100) {
            anchorPoints.fastClear();
        }
    }
    return result;
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
	Array<Point3> fruitLocations = Array<Point3>();

	if (m_typesIndex == 0) {
		makeLTree("firstTree", fruitLocations);
	}
	else {
		shared_ptr<Tree> skeleton = makeSCTreeSkeleton(m_spaceAnchorCount, [this](float y) {return App::envelopePerimeter(y);}, m_spaceHeight, m_spaceRadius, m_spaceKillDistance, m_spaceTreeDistance, m_spaceAttractionRadius, Point3(0,0,0));
		skeletonToMesh(m_spaceCirclePoints, m_spaceBranchRadius, m_spaceRadiusGrowth, "firstTree", skeleton, fruitLocations);
	}

	FruitDimensions fDims = fruitDims[m_fruitsIndex];

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
	writer.printf("filename = \"models/" + m_fruits[m_fruitsIndex] + "/" + fDims.filename + "\";");
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

    for (int i = 0; i < m_numRows; ++i) {
		float xOffset = m_maxRecursionDepth * i + rand.uniform();
		
		for (int j = 0; j < m_numTrees; ++j) {
			int zVar = rand.integer(0, 10);
			float zOffset = m_maxRecursionDepth * j + rand.uniform();

			writer.writeNewlines(2);
			writer.printf("tree%d%d = VisibleEntity { model = \"treeModel\";", i, j);
			writer.writeNewline();
			writer.printf("frame = CFrame::fromXYZYPRDegrees(%f, 0, %f); };", xOffset, zOffset);
			writer.writeNewlines(2);

			CoordinateFrame frame = CoordinateFrame(Point3(xOffset, 0.0f, zOffset));

			for (int k = 0; k < fruitLocations.length(); ++k) {
				int placeFruit = rand.integer(-2, 1);
				if (placeFruit == 1) {
					Point3 location = frame.pointToWorldSpace(fruitLocations[k]);
					writer.printf("fruit%d%d%d = VisibleEntity { model = \"fruitModel\";", i, j, k);
					writer.writeNewline();
					writer.printf("frame = CFrame::fromXYZYPRDegrees(%f, %f, %f, %d, %d, %d); };", location.x, location.y - fDims.yOffset, location.z, rand.integer(0,360), rand.integer(-20,20), 0);
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

void App::generateAnchorPoints(Array<Point3>& anchorPoints, int count, float height, float radius, std::function<float(float)> radiusCurve) {
    Random& rng = Random::threadCommon();
    for(int i = 0; i < count; ++i) {
        Point3 newPoint;
        do {
            newPoint = Point3((2.0f * rng.uniform()) - 1.0f, (2.0f * rng.uniform()) - 1.0f, (2.0f * rng.uniform()) - 1.0f);
        } while(newPoint.xz().length() >= radiusCurve(newPoint.y));
        anchorPoints.push(Point3(newPoint.x * radius, newPoint.y * height, newPoint.z * radius));
    }
}

void App::skeletonToMesh(int circlePoints, float initRadius, float radiusGrowth, String filename, shared_ptr<Tree>& skeleton, Array<Point3>& fruitLocations) {
    Mesh treeMesh = Mesh(filename);
    Mesh leafMesh = Mesh("leaf");

    Array<shared_ptr<Tree>> stack;
    Array<shared_ptr<Tree>> bottomUpIter;
    stack.push(skeleton);
    bottomUpIter.push(skeleton);

    std::map<shared_ptr<Tree>, int> indexMap;
    std::map<shared_ptr<Tree>, float> radiusMap;

    int curIndex = 0;
    while(stack.size() > 0) {
        shared_ptr<Tree> currentNode = stack.pop();
        Array<shared_ptr<Tree>> children = *currentNode->getChildren();


        indexMap[currentNode] = curIndex;
        curIndex += circlePoints;

        bottomUpIter.append(children);

        stack.append(children);
    }

    //Determine the radii for each branch. Has to use the bottom up iterator because
    //  the radius of each branch depends on the radii of its children.
    bottomUpIter.reverse();
    for(shared_ptr<Tree> treeNode : bottomUpIter) {
        if(treeNode->numChildren() == 0) {
            radiusMap[treeNode] = initRadius;
        } else {
            float tempSum = 0;
            for(shared_ptr<Tree> child : *treeNode->getChildren()) {
                tempSum += pow(radiusMap[child], radiusGrowth);
            }
            radiusMap[treeNode] = pow(tempSum, (1.0f/radiusGrowth));
        }
    }


    //Initial circle
    for(int i = 0; i < circlePoints; ++i) {
       	float angle = (i * 2.0f * pif()) / circlePoints;
        float radius = radiusMap[skeleton];
        Vector3 vec = Vector3(cos(angle) * radius, 0.0f, sin(angle) * radius);
	    treeMesh.addVertex(vec.x, vec.y, vec.z);
    }

    stack.fastClear();
    stack.append(*skeleton->getChildren());
    while(stack.size() > 0) {
        shared_ptr<Tree> currentNode = stack.pop();

        Point3 translation = currentNode->getContents().translation;
        Vector3 x, y, z;
        y = (currentNode->getParent()->getContents().translation - translation).direction();
        y.getTangents(x,z);
        CoordinateFrame nextFrame = CoordinateFrame(Matrix3::fromColumns(x, y, z), translation);
        CoordinateFrame leafFrame = CoordinateFrame(Matrix3::fromColumns(-x, -y, -z), translation);

        addCylindricSection(treeMesh, indexMap[currentNode->getParent()], indexMap[currentNode], circlePoints, nextFrame, radiusMap[currentNode]);
        if(currentNode->getChildren()->length() == 0) {
            float leafSize = 0.5f;
            addLeaf(leafMesh, leafSize, nextFrame);

			addFruits(fruitLocations, leafFrame);
        }
        stack.append(*currentNode->getChildren());
    }

    treeMesh.addMesh(leafMesh);
    treeMesh.toOBJ();
}