/** \file App.cpp */
#include "App.h"
#include "Mesh.h"
#include "Tree.h"
#include <cmath>
#include <map>

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
   // makeTree();

    makeTreeSkeleton(5, [this](float y) {return App::envelopePerimeter(y);}, 10.0f, 10.0f, 2.0f, 0.5f, Point3(0,0,0));

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

	//Code to add cylinder generation GUI.
    GuiPane* treePane = debugPane->addPane("Treez");

    treePane->setNewChildSize(240);
    treePane->addButton("Generate tree"); // TODO:: ADD BUTTON FUNCTIONALITY

    debugWindow->pack();
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}


void App::makeTree() {
    Mesh tree = Mesh("tree.OFF");
    Mesh leafMesh = Mesh("leaf.OFF");
    float length = 1.0f;
    makeBranch(tree, leafMesh, CoordinateFrame() * CoordinateFrame::fromXYZYPRDegrees(0,0,0,0,0,0), length, [this](float t) {return App::spineCurve(t);}, [this](float t, int depth) {return App::branchRadius(t, depth);}, 7, 3, 1);
    tree.toOFF();
    leafMesh.toOFF();
}


void App::makeBranch(Mesh& mesh, Mesh& leafMesh, const CoordinateFrame& initial, float& length, std::function<Vector3(float)> spineCurve, std::function<float(float, int)> branchRadius, int recursionDepth, int circlePoints, int branchSections) const {
    if (recursionDepth != 0) {
        int index = mesh.numVertices();
	    float sectionHeight;
	    float sectionRadius;
        Point3 branchEnd = initial.pointToWorldSpace(Point3(0,length,0));
        Point3 branchMid = initial.pointToWorldSpace(Point3(0,length*(6.0f/10.0f),0));
	    
        // Add vertices of intial circle at the top of the branch we are currently making to mesh
	    for(int i = 0; i < circlePoints; ++i) {
	    	float angle = (i * 2.0f * pif()) / circlePoints;
            sectionRadius = branchRadius(length, recursionDepth);
            Vector3 vec = Vector3(cos(angle) * sectionRadius, length, sin(angle) * sectionRadius);
            vec = initial.pointToWorldSpace(vec);
	    	mesh.addVertex(vec.x, vec.y, vec.z);  
	    }
	    
        // Add vertices of circles underneath the initial circle to mesh
	    for(int i = 0; i < branchSections; ++i) {
	    	sectionRadius = branchRadius(length - (i * length / float(branchSections)), recursionDepth);
	    	// TODO: pass a coordinate frame that is returned by space curve function (instead of initial)
            addCylindricSection(mesh, circlePoints, initial, sectionRadius);
	    }

        float newLength1 = 3*length / 5.0f;
        float newLength2 = 4*length / 5.0f;


        CoordinateFrame branch1 = initial * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, 35.0f, 0.0f);
        branch1.translation = branchEnd;
        CoordinateFrame branch2 = initial * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 90.0f, 0.0f, 0.0f);
        branch2 = branch2 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, 25.0f, 0.0f);        
        branch2.translation = branchEnd;
        CoordinateFrame branch3 = initial * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 90.0f, 0.0f, 0.0f);
        branch3 = branch3 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -25.0f, 0.0f);
        branch3.translation = branchEnd;
        CoordinateFrame branch4 = initial * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -35.0f, 0.0f);
        //branch4 = branch4 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -60.0f, 0.0f);
        branch4.translation = branchEnd;    
        CoordinateFrame branch5 = initial * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
        //branch3 = branch5 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -25.0f, 0.0f);
        branch5.translation = branchEnd;

        if(recursionDepth < 10){
            makeBranch(mesh, leafMesh, branch1, newLength1, spineCurve, branchRadius, recursionDepth - 1, circlePoints, branchSections);
            makeBranch(mesh, leafMesh, branch2, newLength1, spineCurve, branchRadius, recursionDepth - 1, circlePoints, branchSections);
            makeBranch(mesh, leafMesh, branch3, newLength1, spineCurve, branchRadius, recursionDepth - 1, circlePoints, branchSections);
            makeBranch(mesh, leafMesh, branch4, newLength1, spineCurve, branchRadius, recursionDepth - 1, circlePoints, branchSections);
        }
        makeBranch(mesh, leafMesh, branch5, newLength2, spineCurve, branchRadius, recursionDepth - 1, circlePoints, branchSections);

        if(recursionDepth == 1){
            CoordinateFrame leaf = initial;
            leaf.translation = branchEnd;
            addLeaves(leafMesh, length, leaf);
        }
    }
}


void App::addLeaves(Mesh& leafMesh, float& length, const CoordinateFrame& initial) const{
    int index = leafMesh.numVertices();
    float leafSize = length*4.0f;
    Vector3 vec1 = Vector3(leafSize / 2.0f, 0.0f, 0.0f);
    vec1 = initial.pointToWorldSpace(vec1);
    Vector3 vec2 = Vector3(-leafSize / 2.0f, 0.0f, 0.0f);
    vec2 = initial.pointToWorldSpace(vec2);
    Vector3 vec3 = Vector3(0.0f, leafSize, 0.0f);
    vec3 = initial.pointToWorldSpace(vec3);
    leafMesh.addVertex(vec1);
    leafMesh.addVertex(vec2);
    leafMesh.addVertex(vec3);
    leafMesh.addFace(index, index+1, index+2);
    leafMesh.addFace(index+2, index+1, index);

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
		mesh.addFace(offset + ((i + 1) % pts), offset + i + (pts), offset + i);
		mesh.addFace(offset + ((i + 1) % pts), offset + ((i + 1) % pts) + (pts), offset + i + (pts));
	}
}


Vector3 App::spineCurve(float t) {
    return Vector3(0.0f, t, 0.0f);
}


float App::branchRadius(float t, int recursionDepth) {
    return ((float)recursionDepth/100.0f);
}

float App::envelopePerimeter(float y) {
    if(y < 0.5f) {
        return 0.0f;
    } else {
        return 2.0f - (2.0f * y);
    }
}

Tree App::makeTreeSkeleton(int anchorPoints, std::function<float(float)> envelopePerimeter, float height, float radius, float killDistance, float nodeDistance, Point3 initTreeNode) {

    Tree result = Tree(initTreeNode);

    Tree testNode = Tree(initTreeNode);

    Array<Point3> anchorPointsList = generateAnchorPoints(anchorPoints, height, radius, envelopePerimeter);
    
    while(anchorPointsList.size() > 0) {
        Array<Tree*> selectedNodes = Array<Tree*>();
        
        //For each anchor, select the closest tree node.
        for(int i = 0; i < anchorPointsList.size(); ++i) {
            Point3 currentAnchor = anchorPointsList[i];
            Array<Tree*> stack = Array<Tree*>(&result);
            Tree* closestNode = &result;

            while(stack.size() > 0) {
                Tree* challenger = stack.pop();
                if((challenger->getContents() - currentAnchor).magnitude() < (closestNode->getContents() - currentAnchor).magnitude()) {
                    closestNode = challenger;
                }
                Array<Tree*>* tempArray = challenger->getChildren();
                stack.append(*challenger->getChildren());
            }
            selectedNodes.push(closestNode);
        }

        //Assign the directions for new tree nodes.
        std::map<Tree*, Vector3> nodeDirections;
        for(int i = 0; i < anchorPointsList.size(); ++i) {
            Tree* tempNode = selectedNodes[i];
            if(nodeDirections.find(tempNode) == nodeDirections.end()) {
                nodeDirections[tempNode] = (anchorPointsList[i] - tempNode->getContents()).direction();
            } else {
                nodeDirections[tempNode] = (nodeDirections[tempNode] + (anchorPointsList[i] - tempNode->getContents()).direction()).direction();
            }
        }

        //Iterate over the selected nodes and spawn new tree nodes.
        for(auto it = nodeDirections.begin(); it != nodeDirections.end(); ++it) {
            Tree* parent = it->first;
            Point3 pos1 = parent->getContents();
            Point3 pos2 = nodeDistance * nodeDirections[parent];
            Tree child = Tree(pos1 + pos2, parent);
            Array<Tree*>* children = parent->getChildren();
            children->push(&child);
        }

        //Kill anchors too close to the tree.
        for(int i = 0; i < anchorPointsList.size(); ++i) {
            Point3 currentAnchor = anchorPointsList[i];

            Array<Tree*> stack = Array<Tree*>(&result);

            while(stack.size() > 0) {
                Tree* challenger = stack.pop();
                if((challenger->getContents() - currentAnchor).magnitude() < killDistance * nodeDistance) {
                    anchorPointsList.remove(i);
                    break;
                }
                stack.append(*challenger->getChildren());
            }
        }
    }
    return result;
}

Array<Point3> App::generateAnchorPoints(int count, float height, float radius, std::function<float(float)> radiusCurve) {
    Random& rng = Random::threadCommon();
    Array<Point3> result = Array<Point3>();
    for(int i = 0; i < count; ++i) {
        Point3 newPoint;
        do {
            newPoint = Point3(rng.uniform(), rng.uniform(), rng.uniform());
        } while(newPoint.xz().length() >= radiusCurve(newPoint.y));
        result.push(Point3(newPoint.x * radius, newPoint.y * height, newPoint.z * radius));
    }
    return result;
}

void App::skeletonToMesh(float initRadius, float radiusGrowth, String filename, Tree skeleton) {
    
}
