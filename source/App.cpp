/** \file App.cpp */
#include "App.h"
#include "Mesh.h"
#include <cmath>

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
        "G3D Triangle" // Load something simple
        //developerWindow->sceneEditorWindow->selectedSceneName()  // Load the first scene encountered 
        );

    makeTree();
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
    float length = 5.0f;
    makeBranch(tree, CoordinateFrame(), length, [this](float t) {return App::spineCurve(t);}, [this](float t, int depth) {return App::branchRadius(t, depth);}, 2, 50, 50);
    tree.toOFF();
}


void App::makeBranch(Mesh& mesh, const CoordinateFrame& initial, float& length, std::function<Vector3(float)> spineCurve, std::function<float(float, int)> branchRadius, int recursionDepth, int circlePoints, int branchSections) const {
    if (recursionDepth != 0) {
        int index = mesh.numVertices();
	    float sectionHeight;
	    float sectionRadius;
        CoordinateFrame branchEnd = CoordinateFrame::fromXYZYPRDegrees(0, length, 0, 0, 0, 0) * initial;

	    //Adding verticies for outer circle of the top lip.
	    for(int i = 0; i < circlePoints; ++i) {
	    	float angle = (i * 2.0f * pif()) / circlePoints;
            sectionRadius = branchRadius(length, recursionDepth);
            Vector3 vec = Vector3(cos(angle) * sectionRadius, 0.0f, sin(angle) * sectionRadius);
            vec = branchEnd.pointToWorldSpace(vec);
	    	mesh.addVertex(vec.x, vec.y, vec.z);  
	    }
	    //Adding outer and inner vertices/faces for the sides of the glass.
	    for(int i = 0; i < branchSections; ++i) {
	    	sectionRadius = branchRadius(length - (i * length / float(branchSections)), recursionDepth);
	    	addCylindricSection(mesh, circlePoints, CoordinateFrame::fromXYZYPRDegrees(0, length - (i * length / float(branchSections)), 0, 0, 0, 0) * initial , sectionRadius);
	    }

        float newLength = length / 2.0f;
        makeBranch(mesh, branchEnd * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, 30.0f, 0.0f), newLength, spineCurve, branchRadius, recursionDepth - 1, circlePoints, branchSections);
        makeBranch(mesh, branchEnd * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -30.0f, 0.0f), newLength, spineCurve, branchRadius, recursionDepth - 1, circlePoints, branchSections);
    }    
}


void App::addCylindricSection(Mesh& mesh, const int& pts, const CoordinateFrame& origin, const float& radius) const {
	int index = mesh.numVertices();

    for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
        Vector3 vec = Vector3(cos(angle) * radius, 0.0f, sin(angle) * radius);
        vec = origin.pointToWorldSpace(vec);

        //TODO:: add mesh method to take vec directly
		mesh.addVertex(vec.x, vec.y, vec.z);
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
    return 1.0f + recursionDepth;
}