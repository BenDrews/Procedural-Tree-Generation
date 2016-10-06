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
        "G3D Whiteroom" // Load something simple
        //developerWindow->sceneEditorWindow->selectedSceneName()  // Load the first scene encountered 
        );
}


void App::makeGUI() {
    // Initialize the developer HUD
    createDeveloperHUD();

    debugWindow->setVisible(true);
    developerWindow->videoRecordDialog->setEnabled(true);

	//Code to add cylinder generation GUI.
    GuiPane* cylinderPane = debugPane->addPane("Cylinder");

	cylinderPane->setNewChildSize(240);
 
	cylinderPane->beginRow(); {
		cylinderPane->addNumberBox("Radius", &m_cylinderRadius, "m", GuiTheme::LINEAR_SLIDER, 1.0f, 20.0f);
		cylinderPane->addNumberBox("Height", &m_cylinderHeight, "m", GuiTheme::LINEAR_SLIDER, 1.0f, 20.0f);
		cylinderPane->addButton("Generate", [this](){
			
		message("Generating Model");
		makeCylinder();
		
		ArticulatedModel::clearCache();
		loadScene(developerWindow->sceneEditorWindow->selectedSceneName());  
}	);

	//Code to add heightfield generation GUI.
    GuiPane* heightfieldPane = debugPane->addPane("Heightfield");

	heightfieldPane->setNewChildSize(240);
 
	heightfieldPane->beginRow(); {
		heightfieldPane->addNumberBox("yScale", &m_heightfieldYScale, "m", GuiTheme::LOG_SLIDER, 0.01f, 100.0f);
		heightfieldPane->addNumberBox("xzScale", &m_heightfieldXZScale, "m", GuiTheme::LOG_SLIDER, 0.01f, 100.0f);
		heightfieldPane->addTextBox("Input Image", &m_heightfieldSource)->setWidth(210);
		heightfieldPane->addButton("...", [this]() {
			FileDialog::getFilename(m_heightfieldSource, "png", false);
			})->setWidth(30);
		} heightfieldPane->endRow();
    
		heightfieldPane->addButton("Generate", [this](){
		shared_ptr<Image> image;
		try {
			image = Image::fromFile(m_heightfieldSource);
        
			message("Generating Model");
			fromHeightfield(image, m_heightfieldYScale, m_heightfieldXZScale);
        
			ArticulatedModel::clearCache();
			loadScene(developerWindow->sceneEditorWindow->selectedSceneName());  
		 } catch (...) {
			 msgBox("Unable to load the image.", m_heightfieldSource);
		}
}	);

	//Code to add glass generation GUI.
	GuiPane* glassPane = debugPane->addPane("Drinking Glass");

	glassPane->setNewChildSize(240);
 
	glassPane->beginRow(); {
		glassPane->addNumberBox("Circle Points", &m_glassPts, "", GuiTheme::LINEAR_SLIDER, 10, 100);
		glassPane->addNumberBox("Sections", &m_glassSec, "", GuiTheme::LINEAR_SLIDER, 10, 100);
		glassPane->addNumberBox("Height", &m_glassHeight, "m", GuiTheme::LINEAR_SLIDER, 10.0f, 100.0f);
	}glassPane->endRow();

		glassPane->addButton("Generate", [this](){
		message("Generating Model");
		makeGlass(m_glassPts, m_glassSec, m_glassHeight, [this](float f) {return App::glassCallback(f);});
        
		ArticulatedModel::clearCache();
		loadScene(developerWindow->sceneEditorWindow->selectedSceneName());  
}	);

	//Code to add heightfield ring generation GUI.
	GuiPane* heightfieldRingPane = debugPane->addPane("Heightfield Ring");

	heightfieldRingPane->setNewChildSize(240);
 
	heightfieldRingPane->beginRow(); {
		heightfieldRingPane->addNumberBox("Radius", &m_heightfieldRingRadius, "m", GuiTheme::LINEAR_SLIDER, 10.0f, 100.0f);
		heightfieldRingPane->addNumberBox("rScale", &m_heightfieldRingRScale, "m", GuiTheme::LINEAR_SLIDER, 10.0f, 100.0f);
		heightfieldRingPane->addNumberBox("zScale", &m_heightfieldRingZScale, "m", GuiTheme::LINEAR_SLIDER, 0.01f, 1.0f);
		heightfieldRingPane->addTextBox("Input Image", &m_heightfieldRingSource)->setWidth(210);
		heightfieldRingPane->addButton("...", [this]() {
			FileDialog::getFilename(m_heightfieldRingSource, "png", false);
			})->setWidth(30);
		} heightfieldRingPane->endRow();
		heightfieldRingPane->addButton("Generate", [this](){
		shared_ptr<Image> image;
		try {
			image = Image::fromFile(m_heightfieldRingSource);

			message("Generating Model");
			fromHeightfieldRing(image, m_heightfieldRingRScale, m_heightfieldRingZScale, m_heightfieldRingRadius);

			ArticulatedModel::clearCache();
			loadScene(developerWindow->sceneEditorWindow->selectedSceneName());        
		 } catch (...) {
			 msgBox("Unable to load the image.", m_heightfieldSource);
			 ArticulatedModel::clearCache();
		}
}	);
    debugWindow->pack();
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}
}

/** Function to  **/
void App::makeCylinder() {
	Mesh mesh = Mesh("model/cylinder.off");
	const int circlePts = 50;
	for(int i = 0; i < circlePts; ++i) {
		float angle = (i * 2.0f * pif()) / circlePts;
		mesh.addVertex(m_cylinderRadius * cos(angle), 0.5f * m_cylinderHeight, m_cylinderRadius * sin(angle));
	}
	for(int i = 0; i < circlePts; ++i) {
		float angle = (i * 2.0f * pif()) / circlePts;
		mesh.addVertex(m_cylinderRadius * cos(angle), -0.5f * m_cylinderHeight, m_cylinderRadius * sin(angle));
	}
	for(int i = 0; i < circlePts; ++i) {
		mesh.addFace((i + 1) % circlePts, i + circlePts, i);
		mesh.addFace((i + 1) % circlePts, ((i + 1) % circlePts) + circlePts, i + circlePts);
	}
	for(int i = 0; i < circlePts - 2; ++i) {
		mesh.addFace(i + 2, i + 1, 0);
		mesh.addFace(i + circlePts, i + circlePts + 1, circlePts + 2);
	}
	mesh.toOFF();
}

void App::fromHeightfield(const shared_ptr<Image> img, const float& yScale, const float& xzScale) {
	int height = img->height();
	int width = img->width();
	
	Mesh mesh = Mesh("model/heightfield.off");

	Color1 color;
	float depth;
	for(int32 x = 0; x < width; ++x) {
		for(int32 z = 0; z < height; ++z) {
		img->get(Point2int32(x,z),color);
		depth = color.value * yScale;
		mesh.addVertex(x * xzScale, depth, z * xzScale);
		}
	}
	int offset;
	for(int32 x = 1; x < width; ++x) {
		for(int32 z = 1; z < height; ++z) {
        offset = (x - 1) * height;
		mesh.addFace(offset + (z - 1), offset + z, offset + height + (z - 1));
		mesh.addFace(offset + height + (z - 1), offset + z, offset + height + z);
		}
	}
	mesh.toOFF();
}

float App::glassCallback(float f) {
	float result;
	if (f > 0.4f) {
		result = max(abs(5.0f * atan(10.0f * (f-0.4f))), 1.0f);
	} else if (f < 0.1) {
		result = max(5.0f * (pow(2.0f, -1.0f * (20.0f * f))), 1.0f);
	} else {
		result = 1.0f;
	}
	return result;
}

void App::makeGlass(const int& pts, const int& sections, const float& height, std::function<float(float)> callback) {
	shared_ptr<Mesh> mesh(new Mesh("model/glass.off"));

	float sectionHeight;
	float sectionRadius;

	//Adding verticies for outer circle of the top lip.
	for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
		sectionHeight = height; 
		sectionRadius = callback(1.0f);
		mesh->addVertex(sectionRadius * cos(angle), sectionHeight, sectionRadius * sin(angle));  
	}
	//Adding verticies for the inner cirlce of the top lip.
	for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
		sectionHeight = height; 
		sectionRadius = callback(1.0f);
		mesh->addVertex((sectionRadius - 1.0f) * cos(angle), sectionHeight, (sectionRadius - 1.0f) * sin(angle));  
	}
	//Adding faces along the top lip.
	for(int i = 0; i < pts; ++i) {
		mesh->addFace(i, i + pts, (i + 1) % pts);
		mesh->addFace(((i + 1) % pts), i + pts, ((i + 1) % pts) + pts);
	}
	//Adding outer and inner vertices/faces for the sides of the glass.
	for(int i = 0; i < sections; ++i) {
		sectionHeight = height - (i * height / sections);
		sectionRadius = callback(1.0f - (float(i) / float(sections)));
		addCylindricSection(mesh, i, pts, sectionHeight, sectionRadius);
	}
	//Add bottom faces.
	for(int i = 1; i < pts - 1; ++i) {
		int offset = sections * pts * 2;
		mesh->addFace(offset, offset + i, offset + i + 1);
	}
	mesh->toOFF();
}


void App::addCylindricSection(const shared_ptr<Mesh> mesh, const int& index, const int& pts, const float& height, const float& radius) {
	for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
		mesh->addVertex(radius * cos(angle), height, radius * sin(angle));
	}
	for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
		mesh->addVertex((radius - 1.0f) * cos(angle), height, (radius - 1.0f) * sin(angle));
	}
	int offset = index * pts * 2;
	for(int i = 0; i < pts; ++i) {
		mesh->addFace(offset + ((i + 1) % pts), offset + i + (2 * pts), offset + i);
		mesh->addFace(offset + ((i + 1) % pts), offset + ((i + 1) % pts) + (2 * pts), offset + i + (2 * pts));
	}
	offset += pts;
	for(int i = 0; i < pts; ++i) {
		mesh->addFace(offset + ((i + 1) % pts), offset + i, offset + i + (2 * pts));
		mesh->addFace(offset + ((i + 1) % pts) + (2 * pts), offset + ((i + 1) % pts), offset + i + (2 * pts));
	}
}

//Generates a scene from a heightfield along the inside of a ring.
void App::fromHeightfieldRing(const shared_ptr<Image> img, const float& rScale, const float& zScale, const float& radius) {
	int height = img->height();
	int width = img->width();
	
	Mesh mesh = Mesh("model/heightfieldRing.off");

	Color1 color;
	float depth;
	float angle;
	for(int32 x = 0; x < width; ++x) {
		for(int32 z = 0; z < height; ++z) {
		if(z == 0 || z == height - 1) {
			depth = 0;
		} else {
			img->get(Point2int32(x,z),color);
			depth = color.value * rScale;
		}
		angle = x * pif() / (width - 2);
		mesh.addVertex( cos(angle) * (radius - depth), sin(angle) * (radius - depth), z * zScale);
		}
	}
	for(int32 x = 0; x < width; ++x) {
		for(int32 z = 0; z < height; ++z) {
			if(z == 0 || z == height - 1) {
				depth = 0;
			} else {
				img->get(Point2int32(x,z), color);
				depth = color.value * rScale;
			}
			angle = -1.0f * x * pif() / (width - 2);
			mesh.addVertex( cos(angle) * (radius - depth), sin(angle) * (radius - depth), z * zScale);
		}
	}

	//Generate the faces, the offset is used for the mirrored half of the ring.
	int offset = height * width;
	for(int32 x = 1; x < width; ++x) {
		for(int32 z = 1; z < height; ++z) {
			mesh.addFace((x - 1) * height + (z - 1), (x - 1) * height + z, (x * height) + (z - 1));
			mesh.addFace((x * height) + (z - 1), (x - 1) * height + z, (x * height) + z);
			mesh.addFace(offset + (x - 1) * height + z, offset + (x - 1) * height + (z - 1), offset + (x * height) + (z - 1));
			mesh.addFace(offset + (x * height) + z, offset + (x - 1) * height + z, offset + (x * height) + (z - 1));
		}
	}
	mesh.toOFF();
}

void App::message(const String& msg) const {
    renderDevice->clear();
    renderDevice->push2D();
    debugFont->draw2D(renderDevice, msg, renderDevice->viewport().center(), 12,
        Color3::white(), Color4::clear(), GFont::XALIGN_CENTER, GFont::YALIGN_CENTER);
    renderDevice->pop2D();

    // Force update so that we can see the message
    renderDevice->swapBuffers();
}