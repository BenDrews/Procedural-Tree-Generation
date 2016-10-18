/** \file App.cpp */
#include "App.h"
#include "PathTracer.h"

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
    settings.window.width               = 1280; settings.window.height       = 720;
    settings.window.fullScreen          = false;
    settings.window.resizable           = ! settings.window.fullScreen;
    settings.window.framed              = ! settings.window.fullScreen;
    settings.dataDir                    = FileSystem::currentDirectory();
    settings.screenshotDirectory        = "../journal/";

    settings.renderer.deferredShading = true;
    settings.renderer.orderIndependentTransparency = true;

    return App(settings).run();
}


App::App(const GApp::Settings& settings) : GApp(settings) {}


void App::onInit() {
    GApp::onInit();
    setFrameDuration(1.0f / 30.0f);

    showRenderingStats      = true;

    m_pathTracer = PathTracer::create();

    makeGUI();
    developerWindow->videoRecordDialog->setCaptureGui(false);
    developerWindow->cameraControlWindow->moveTo(Point2(developerWindow->cameraControlWindow->rect().x0(), 0));
    loadScene(
        "G3D Sponza (White)"
        //"G3D Triangle"
        //"G3D Cornell Box (Empty CO)"
 //       "G3D Cornell Box"
//        "Glossy floor"
        );
}


void App::makeGUI() {
    // Initialize the developer HUD (using the existing scene)
    createDeveloperHUD();
    debugWindow->setVisible(true);
    developerWindow->videoRecordDialog->setEnabled(true);

    m_resolution = debugPane->addDropDownList("Resolution", Array<String>({"1 x 1", "32 x 32", "64 x 40", "320 x 200", "640 x 400", "400 x 640", "1280 x 720", "720 x 1280", "1920 x 1080"}));
    debugPane->setNewChildSize(500, -1, 300);
    debugPane->addNumberBox("Rays per pixel",        &m_options.raysPerPixel, "", GuiTheme::LOG_SLIDER, 1, 8192 * 2);
    debugPane->addNumberBox("Max scattering events", &m_options.maxScatteringEvents, "", GuiTheme::LINEAR_SLIDER, 1, 10);
    debugPane->addButton("Render", this, &App::onRender);
    debugPane->addButton("Render Convergence Images", this, &App::onRenderConvergence);

#	ifdef G3D_DEBUG
		m_resolution->setSelectedIndex(0);
#	else
	    m_resolution->setSelectedIndex(3);
#	endif

    debugWindow->pack();
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}


Vector2int32 App::resolution() const {
    TextInput ti(TextInput::FROM_STRING, m_resolution->selectedValue());
    const int width = ti.readInteger();
    ti.readSymbol("x");
    const int height = ti.readInteger();
    return Vector2int32(width, height);
}


void App::onRender() {
    const Vector2int32 res = resolution();

    drawMessage("Rendering...");

    const shared_ptr<PathTracer>& PathTracer = PathTracer::create();

    m_pathTracer->setScene(scene());
    const shared_ptr<Image>& image = Image::create(res.x, res.y, ImageFormat::RGB32F());

    Stopwatch timer;
    timer.tick();
    m_pathTracer->traceImage(image, activeCamera(), m_options, [this](const String& msg, float pct) {
        debugPrintf("%d%% (%s)\n", iRound(100.0f * pct), msg.c_str());
    });
    timer.tock();

    const shared_ptr<Texture>& src = Texture::fromImage("Source", image);
    shared_ptr<Texture> dst;
    m_film->exposeAndRender(renderDevice, activeCamera()->filmSettings(), src, 0, 0, dst);    
    show(dst, format("%ds @ \n", iRound(timer.elapsedTime())) + System::currentTimeString());
}


void App::onRenderConvergence() {
    const Vector2int32 res = resolution();

    drawMessage("Rendering...");

    const shared_ptr<PathTracer>& PathTracer = PathTracer::create();

    m_pathTracer->setScene(scene());
    const shared_ptr<Image>& image = Image::create(res.x, res.y, ImageFormat::RGB32F());

    const PathTracer::Options oldOptions = m_options;

    for (int c = 1; c < oldOptions.raysPerPixel; c = iCeil(c * 1.25f)) {
        m_options.raysPerPixel = c;
        Stopwatch timer;
        timer.tick();
        m_pathTracer->traceImage(image, activeCamera(), m_options, [this](const String& msg, float pct) {
            debugPrintf("%d%% (%s)\n", iRound(100.0f * pct), msg.c_str());
        });
        timer.tock();

        const shared_ptr<Texture>& src = Texture::fromImage("Source", image);
        shared_ptr<Texture> dst;
        m_film->exposeAndRender(renderDevice, activeCamera()->filmSettings(), src, 0, 0, dst);
        const shared_ptr<Image>& result = dst->toImage(ImageFormat::RGB8());
        result->save(format("%s-convergence-%02d-rays-per-pixel.png", FilePath::makeLegalFilename(scene()->name()).c_str(), c));
        logPrintf("Time to render: %ds for %d rays/pixel\\n", iRound(timer.elapsedTime()), c);
    }

    m_options = oldOptions;
}
