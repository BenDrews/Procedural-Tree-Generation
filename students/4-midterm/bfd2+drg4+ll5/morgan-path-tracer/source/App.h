/**
  \file App.h
 */
#pragma once
#include <G3D/G3DAll.h>
#include "PathTracer.h"

/** \brief Application framework. */
class App : public GApp {
protected:

    GuiDropDownList*        m_resolution;
    PathTracer::Options     m_options;
    shared_ptr<PathTracer>  m_pathTracer;

    /** Called from onInit */
    void makeGUI();

    Vector2int32 resolution() const;

public:
    
    App(const GApp::Settings& settings = GApp::Settings());

    virtual void onInit() override;

    void onRender();
    void onRenderConvergence();
};
