/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3DAll.h>
#include "Mesh.h"

/** \brief Application framework. */
class App : public GApp {
protected:
	/** Set by Gui, ready by makeCylinder**/
	float m_cylinderRadius;
	float m_cylinderHeight;
	String m_heightfieldSource;
	float m_heightfieldYScale;
	float m_heightfieldXZScale;
	String m_heightfieldRingSource;
	float m_heightfieldRingRadius;
	float m_heightfieldRingRScale;
	float m_heightfieldRingZScale;
	int m_glassPts;
	int m_glassSec;
	float m_glassHeight;
    /** Called from onInit */
    void makeGUI();
	void makeCylinder();
	void fromHeightfield(const shared_ptr<Image> img, const float& yScale, const float& xzScale);
	void fromHeightfieldRing(const shared_ptr<Image> img, const float& rScale, const float& zScale, const float& radius);
	void makeGlass(const int& pts, const int& sections, const float& height, std::function<float(float)> callback);
	void addCylindricSection(const shared_ptr<Mesh> mesh, const int& index, const int& pts, const float& height, const float& radius);
	float glassCallback(float f);
	void App::message(const String& msg) const;
public:
    
    App(const GApp::Settings& settings = GApp::Settings());

    virtual void onInit() override;
};
