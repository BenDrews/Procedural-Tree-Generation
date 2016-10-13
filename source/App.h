/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3DAll.h>
#include "Mesh.h"
#include "BranchDimensions.h"
//#include "Tree.h"

/** \brief Application framework. */
class App : public GApp {
protected:
    /** Called from onInit */
    void makeGUI();

    void makeTree();
    void makeBranch(Mesh& mesh, Mesh& leafMesh, const CoordinateFrame& initial, float& length, std::function<Vector3(float)> spineCurve, std::function<float(float, int)> branchRadius, 
        std::function<void(Array<BranchDimensions>&, float, const CoordinateFrame&, Point3&, int, int)> phenotype, int maxRecursionDepth, int currentRecursionDepth, int circlePoints, int branchSections) const;
	void addCylindricSection(Mesh& mesh, const int& pts, const CoordinateFrame& origin, const float& radius) const;
    Vector3 spineCurve(float t);
    float branchRadius(float t, int recursionDepth);
    void App::spikyTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);
    float App::envelopePerimeter(float y);

    //Tree makeTreeSkeleton(int anchorPoints, std::function<float(float)> envelopePerimeter, float height, float radius, float killDistance, float nodeDistance, Point3 initTreeNode);
    //Array<Point3> generateAnchorPoints(int count, float height, float radius, std::function<float(float)> radiusCurve);
    //void skeletonToMesh(float initRadius, float radiusGrowth, String filename, Tree Skeleton);
    void addLeaves(Mesh& leafMesh, float& length, const CoordinateFrame& initial) const;

public:
    
    App(const GApp::Settings& settings = GApp::Settings());

    virtual void onInit() override;
};
