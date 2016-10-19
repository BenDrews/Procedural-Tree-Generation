/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3DAll.h>
#include "Mesh.h"
#include "Tree.h"
#include "BranchDimensions.h"
#include "FruitDimensions.h"

/** \brief Application framework. */
class App : public GApp {
protected:
	int m_maxRecursionDepth = 4;
	float m_initialHeight = 1.0f;
	int m_circlePts = 3;
	int m_branchSections = 1;
	const Array<String> m_phenotypes = Array<String>("Normal", "Random", "Bush", "Pine");
	int m_phenotypesIndex = 0;

	int m_numRows = 2;
	int m_numTrees = 3;
	const Array<String> m_types = Array<String>("L-System", "Space Colonization");
	int m_typesIndex = 0;
	const Array<String> m_fruits = Array<String>("Apple", "Lemon", "Pear", "Banana", "Money", "Teapot");
	int m_fruitsIndex = 0;
	Array<FruitDimensions> fruitDims;

    int m_spaceAnchorCount = 1000;
    float m_spaceHeight = 20.0f;
    float m_spaceRadius = 10.0f;
    int m_spaceCirclePoints = 10;
    float m_spaceTreeDistance = 0.6f;
    float m_spaceKillDistance = 2.0f;
    float m_spaceBranchRadius = 0.01f;
    float m_spaceRadiusGrowth = 2.0f;
    float m_spaceAttractionRadius = 100.0f;
        

    /** Called from onInit */
    void makeGUI();

    void makeLTree(String filename, Array<Point3>& fruitLocations);
    shared_ptr<Tree> makeLTreeSkeleton(const CoordinateFrame& initial, std::function<void(Array<BranchDimensions>&, float, const CoordinateFrame&, Point3&, int, int)> phenotype, std::function<Vector3(float)> spineCurve, float length, int maxRecursionDepth, int currentRecursionDepth, shared_ptr<Tree> parent = nullptr);
    void buildTree(Mesh& mesh, Mesh& leafMesh, const shared_ptr<Tree> tree, std::function<Vector3(float)> spineCurve, std::function<float(float, shared_ptr<Tree>)> branchRadius, Array<Point3>& fruitLocations, int circlePoints, int branchSections, float initialLength);
    float distanceAlongBranch = 0.0f;
    void addLeaf(Mesh& leafMesh, float& length, const CoordinateFrame& initial) const;
    void addLeaves(CoordinateFrame& initial, float length, Mesh& leafMesh, Array<Point3>& fruitLocations);

	void addFruits(Array<Point3>& fruitLocations, const CoordinateFrame& fruitFrame);
    void addCylindricSection(Mesh& mesh, const int& pts, const CoordinateFrame& origin, const float& radius) const;
    
    Vector3 straight(float t);
	Vector3 curvy(float t);
    
    float branchRadius(float t, shared_ptr<Tree> tree);
    
    shared_ptr<Tree> makeSCTreeSkeleton(int anchorPoints, std::function<float(float)> envelopePerimeter, float height, float radius, float killDistance, float nodeDistance, float attractionRadius, Point3 initTreeNode);
    void generateAnchorPoints(Array<Point3>& anchorPoints, int count, float height, float radius, std::function<float(float)> radiusCurve);
    void skeletonToMesh(int circlePoints, float initRadius, float radiusGrowth, String filename, shared_ptr<Tree>& skeleton, Array<Point3>& fruitLocations);
    void App::randomTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);
    void App::normalTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);
    void App::bushTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);
    void App::pineTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);
	void App::generateOrchard();
    void addCylindricSection(Mesh& mesh, const int& parentIndex, const int& currentIndex, const int& pts, const CoordinateFrame& origin, const float& radius) const;

    float App::envelopePerimeter(float y);
  
public:
    App(const GApp::Settings& settings = GApp::Settings());
    virtual void onInit() override;
};