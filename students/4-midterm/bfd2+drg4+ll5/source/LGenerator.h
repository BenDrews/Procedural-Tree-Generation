/**
  \file LGenerator.h
 */
#pragma once
#include <G3D/G3DAll.h>
#include "App.h"
#include "Mesh.h"
#include "BranchDimensions.h"
#include "Tree.h"
#include "FruitDimensions.h"
#include "SCGenerator.h"
#include <cmath>
#include <map>
#include <tuple>
#include <stdlib.h>

/** \brief Logic to generate trees using an L-System */
class LGenerator {
protected:   
    
    float distanceAlongBranch = 0.0f;
    void addLeaf(Mesh& leafMesh, float& length, const CoordinateFrame& initial, String& leafName) const;
    void addLeaves(CoordinateFrame& initial, float length, Mesh& leafMesh, Array<Point3>& fruitLocations, bool fall);

	void addFruits(Array<Point3>& fruitLocations, const CoordinateFrame& fruitFrame);
    void addCylindricSection(Mesh& mesh, const int& pts, const CoordinateFrame& origin, const float& radius, String& barkName) const;
    
public:
    shared_ptr<Tree> makeLTreeSkeleton(const CoordinateFrame& initial, std::function<void(Array<BranchDimensions>&, float, const CoordinateFrame&, Point3&, int, int)> phenotype, std::function<Vector3(float)> spineCurve, float length, int maxRecursionDepth, int currentRecursionDepth, shared_ptr<Tree> parent = nullptr);
    void skeletonToMeshL(Mesh& mesh, Mesh& leafMesh, const shared_ptr<Tree> tree, std::function<Vector3(float)> spineCurve, std::function<float(float, shared_ptr<Tree>)> branchRadius, Array<Point3>& fruitLocations, int circlePoints, int branchSections, float initialLength, bool fall, String& bark);
  
    static Vector3 straight(float t);
	static Vector3 curvy(float t);
    static Vector3 corkscrew(float t);
    static Vector3 gentleCurve(float t);
    
    static float branchRadius(float t, shared_ptr<Tree> tree);
    static void randomTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);
    static void normalTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);
    static void bushTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);
    static void pineTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth);

};