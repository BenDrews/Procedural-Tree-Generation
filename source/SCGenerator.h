/**
  \file SCGenerator.h
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

/** \brief Application framework. */
class SCGenerator {
protected:   
    void generateAnchorPoints(Array<Point3>& anchorPoints, int count, float height, float radius, std::function<float(float)> radiusCurve);
    void addIntermediateAnchors(Array<Point3>& anchorPoints, int& anchorPointsCount, float height, float radius, float discountRate, std::function<float(float)> envelopePerimeter);
    void findClosestNodes(Array<Point3>& anchorPoints, Array<shared_ptr<Tree>>& growingNodes, std::map<shared_ptr<Tree>, Vector3>& growthDirections, shared_ptr<Tree>& result, float attractionRadius);
    void growTreeNodes(Array<Point3>& anchorPoints, std::map<shared_ptr<Tree>, Vector3>& growthDirections, float nodeDistance);
    void killAnchors(Array<Point3>& anchorPoints, shared_ptr<Tree>& result, float nodeDistance, float killDistance);

    void addCylindricSection(Mesh& mesh, const int& parentIndex, const int& currentIndex, const int& pts, const CoordinateFrame& origin, const float& radius) const;
    void addLeaf(Mesh& leafMesh, float& length, const CoordinateFrame& initial) const;
	void addFruits(Array<Point3>& fruitLocations, const CoordinateFrame& fruitFrame);
public:
    shared_ptr<Tree> makeSCTreeSkeleton(int anchorPoints, std::function<float(float)> envelopePerimeter, float height, float radius, float killDistance, float nodeDistance, float attractionRadius, float discountRate, Point3 initTreeNode);
    
    void skeletonToMeshSC(int circlePoints, float initRadius, float radiusGrowth, float leafiness, String filename, shared_ptr<Tree>& skeleton, Array<Point3>& fruitLocations);
    
    static float envelopePerimeter(float y);
};