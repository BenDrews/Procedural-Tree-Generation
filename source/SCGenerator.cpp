/** \file SCGenerator.cpp */
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

/**Main methods for the space colonization algorithm **/
shared_ptr<Tree> SCGenerator::makeSCTreeSkeleton(int anchorPointsCount, std::function<float(float)> envelopePerimeter, float height, float radius, float killDistance, float nodeDistance, float attractionRadius, float discountRate, Point3 initTreeNode) {

    shared_ptr<Tree> result = Tree::create(initTreeNode, nullptr);
    
	//Points in space
    Array<Point3> anchorPoints;
    generateAnchorPoints(anchorPoints, anchorPointsCount, height, radius, envelopePerimeter);
    Array<Point3> newAnchors;

    //Tree nodes about to grow
    Array<shared_ptr<Tree>> growingNodes;

    //Directions in which the growing nodes will grow
    std::map<shared_ptr<Tree>, Vector3> growthDirections;

    //Keep generating more tree nodes until all anchor points are killed.
    while(anchorPoints.size() > 0) {
        debugPrintf("Anchor nodes remaining: %d\n", anchorPoints.size());
        growingNodes.fastClear();
        growthDirections.clear();
        newAnchors.fastClear();
        
       addIntermediateAnchors(anchorPoints, anchorPointsCount, height, radius, discountRate, envelopePerimeter);

       findClosestNodes(anchorPoints, growingNodes, growthDirections, result, attractionRadius);

       growTreeNodes(anchorPoints, growthDirections, nodeDistance);

       killAnchors(anchorPoints, result, nodeDistance, killDistance);   
    }
    return result;
}

/** Generates a set of anchor points within the specified envelope **/
void SCGenerator::generateAnchorPoints(Array<Point3>& anchorPoints, int count, float height, float radius, std::function<float(float)> radiusCurve) {
    Random& rng = Random::threadCommon();
    for (int i = 0; i < count; ++i) {
        Point3 newPoint;
        do {
            newPoint = Point3((2.0f * rng.uniform()) - 1.0f, (2.0f * rng.uniform()) - 1.0f, (2.0f * rng.uniform()) - 1.0f);
        } while (newPoint.xz().length() >= radiusCurve(newPoint.y));
        anchorPoints.push(Point3(newPoint.x * radius, newPoint.y * height, newPoint.z * radius));
    }
}

/**Generate a new set of anchors within the envelope and add them to the set of anchors.**/
void SCGenerator::addIntermediateAnchors(Array<Point3>& anchorPoints, int& anchorPointsCount, float height, float radius, float discountRate, std::function<float(float)> envelopePerimeter) {
        
        Array<Point3> newAnchors;
        generateAnchorPoints(newAnchors, anchorPointsCount, height, radius, envelopePerimeter);
        anchorPoints.append(newAnchors);

        //Reduce the number of anchors added for the next iteration.
        anchorPointsCount *= discountRate;
}

/**For each anchor, select the closest tree node. **/
void SCGenerator::findClosestNodes(Array<Point3>& anchorPoints, Array<shared_ptr<Tree>>& growingNodes, std::map<shared_ptr<Tree>, Vector3>& growthDirections, shared_ptr<Tree>& result, float attractionRadius) {
    
        for(const Point3& currentAnchor : anchorPoints) {

            //Create a stack to traverse the tree with.
            Array<shared_ptr<Tree>> stack;
            stack.push(result);  
            shared_ptr<Tree> closestNode = nullptr;
            shared_ptr<Array<shared_ptr<Tree>>> children;

            while(stack.size() > 0) {

                //For each tree node, check to see if that node is within the current anchor's radius of attraction.
                //If it is, also check to see if it is closer than the current closest node.
                shared_ptr<Tree> challenger = stack.pop();
                float distanceToNode = (challenger->getContents().translation - currentAnchor).magnitude();
                if(distanceToNode < attractionRadius && (isNull(closestNode) || distanceToNode < (closestNode->getContents().translation - currentAnchor).magnitude())) {
                    closestNode = challenger;
                }
                children = challenger->getChildren();

                stack.append(*children);
            }

            //If a node was found within the anchor's radius of attraction, add that node to the set of growing nodes.
            if(!isNull(closestNode)) {
                growingNodes.push(closestNode);

                //Check to see if the closest node is in growth directions already. If not add it.
                //If it is, add the new direction to the result.
                if(growthDirections.find(closestNode) == growthDirections.end()) {
                    growthDirections[closestNode] = (currentAnchor - closestNode->getContents().translation).direction();
                } else {
                    growthDirections[closestNode] = (growthDirections[closestNode] + (currentAnchor - closestNode->getContents().translation).direction());
                }
            }
        }

        for(auto it = growthDirections.begin(); it != growthDirections.end(); ++it) {
            growthDirections[it->first] = it->second.direction();
        }
}

/** Grow selected nodes in the correct directions **/
void SCGenerator::growTreeNodes(Array<Point3>& anchorPoints, std::map<shared_ptr<Tree>, Vector3>& growthDirections, float nodeDistance) {
    bool noNewNodes = true;
        for(auto it = growthDirections.begin(); it != growthDirections.end(); ++it) {
            shared_ptr<Tree> parent = it->first;
        
            Point3 newPos = parent->getContents().translation + (nodeDistance * growthDirections[parent]);
            debugPrintf("New tree node at: %f, %f, %f\n", newPos.x, newPos.y, newPos.z);
            shared_ptr<Tree> child = Tree::create(newPos, parent);
            shared_ptr<Array<shared_ptr<Tree>>> children = parent->getChildren();
            bool isNewNode = true;
            for(int i = 0; i < children->length(); ++i) {
                if(children->operator[](i)->getContents().fuzzyEq(newPos)) {
                    isNewNode = false;
                }
            }
            if(isNewNode) {
               children->push(child);
               noNewNodes = false;
            }
        }
        if(noNewNodes) {
            anchorPoints.clear();
        }
}

/**Kill anchors too close to the tree.**/
void SCGenerator::killAnchors(Array<Point3>& anchorPoints, shared_ptr<Tree>& result, float nodeDistance, float killDistance) {
 
    for(Point3 currentAnchor : anchorPoints) {

        //Create a stack to traverse the tree with.
        Array<shared_ptr<Tree>> stack;
        stack.push(result); 

        while(stack.size() > 0) {

            shared_ptr<Tree> challenger = stack.pop();
            if((challenger->getContents().translation - currentAnchor).magnitude() < killDistance * nodeDistance) {
                if(anchorPoints.findIndex(currentAnchor) > -1) {
                    anchorPoints.remove(anchorPoints.findIndex(currentAnchor));
                }
                break;
            }
            stack.append(*challenger->getChildren());
        }
    }
}

//Version that takes a hard coded index
void SCGenerator::addCylindricSection(Mesh& mesh, const int& parentIndex, const int& currentIndex, const int& pts, const CoordinateFrame& origin, const float& radius) const {

    for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
        Vector3 vec = Vector3(cos(angle) * radius, 0.0f, sin(angle) * radius);
        vec = origin.pointToWorldSpace(vec);

		mesh.addVertex(vec);
	}

	for(int i = 0; i < pts; ++i) {
		mesh.addFace(parentIndex + ((i + 1) % pts), parentIndex + i, currentIndex + i, 2, 3, 1, "bark");
		mesh.addFace(currentIndex + ((i + 1) % pts), parentIndex + ((i + 1) % pts), currentIndex + i, 2, 4, 3, "bark");
	}
}

/*
    Functions to define different envelope perimeters
*/


float SCGenerator::sphericalEnvelope(float y) {
    if(y < 0.20f) {
        return 0.0f;
    } else {
        return sin((5*pif() * y / 4) - pif()/4.0f);
    }
}

float SCGenerator::cylindricEnvelope(float y) {
    if(y < 0.20f) {
        return 0.0f;
    } else {
        return 1.0f;
    }
}

float SCGenerator::conicalEnvelope(float y) {
    if(y < 0.2f) {
        return 0.0f;
    } else {
        return -0.8f * y + 0.8f;
    }
}

float SCGenerator::bulbEnvelope(float y) {
    if(y < 0.20f) {
        return 0.0f;
    } else if (y >= 0.20f && y < 0.4f) {
        return sin((5*pif() * y / 2) - pif() / 2);
    } else if( y >= 0.4f) {
        return (-2.7777f * pow(y - 0.4f, 2.0f)) + 1.0f;
    }
}

//Lay a mesh over a Space Colonization tree skeleton
void SCGenerator::skeletonToMeshSC(int circlePoints, float initRadius, float radiusGrowth, float leafiness, String filename, shared_ptr<Tree>& skeleton, Array<Point3>& fruitLocations) {

    Mesh treeMesh = Mesh(filename);
    Mesh leafMesh = Mesh("leaf");

    //Stack will be used for forwared traversal of the tree
    //Bottom Up Iter will be used for leaf to base traversal of the tree
    Array<shared_ptr<Tree>> stack;
    Array<shared_ptr<Tree>> bottomUpIter;
    stack.push(skeleton);
    bottomUpIter.push(skeleton);

    //Maps will be used to store information about each branch section
    std::map<shared_ptr<Tree>, int> indexMap;
    std::map<shared_ptr<Tree>, float> radiusMap;

    //Current number of vertices in the mesh
    int curIndex = 0;
    while(stack.size() > 0) {
        shared_ptr<Tree> currentNode = stack.pop();
        Array<shared_ptr<Tree>> children = *currentNode->getChildren();

        indexMap[currentNode] = curIndex;

        //If the branch has multiple children, a circle is made for each so that branches don't twist out of each other in impossible angles
        if(children.length() > 1) {
            curIndex += circlePoints * children.length();   
        } 

        //Each branch section adds one circle to the mesh
        curIndex += circlePoints;

        bottomUpIter.append(children);
        stack.append(children);
    }

    //Determine the radii for each branch. Has to use the bottom up iterator because
    //  the radius of each branch depends on the radii of its children.
    bottomUpIter.reverse();
    for(shared_ptr<Tree> treeNode : bottomUpIter) {
        if(treeNode->numChildren() == 0) {
            radiusMap[treeNode] = initRadius;
        } else {
            float tempSum = 0;
            for(shared_ptr<Tree> child : *treeNode->getChildren()) {
                tempSum += pow(radiusMap[child], radiusGrowth);
            }
            radiusMap[treeNode] = pow(tempSum, (1.0f/radiusGrowth));
        }
    }


    //Initial circle
    for(int i = 0; i < circlePoints; ++i) {
       	float angle = (i * 2.0f * pif()) / circlePoints;
        float radius = radiusMap[skeleton];
        Vector3 vec = Vector3(cos(angle) * radius, 0.0f, sin(angle) * radius);
	    treeMesh.addVertex(vec.x, vec.y, vec.z);
    }

    stack.fastClear();
    stack.append(*skeleton->getChildren());
    while(stack.size() > 0) {
        shared_ptr<Tree> currentNode = stack.pop();

        //Find the frame of the new branch section
        Point3 translation = currentNode->getContents().translation;
        Vector3 x, y, z;
        y = (currentNode->getParent()->getContents().translation - translation).direction();
        y.getTangents(x,z);
        CoordinateFrame nextFrame = CoordinateFrame(Matrix3::fromColumns(x, y, z), translation);

        //If the branchs parent had multiple children, bind to the correctly oriented circle.
        //Else bind to the parents circle.
        if(currentNode->getParent()->numChildren() > 1) {
            addCylindricSection(treeMesh, indexMap[currentNode->getParent()] + (circlePoints * (1 + currentNode->getParent()->getChildren()->findIndex(currentNode))), indexMap[currentNode], circlePoints, nextFrame, radiusMap[currentNode]);    
        } else {
            addCylindricSection(treeMesh, indexMap[currentNode->getParent()], indexMap[currentNode], circlePoints, nextFrame, radiusMap[currentNode]);
        }

        //If the current node has multiple children, make a circle correctly oriented for each new child.
        //Bind that circle to this node's circle
        if(currentNode->numChildren() > 1) {
              
              for(int i = 0; i < currentNode->numChildren(); ++i) {

                //Find the frame pointing at the child.
                shared_ptr<Tree> childNode = currentNode->getChildren()->operator[](i);
                Vector3 childX, childY, childZ;
                childY = (translation - childNode->getContents().translation).direction();
                childY.getTangents(childX, childZ);
                CoordinateFrame branchStartFrame = CoordinateFrame(Matrix3::fromColumns(childX, childY, childZ), translation);
                addCylindricSection(treeMesh, indexMap[currentNode], indexMap[currentNode] + circlePoints * (1 + i), circlePoints, branchStartFrame, radiusMap[currentNode]);
              }
        }

        //Determine if this branch is section is thin enough to add leaves.
        if(radiusMap[currentNode] < leafiness * initRadius) {
            float leafSize = 0.5f;
            Random& rng = Random::threadCommon();
            CoordinateFrame leafFrame = CoordinateFrame(Matrix3::fromColumns(-x, -y, -z), translation);
            addLeaf(leafMesh, leafSize, leafFrame * CoordinateFrame::fromXYZYPRDegrees(0,0,0, rng.uniform(0.0f, 360.0f), rng.uniform(0.0f, 360.0f), rng.uniform(45.0f, 135.0f)));
			addFruits(fruitLocations, leafFrame);
        }

        stack.append(*currentNode->getChildren());
    }

    treeMesh.addMesh(leafMesh);
    treeMesh.toOBJ();
}

/**Function to add a leaf at a coordinate frame.**/
void SCGenerator::addLeaf(Mesh& leafMesh, float& leafSize, const CoordinateFrame& leafFrame) const {
    int index = leafMesh.numVertices();
    Vector3 vec1 = Vector3(leafSize / 2.0f, leafSize, 0.0f);
    vec1 = leafFrame.pointToWorldSpace(vec1);
    Vector3 vec2 = Vector3(-leafSize / 2.0f, leafSize, 0.0f);
    vec2 = leafFrame.pointToWorldSpace(vec2);
    Vector3 vec3 = Vector3(0.0f, 0.0f, 0.0f);
    vec3 = leafFrame.pointToWorldSpace(vec3);
    leafMesh.addVertex(vec1);
    leafMesh.addVertex(vec2);
    leafMesh.addVertex(vec3);
    leafMesh.addFace(index, index+1, index+2, 4, 3, 5, "leaf");
    leafMesh.addFace(index+2, index+1, index, 5, 3, 4, "leaf");
}

/**Function to add a fruit at a coordinate frame.**/
void SCGenerator::addFruits(Array<Point3>& fruitLocations, const CoordinateFrame& frame) {
	fruitLocations.push( frame.pointToWorldSpace(Vector3(0.0f, 0.0f, 0.0f)) );
}
