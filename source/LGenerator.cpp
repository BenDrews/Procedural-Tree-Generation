/** \file SCGenerator.cpp */
#include "App.h"
#include "Mesh.h"
#include "BranchDimensions.h"
#include "Tree.h"
#include "FruitDimensions.h"
#include "LGenerator.h"
#include <cmath>
#include <map>
#include <tuple>
#include <stdlib.h>

/**Main method for L-System tree generation **/
shared_ptr<Tree> LGenerator::makeLTreeSkeleton(const CoordinateFrame& initial, std::function<void(Array<BranchDimensions>&, float, const CoordinateFrame&, Point3&, int, int)> phenotype, std::function<Vector3(float)> spineCurve, float length, int maxRecursionDepth, int currentRecursionDepth, shared_ptr<Tree> parent){
    shared_ptr<Tree> tree = Tree::create(initial, parent);
    Point3 branchEnd = initial.pointToWorldSpace(length * spineCurve(19.0f/20.0f));

    // callback function to decide how to recurse, populates an array of BranchDimensions which contain coordinate frames and lengths for the next branches
    Array<BranchDimensions> nextBranches;
    phenotype(nextBranches, length, initial, branchEnd, maxRecursionDepth, currentRecursionDepth);

    if (currentRecursionDepth > 0) {
        for (int i = 0; i < nextBranches.length(); ++i) {
            BranchDimensions nextBranch = nextBranches[i];
            CoordinateFrame branch = nextBranch.frame;
            float newLength = nextBranch.length;
            
            shared_ptr<Tree> childTree = makeLTreeSkeleton(branch, phenotype, spineCurve, newLength, maxRecursionDepth, currentRecursionDepth - 1, tree);
            if(notNull(childTree)){
                tree->addChild(childTree);
            }
        }
        return tree;
    } else {
        return nullptr;
    }
}

/**Converts an L-System tree skeleton to a mesh. **/
void LGenerator::skeletonToMeshL(Mesh& mesh, Mesh& leafMesh, const shared_ptr<Tree> tree, std::function<Vector3(float)> spineCurve, std::function<float(float, shared_ptr<Tree>)> branchRadius, Array<Point3>& fruitLocations, int circlePoints, int branchSections, float initialLength){
    float distanceAlongBranch = 0.0f;
    shared_ptr<Array<shared_ptr<Tree>>> children = tree->getChildren();
    float sectionRadius;
    float length = initialLength;
    CoordinateFrame initial = tree->getContents();


    if(children->size() != 0){
       length = (initial.translation - children->operator[](0)->getContents().translation).magnitude();
    }

    
    // Add vertices of intial circle at the bottom of the branch we are currently making to mesh
	for(int i = 0; i < circlePoints; ++i) {
		float angle = (i * 2.0f * pif()) / circlePoints;
        sectionRadius = branchRadius(distanceAlongBranch, tree);
        Vector3 vec = Vector3(cos(angle) * sectionRadius, 0.0f, sin(angle) * sectionRadius);
        vec = initial.pointToWorldSpace(vec);
		mesh.addVertex(vec.x, vec.y, vec.z);  
	}

     // Add vertices of circles on top of the initial circle to mesh
	for(int i = 1; i <= branchSections; ++i) {
        distanceAlongBranch =  (float)(i) / float(branchSections);
		sectionRadius = branchRadius(distanceAlongBranch, tree);
		// TODO:: pass a coordinate frame that is returned by space curve function (instead of initial)
        CoordinateFrame section = initial;
        section.translation = initial.pointToWorldSpace(length * spineCurve(distanceAlongBranch));
        addCylindricSection(mesh, circlePoints, section, sectionRadius);
	}

    if(children->size() == 0){
        addLeaves(initial, length, leafMesh, fruitLocations);
        
    }else{
        float newLength = ((3.0f/5.0f)*length);
            
        for(int i = 0; i < children->size(); ++i){
            skeletonToMeshL(mesh, leafMesh, children->operator[](i), spineCurve, branchRadius, fruitLocations, circlePoints, branchSections, newLength);
        }
    
    }
}

/**Adds leaves and fruits to a mesh**/
void LGenerator::addLeaves(CoordinateFrame& initial, float length, Mesh& leafMesh, Array<Point3>& fruitLocations){
    Random& rand = Random::threadCommon();

    //Add one leaf on the end of the branch
    CoordinateFrame leafFrame = initial;
    leafFrame.translation = initial.pointToWorldSpace(Point3(0.0f, (19.0f/20.0f)*length, 0.0f));
    float leafSize = 0.15f;
    addLeaf(leafMesh, leafSize, leafFrame);
	addFruits(fruitLocations, leafFrame);
    
    //Add random leaves along branch
    int leafNumber = 5;
    float rPitch;
    float rYaw;
    float rDisplacement;
    for(int i = 0; i < 5; ++i){
        leafFrame = initial;
        rDisplacement = ( (float)rand.integer(0,100) / 100.0f ) * length;
        rYaw = rand.integer(0, 360);
        rPitch = rand.integer(0, 40) + 30.0f;  
        leafFrame.translation = initial.pointToWorldSpace(Point3(0.0f, rDisplacement, 0.0f));
        leafFrame = leafFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, rYaw);
        leafFrame = leafFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, rPitch);
        addLeaf(leafMesh, leafSize, leafFrame);
    }

}

/**Adds a single leaf to a mesh **/
void LGenerator::addLeaf(Mesh& leafMesh, float& leafSize, const CoordinateFrame& leafFrame) const {
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

/**Adds a single fruit to a mesh**/
void LGenerator::addFruits(Array<Point3>& fruitLocations, const CoordinateFrame& frame) {
	fruitLocations.push( frame.pointToWorldSpace(Vector3(0.0f, 0.0f, 0.0f)) );
}

/**Creates a new branch section in a mesh**/
void LGenerator::addCylindricSection(Mesh& mesh, const int& pts, const CoordinateFrame& origin, const float& radius) const {
	int index = mesh.numVertices();

    for(int i = 0; i < pts; ++i) {
		float angle = (i * 2.0f * pif()) / pts;
        Vector3 vec = Vector3(cos(angle) * radius, 0.0f, sin(angle) * radius);
        vec = origin.pointToWorldSpace(vec);

		mesh.addVertex(vec);
	}
	int offset = index - pts;
	for(int i = 0; i < pts; ++i) {
		mesh.addFace( offset + i, offset + i + (pts),offset + ((i + 1) % pts), 2, 3, 1, "bark");
		mesh.addFace(offset + i + (pts), offset + ((i + 1) % pts) + (pts), offset + ((i + 1) % pts), 2, 4, 3, "bark");
	}
}

/**
	Callback functions for the curve of the tree
*/
Vector3 LGenerator::straight(float t) {
    return Vector3(0.0f, t, 0.0f);
}

Vector3 LGenerator::curvy(float t) { 
	return Vector3(-sqrt(t), t, 0.0f);
}

Vector3 LGenerator::corkscrew(float t) {
    return Vector3(sin(pif() * t * 4), t, cos(pif() * t * 4) - 1);
}


/**
	Callback functions for the radii of the branches
*/
float LGenerator::branchRadius(float t, shared_ptr<Tree> tree) {
    shared_ptr<Array<shared_ptr<Tree>>> children = tree->getChildren();
    //float end = branchRadius(0.0f, branchingNumber, recursionDepth - 1);
    if(children->size() == 0){
        return 0.005f;
    }
    float squareSum = 0.0f;

    float largestRadius = 0.0f;
    for(int i = 0; i < children->size(); ++i){
        float childRadius = branchRadius(0.0f, children->operator[](i));
        squareSum += pow(childRadius, 2);
        
        if(largestRadius <= childRadius){
            largestRadius = childRadius;
        }
        
    }
    float end = largestRadius;
    float base = sqrt(squareSum);


    debugAssert(t >= 0.0f && t <= 1.0f);
    //a, b, and c from the quadratic eq. ax^2 + bx + c = y
    float a = (base - end);
    float b = 2.0f*(end - base);
    float c = base;

    return (t*t*a) + (b*t) + c;
}


/**
	Callback functions for the phenotype of the tree
*/
void LGenerator::randomTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth) {  

    float newBranchLength = 3.0f * initialLength / 5.0f;
    float rPitch = rand() % 40 + 15.0f; //This will give a random pitch within the range 25 - 45
    float rYaw;
    CoordinateFrame branch1 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, rPitch, 0.0f);
    branch1.translation = branchEnd;
    
    rPitch = rand() % 40 + 15.0f; 
    rYaw = rand() % 40 + 70.0f; //This will give a random yaw fom 70 to 110
    CoordinateFrame branch2 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, rYaw, 0.0f, 0.0f);
    branch2 = branch2 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, rPitch, 0.0f);
    branch2.translation = branchEnd;
    
    rPitch = rand() % 40 + 15.0f; 
    rYaw = rand() % 40 + 70.0f;
    CoordinateFrame branch3 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, rYaw, 0.0f, 0.0f);
    branch3 = branch3 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -rPitch, 0.0f);
    branch3.translation = branchEnd;
    
    rPitch = rand() % 40 + 15.0f; 
    CoordinateFrame branch4 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -rPitch, 0.0f);
    branch4.translation = branchEnd;

    
    //These random values will determine how many child branches will be added (Fix realistic Radii with random child branches)
    bool makeBranch1 = (rand() % 100) > 30;
    bool makeBranch2 = (rand() % 100) > 30;
    bool makeBranch3 = (rand() % 100) > 30;
    bool makeBranch4 = (rand() % 100) > 30;

    if (makeBranch1) {
        nextBranches.push(BranchDimensions(branch1, newBranchLength));
    }
    if (makeBranch2) { 
        nextBranches.push(BranchDimensions(branch2, newBranchLength));
    }
    if (makeBranch3) {
        nextBranches.push(BranchDimensions(branch3, newBranchLength));
    }
    if (makeBranch4) {
        nextBranches.push(BranchDimensions(branch4, newBranchLength));
    }
    debugAssert(nextBranches.size() > 0);
}
    
void LGenerator::normalTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth) {  

    float newBranchLength = 3.0f * initialLength / 5.0f;


    CoordinateFrame branch1 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, 35.0f, 0.0f);
    branch1.translation = branchEnd;

    CoordinateFrame branch2 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -35.0f, 0.0f);
    branch2.translation = branchEnd;

    CoordinateFrame branch3 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 90.0f, 0.0f, 0.0f);
    branch3 = branch3 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, 35.0f, 0.0f);
    branch3.translation = branchEnd;
    
    CoordinateFrame branch4 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 90.0f, 0.0f, 0.0f);
    branch4 = branch4 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -35.0f, 0.0f);
    branch4.translation = branchEnd;
    
    nextBranches.push(BranchDimensions(branch1, newBranchLength));
    nextBranches.push(BranchDimensions(branch2, newBranchLength));
    nextBranches.push(BranchDimensions(branch3, newBranchLength));
    nextBranches.push(BranchDimensions(branch4, newBranchLength));
 
    debugAssert(nextBranches.size() > 0);
}

void LGenerator::bushTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth) {  

    float trunkLength = 4.0f/5.0f * initialLength;
    float newBranchLength = 3.0f * initialLength / 5.0f;

    float levelPitch = 60.0f*((float)(currentRecursionDepth) / (float)(maxRecursionDepth));
    float pitch = rand()%20 + levelPitch;

    float yaw = rand()%20 - 10;

    CoordinateFrame branch1 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch1 = branch1 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, pitch, 0.0f);
    branch1.translation = branchEnd;

    yaw = rand()%20 - 10;
    pitch = rand()%20 + levelPitch;

    CoordinateFrame branch2 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch2 = branch2 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -pitch, 0.0f);
    branch2.translation = branchEnd;

    yaw = rand()%20 + 80;
    pitch = rand()%20 + levelPitch;

    CoordinateFrame branch3 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch3 = branch3 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, pitch, 0.0f);
    branch3.translation = branchEnd;
    
    yaw = rand()%20 + 80;
    pitch = rand()%20 + levelPitch;

    CoordinateFrame branch4 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch4 = branch4 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -pitch, 0.0f);
    branch4.translation = branchEnd;
    
    CoordinateFrame branch5 = initialFrame;
    branch5.translation = branchEnd;

    //These random values will determine how many child branches will be added (Fix realistic Radii with random child branches)
    bool makeBranch1 = (rand() % 100) > 20;
    bool makeBranch2 = (rand() % 100) > 20;
    bool makeBranch3 = (rand() % 100) > 20;
    bool makeBranch4 = (rand() % 100) > 20;

    if(makeBranch1){
        nextBranches.push(BranchDimensions(branch1, newBranchLength));
    }
    if(makeBranch2){
        nextBranches.push(BranchDimensions(branch2, newBranchLength));
    }
    if(makeBranch3){
        nextBranches.push(BranchDimensions(branch3, newBranchLength));
    }
    if(makeBranch4){
        nextBranches.push(BranchDimensions(branch4, newBranchLength));
    }
    nextBranches.push(BranchDimensions(branch5, trunkLength));

    debugAssert(nextBranches.size() > 0);
}

void LGenerator::pineTree(Array<BranchDimensions>& nextBranches, const float initialLength, const CoordinateFrame& initialFrame, const Point3& branchEnd, const int maxRecursionDepth, const int currentRecursionDepth) {  

    float trunkLength = 4.0f/5.0f * initialLength;
    float newBranchLength = 2.0f * initialLength / 5.0f;

    float levelPitch = 100 + 60.0f*((float)(maxRecursionDepth - currentRecursionDepth) / (float)(maxRecursionDepth));

    float pitch = rand()%20 + levelPitch;
    float yaw = rand()%20 - 10;
    Point3 randomPos = initialFrame.pointToWorldSpace(Point3(0.0f,(((float)(rand()%10) / 10.0f)*initialLength),0.0f));

    CoordinateFrame branch1 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch1 = branch1 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, pitch, 0.0f);
    branch1.translation = randomPos;

    yaw = rand()%20 - 10;
    pitch = rand()%20 + levelPitch;
    randomPos = initialFrame.pointToWorldSpace(Point3(0.0f,(((float)(rand()%10) / 10.0f)*initialLength),0.0f));

    CoordinateFrame branch2 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch2 = branch2 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -pitch, 0.0f);
    branch2.translation = randomPos;

    yaw = rand()%20 + 80;
    pitch = rand()%20 + levelPitch;
    randomPos = initialFrame.pointToWorldSpace(Point3(0.0f,(((float)(rand()%10) / 10.0f)*initialLength),0.0f));

    CoordinateFrame branch3 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch3 = branch3 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, pitch, 0.0f);
    branch3.translation = randomPos;
    
    yaw = rand()%20 + 80;
    pitch = rand()%20 + levelPitch;
    randomPos = initialFrame.pointToWorldSpace(Point3(0.0f,(((float)(rand()%10) / 10.0f)*initialLength),0.0f));

    CoordinateFrame branch4 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch4 = branch4 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -pitch, 0.0f);
    branch4.translation = randomPos;

    yaw = rand()%20 + 45;
    pitch = rand()%20 + levelPitch;
    randomPos = initialFrame.pointToWorldSpace(Point3(0.0f,(((float)(rand()%10) / 10.0f)*initialLength),0.0f));

    CoordinateFrame branch5 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch5 = branch5 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, pitch, 0.0f);
    branch5.translation = randomPos;
    
    yaw = rand()%20 + 45;
    pitch = rand()%20 + levelPitch;
    randomPos = initialFrame.pointToWorldSpace(Point3(0.0f,(((float)(rand()%10) / 10.0f)*initialLength),0.0f));

    CoordinateFrame branch6 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch5 = branch5 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -pitch, 0.0f);
    branch5.translation = randomPos;

    yaw = rand()%20 + 135;
    pitch = rand()%20 + levelPitch;
    randomPos = initialFrame.pointToWorldSpace(Point3(0.0f,(((float)(rand()%10) / 10.0f)*initialLength),0.0f));

    CoordinateFrame branch7 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch7 = branch7 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, pitch, 0.0f);
    branch7.translation = randomPos;
    
    yaw = rand()%20 + 135;
    pitch = rand()%20 + levelPitch;
    randomPos = initialFrame.pointToWorldSpace(Point3(0.0f,(((float)(rand()%10) / 10.0f)*initialLength),0.0f));

    CoordinateFrame branch8 = initialFrame * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, yaw, 0.0f, 0.0f);
    branch8 = branch8 * CoordinateFrame::fromXYZYPRDegrees(0.0f, 0.0f, 0.0f, 0.0f, -pitch, 0.0f);
    branch8.translation = randomPos;
    
    CoordinateFrame trunk = initialFrame;
    trunk.translation = branchEnd;

    //These random values will determine how many child branches will be added (Fix realistic Radii with random child branches)
    bool makeBranch1 = true;//(rand() % 100) > 20;
    bool makeBranch2 = true;//(rand() % 100) > 20;
    bool makeBranch3 = true;//(rand() % 100) > 20;
    bool makeBranch4 = true;//(rand() % 100) > 20;
    bool makeBranch5 = true;//(rand() % 100) > 20;
    bool makeBranch6 = true;//(rand() % 100) > 20;
    bool makeBranch7 = true;//(rand() % 100) > 20;
    bool makeBranch8 = true;//(rand() % 100) > 20;

    //if(currentRecursionDepth == maxRecursionDepth){
    //    makeBranch1 = false;
    //    makeBranch2 = false;
    //    makeBranch3 = false;
    //    makeBranch4 = false;
    //    makeBranch5 = false;
    //    makeBranch6 = false;
    //    makeBranch7 = false;
    //    makeBranch8 = false;
    //}

    nextBranches.push(BranchDimensions(trunk, trunkLength));

    if(makeBranch1){
        nextBranches.push(BranchDimensions(branch1, newBranchLength));
    }
    if(makeBranch2){
        nextBranches.push(BranchDimensions(branch2, newBranchLength));
    }
    if(makeBranch3){
        nextBranches.push(BranchDimensions(branch3, newBranchLength));
    }
    if(makeBranch4){
        nextBranches.push(BranchDimensions(branch4, newBranchLength));
    }
    if(makeBranch5){
        nextBranches.push(BranchDimensions(branch5, newBranchLength));
    }
    if(makeBranch6){
        nextBranches.push(BranchDimensions(branch6, newBranchLength));
    }
    if(makeBranch7){
        nextBranches.push(BranchDimensions(branch7, newBranchLength));
    }
    if(makeBranch8){
        nextBranches.push(BranchDimensions(branch8, newBranchLength));
    }
    debugAssert(nextBranches.size() > 0);
}
