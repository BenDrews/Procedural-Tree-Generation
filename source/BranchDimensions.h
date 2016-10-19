#pragma once
#include <G3D/G3DAll.h>
/**
  \file BranchDimensions.h

  Data structure used to hold information about the size of a branch.
 */

class BranchDimensions {
	public:
        CoordinateFrame frame;
        float length;
        BranchDimensions();
        BranchDimensions(CoordinateFrame frame, float length);
};