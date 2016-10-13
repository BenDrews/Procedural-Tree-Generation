#pragma once
#include <G3D/G3DAll.h>
#include <vector>
/**
  \file BranchDimensions.h

  Data structure used to hold information about the size of a branch.
 */

class BranchDimensions {
	protected:
	public:
        CoordinateFrame frame;
        float length;
        BranchDimensions();
        BranchDimensions(CoordinateFrame frame, float length);
        //CoordinateFrame frame();
        //float length();
};