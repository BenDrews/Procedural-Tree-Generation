#pragma once
#include <G3D/G3DAll.h>
/**
  \file FruitDimensions.h

  Data structure used to hold information about the filename, size, and placement of a fruit.
 */

class FruitDimensions {
	public:
        String filename;
        float scale;
		float yOffset;
		float roll;

        FruitDimensions();
        FruitDimensions(String filename, float scale, float yOffset, float roll);
};