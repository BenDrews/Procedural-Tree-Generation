/** \file BranchDimensions.cpp */
#include "BranchDimensions.h"

BranchDimensions::BranchDimensions() {
    frame = CoordinateFrame();
    length = 0.0f;
}

BranchDimensions::BranchDimensions(CoordinateFrame f, float len) {
    frame = f;
    length = len;
};