/** \file BranchDimensions.cpp */
#include "BranchDimensions.h"

BranchDimensions::BranchDimensions() {
    frame = CoordinateFrame();
    length = 0.0f;
}

BranchDimensions::BranchDimensions(CoordinateFrame f, float len) {
    frame = f;
    length = len;
}

//CoordinateFrame BranchDimensions::frame() {
//    return frame;
//}
//
//float BranchDimensions::length() {
//    return length;
//}