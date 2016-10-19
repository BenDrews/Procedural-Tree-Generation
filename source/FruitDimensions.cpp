/** \file FruitDimensions.cpp */
#include "FruitDimensions.h"


FruitDimensions::FruitDimensions() {
	filename = "";
	scale = 0.0f;
	yOffset = 0.0f;
}

FruitDimensions::FruitDimensions(String fn, float s, float yOff) {
	filename = fn;
	scale = s;
	yOffset = yOff;	
};