/** \file FruitDimensions.cpp */
#include "FruitDimensions.h"


FruitDimensions::FruitDimensions() {
	filename = "";
	scale = 0.0f;
	yOffset = 0.0f;
	roll = 0.0f;
}

FruitDimensions::FruitDimensions(String fn, float s, float yOff, float r) {
	filename = fn;
	scale = s;
	yOffset = yOff;
	roll = r;
};