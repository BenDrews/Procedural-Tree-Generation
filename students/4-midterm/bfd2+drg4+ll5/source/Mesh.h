#pragma once
#include <G3D/G3DAll.h>
#include <vector>
/**
  \file Mesh.h

  Data structure used to hold information about a mesh as it is constructed.
 */

class Mesh {
	protected:
		String filename;
		std::vector<Point3> vertices;
		std::vector<SmallArray<int, 3>> faces;
	public:
		void addVertex(float x, float y, float z);
		void addFace(int f1, int f2, int f3);
		Mesh(const String& fn);
		void toOFF();
        int numVertices();
};