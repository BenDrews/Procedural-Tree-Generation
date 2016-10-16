#pragma once
#include <G3D/G3DAll.h>
#include <vector>
/**
  \file Mesh.h

  Data structure used to hold information about a mesh as it is constructed.
 */
class Face {
    public:
        Face(const String texture);
        void increment(int pts);
        String texture;
        SmallArray<int,3> points;
        SmallArray<int,3> txtPoints;

};

class Mesh {
	protected:
		String filename;
		std::vector<Point3> vertices;
		std::vector<Face> faces;
	public:
		void addVertex(float x, float y, float z);
        void addVertex(Point3 vertex);
		void addFace(int f1, int f2, int f3, int t1, int t2, int t3, String type);
        void addFace(Face);
        void addMesh(Mesh other);
		Mesh(const String& fn);
		void toOFF();
        void toOBJ();
        int numVertices();
};

