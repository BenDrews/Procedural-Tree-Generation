/** \file Mesh.cpp */
#include "Mesh.h"
#include <vector>

typedef SmallArray<int, 3> Poly;

Mesh::Mesh(const String& fn) {
	filename = String(fn);
	vertices = std::vector<Point3>();
	faces = std::vector<Poly>();
}

void Mesh::addVertex(float x, float y, float z) {
	vertices.push_back(Point3(x, y, z));
}

void Mesh::addFace(int f1, int f2, int f3) {
	Poly face = Poly::SmallArray();
	face.push(f1);
	face.push(f2);
	face.push(f3);
	faces.push_back(face);
}

void Mesh::toOFF() { 
	int vertexCount = vertices.size();
	int faceCount = faces.size();
	TextOutput tOut = TextOutput(filename);
	tOut.printf("OFF\n%d %d 0\n", vertexCount, faceCount);
	for(int i = 0; i < vertexCount; ++i) {
		Point3 vertex = vertices.at(i);
		tOut.printf("%f %f %f\n", vertex.x, vertex.y, vertex.z);
	}
	for(int i = 0; i < faceCount; ++i) {
		Poly face = faces.at(i);
		tOut.printf("3 %d %d %d\n", face[0], face[1], face[2]);
	}
	tOut.commit();
}

int Mesh::numVertices() {
    return vertices.size();
}