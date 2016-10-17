/** \file Mesh.cpp */
#include "Mesh.h"
#include <vector>

//typedef SmallArray<int, 3> Poly;

Mesh::Mesh(const String& fn) {
	filename = String(fn);
	vertices = std::vector<Point3>();
	faces = std::vector<Face>();
}

void Mesh::addVertex(float x, float y, float z) {
	vertices.push_back(Point3(x, y, z));
}

void Mesh::addVertex(Point3 vertex) {
	vertices.push_back(vertex);
}

void Mesh::addMesh(Mesh other){
    int initVertices = numVertices();
    for (int i = 0; i < other.vertices.size(); ++i) {
        addVertex(other.vertices[i]);
    }
    debugAssert(initVertices < numVertices());
    for (int i = 0; i < other.faces.size(); ++i) {
        Face face = other.faces[i];
        face.increment(initVertices);
        addFace(face);
    }
}

void Mesh::addFace(int f1, int f2, int f3, int t1, int t2, int t3, String type) {
	Face face(type);
	face.points.push(f1);
	face.points.push(f2);
	face.points.push(f3);
    face.txtPoints.push(t1);
    face.txtPoints.push(t2);
    face.txtPoints.push(t3);

	faces.push_back(face);
}

void Mesh::addFace(Face face){
    faces.push_back(face);
}

void Mesh::toOFF() { 
	int vertexCount = vertices.size();
	int faceCount = faces.size();
    String offFilename = filename + (String)".OFF";
	TextOutput tOut = TextOutput(offFilename);
	tOut.printf("OFF\n%d %d 0\n", vertexCount, faceCount);

	for (int i = 0; i < vertexCount; ++i) {
		Point3 vertex = vertices.at(i);
		tOut.printf("%f %f %f\n", vertex.x, vertex.y, vertex.z);
	}
	for (int i = 0; i < faceCount; ++i) {
		Face face = faces.at(i);
		tOut.printf("3 %d %d %d\n", face.points[0], face.points[1], face.points[2]);
	}
	tOut.commit();
}

void Mesh::toOBJ() {
    int vertexCount = vertices.size();
	int faceCount = faces.size();
    String objFilename = filename + (String)".OBJ";
	TextOutput tOut = TextOutput(objFilename);
	
    tOut.printf("mtllib " + filename + ".mtl\n");

	for (int i = 0; i < vertexCount; ++i) {
		Point3 vertex = vertices.at(i);
		tOut.printf("v %f %f %f\n", vertex.x, vertex.y, vertex.z);
	}
    
    // Add texture coordinates
    tOut.printf("vt 0.000000 0.000000\nvt 1.000000 0.000000\nvt 0.000000 1.000000\nvt 1.000000 1.000000\nvt 0.500000 0.000000\n");

	for (int i = 0; i < faceCount; ++i) {
		Face face = faces.at(i);
        String mtl = face.texture;
		tOut.printf("usemtl " + mtl + "\nf %d/%d %d/%d %d/%d\n", face.points[0] + 1, face.txtPoints[0], face.points[1] + 1, face.txtPoints[1], face.points[2] + 1, face.txtPoints[2]);
	}
	tOut.commit();
}

int Mesh::numVertices() {
    return vertices.size();
}


Face::Face(const String nTexture){
    texture = nTexture;
}

void Face::increment(int pts){
    for (int i = 0; i < points.size(); ++i) {
        points[i] += pts;
    }
}