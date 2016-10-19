#include "Tree.h"
/**
  \file Tree.cpp

   Simple implementation of a tree data structure.
 */

shared_ptr<Tree> Tree::create(CoordinateFrame contents, shared_ptr<Tree> parent) {
    return std::make_shared<Tree>(contents, parent);
}

Tree::Tree() {
    mChildren = std::make_shared<Array<shared_ptr<Tree>>>();
    mParent = nullptr;
    mContents = CoordinateFrame();
}

Tree::Tree(CoordinateFrame contents, shared_ptr<Tree> parent) {
    mChildren = std::make_shared<Array<shared_ptr<Tree>>>();
    mParent = parent;
    mContents = contents;
}

void Tree::setContents(CoordinateFrame contents) {
    mContents = contents;
}

CoordinateFrame Tree::getContents() const{
    return mContents;
}

void Tree::addChild(shared_ptr<Tree> child) {
    mChildren->push(child);
}

shared_ptr<Array<shared_ptr<Tree>>> Tree::getChildren() {
    return mChildren;
}

void Tree::setParent(shared_ptr<Tree> parent) {
    mParent = parent;
}

shared_ptr<Tree> Tree::getParent() {
    return mParent;
}

int Tree::numChildren(){ 
    return mChildren->size();
}

bool Tree::operator<(const Tree &other) const {
    return (mContents.translation.x) < (other.getContents().translation.x);
}