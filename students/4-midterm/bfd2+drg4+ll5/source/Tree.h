#pragma once
#include <G3D/G3DAll.h>
/**
  \file Tree.h

   Simple implementation of a tree data structure.
 */

class Tree {
	protected:
        Point3 mContents;
        shared_ptr<Array<shared_ptr<Tree>>> mChildren;
        shared_ptr<Tree> mParent;
	public:
        Tree();
        Tree(Point3 contents, shared_ptr<Tree> parent);
        
        static shared_ptr<Tree> create(Point3 contents, shared_ptr<Tree> parent);
        void setContents(Point3 contents);
        Point3 getContents() const;

        void addChild(shared_ptr<Tree> child);
        shared_ptr<Array<shared_ptr<Tree>>> getChildren();

        void setParent(shared_ptr<Tree> parent);
        shared_ptr<Tree> getParent();
        
        int numChildren();
        
        bool operator<(const Tree &other) const;
};


