#pragma once

#include "Mesh.hpp"
#include <Eigen/Core>
#include "MeshMatrix.hpp"

using namespace Eigen;
using namespace Jaal;

struct JSurfaceVectorField {
    JSurfaceVectorField()
    {
        numVecPerFace = 1;
    }

    void addConstraint( const JFacePtr &p, const Vec3D &v);

    Vec3F  getRandomVector( const JFacePtr &f);
    Eigen::VectorXd genRandomVector(const Eigen::VectorXd& b1,
                                    const Eigen::VectorXd& b2, int n);

    // Read constraints from the file. The first file constaints the
    // faceIDs and the second file contains the vectors on the faces..
    // Unfortunately, libigl separates this in two files :)...
    int readConstraints( const string &f1, const string &f2);

    JMeshPtr getVectorField()
    {
        return vecField;
    }

    JMeshPtr getConstrainedField();

    JNodeSequence getPositiveSingularNodes()
    {
        return positiveSingularNodes;
    }

    JNodeSequence getNegativeSingularNodes()
    {
        return negativeSingularNodes;
    }

    JFaceSequence getConstrainedFaces();

    int getNumSingularities() const
    {
        return negativeSingularNodes.size() + positiveSingularNodes.size();
    }

    int getNumConstraints() const
    {
        return constrainedFaces.size();
    }

    void setNumVectorsPerFace(int n)
    {
        numVecPerFace = n;
    }
    int  getNumVectorsPerFace() const
    {
        return numVecPerFace;
    }

    JMeshPtr mesh;
    JMeshPtr vecField;

    int      numVecPerFace;
    double   veclen;
    JNodeSequence positiveSingularNodes, negativeSingularNodes;
    JFaceSequence constrainedFaces;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

// Face barycenters
    Eigen::MatrixXd B;

//  Basis vector on each face of the triangle ...
    Eigen::MatrixXd B1, B2, B3;

// Vector of constrained faces
    Eigen::VectorXi b;

// Matrix of constraints
    Eigen::MatrixXd bc;

// Scale for visualizing the fields
    double global_scale;

// Random length factor
    double rand_factor = 5;

    int  readVecField();
    int  saveConstraints();
    void getConstraints();
    void addEdges(const JMeshPtr &m, const Eigen::MatrixXd &tail, const Eigen::MatrixXd &head);
    bool isPlanar(const JFacePtr &f, const Vec3D &v);
};

