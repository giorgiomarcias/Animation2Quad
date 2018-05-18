/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef SMOOTHERWRAPPER_H
#define SMOOTHERWRAPPER_H

#include <vector>
#include <list>
#include <utility>
#include <Eigen/Dense>
#include <igl/copyleft/comiso/nrosy.h>

using namespace std;
using namespace igl::copyleft::comiso;

template < typename TriMeshType >
class SmootherWrapper{
public:
    // define types
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::ScalarType ScalarType;

    /**
     * @brief SmootherWrapper
     */
    SmootherWrapper() {
        smoother = 0;
    }

    /**
     * @brief ~SmootherWrapper
     */
    ~SmootherWrapper() {
        reset();
    }

    /**
     * @brief reset
     */
    inline void reset() {
        if (smoother)
            delete smoother;
        smoother = 0;
    }

    /**
     * @brief initializeSmoother
     * @param mesh
     */
    void initializeSmoother(TriMeshType &mesh) {
        // reset
        reset();

        // create eigen matrix of vertices
        Eigen::MatrixXd v(mesh.VN(), 3);

        // copy vertices
        for (int i = 0; i < mesh.VN(); i++)
            for (int j = 0; j < 3; j++)
                v(i,j) = mesh.vert[i].cP()[j];

        // create eigen matrix of faces
        Eigen::MatrixXi f(mesh.FN(), 3);

        // copy faces
        VertexType *v0 = &mesh.vert[0];
        for (int i = 0; i < mesh.FN(); i++)
            for (int j = 0; j < 3; j++) {
                f(i,j) = (int)(mesh.face[i].V(j) - v0);
                assert(f(i,j) >= 0 && f(i,j) < mesh.VN());
            }

        // create the smoother
        smoother = new NRosyField(v, f);
    }

    /**
     * @brief setSoftAlpha
     * @param alpha
     */
    void setSoftAlpha(double alpha) {
        assert(smoother);

        // set smoother alpha
        smoother->setSoftAlpha(alpha);
    }

    /**
     * @brief setConstraintSoft
     * @param fid
     * @param w
     * @param v
     */
    void setConstraintSoft(const int fid, const double w, const CoordType &v) {
        assert(smoother);

        // create eigen vector
        Eigen::Vector3d c;

        // copy coordinates
        for (int i = 0; i < 3; i++)
            c(i) = v[i];

        // set smoother soft constraint
        smoother->setConstraintSoft(fid, w, c);
    }

    /**
     * @brief setConstraintHard
     * @param fid
     * @param v
     */
    void setConstraintHard(const int fid, const CoordType &v) {
        assert(smoother);

        // create eigen vector
        Eigen::Vector3d c;

        // copy coordinates
        for (int i = 0; i < 3; i++)
            c(i) = v[i];

        // set smoother hard constraint
        smoother->setConstraintHard(fid, c);
    }

    /**
     * @brief resetConstraints
     */
    void resetConstraints() {
        assert(smoother);

        // reset constraints
        smoother->resetConstraints();
    }

    /**
     * @brief solve
     * @param N
     */
    void solve(const int N = 4) {
        assert(smoother);

        // solve
        smoother->solve(N);
    }

    /**
     * @brief getFieldPerFace
     * @param fields
     */
    void getFieldPerFace(vector<CoordType> &fields) {
        assert(smoother);

        // get fields
        Eigen::MatrixXd fs = smoother->getFieldPerFace();

        // reset output vector of fields
        fields.clear();

        // resize output vector of fields
        fields.resize(fs.rows());

        // copy fields in output vector
        for (int i = 0; i < fs.rows(); i++)
            for (int j = 0; j < 3; j++)
                fields[i][j] = fs(i,j);
    }

    /**
     * @brief getFFieldPerFace
     * @param ffields
     */
    void getFFieldPerFace(vector< pair<CoordType,CoordType> > &ffields) {
        assert(smoother);

        // get fields
        Eigen::MatrixXd ffs = smoother->getFFieldPerFace();

        // reset output vector of fields
        ffields.clear();

        // resize output vector of fields
        ffields.resize(ffs.rows());

        // copy fields in output vector
        for (int i = 0; i < ffs.rows(); i++)
            for (int j = 0; j < 3; j++) {
                ffields[i].first[j] = ffs(i,j);
                ffields[i].second[j] = ffs(i,3+j);
            }
    }

    /**
     * @brief findCones
     * @param N
     */
    void findCones(int N = 4) {
        assert(smoother);

        // compute singularity indexes
        smoother->findCones(N);
    }

    /**
     * @brief getSingularityIndexPerVertex
     * @param sings
     */
    void getSingularityIndexPerVertex(vector<ScalarType> &sings) {
        assert(smoother);

        // get fields
        Eigen::VectorXd s = smoother->getSingularityIndexPerVertex();

        // reset output vector of singularities
        sings.clear();

        // resize output vector of singularities
        sings.resize(s.rows());

        // copy fields in output vector
        for (int i = 0; i < s.rows(); i++)
            sings[i] = s(i);
    }

    /**
     * @brief getSingularitiesIndexPerVertexList
     * @param sIndexes
     * @param t
     */
    void getSingularitiesIndexPerVertexList(list<int> &sIndexes, const ScalarType t = ScalarType(0)) {
        assert(smoother);

        // get singularities vector
        vector<ScalarType> sings;
        getSingularityIndexPerVertex(sings);

        // reset output list of indexes
        sIndexes.clear();

        // search for indexes with singularty greater than or equal to t
        for (unsigned long i = 0; i < sings.size(); i++)
            if (sings[i] >= t)
                sIndexes.push_back(i);
    }

private:
    NRosyField *smoother;
};

#endif // SMOOTHERWRAPPER_H
