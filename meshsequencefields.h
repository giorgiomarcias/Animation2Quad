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

#ifndef MESH_SEQUENCEFIELDS_H
#define MESH_SEQUENCEFIELDS_H

#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include <string>
#include <sstream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
using namespace std;

#include <vcg/space/point3.h>
#include <vcg/space/point2.h>
//#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>

#include "compute.h"

template < typename ScalarType >
class MeshSequenceFields
{
public:
    /// type defines
    typedef vcg::Point3<ScalarType> FieldType;
    typedef ScalarType              ValueType;

    // type of stretch to consider when averaging
    typedef enum{
        SINGVAL1,   // s1
        SINGVAL2,   // s2
        MAXWARP,    // max(|s1|, |1/s2|)
        RATIO       // |s1 / s2|
    } StretchType;

    /**
     * @brief MeshSequenceFields
     */
    MeshSequenceFields() {
        numOfMeshes = 0;
        toLog = 0;
        minAbsSequenceSingularValue1 = minAbsSequenceSingularValue2 = minAbsSequenceStretch1 = minAbsSequenceStretch2
                = numeric_limits<ValueType>::max();
        maxAbsSequenceSingularValue1 = maxAbsSequenceSingularValue2 = maxAbsSequenceStretch1 = maxAbsSequenceStretch2
                = -numeric_limits<ValueType>::max();
        singularValue2LowerBound = std::numeric_limits<ValueType>::epsilon();
        singularValue1UpperBound = numeric_limits<ValueType>::max();
    }

    /**
     * @brief clearAt
     * @param meshIndex
     */
    void clearAt(unsigned long meshIndex) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        // clear vectors at meshIndex
        fields[meshIndex].clear();
        fields2D[meshIndex].clear();
#ifdef NOT_NORMALIZE
        fields2Dmap[meshIndex].clear();
#endif
        singularValues[meshIndex].clear();
        minAbsSingularValues[meshIndex] = pair<ValueType,ValueType>(numeric_limits<ValueType>::max(),
                                                                    numeric_limits<ValueType>::max());
        maxAbsSingularValues[meshIndex] = pair<ValueType,ValueType>(-numeric_limits<ValueType>::max(),
                                                                    -numeric_limits<ValueType>::max());
        stretch[meshIndex].clear();
        minStretch[meshIndex] = pair<ValueType,ValueType>(numeric_limits<ValueType>::max(),
                                                          numeric_limits<ValueType>::max());
        maxStretch[meshIndex] = pair<ValueType,ValueType>(-numeric_limits<ValueType>::max(),
                                                          -numeric_limits<ValueType>::max());
    }

    /**
     * @brief clear
     */
    void clear() {
        // clear all vectors
        fields.clear();
        fields2D.clear();
#ifdef NOT_NORMALIZE
        fields2Dmap.clear();
#endif
        singularValues.clear();
        minAbsSingularValues.clear();
        maxAbsSingularValues.clear();
        stretch.clear();
        minStretch.clear();
        maxStretch.clear();
        numOfMeshes = 0;
        minAbsSequenceSingularValue1 = minAbsSequenceSingularValue2 = minAbsSequenceStretch1 = minAbsSequenceStretch2
                = numeric_limits<ValueType>::max();
        maxAbsSequenceSingularValue1 = maxAbsSequenceSingularValue2 = maxAbsSequenceStretch1 = maxAbsSequenceStretch2
                = -numeric_limits<ValueType>::max();
    }

    /**
     * @brief resizeFieldsAt
     * @param meshIndex
     * @param num
     */
    void resizeFieldsAt(unsigned long meshIndex, unsigned long num) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        // resize vector
        fields[meshIndex].resize(num);
        fields2D[meshIndex].resize(num);
#ifdef NOT_NORMALIZE
        fields2Dmap[meshIndex].resize(num);
#endif
        singularValues[meshIndex].resize(num);
        minAbsSingularValues[meshIndex] = pair<ValueType,ValueType>(numeric_limits<ValueType>::max(),
                                                                    numeric_limits<ValueType>::max());
        maxAbsSingularValues[meshIndex] = pair<ValueType,ValueType>(-numeric_limits<ValueType>::max(),
                                                                    -numeric_limits<ValueType>::max());
        stretch[meshIndex].resize(num);
        minStretch[meshIndex] = pair<ValueType,ValueType>(numeric_limits<ValueType>::max(),
                                                          numeric_limits<ValueType>::max());
        maxStretch[meshIndex] = pair<ValueType,ValueType>(-numeric_limits<ValueType>::max(),
                                                          -numeric_limits<ValueType>::max());
    }

    /**
     * @brief resizeNumOfMeshes
     * @param num
     */
    void resizeNumOfMeshes(unsigned long num) {
        // resize vectors
        fields.resize(num);
        fields2D.resize(num);
#ifdef NOT_NORMALIZE
        fields2Dmap.resize(num);
#endif
        singularValues.resize(num);
        minAbsSingularValues.resize(num);
        maxAbsSingularValues.resize(num);
        stretch.resize(num);
        minStretch.resize(num);
        maxStretch.resize(num);
        numOfMeshes = num;
    }

    /**
     * @brief getNumberOfMeshes
     * @return
     */
    inline unsigned long getNumberOfMeshes() {
        return numOfMeshes;
    }

    /**
     * @brief getNumberOfFieldsAt
     * @param meshIndex
     * @return
     */
    inline unsigned long getNumberOfFieldsAt(unsigned long meshIndex) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        return fields[meshIndex].size();
    }

    /**
     * @brief getFirstFieldAt
     * @param meshIndex
     * @param faceIndex
     * @param field
     */
    inline void getFirstFieldAt(unsigned long meshIndex, unsigned long faceIndex, FieldType & field) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        // check parameter
        assert(faceIndex < fields[meshIndex].size());

        // copy field
        field = fields[meshIndex][faceIndex].first;
    }

    /**
     * @brief getSecondFieldAt
     * @param meshIndex
     * @param faceIndex
     * @param field
     */
    inline void getSecondFieldAt(unsigned long meshIndex, unsigned long faceIndex, FieldType & field) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        // check parameter
        assert(faceIndex < fields[meshIndex].size());

        // copy field
        field = fields[meshIndex][faceIndex].second;
    }

    /**
     * @brief getStretchMeasureAt
     * @param meshIndex
     * @param faceIndex
     * @param stretchType
     * @return
     */
    inline ValueType getStretchMeasureAt(unsigned long meshIndex, unsigned long faceIndex, const StretchType stretchType) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        // check parameter
        assert(faceIndex < stretch[meshIndex].size());

        switch(stretchType) {
        case SINGVAL1:
            return singularValues[meshIndex][faceIndex].first;

        case SINGVAL2:
            return singularValues[meshIndex][faceIndex].second;

        case MAXWARP:
            return stretch[meshIndex][faceIndex].first;

        case RATIO:
            return stretch[meshIndex][faceIndex].second;

        default:
            assert(false);
            return stretch[meshIndex][faceIndex].first;
        }
    }

    /**
     * @brief getMinAbsStretchMeasureAt
     * @param meshIndex
     * @param stretchType
     * @return
     */
    inline ValueType getMinAbsStretchMeasureAt(unsigned long meshIndex, const StretchType stretchType) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        switch(stretchType) {
        case SINGVAL1:
            return minAbsSingularValues[meshIndex].first;

        case SINGVAL2:
            return minAbsSingularValues[meshIndex].second;

        case MAXWARP:
            return minStretch[meshIndex].first;

        case RATIO:
            return minStretch[meshIndex].second;

        default:
            assert(false);
            return -numeric_limits<ValueType>::max();
        }
    }

    /**
     * @brief getMaxAbsStretchMeasureAt
     * @param meshIndex
     * @param stretchType
     * @return
     */
    inline ValueType getMaxAbsStretchMeasureAt(unsigned long meshIndex, const StretchType stretchType) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        switch(stretchType) {
        case SINGVAL1:
            return maxAbsSingularValues[meshIndex].first;

        case SINGVAL2:
            return maxAbsSingularValues[meshIndex].second;

        case MAXWARP:
            return maxStretch[meshIndex].first;

        case RATIO:
            return maxStretch[meshIndex].second;

        default:
            assert(false);
            return numeric_limits<ValueType>::max();
        }
    }

    /**
     * @brief getMinAbsSequenceStretchMeasure
     * @param stretchType
     * @return
     */
    inline ValueType getMinAbsSequenceStretchMeasure(const StretchType stretchType) {
        switch(stretchType) {
        case SINGVAL1:
            return minAbsSequenceSingularValue1;

        case SINGVAL2:
            return minAbsSequenceSingularValue2;

        case MAXWARP:
            return minAbsSequenceStretch1;

        case RATIO:
            return minAbsSequenceStretch2;

        default:
            assert(false);
            return -numeric_limits<ValueType>::max();
        }
    }

    /**
     * @brief getMaxAbsSequenceStretchMeasure
     * @param stretchType
     * @return
     */
    inline ValueType getMaxAbsSequenceStretchMeasure(const StretchType stretchType) {
        switch(stretchType) {
        case SINGVAL1:
            return maxAbsSequenceSingularValue1;

        case SINGVAL2:
            return maxAbsSequenceSingularValue2;

        case MAXWARP:
            return maxAbsSequenceStretch1;

        case RATIO:
            return maxAbsSequenceStretch2;

        default:
            assert(false);
            return numeric_limits<ValueType>::max();
        }
    }

    /**
     * @brief setAbsSingularValue2LowerBound
     * @param lowerBound
     */
    inline void setAbsSingularValue2LowerBound(const ValueType lowerBound) {
        assert(lowerBound >= ValueType(0));
        singularValue2LowerBound = lowerBound;
    }

    /**
     * @brief setAbsSingularValue1UpperBound
     * @param upperBound
     */
    inline void setAbsSingularValue1UpperBound(const ValueType upperBound) {
        assert(upperBound >= ValueType(0));
        singularValue1UpperBound = upperBound;
    }

    /**
     * @brief getAbsStretchMeasureLowerBound
     * @param stretchType
     * @return
     */
    inline ValueType getAbsStretchMeasureLowerBound(const StretchType stretchType) {
        switch(stretchType) {
        case SINGVAL1:
        case SINGVAL2:
            return singularValue2LowerBound;

        case MAXWARP:
            return std::max(singularValue2LowerBound, ValueType(1) / singularValue2LowerBound);

        case RATIO:
            return ValueType(0);

        default:
            assert(false);
            return ValueType(0);
        }
    }

    /**
     * @brief getAbsStretchMeasureUpperBound
     * @param stretchType
     * @return
     */
    inline ValueType getAbsStretchMeasureUpperBound(const StretchType stretchType) {
        switch(stretchType) {
        case SINGVAL1:
        case SINGVAL2:
            return singularValue1UpperBound;

        case MAXWARP:
            return std::max(singularValue1UpperBound, ValueType(1) / singularValue2LowerBound);

        case RATIO:
            return singularValue1UpperBound / singularValue2LowerBound - ValueType(1);

        default:
            assert(false);
            return numeric_limits<ValueType>::max();
        }
    }

    /**
     * @brief getClampedAbsStretchMeasure
     * @param value
     * @param stretchType
     * @return
     */
    inline ValueType getClampedAbsStretchMeasure(ValueType value, const StretchType stretchType) {
        // check parameter
        // debug:
        assert(value >= ValueType(0));
        if (value < ValueType(0))
            value = -value;

        // check if less than lower bound
        if (value < getAbsStretchMeasureLowerBound(stretchType))
            return getAbsStretchMeasureLowerBound(stretchType);

        // check if greater than upper bound
        if (value > getAbsStretchMeasureUpperBound(stretchType))
            return getAbsStretchMeasureUpperBound(stretchType);

        // otherwise it is into interval
        return value;
    }

    /**
     * @brief getClampedAbsStretchMeasureAt
     * @param meshIndex
     * @param faceIndex
     * @param stretchType
     * @return
     */
    inline ValueType getClampedAbsStretchMeasureAt(unsigned long meshIndex, unsigned long faceIndex, const StretchType stretchType) {
        // check parameter
        assert(meshIndex < numOfMeshes);

        // check parameter
        assert(faceIndex < stretch[meshIndex].size());

        // temporary variables
        ValueType s1, s2, s;
        switch(stretchType) {
        case SINGVAL1:
            return getClampedAbsStretchMeasure(singularValues[meshIndex][faceIndex].first, stretchType);

        case SINGVAL2:
            return getClampedAbsStretchMeasure(singularValues[meshIndex][faceIndex].second, stretchType);

        case MAXWARP:
            s1 = getClampedAbsStretchMeasure(singularValues[meshIndex][faceIndex].first, SINGVAL1);
            s2 = getClampedAbsStretchMeasure(singularValues[meshIndex][faceIndex].second, SINGVAL2);
            s = std::max(std::fabs(s1), std::fabs(ValueType(1) / s2));
            return getClampedAbsStretchMeasure(s, stretchType);

        case RATIO:
            s1 = getClampedAbsStretchMeasure(singularValues[meshIndex][faceIndex].first, SINGVAL1);
            s2 = getClampedAbsStretchMeasure(singularValues[meshIndex][faceIndex].second, SINGVAL2);
            assert(s1 >= s2);
            s = std::fabs(s1 / s2) - ValueType(1);
            return getClampedAbsStretchMeasure(s, stretchType);

        default:
            assert(false);
            return getAbsStretchMeasureLowerBound(stretchType);
        }
    }

    /**
     * @brief computeMeshFieldsAt
     * @param index
     * @param restMesh
     * @param referenceMesh
     */
    template < typename TriMeshType >
    void computeMeshFieldsAt(unsigned long index,
                             const TriMeshType &restMesh,
                             const TriMeshType &referenceMesh) {
        // check parameter
        assert(index < numOfMeshes);
        assert(restMesh.FN() == referenceMesh.FN());

        // check the size of both triangle meshes
        if (restMesh.FN() != referenceMesh.FN()) {
            if (toLog) {
                ostringstream numtostr;
                numtostr << index;
                toLog(string("Error. In ") + numtostr.str() + string(" rest mesh and reference mesh have not the same number of faces.\n"));
            }
            return;
        }

        // resize vector
        resizeFieldsAt(index, (unsigned long)restMesh.FN());

        // create matrices A, U, V, P
        Eigen::Matrix<ValueType,2,2> A, U, V, P;
        // create vector S
        Eigen::Matrix<ValueType,2,1> S;

        // temporary data
        FieldType triangle1[3], triangle2[3], base1[3], base2[3];
        FieldType field, origin;
        vcg::Point2<ScalarType> triangle1_2d[2], triangle2_2d[2];
        ValueType den1, den2;

        // for each triangle
        for (int i = 0; i < restMesh.FN(); i++) {
            // create triangles
            for (int j = 0; j < 3; j++) {
                triangle1[j] = restMesh.face[i].cP(j);
                triangle2[j] = referenceMesh.face[i].cP(j);
            }

            // centrate
            origin = triangle1[0];
            for (int j = 0; j < 3; j++)
                triangle1[j] -= origin;
            origin = triangle2[0];
            for (int j = 0; j < 3; j++)
                triangle2[j] -= origin;

            // new reference base for t1
            den1 = (triangle1[1] - triangle1[0]).Norm();
            base1[0] = (triangle1[1] - triangle1[0]) / den1;
            base1[2] = restMesh.face[i].cN();
            base1[1] = base1[2] ^ base1[0];

            // express t1 in new reference
            triangle1_2d[0] = vcg::Point2<ScalarType>(den1 ,0);
            triangle1_2d[1] = vcg::Point2<ScalarType>((triangle1[2] - triangle1[0]) * base1[0], (triangle1[2] - triangle1[0]) * base1[1]);

            // new reference base for t2
            den2 = (triangle2[1] - triangle2[0]).Norm();
            base2[0] = (triangle2[1] - triangle2[0]) / den2;
            base2[2] = referenceMesh.face[i].cN();
            base2[1] = base2[2] ^ base2[0];

            // express t2 in new reference
            triangle2_2d[0] = vcg::Point2<ScalarType>(den2, 0);
            triangle2_2d[1] = vcg::Point2<ScalarType>((triangle2[2] - triangle2[0]) * base2[0], (triangle2[2] - triangle2[0]) * base2[1]);

            // compute map A (in 2D)
            getTriangle2TriangleMap(triangle1_2d, triangle2_2d, A);

            // compute P
            U <<    A(0,0)-A(1,1), A(0,1)+A(1,0),
                    A(1,0)+A(0,1), A(1,1)-A(0,0);
            U = U / 2;
            ScalarType det = std::sqrt(std::fabs(U.determinant())); // it should be greater than 0
            U = U / (det > 0.0 ? det : 1.0);
            P = U.transpose() * A;

            // compute svd (in 2D) of P
            getSVD(P, U, V, S);

            // copy 2D fields
            fields2D[index][i] = U.col(0);
#ifdef NOT_NORMALIZE
            // non-normalized:
            fields2Dmap[index][i] = A;
#endif

            // remap first field to global reference
            field = base1[0] * U(0,0) + base1[1] * U(1,0);

            // store first field
            fields[index][i].first = field;

            // store first singular value
            singularValues[index][i].first = S(0);

            // update min and max
            if (std::fabs(S(0)) < minAbsSingularValues[index].first)
                minAbsSingularValues[index].first = std::fabs(S(0));
            if (std::fabs(S(0)) > maxAbsSingularValues[index].first)
                maxAbsSingularValues[index].first = std::fabs(S(0));
            if (std::fabs(S(0)) < minAbsSequenceSingularValue1)
                minAbsSequenceSingularValue1 = std::fabs(S(0));
            if (std::fabs(S(0)) > maxAbsSequenceSingularValue1)
                maxAbsSequenceSingularValue1 = std::fabs(S(0));

            // remap second field to global reference
            field = base1[0] * U(0,1) + base1[1] * U(1,1);

            // store second field
            fields[index][i].second = field;
            singularValues[index][i].second = S(1);

            // update min and max
            if (std::fabs(S(1)) < minAbsSingularValues[index].second)
                minAbsSingularValues[index].second = std::fabs(S(1));
            if (std::fabs(S(1)) > maxAbsSingularValues[index].second)
                maxAbsSingularValues[index].second = std::fabs(S(1));
            if (std::fabs(S(1)) < minAbsSequenceSingularValue2)
                minAbsSequenceSingularValue2 = std::fabs(S(1));
            if (std::fabs(S(1)) > maxAbsSequenceSingularValue2)
                maxAbsSequenceSingularValue2 = std::fabs(S(1));

            // store first type of stretch
            stretch[index][i].first = std::max(std::fabs(singularValues[index][i].first),
                                               singularValues[index][i].second != ValueType(0) ?
                                                    std::fabs((float)1.0 / singularValues[index][i].second) :
                                                    ValueType(0));
            // store second type of stretch
            stretch[index][i].second = std::fabs(singularValues[index][i].first /
                                                    (singularValues[index][i].second != ValueType(0) ?
                                                        singularValues[index][i].second :
                                                        ValueType(1))) - ValueType(1);

            // update min first type of stretch
            if (stretch[index][i].first < minStretch[index].first)
                minStretch[index].first = stretch[index][i].first;
            if (stretch[index][i].first < minAbsSequenceStretch1)
                minAbsSequenceStretch1 = stretch[index][i].first;

            // update min second type of stretch
            if (stretch[index][i].second < minStretch[index].second)
                minStretch[index].second = stretch[index][i].second;
            if (stretch[index][i].second < minAbsSequenceStretch2)
                minAbsSequenceStretch2 = stretch[index][i].second;

            // update max first type of stretch
            if (stretch[index][i].first > maxStretch[index].first)
                maxStretch[index].first = stretch[index][i].first;
            if (stretch[index][i].first > maxAbsSequenceStretch1)
                maxAbsSequenceStretch1 = stretch[index][i].first;

            // update min second type of stretch
            if (stretch[index][i].second > maxStretch[index].second)
                maxStretch[index].second = stretch[index][i].second;
            if (stretch[index][i].second > maxAbsSequenceStretch2)
                maxAbsSequenceStretch2 = stretch[index][i].second;
        }

        // log
        if (toLog) {
            ostringstream numtostr;
            numtostr << index;
            toLog(string("Fields for reference mesh ") + numtostr.str() + string(" computed.\n"));
        }
    }

    /**
     * @brief computeAverageFieldsAt
     * @param restIndex
     * @param restMesh
     * @param stretchType
     * @param threshold
     */
    template < typename TriMeshType >
    void computeAverageFieldsAt(unsigned long restIndex, const TriMeshType &restMesh,
                                const StretchType stretchType, const ValueType threshold = ValueType(0)) {
        // check parameters
        assert(restIndex < numOfMeshes && restMesh.FN() > 0);
        assert(threshold >= getAbsStretchMeasureLowerBound(stretchType)
               && threshold <= getAbsStretchMeasureUpperBound(stretchType));

        // clear rest mesh fields
        clearAt(restIndex);

        // resize rest mesh fields vector
        resizeFieldsAt(restIndex, restMesh.FN());

        // temporary variables
//        vector< FieldType > tangentsVector(numOfMeshes - 1);
//        vector< ValueType > weigthsVector(numOfMeshes - 1);
//        vector< FieldType > normalsVector(numOfMeshes - 1);
//        FieldType baseNorm;
//        FieldType baseDir;
        Eigen::Matrix<ValueType,Eigen::Dynamic,2> V;
        Eigen::Matrix<ValueType,Eigen::Dynamic,1> w;
        list< Eigen::Matrix<ValueType,2,1> > VList;
#ifdef NOT_NORMALIZE
        list< Eigen::Matrix<ValueType,2,2> > VListmap;
#endif
        list<ValueType> wList;
        Eigen::Matrix<ValueType,2,1> R;
        FieldType triangle[3], base[3];
        FieldType origin;
        ValueType den, den1, den2, tmp;
        int i;
//        ValueType min = getMinAbsSequenceStretchMeasure(stretchType);
//        min = getClampedAbsStretchMeasure(min, stretchType);
//        ValueType max = getMaxAbsSequenceStretchMeasure(stretchType);
//        max = getClampedAbsStretchMeasure(max, stretchType);
        ValueType min = getAbsStretchMeasureLowerBound(stretchType);
        ValueType max = getAbsStretchMeasureUpperBound(stretchType);

        // log
        if (toLog) {
            ostringstream numtostr;
            numtostr << restIndex;
            toLog(string("Averaging fields for rest mesh ") + numtostr.str() + string(" ..."));
        }

        // for each face
        for (unsigned long findex = 0; findex < (unsigned long)restMesh.FN(); findex++) {
//            // construct temporary variables
//            baseNorm = restMesh.face[findex].cN();
//            baseDir = (restMesh.face[findex].cP(1) - restMesh.face[findex].cP(0)).Normalize();

            // reset variables
            singularValues[restIndex][findex] = pair<ValueType,ValueType>(ValueType(0), ValueType(0));
            stretch[restIndex][findex] = pair<ValueType,ValueType>(ValueType(0), ValueType(0));
            den1 = ValueType(0);
            den2 = ValueType(0);
            VList.clear();
            wList.clear();
#ifdef NOT_NORMALIZE
            VListmap.clear();
#endif

            // for all other meshes
            for (unsigned long mindex = 0; mindex < numOfMeshes; mindex++) {
                // of course, skip current mesh index (rest pose)
                if (mindex == restIndex)
                    continue;

                // check same number of fields
                assert(findex < fields[mindex].size());

//                // update temporary vectors
//                tangentsVector[i] = fields[mindex][findex].first;
//                weigthsVector[i] = (getClampedAbsStretchMeasureAt(mindex, findex, stretchType) - min) / (max - min);

                // threshold
                if (getClampedAbsStretchMeasureAt(mindex, findex, stretchType) >= threshold) {
                    // copy fields
                    VList.push_back(fields2D[mindex][findex]);
#ifdef NOT_NORMALIZE
                    VListmap.push_back(fields2Dmap[mindex][findex]);
#endif

                    // copy stretch as weight (normalized)
                    wList.push_back((getClampedAbsStretchMeasureAt(mindex, findex, stretchType) - min) / (max - min));

//                    normalsVector[i] = baseNorm;

                    // sum singular values
                    singularValues[restIndex][findex].first += getClampedAbsStretchMeasureAt(mindex, findex, SINGVAL1);
                    singularValues[restIndex][findex].second += getClampedAbsStretchMeasureAt(mindex, findex, SINGVAL2);

                    // (weighting) sum stretch values of first type
                    tmp = getClampedAbsStretchMeasureAt(mindex, findex, MAXWARP);
                    den1 += tmp;
                    tmp *= tmp;
                    stretch[restIndex][findex].first += tmp;

                    // (weighting) sum stretch values of second type
                    tmp = getClampedAbsStretchMeasureAt(mindex, findex, RATIO);
                    den2 += tmp;
                    tmp *= tmp;
                    stretch[restIndex][findex].second += tmp;
                }
            }

            // copy 2D fields
            if (VList.size() > 0) {
                i = 0;
                V.resize(VList.size(),2);
#ifdef NOT_NORMALIZE
                for (typename list< Eigen::Matrix<ValueType,2,2> >::iterator it = VListmap.begin(); it != VListmap.end(); it++) {
                    V(i,0) = (*it)(0,0);
                    V(i,1) = (*it)(0,1);
                    i++;
                }
#else
                for (typename list< Eigen::Matrix<ValueType,2,1> >::iterator it = VList.begin(); it != VList.end(); it++) {
                    V(i,0) = (*it)(0);
                    V(i,1) = (*it)(1);
                    i++;
                }
#endif

                // copy weights
                i = 0;
                w.resize(wList.size());
                for (typename list<ValueType>::iterator it = wList.begin(); it != wList.end(); it++)
                    w(i++) = *it;

                // average vectors (in 2D)
#ifdef NOT_NORMALIZE
                R.setZero();
                for (int j = 0; j < V.rows(); j++)
                  R += w(j) * V.row(j).transpose();
                R /= w.sum();
#else
                R = averageVectors(V, w, 4);
#endif

                // copy field in 2D
                fields2D[restIndex][findex](0) = R(0);
                fields2D[restIndex][findex](1) = R(1);

                // copy triangle
                for (int j = 0; j < 3; j++)
                    triangle[j] = restMesh.face[findex].cP(j);

                // centrate
                origin = triangle[0];
                for (int j = 0; j < 3; j++)
                    triangle[j] -= origin;

                // new reference base for t1
                den = (triangle[1] - triangle[0]).Norm();
                base[0] = (triangle[1] - triangle[0]) / den;
                base[2] = restMesh.face[findex].cN();
                base[1] = base[2] ^ base[0];

                // remap average field to global reference
                fields[restIndex][findex].first = base[0] * R(0) + base[1] * R(1);

                // average/interpolate first field
//                fields[restIndex][findex].first = vcg::tri::CrossField<TriMeshType>::InterpolateCrossField(tangentsVector,
//                                                                                                           weigthsVector,
//                                                                                                           normalsVector,
//                                                                                                           baseNorm,
//                                                                                                           baseDir);

#ifdef NOT_NORMALIZE
                i = 0;
                for (typename list< Eigen::Matrix<ValueType,2,2> >::iterator it = VListmap.begin(); it != VListmap.end(); it++) {
                    V(i,0) = (*it)(1,0);
                    V(i,1) = (*it)(1,1);
                    i++;
                }
                // average vectors (in 2D)
                R.setZero();
                for (int j = 0; j < V.rows(); j++)
                  R += w(j) * V.row(j).transpose();
                R /= w.sum();
                // remap average field to global reference
                fields[restIndex][findex].second = base[0] * R(0) + base[1] * R(1);
#else
                // compute second field as perpendicolar to the first one and the face normal
                fields[restIndex][findex].second = restMesh.face[findex].cN() ^ fields[restIndex][findex].first;
#endif

                // average singular values
                singularValues[restIndex][findex].first /= VList.size();
                singularValues[restIndex][findex].second /= VList.size();

                // (weighting) average stretch values
                stretch[restIndex][findex].first /= den1;
                stretch[restIndex][findex].second /= den2;

            } else {
                // otherwise, search for maximum

                // reset variables
                singularValues[restIndex][findex] = pair<ValueType,ValueType>(ValueType(0), ValueType(0));
                stretch[restIndex][findex] = pair<ValueType,ValueType>(ValueType(0), ValueType(0));

                // temporary variables
                unsigned long indexOfMaxStretch = 0;
                ValueType tmpVal = ValueType(0);

                // for all other meshes
                for (unsigned long mindex = 0; mindex < numOfMeshes; mindex++) {
                    // of course, skip current mesh index (rest pose)
                    if (mindex == restIndex)
                        continue;

                    // check same number of fields
                    assert(findex < fields[mindex].size());

                    // search for maximal stretch value
                    if (getClampedAbsStretchMeasureAt(mindex, findex, stretchType) > tmpVal) {
                        tmpVal = getClampedAbsStretchMeasureAt(mindex, findex, stretchType);
                        indexOfMaxStretch = mindex;
                    }
                }

                // copy fields
                fields[restIndex][findex] = fields[indexOfMaxStretch][findex];
                fields2D[restIndex][findex] = fields2D[indexOfMaxStretch][findex];

                // copy singular values
                singularValues[restIndex][findex].first = getClampedAbsStretchMeasureAt(indexOfMaxStretch, findex, SINGVAL1);
                singularValues[restIndex][findex].second = getClampedAbsStretchMeasureAt(indexOfMaxStretch, findex, SINGVAL2);

                // copy stretches
                stretch[restIndex][findex].first = getClampedAbsStretchMeasureAt(indexOfMaxStretch, findex, MAXWARP);
                stretch[restIndex][findex].second = getClampedAbsStretchMeasureAt(indexOfMaxStretch, findex, RATIO);
            }

            // update min of first singular value
            if (std::fabs(singularValues[restIndex][findex].first) < minAbsSingularValues[restIndex].first)
                minAbsSingularValues[restIndex].first = std::fabs(singularValues[restIndex][findex].first);
            if (std::fabs(singularValues[restIndex][findex].first) < minAbsSequenceSingularValue1)
                minAbsSequenceSingularValue1 = std::fabs(singularValues[restIndex][findex].first);

            // update max of first singular value
            if (std::fabs(singularValues[restIndex][findex].first) > maxAbsSingularValues[restIndex].first)
                maxAbsSingularValues[restIndex].first = std::fabs(singularValues[restIndex][findex].first);
            if (std::fabs(singularValues[restIndex][findex].first) > maxAbsSequenceSingularValue1)
                maxAbsSequenceSingularValue1 = std::fabs(singularValues[restIndex][findex].first);

            // update min of second singular value
            if (std::fabs(singularValues[restIndex][findex].second) < minAbsSingularValues[restIndex].second)
                minAbsSingularValues[restIndex].second = std::fabs(singularValues[restIndex][findex].second);
            if (std::fabs(singularValues[restIndex][findex].second) < minAbsSequenceSingularValue2)
                minAbsSequenceSingularValue2 = std::fabs(singularValues[restIndex][findex].second);

            // update max of second singular value
            if (std::fabs(singularValues[restIndex][findex].second) > maxAbsSingularValues[restIndex].second)
                maxAbsSingularValues[restIndex].second = std::fabs(singularValues[restIndex][findex].second);
            if (std::fabs(singularValues[restIndex][findex].second) > maxAbsSequenceSingularValue2)
                maxAbsSequenceSingularValue2 = std::fabs(singularValues[restIndex][findex].second);

            // update min first type of stretch
            if (stretch[restIndex][findex].first < minStretch[restIndex].first)
                minStretch[restIndex].first = stretch[restIndex][findex].first;
            if (stretch[restIndex][findex].first < minAbsSequenceStretch1)
                minAbsSequenceStretch1 = stretch[restIndex][findex].first;

            // update min second type of stretch
            if (stretch[restIndex][findex].second < minStretch[restIndex].second)
                minStretch[restIndex].second = stretch[restIndex][findex].second;
            if (stretch[restIndex][findex].second < minAbsSequenceStretch2)
                minAbsSequenceStretch2 = stretch[restIndex][findex].second;

            // update max first type of stretch
            if (stretch[restIndex][findex].first > maxStretch[restIndex].first)
                maxStretch[restIndex].first = stretch[restIndex][findex].first;
            if (stretch[restIndex][findex].first > maxAbsSequenceStretch1)
                maxAbsSequenceStretch1 = stretch[restIndex][findex].first;

            // update min second type of stretch
            if (stretch[restIndex][findex].second > maxStretch[restIndex].second)
                maxStretch[restIndex].second = stretch[restIndex][findex].second;
            if (stretch[restIndex][findex].second > maxAbsSequenceStretch2)
                maxAbsSequenceStretch2 = stretch[restIndex][findex].second;
        }

        // log
        if (toLog)
            toLog(string(" done.\n"));
    }

    /**
     * @brief computeMaximalFieldsAt
     * @param restIndex
     * @param restMesh
     * @param stretchType
     */
    template < typename TriMeshType >
    void computeMaximalFieldsAt(unsigned long restIndex, const TriMeshType &restMesh, const StretchType stretchType) {
        // check parameters
        assert(restIndex < numOfMeshes && restMesh.FN() > 0);

        // clear rest mesh fields
        clearAt(restIndex);

        // resize rest mesh fields vector
        resizeFieldsAt(restIndex, restMesh.FN());

        // log
        if (toLog) {
            ostringstream numtostr;
            numtostr << restIndex;
            toLog(string("Search for maximal fields for rest mesh ") + numtostr.str() + string(" over all other meshes ..."));
        }

        // for each face
        for (unsigned long findex = 0; findex < (unsigned long)restMesh.FN(); findex++) {
            // reset variables
            singularValues[restIndex][findex] = pair<ValueType,ValueType>(ValueType(0), ValueType(0));
            stretch[restIndex][findex] = pair<ValueType,ValueType>(ValueType(0), ValueType(0));

            // temporary variables
            unsigned long indexOfMaxStretch = 0;
            ValueType tmpVal = ValueType(0);

            // for all other meshes
            for (unsigned long mindex = 0; mindex < numOfMeshes; mindex++) {
                // of course, skip current mesh index (rest pose)
                if (mindex == restIndex)
                    continue;

                // check same number of fields
                assert(findex < fields[mindex].size());

                // search for maximal stretch value
                if (getClampedAbsStretchMeasureAt(mindex, findex, stretchType) > tmpVal) {
                    tmpVal = getClampedAbsStretchMeasureAt(mindex, findex, stretchType);
                    indexOfMaxStretch = mindex;
                }
            }

            // copy fields
            fields[restIndex][findex] = fields[indexOfMaxStretch][findex];
            fields2D[restIndex][findex] = fields2D[indexOfMaxStretch][findex];

            // copy singular values
            singularValues[restIndex][findex].first = getClampedAbsStretchMeasureAt(indexOfMaxStretch, findex, SINGVAL1);
            singularValues[restIndex][findex].second = getClampedAbsStretchMeasureAt(indexOfMaxStretch, findex, SINGVAL2);

            // copy stretches
            stretch[restIndex][findex].first = getClampedAbsStretchMeasureAt(indexOfMaxStretch, findex, MAXWARP);
            stretch[restIndex][findex].second = getClampedAbsStretchMeasureAt(indexOfMaxStretch, findex, RATIO);

            // update min of first singular value
            if (std::fabs(singularValues[restIndex][findex].first) < minAbsSingularValues[restIndex].first)
                minAbsSingularValues[restIndex].first = std::fabs(singularValues[restIndex][findex].first);
            if (std::fabs(singularValues[restIndex][findex].first) < minAbsSequenceSingularValue1)
                minAbsSequenceSingularValue1 = std::fabs(singularValues[restIndex][findex].first);

            // update max of first singular value
            if (std::fabs(singularValues[restIndex][findex].first) > maxAbsSingularValues[restIndex].first)
                maxAbsSingularValues[restIndex].first = std::fabs(singularValues[restIndex][findex].first);
            if (std::fabs(singularValues[restIndex][findex].first) > maxAbsSequenceSingularValue1)
                maxAbsSequenceSingularValue1 = std::fabs(singularValues[restIndex][findex].first);

            // update min of second singular value
            if (std::fabs(singularValues[restIndex][findex].second) < minAbsSingularValues[restIndex].second)
                minAbsSingularValues[restIndex].second = std::fabs(singularValues[restIndex][findex].second);
            if (std::fabs(singularValues[restIndex][findex].second) < minAbsSequenceSingularValue2)
                minAbsSequenceSingularValue2 = std::fabs(singularValues[restIndex][findex].second);

            // update max of second singular value
            if (std::fabs(singularValues[restIndex][findex].second) > maxAbsSingularValues[restIndex].second)
                maxAbsSingularValues[restIndex].second = std::fabs(singularValues[restIndex][findex].second);
            if (std::fabs(singularValues[restIndex][findex].second) > maxAbsSequenceSingularValue2)
                maxAbsSequenceSingularValue2 = std::fabs(singularValues[restIndex][findex].second);

            // update min first type of stretch
            if (stretch[restIndex][findex].first < minStretch[restIndex].first)
                minStretch[restIndex].first = stretch[restIndex][findex].first;
            if (stretch[restIndex][findex].first < minAbsSequenceStretch1)
                minAbsSequenceStretch1 = stretch[restIndex][findex].first;

            // update min second type of stretch
            if (stretch[restIndex][findex].second < minStretch[restIndex].second)
                minStretch[restIndex].second = stretch[restIndex][findex].second;
            if (stretch[restIndex][findex].second < minAbsSequenceStretch2)
                minAbsSequenceStretch2 = stretch[restIndex][findex].second;

            // update max first type of stretch
            if (stretch[restIndex][findex].first > maxStretch[restIndex].first)
                maxStretch[restIndex].first = stretch[restIndex][findex].first;
            if (stretch[restIndex][findex].first > maxAbsSequenceStretch1)
                maxAbsSequenceStretch1 = stretch[restIndex][findex].first;

            // update min second type of stretch
            if (stretch[restIndex][findex].second > maxStretch[restIndex].second)
                maxStretch[restIndex].second = stretch[restIndex][findex].second;
            if (stretch[restIndex][findex].second > maxAbsSequenceStretch2)
                maxAbsSequenceStretch2 = stretch[restIndex][findex].second;
        }

        // log
        if (toLog)
            toLog(string(" computed.\n"));
    }

    /**
     * @brief storeFields
     * @param file
     */
    void storeFields(string file) {
        // check state
        if (numOfMeshes == 0) {
            if (toLog)
                toLog("Error. Nothing to store.\n");
            return;
        }

        // use a file stream to write on files
        ofstream output(file.c_str(), ios_base::out | ios_base::app | ios_base::ate);

        // check stream state
        if (!output.is_open()) {
            if (toLog)
                toLog("Error. Not able to open file '" + file + "' for writing.\n");
            return;
        }

        // store bounds
        output << "Singular values lower bound: " << singularValue2LowerBound << endl;
        output << "Singular values upper bound: " << singularValue1UpperBound << endl << endl << endl;

        // store the number of meshes
        output << "Number of meshes: " << numOfMeshes << endl << endl;

        // for each mesh
        for (unsigned long i = 0; i < numOfMeshes; i++) {
            // store the number of fields
            output << "Number of faces: " << fields[i].size() << endl;

            // store header
            output << "1st field | 1st singular value | 2nd field | 2nd singular value | 1st type of stretch | 2nd type of stretch" << endl;

            // for each face
            for (unsigned long j = 0; j < fields[i].size(); j++) {
                // store 1st field and singular value
                output << "( " << fields[i][j].first[0] << " , "
                       << fields[i][j].first[1] << " , "
                       << fields[i][j].first[2] << " ) | "
                       << singularValues[i][j].first << " | ";
                // store 2nd field and singular value
                output << "( " << fields[i][j].second[0] << " , "
                       << fields[i][j].second[1] << " , "
                       << fields[i][j].second[2] << " ) | "
                       << singularValues[i][j].second << " | ";
                // store 1st and 2nd type of stretch
                output << stretch[i][j].first << " | "
                       << stretch[i][j].second << endl;
            }

            // leave one raw
            output << endl;
        }

        // close the file
        output.close();
    }

    /**
     * @brief loadFields
     * @param file
     * @param pos
     */
    void loadFields(string file, fstream::pos_type pos) {
        // first clear
        clear();

        // use a file stream to read from files
        ifstream input(file.c_str());

        // check stream state
        if (!input.is_open()) {
            if (toLog)
                toLog("Error. Not able to open file '" + file + "' for reading.\n");
            return;
        }

        // set position on file
        input.seekg(pos);

        // temporary variables
        unsigned long num = 0;
        ValueType val = ValueType(0);
        string tmp;

        // load bounds
        tmp = "Singular values lower bound: ";
        input.ignore(tmp.length());
        input >> val;
        assert(input.good());
        setAbsSingularValue2LowerBound(val);
        input.ignore();
        tmp = "Singular values upper bound: ";
        input.ignore(tmp.length());
        input >> val;
        assert(input.good());
        setAbsSingularValue1UpperBound(val);
        input.ignore();
        input.ignore();
        input.ignore();

        // load the number of meshes
        tmp = "Number of meshes: ";
        input.ignore(tmp.length());
        input >> num;
        assert(input.good());
        input.ignore();

        // resize vectors
        resizeNumOfMeshes(num);

        // log
        ostringstream text;
        text << "Loading fields for " << num << " meshes." << endl;
        toLog(text.str());

        // for each mesh
        for (unsigned long i = 0; i < numOfMeshes; i++) {
            // load the number of fields
            tmp = "Number of faces: ";
            input.ignore(tmp.length());
            input >> num;
            assert(input.good());
            input.ignore();

            // resize vector
            resizeFieldsAt(i, num);

            // ignore header
            tmp = "1st field | 1st singular value | 2nd field | 2nd singular value | 1st type of stretch | 2nd type of stretch";
            input.ignore(tmp.length());
            input.ignore();

            // log
            ostringstream text;
            text << "Loading fields for mesh " << i << " ..." << endl;
            toLog(text.str());

            // for each face
            for (unsigned long j = 0; j < num; j++) {
                // load 1st field and singular value
                tmp = "( ";
                input.ignore(tmp.length());
                input >> fields[i][j].first[0];
                assert(input.good());
                tmp = " , ";
                input.ignore(tmp.length());
                input >> fields[i][j].first[1];
                assert(input.good());
                tmp = " , ";
                input.ignore(tmp.length());
                input >> fields[i][j].first[2];
                assert(input.good());
                tmp = " ) | ";
                input.ignore(tmp.length());
                input >> singularValues[i][j].first;
                assert(input.good());

                // load 2nd field and singular value
                tmp = " | ( ";
                input.ignore(tmp.length());
                input >> fields[i][j].second[0];
                assert(input.good());
                tmp = " , ";
                input.ignore(tmp.length());
                input >> fields[i][j].second[1];
                assert(input.good());
                tmp = " , ";
                input.ignore(tmp.length());
                input >> fields[i][j].second[2];
                assert(input.good());
                tmp = " ) | ";
                input.ignore(tmp.length());
                input >> singularValues[i][j].second;
                assert(input.good());

                // load 1st and 2nd type of stretch
                tmp = " | ";
                input.ignore(tmp.length());
                input >> stretch[i][j].first;
                assert(input.good());
                tmp = " | ";
                input.ignore(tmp.length());
                input >> stretch[i][j].second;
                assert(input.good());

                // ignore new line char
                input.ignore();

                // update min of first singular value
                if (std::fabs(singularValues[i][j].first) < minAbsSingularValues[i].first)
                    minAbsSingularValues[i].first = std::fabs(singularValues[i][j].first);
                if (std::fabs(singularValues[i][j].first) < minAbsSequenceSingularValue1)
                    minAbsSequenceSingularValue1 = std::fabs(singularValues[i][j].first);

                // update max of first singular value
                if (std::fabs(singularValues[i][j].first) > maxAbsSingularValues[i].first)
                    maxAbsSingularValues[i].first = std::fabs(singularValues[i][j].first);
                if (std::fabs(singularValues[i][j].first) > maxAbsSequenceSingularValue1)
                    maxAbsSequenceSingularValue1 = std::fabs(singularValues[i][j].first);

                // update min of second singular value
                if (std::fabs(singularValues[i][j].second) < minAbsSingularValues[i].second)
                    minAbsSingularValues[i].second = std::fabs(singularValues[i][j].second);
                if (std::fabs(singularValues[i][j].second) < minAbsSequenceSingularValue2)
                    minAbsSequenceSingularValue2 = std::fabs(singularValues[i][j].second);

                // update max of second singular value
                if (std::fabs(singularValues[i][j].second) > maxAbsSingularValues[i].second)
                    maxAbsSingularValues[i].second = std::fabs(singularValues[i][j].second);
                if (std::fabs(singularValues[i][j].second) > maxAbsSequenceSingularValue2)
                    maxAbsSequenceSingularValue2 = std::fabs(singularValues[i][j].second);

                // update min first type of stretch
                if (stretch[i][j].first < minStretch[i].first)
                    minStretch[i].first = stretch[i][j].first;
                if (stretch[i][j].first < minAbsSequenceStretch1)
                    minAbsSequenceStretch1 = stretch[i][j].first;

                // update min second type of stretch
                if (stretch[i][j].second < minStretch[i].second)
                    minStretch[i].second = stretch[i][j].second;
                if (stretch[i][j].second < minAbsSequenceStretch2)
                    minAbsSequenceStretch2 = stretch[i][j].second;

                // update max first type of stretch
                if (stretch[i][j].first > maxStretch[i].first)
                    maxStretch[i].first = stretch[i][j].first;
                if (stretch[i][j].first > maxAbsSequenceStretch1)
                    maxAbsSequenceStretch1 = stretch[i][j].first;

                // update min second type of stretch
                if (stretch[i][j].second > maxStretch[i].second)
                    maxStretch[i].second = stretch[i][j].second;
                if (stretch[i][j].second > maxAbsSequenceStretch2)
                    maxAbsSequenceStretch2 = stretch[i][j].second;
            }

            // leave one raw
            input.ignore();
        }

        // log
        toLog("Fields loaded.\n");

        // close the file
        input.close();
    }

    /**
     * @brief setLogFunction
     */
    inline void setLogFunction(void (*logFunction)(string)) {
        toLog = logFunction;
    }

private:
    // fields
    vector< vector< pair<FieldType,FieldType> > > fields;
    vector< vector< Eigen::Matrix<ValueType,2,1> > > fields2D;
#ifdef NOT_NORMALIZE
    vector< vector< Eigen::Matrix<ValueType,2,2> > > fields2Dmap;
#endif
    vector< vector< pair<ValueType,ValueType> > > singularValues;
    vector< pair<ValueType,ValueType> > minAbsSingularValues;
    vector< pair<ValueType,ValueType> > maxAbsSingularValues;
    vector< vector< pair<ValueType,ValueType> > > stretch;
    vector< pair<ValueType,ValueType> > minStretch;
    vector< pair<ValueType,ValueType> > maxStretch;
    ValueType minAbsSequenceSingularValue1;
    ValueType minAbsSequenceSingularValue2;
    ValueType minAbsSequenceStretch1;
    ValueType minAbsSequenceStretch2;
    ValueType maxAbsSequenceSingularValue1;
    ValueType maxAbsSequenceSingularValue2;
    ValueType maxAbsSequenceStretch1;
    ValueType maxAbsSequenceStretch2;
    ValueType singularValue1UpperBound;
    ValueType singularValue2LowerBound;

    // state variables
    unsigned long numOfMeshes;

    // pointer to log function
    void (*toLog)(string text);
};

#endif // MESH_SEQUENCEFIELDS_H
