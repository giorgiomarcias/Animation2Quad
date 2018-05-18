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

#ifndef QMESHSEQUENCEFIELDSWRAPPER_H
#define QMESHSEQUENCEFIELDSWRAPPER_H

#include <QObject>
#include <QVector>
#include <QReadWriteLock>
#include <vcg/complex/complex.h>
#include "meshsequencefields.h"

typedef float ScalarType;

class QMeshSequenceFieldsWrapper : public QObject, private MeshSequenceFields<ScalarType>
{
    Q_OBJECT
public:
    // turn StretchType again into a public enum
    using MeshSequenceFields<ScalarType>::StretchType;
    using MeshSequenceFields<ScalarType>::SINGVAL1;
    using MeshSequenceFields<ScalarType>::SINGVAL2;
    using MeshSequenceFields<ScalarType>::MAXWARP;
    using MeshSequenceFields<ScalarType>::RATIO;

    /**
     * @brief QMeshSequenceFieldsWrapper
     * @param parent
     */
    explicit QMeshSequenceFieldsWrapper(QObject *parent = 0) : QObject(parent) {
        numOfMeshes = 0;
    }

    /**
     * @brief setLogFunction
     */
    inline void setLogFunction(void (*logFunction)(string)) {
        // lock mutex
        resetMutex.lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::setLogFunction(logFunction);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief clearAt
     * @param meshIndex
     */
    void clearAt(unsigned long meshIndex) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(meshIndex < numOfMeshes);

        // lock single mutex
        mutexes[meshIndex]->lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::clearAt(meshIndex);

        // unlock single mutex
        mutexes[meshIndex]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief clear
     */
    void clear() {
        // lock mutex
        resetMutex.lockForWrite();

        // reset mutexes
        for (unsigned long k = 0; k < numOfMeshes; k++)
            delete mutexes[k];

        // resize mutexes vector
        mutexes.clear();

        // reset number of meshes
        numOfMeshes = 0;

        // do work
        MeshSequenceFields<ScalarType>::clear();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief getNumberOfMeshes
     * @return
     */
    inline unsigned long getNumberOfMeshes() {
        // lock mutex
        resetMutex.lockForRead();

        unsigned long ret = numOfMeshes;

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief getNumberOfFieldsAt
     * @param meshIndex
     * @return
     */
    inline unsigned long getNumberOfFieldsAt(unsigned long meshIndex) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(meshIndex < numOfMeshes);

        // lock single mutex
        mutexes[meshIndex]->lockForRead();

        // do work
        unsigned long ret = MeshSequenceFields<ScalarType>::getNumberOfFieldsAt(meshIndex);

        // unlock single mutex
        mutexes[meshIndex]->unlock();

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief resizeFieldsAt
     * @param meshIndex
     * @param num
     */
    void resizeFieldsAt(unsigned long meshIndex, unsigned long num) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(meshIndex < numOfMeshes);

        // lock single mutex
        mutexes[meshIndex]->lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::resizeFieldsAt(meshIndex, num);

        // unlock single mutex
        mutexes[meshIndex]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief resizeNumOfMeshes
     * @param num
     */
    void resizeNumOfMeshes(unsigned long num) {
        // first clear
        clear();

        // lock mutex
        resetMutex.lockForWrite();

        // reset number of meshes
        numOfMeshes = num;

        // resize mutexes vector
        mutexes.resize(numOfMeshes);

        // reset mutexes
        for (unsigned long k = 0; k < numOfMeshes; k++)
            mutexes[k] = new QReadWriteLock;

        // do work
        MeshSequenceFields<ScalarType>::resizeNumOfMeshes(num);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief getFirstFieldAt
     * @param meshIndex
     * @param faceIndex
     * @param field
     */
    inline void getFirstFieldAt(unsigned long meshIndex, unsigned long faceIndex, FieldType & field) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(meshIndex < numOfMeshes);

        // lock single mutex
        mutexes[meshIndex]->lockForRead();

        // do work
        MeshSequenceFields<ScalarType>::getFirstFieldAt(meshIndex, faceIndex, field);

        // unlock single mutex
        mutexes[meshIndex]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief getSecondFieldAt
     * @param meshIndex
     * @param faceIndex
     * @param field
     */
    inline void getSecondFieldAt(unsigned long meshIndex, unsigned long faceIndex, FieldType & field) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(meshIndex < numOfMeshes);

        // lock single mutex
        mutexes[meshIndex]->lockForRead();

        // do work
        MeshSequenceFields<ScalarType>::getSecondFieldAt(meshIndex, faceIndex, field);

        // unlock single mutex
        mutexes[meshIndex]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief getStretchMeasureAt
     * @param meshIndex
     * @param faceIndex
     * @param stretchType
     * @return
     */
    inline ValueType getStretchMeasureAt(unsigned long meshIndex, unsigned long faceIndex, const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(meshIndex < numOfMeshes);

        // lock single mutex
        mutexes[meshIndex]->lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getStretchMeasureAt(meshIndex, faceIndex, stretchType);

        // unlock single mutex
        mutexes[meshIndex]->unlock();

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief getMinAbsStretchMeasureAt
     * @param meshIndex
     * @param stretchType
     * @return
     */
    inline ValueType getMinAbsStretchMeasureAt(unsigned long meshIndex, const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(meshIndex < numOfMeshes);

        // lock single mutex
        mutexes[meshIndex]->lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getMinAbsStretchMeasureAt(meshIndex, stretchType);

        // unlock single mutex
        mutexes[meshIndex]->unlock();

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief getMaxAbsStretchMeasureAt
     * @param meshIndex
     * @param stretchType
     * @return
     */
    inline ValueType getMaxAbsStretchMeasureAt(unsigned long meshIndex, const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(meshIndex < numOfMeshes);

        // lock single mutex
        mutexes[meshIndex]->lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getMaxAbsStretchMeasureAt(meshIndex, stretchType);

        // unlock single mutex
        mutexes[meshIndex]->unlock();

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief getMinAbsSequenceStretchMeasure
     * @param stretchType
     * @return
     */
    inline ValueType getMinAbsSequenceStretchMeasure(const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getMinAbsSequenceStretchMeasure(stretchType);

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief getMaxAbsSequenceStretchMeasure
     * @param stretchType
     * @return
     */
    inline ValueType getMaxAbsSequenceStretchMeasure(const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getMaxAbsSequenceStretchMeasure(stretchType);

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief setAbsSingularValue2LowerBound
     * @param lowerBound
     */
    inline void setAbsSingularValue2LowerBound(const ValueType lowerBound) {
        // lock mutex
        resetMutex.lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::setAbsSingularValue2LowerBound(lowerBound);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief setAbsSingularValue1UpperBound
     * @param upperBound
     */
    inline void setAbsSingularValue1UpperBound(const ValueType upperBound) {
        // lock mutex
        resetMutex.lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::setAbsSingularValue1UpperBound(upperBound);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief getAbsStretchMeasureLowerBound
     * @param stretchType
     * @return
     */
    inline ValueType getAbsStretchMeasureLowerBound(const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getAbsStretchMeasureLowerBound(stretchType);

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief getAbsStretchMeasureUpperBound
     * @param stretchType
     * @return
     */
    inline ValueType getAbsStretchMeasureUpperBound(const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getAbsStretchMeasureUpperBound(stretchType);

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief getClampedAbsStretchMeasure
     * @param value
     * @param stretchType
     * @return
     */
    inline ValueType getClampedAbsStretchMeasure(ValueType value, const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getClampedAbsStretchMeasure(value, stretchType);

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief getClampedAbsStretchMeasureAt
     * @param meshIndex
     * @param faceIndex
     * @param stretchType
     * @return
     */
    inline ValueType getClampedAbsStretchMeasureAt(unsigned long meshIndex, unsigned long faceIndex, const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForRead();

        // do work
        ScalarType ret = MeshSequenceFields<ScalarType>::getClampedAbsStretchMeasureAt(meshIndex, faceIndex, stretchType);

        // unlock mutex
        resetMutex.unlock();

        return ret;
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
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(index < numOfMeshes);

        // lock single mutex
        mutexes[index]->lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::computeMeshFieldsAt(index, restMesh, referenceMesh);

        // unlock single mutex
        mutexes[index]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief computeAverageFieldsAt
     * @param restIndex
     * @param restMesh
     */
    template < typename TriMeshType >
    void computeAverageFieldsAt(unsigned long restIndex, const TriMeshType &restMesh,
                                const StretchType stretchType, const ValueType threshold = ValueType(0)) {
        // lock mutex
        resetMutex.lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::computeAverageFieldsAt(restIndex, restMesh, stretchType, threshold);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief computeMaximalFieldsAt
     * @param restIndex
     * @param restMesh
     * @param stretchType
     */
    template < typename TriMeshType >
    void computeMaximalFieldsAt(unsigned long restIndex, const TriMeshType &restMesh, const StretchType stretchType) {
        // lock mutex
        resetMutex.lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::computeMaximalFieldsAt(restIndex, restMesh, stretchType);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief storeFields
     * @param file
     */
    void storeFields(string file) {
        // lock mutex
        resetMutex.lockForWrite();

        // do work
        MeshSequenceFields<ScalarType>::storeFields(file);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief loadFields
     * @param file
     * @param pos
     */
    void loadFields(string file, fstream::pos_type pos) {
        // first clear
        clear();

        // lock mutex
        resetMutex.lockForRead();

        // delete mutexes
        mutexes.clear();

        // do work
        MeshSequenceFields<ScalarType>::loadFields(file, pos);

        // reset number of meshes
        numOfMeshes = MeshSequenceFields<ScalarType>::getNumberOfMeshes();

        // resize mutexes vector
        mutexes.resize(numOfMeshes);

        // reset mutexes
        for (unsigned long k = 0; k < numOfMeshes; k++)
            mutexes[k] = new QReadWriteLock;

        // unlock mutex
        resetMutex.unlock();
    }
    
private:
    QVector<QReadWriteLock *> mutexes;
    QReadWriteLock resetMutex;
    unsigned long numOfMeshes;
};

#endif // QMESHSEQUENCEFIELDSWRAPPER_H
