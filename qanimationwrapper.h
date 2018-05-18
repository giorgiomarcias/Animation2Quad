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

#ifndef QANIMATIONWRAPPER_H
#define QANIMATIONWRAPPER_H

#include <QObject>
#include <QVector>
#include <QReadWriteLock>
#include "animation.h"

/**
 * @brief The QAnimationWrapper class Gives thread-safeness to Animation class.
 */
class QAnimationWrapper : public QObject, private Animation
{
    Q_OBJECT
public:
    /**
     * @brief QAnimationWrapper
     * @param parent
     */
    explicit QAnimationWrapper(QObject *parent = 0) : QObject(parent) {
        numOfMeshes = 0;
    }

    /**
     * @brief setLogFunction
     */
    void setLogFunction(void (*logFunction)(string)) {
        // lock mutex
        resetMutex.lockForWrite();

        // do work
        Animation::setLogFunction(logFunction);

        // unlock mutex
        resetMutex.unlock();
    }
    
    /**
     * @brief isMeshLoadedAt
     * @param i
     * @return
     */
    bool isMeshLoadedAt(unsigned long i) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(i < numOfMeshes);

        // lock single mutex
        mutexes[i]->lockForRead();

        // read
        bool ret = Animation::isMeshLoadedAt(i);

        // unlock single mutex
        mutexes[i]->unlock();

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief loadMeshAt
     * @param i
     */
    void loadMeshAt(const unsigned long i) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(i < numOfMeshes);

        // lock single mutex
        mutexes[i]->lockForWrite();

        // do work
        Animation::loadMeshAt(i);

        // unlock single mutex
        mutexes[i]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief unloadMeshAt
     * @param i
     */
    void unloadMeshAt(const unsigned long i) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(i < numOfMeshes);

        // lock single mutex
        mutexes[i]->lockForWrite();

        // do work
        Animation::unloadMeshAt(i);

        // unlock single mutex
        mutexes[i]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief setNames
     * @param fnames
     */
    void setNames(vector<string> &fnames) {
        // lock mutex
        resetMutex.lockForWrite();

        // do work
        Animation::setNames(fnames);

        // get number of meshes
        numOfMeshes = fnames.size();

        // resize mutexes vector
        mutexes.resize(numOfMeshes);

        // reset mutexes
        for (unsigned long k = 0; k < numOfMeshes; k++)
            mutexes[k] = new QReadWriteLock;

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief setNames
     * @param fnames
     */
    void setNames(list<string> &fnames) {
        // construct vector of names
        vector<string> vnames(fnames.begin(), fnames.end());

        // initialize
        setNames(vnames);
    }

    /**
     * @brief load
     * @param fnames
     */
    void load(vector<string> &fnames) {
        // first initialize
        setNames(fnames);

        // lock mutex
        resetMutex.lockForWrite();

        // do work
        Animation::load(fnames);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief load
     * @param fnames
     */
    void load(list<string> &fnames) {
        // first initialize
        setNames(fnames);

        // lock mutex
        resetMutex.lockForWrite();

        // do work
        Animation::load(fnames);

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief unload
     */
    void unload() {
        // lock mutex
        resetMutex.lockForWrite();

        // reset mutexes
        for (unsigned long k = 0; k < numOfMeshes; k++)
            delete mutexes[k];

        // resize mutexes vector
        mutexes.clear();

        // get number of meshes
        numOfMeshes = 0;

        // do work
        Animation::unload();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief getMeshNameAt
     * @param i
     * @param name
     */
    void getMeshNameAt(const unsigned long i, string & name) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(i < numOfMeshes);

        // lock single mutex
        mutexes[i]->lockForRead();

        // do work
        Animation::getMeshNameAt(i, name);

        // unlock single mutex
        mutexes[i]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief getTriMeshAt
     * @param i
     * @param mesh
     */
    void getTriMeshAt(const unsigned long i, TMesh & mesh) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(i < numOfMeshes);

        // lock single mutex
        mutexes[i]->lockForRead();

        // do work
        Animation::getTriMeshAt(i, mesh);

        // unlock single mutex
        mutexes[i]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    /**
     * @brief getTriMeshAt
     * @param i
     * @return
     * @note Important: call .unlockAt(i) when finished of using both TMesh and PMesh.
     */
    TMesh * getTriMeshAt(const unsigned long i) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(i < numOfMeshes);

        // lock single mutex
        mutexes[i]->lockForRead();

        // do work
        TMesh * ret = Animation::getTriMeshAt(i);

        // unlock mutex
        resetMutex.unlock();

        // return mesh pointer
        return ret;
    }

//    /**
//     * @brief getPolyMeshAt
//     * @param i
//     * @param mesh
//     */
//    void getPolyMeshAt(const unsigned long i, PMesh & mesh) {
//        // lock mutex
//        resetMutex.lockForRead();

//        // check parameter
//        assert(i < numOfMeshes);

//        // lock single mutex
//        mutexes[i]->lockForRead();

//        // do work
//        Animation::getPolyMeshAt(i, mesh);   // ****************** doesn't work: only for tri meshes!!!!!!!

//        // unlock single mutex
//        mutexes[i]->unlock();

//        // unlock mutex
//        resetMutex.unlock();
//    }

//    /**
//     * @brief getPolyMeshAt
//     * @param i
//     * @return
//     * @note Important: call .unlockAt(i) when finished of using both TMesh and PMesh.
//     */
//    PMesh * getPolyMeshAt(const unsigned long i) {
//        // lock mutex
//        resetMutex.lockForRead();

//        // check parameter
//        assert(i < numOfMeshes);

//        // lock single mutex
//        mutexes[i]->lockForRead();

//        // do work
//        PMesh * ret = Animation::getPolyMeshAt(i);

//        // unlock mutex
//        resetMutex.unlock();

//        // return mesh pointer
//        return ret;
//    }

    /**
     * @brief Animation::getNumberOfMeshes
     * @return
     */
    unsigned long getNumberOfMeshes() {
        // lock mutex
        resetMutex.lockForRead();

        // do work
        unsigned long ret = numOfMeshes;

        // unlock mutex
        resetMutex.unlock();

        return ret;
    }

    /**
     * @brief unlockAt
     * @param i
     */
    void unlockAt(unsigned long i) {
        // lock mutex
        resetMutex.lockForRead();

        // check parameter
        assert(i < numOfMeshes);

        // unlock single mutex
        mutexes[i]->unlock();

        // unlock mutex
        resetMutex.unlock();
    }

    // set public
    using Animation::computeAnisotropyTau;
    
private:
    QVector<QReadWriteLock *> mutexes;
    QReadWriteLock resetMutex;
    unsigned long numOfMeshes;
};

#endif // QANIMATIONWRAPPER_H
