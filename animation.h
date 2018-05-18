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

#ifndef ANIMATION_H
#define ANIMATION_H

#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
using namespace std;

#include "meshes.h"

#include <wrap/io_trimesh/import.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/polygon_support.h>
#include <vcg/complex/append.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/curvature_fitting.h>
#include <vcg/space/point2.h>
#include "compute.h"

/**
 * @brief The Animation class Represents a sequence of meshes.
 */
class Animation
{
public:
    /**
     * @brief Animation::Animation Construct an empty mesh sequence.
     */

    Animation() {
        // no meshes
        numOfMeshes = 0;
        toLog = 0;
    }

    /**
     * @brief isMeshLoadedAt
     * @param i
     * @return
     */
    bool isMeshLoadedAt(unsigned long i) {
        // check parameter
        assert(i < numOfMeshes);

        return meshesLoaded[i];
    }

    /**
     * @brief loadMeshAt
     * @param i
     */
    void loadMeshAt(const unsigned long i) {
        // check parameter
        assert(i < numOfMeshes);

        // check if already loaded
        if (meshesLoaded[i])
            return;

        // check filename
        if (filenames[i].empty())
            return;

        // check allocation of mesh
        if (!triMeshes[i])
            triMeshes[i] = new TMesh;

        if (toLog)
            toLog(filenames[i] + ": loading...");

        // import triangle mesh from file
        int err = vcg::tri::io::Importer<TMesh>::Open(*triMeshes[i], filenames[i].c_str());

        // check error(s)
        if (vcg::tri::io::Importer<TMesh>::ErrorCritical(err)) {
            if (toLog)
                toLog(filenames[i] + ": error while loading into triangle mesh. Not imported.\n");
            return;
        }

        // set features and properties
        vcg::tri::Clean<TMesh>::RemoveDegenerateFace(*triMeshes[i]);
        vcg::tri::UpdateTopology<TMesh>::VertexFace(*triMeshes[i]);
        vcg::tri::UpdateTopology<TMesh>::FaceFace(*triMeshes[i]);
        vcg::tri::UpdateBounding<TMesh>::Box(*triMeshes[i]);
        vcg::tri::UpdateNormal<TMesh>::PerVertexNormalizedPerFaceNormalized(*triMeshes[i]);
        vcg::tri::UpdateNormal<TMesh>::PerBitQuadFaceNormalized(*triMeshes[i]);
        vcg::tri::BitQuad<TMesh>::UpdateValencyInFlags(*triMeshes[i]);

//        // check allocation of mesh
//        if (!polyMeshes[i])
//            polyMeshes[i] = new PMesh;

//        // check if it is already a polygon mesh
//        int j = 0;
//        for (; j < triMeshes[i]->FN(); j++)
//            if (triMeshes[i]->face[j].VN() > 3)
//                break;

//        // import polygon mesh from triangle one
//        if (j < triMeshes[i]->FN())
//            vcg::tri::io::Importer<PMesh>::Open(*polyMeshes[i], filenames[i].c_str());
//        else
//            vcg::tri::PolygonSupport<TMesh,PMesh>::ImportFromTriMesh(*polyMeshes[i], *triMeshes[i]);

//        // update features
//        vcg::tri::UpdateTopology<PMesh>::VertexFace(*polyMeshes[i]);
//        vcg::tri::UpdateTopology<PMesh>::FaceFace(*polyMeshes[i]);
//        vcg::tri::UpdateBounding<PMesh>::Box(*polyMeshes[i]);
//        vcg::tri::UpdateNormal<PMesh>::PerVertexNormalizedPerFaceNormalized(*polyMeshes[i]);

        // log
        if (toLog)
            toLog("computing curvature...");

        // compute curvatures
        computeCurvaturePerVertexPerFace(i);

        // set loaded
        meshesLoaded[i] = true;
        if (toLog)
            toLog("loaded.\n");
    }

    /**
     * @brief unloadMeshAt
     * @param i
     */
    void unloadMeshAt(const unsigned long i) {
        // check parameter
        assert(i < numOfMeshes);

        // check if already unloaded
        if (!meshesLoaded[i])
            return;

        // unload by deleting it
        delete triMeshes[i];
        triMeshes[i] = 0;
//        delete polyMeshes[i];
//        polyMeshes[i] = 0;

        // set loaded
        meshesLoaded[i] = false;
        if (toLog)
            toLog(filenames[i] + ": unloaded.\n");
    }

    /**
     * @brief setNames
     * @param fnames
     */
    void setNames(vector<string> &fnames) {
        // first unload
        unload();

        // set num of meshes
        numOfMeshes = fnames.size();

        // resize containers
        filenames.resize(numOfMeshes);
        triMeshes.resize(numOfMeshes);
//        polyMeshes.resize(numOfMeshes);
        meshesLoaded.resize(numOfMeshes);

        // initialize containers
        for (unsigned long i = 0; i < numOfMeshes; i++)
            filenames[i] = fnames[i];
        std::fill(triMeshes.begin(), triMeshes.end(), (TMesh *)0);
//        std::fill(polyMeshes.begin(), polyMeshes.end(), (PMesh *)0);
        std::fill(meshesLoaded.begin(), meshesLoaded.end(), false);
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
     * @brief Animation::load
     * @param fnames
     * @return
     */
    void load(vector<string> &fnames) {
        // first unload all
        unload();

        // initialize
        setNames(fnames);

        // load all meshes
        for (unsigned long i = 0; i < numOfMeshes; i++)
            // load mesh i-th
            loadMeshAt(i);
    }

    /**
     * @brief Animation::load
     * @param fnames
     * @return
     */
    void load(list<string> &fnames) {
        // construct vector of names
        vector<string> vnames(fnames.begin(), fnames.end());

        // load meshes
        load(vnames);
    }

    /**
     * @brief unload
     * @return
     */
    void unload() {
        for (unsigned long i = 0; i < numOfMeshes; i++)
            // unload i-th mesh
            unloadMeshAt(i);

        // resize vectors
        filenames.clear();
        triMeshes.clear();
//        polyMeshes.clear();
        meshesLoaded.clear();

        // reset number of meshes
        numOfMeshes = 0;
    }

    /**
     * @brief getMeshNameAt
     * @param i
     * @param name
     */
    void getMeshNameAt(const unsigned long i, string & name) {
        // check parameter
        assert(i < numOfMeshes);

        name = filenames[i];
    }

    /**
     * @brief getTriMeshAt
     * @param i
     * @param mesh
     */
    void getTriMeshAt(const unsigned long i, TMesh & mesh) {
        // check parameter
        assert(i < numOfMeshes);

        // check if unloaded, then load it
        if (!meshesLoaded[i])
            loadMeshAt(i);

        // copy mesh
        vcg::tri::Append<TMesh,TMesh>::MeshCopy(mesh, *triMeshes[i]);
    }

    /**
     * @brief getTriMeshAt
     * @param i
     * @return
     */
    TMesh * getTriMeshAt(const unsigned long i) {
        // check parameter
        assert(i < numOfMeshes);

        // check if unloaded, then load it
        if (!meshesLoaded[i])
            loadMeshAt(i);

        // return mesh pointer
        return triMeshes[i];
    }

//    /**
//     * @brief getPolyMeshAt
//     * @param i
//     * @param mesh
//     */
//    void getPolyMeshAt(const unsigned long i, PMesh & mesh) {
//        // check parameter
//        assert(i < numOfMeshes);

//        // check if unloaded, then load it
//        if (!meshesLoaded[i])
//            loadMeshAt(i);

//        // copy mesh
//        vcg::tri::Append<PMesh,PMesh>::MeshCopy(mesh, *polyMeshes[i]); // dowsn't work: only for tri meshes
//    }

//    /**
//     * @brief getPolyMeshAt
//     * @param i
//     * @return
//     */
//    PMesh * getPolyMeshAt(const unsigned long i) {
//        // check parameter
//        assert(i < numOfMeshes);

//        // check if unloaded, then load it
//        if (!meshesLoaded[i])
//            loadMeshAt(i);

//        // return mesh pointer
//        return polyMeshes[i];
//    }

    /**
     * @brief Animation::getNumberOfMeshes
     * @return
     */
    unsigned long getNumberOfMeshes() {
        return numOfMeshes;
    }

    /**
     * @brief Animation::setLogFunction
     */
    void setLogFunction(void (*logFunction)(string)) {
        toLog = logFunction;
    }

    /**
     * @brief computeAnisotropyTau
     * @param k1
     * @param k2
     * @return
     */
    static TMesh::ScalarType computeAnisotropyTau(const float k1, float k2) {
        // get face curvatures to compute anisotropy (from MIQ paper)
//        return k1 == 0.0 ? 0.0 : std::fabs(std::fabs(k1) - std::fabs(k2)) / std::fabs(k1);

        // get face curvatures to compute anisotropy (by Enrico)
        return std::fabs(k1) <= 0.1 && std::fabs(k2) <= 0.1 ? 0.0 : std::fabs(std::fabs(k1) - std::fabs(k2))
                                                / std::max(std::fabs(k1),std::fabs(k2));
    }

private:
    // containers
    unsigned long numOfMeshes;
    vector<string> filenames;
    vector<TMesh *> triMeshes;
//    vector<PMesh *> polyMeshes;

    // state variables
    vector<bool> meshesLoaded;

    // pointer to log function
    void (*toLog)(string text);

    /**
     * @brief computeCurvaturePerVertexPerFace
     * @param mesh
     */
    void computeCurvaturePerVertexPerFace(const unsigned long i) {
        assert(i < numOfMeshes);
        assert(triMeshes[i]);
        TMesh &mesh = *triMeshes[i];

        // try to read curvature from file
        bool readCurv = readVertexCurvatureFromFile(i);

        if (!readCurv) {
            // compute curvature per vertex
//            vcg::tri::UpdateCurvatureFitting<TMesh>::computeCurvature(mesh);
            vcg::tri::UpdateCurvature<TMesh>::PrincipalDirectionsPCA(mesh, mesh.bbox.Diag()/20);
//            vcg::tri::UpdateCurvature<TMesh>::PrincipalDirections(mesh);
//            vcg::tri::UpdateCurvature<TMesh>::PrincipalDirectionsNormalCycle(mesh);
        }

        // compute per vertex anisotropy
        for (int i = 0; i < mesh.VN(); i++)
            mesh.vert[i].Q() = computeAnisotropyTau(mesh.vert[i].K1(), mesh.vert[i].K2());

        // temporary variables
        TMesh::CoordType nf, nv, cv3D;
        vcg::Matrix33<TMesh::ScalarType> R;
        vcg::Point2<TMesh::ScalarType> cv2D;
        Eigen::Matrix<TMesh::ScalarType,Eigen::Dynamic,2> V(3, 2);
        Eigen::Matrix<TMesh::ScalarType,Eigen::Dynamic,1> w(3, 1);
        Eigen::Matrix<TMesh::ScalarType,2,1> curv2D;
        TMesh::CoordType triangle[3], base[3];
        TMesh::CoordType origin;
        TMesh::ScalarType den, /*den1,*/ den2, tau;

        // for each face
        for (int i = 0; i < mesh.FN(); i++) {
            // init curvature value
            mesh.face[i].K1() = TMesh::ScalarType(0);
            mesh.face[i].K2() = TMesh::ScalarType(0);
            mesh.face[i].Q() = TMesh::ScalarType(0);

            // get normal at face i
            nf = mesh.face[i].N().Normalize();

            // copy triangle
            for (int j = 0; j < 3; j++)
                triangle[j] = mesh.face[i].P(j);

            // centrate
            origin = triangle[0];
            for (int j = 0; j < 3; j++)
                triangle[j] -= origin;

            // new reference base for triangle i
            den = (triangle[1] - triangle[0]).Norm();
            base[0] = (triangle[1] - triangle[0]) / den;
            base[2] = nf;
            base[1] = base[2] ^ base[0];

//            den1 = TMesh::ScalarType(0);
//            for (int j = 0; j < 3; j++)
//                if (std::max(std::fabs(mesh.face[i].V(j)->K1()), std::fabs(mesh.face[i].V(j)->K2())) > den1)
//                    den1 = std::max(std::fabs(mesh.face[i].V(j)->K1()), std::fabs(mesh.face[i].V(j)->K2()));

            den2 = TMesh::ScalarType(0);
            tau = TMesh::ScalarType(0);
            // rotate curvature dirs at vertex to lie on triangle plane
            for (int j = 0; j < 3; j++) {
                // get normals at vertices
                nv = mesh.face[i].V(j)->N().Normalize();

                // get rotation matrixes for alineating vertex normals with face normal
                R = vcg::RotationMatrix(nv, nf);

                // rotate curvature 1
                cv3D = R * mesh.face[i].V(j)->cPD1();
                cv3D.Normalize();

                // get curvature dir in 2D
                cv2D = vcg::Point2<TMesh::ScalarType>(cv3D * base[0], cv3D * base[1]);

                // copy into eigen container
                V(j,0) = cv2D[0];
                V(j,1) = cv2D[1];

                // set weight
//                w(j) = den1 <= TMesh::ScalarType(0.01) ? TMesh::ScalarType(0) :
//                                std::max(std::fabs(mesh.face[i].V(j)->K1()), std::fabs((mesh.face[i].V(j)->K2())) / den1;
                w(j) = mesh.face[i].V(j)->Q();

                // sum curvature values
                mesh.face[i].K1() += mesh.face[i].V(j)->Q() * mesh.face[i].V(j)->K1();
                mesh.face[i].K2() += mesh.face[i].V(j)->Q() * mesh.face[i].V(j)->K2();
                den2 += mesh.face[i].V(j)->Q();

//                // search for max vertex anisotropy
//                if (mesh.face[i].V(j)->Q() > tau)
//                    tau = mesh.face[i].V(j)->Q();
                // sum anisotropy
                mesh.face[i].Q() += mesh.face[i].V(j)->Q() * mesh.face[i].V(j)->Q();
            }

            // average curvature values
            mesh.face[i].K1() /= den2;
            mesh.face[i].K2() /= den2;

//            // set anisotropy
//            mesh.face[i].Q() = tau;
            // average anisotropy
            mesh.face[i].Q() /= den2;

            // compute average in 2D
            curv2D = averageVectors(V, w, 4);

            // get curvature dir in 3D
            mesh.face[i].PD1() = base[0] * curv2D[0] + base[1] * curv2D[1];

            // compute min curvature dir as cross product of the max one with the normal
            mesh.face[i].PD2() = nf ^ mesh.face[i].PD1();
        }
    }

    /**
     * @brief readVertexCurvatureFromFile
     * @param i
     * @return
     */
    bool readVertexCurvatureFromFile(const unsigned long i) {
        // construct curvature filename as mesh filename without
        string file = filenames[i] + string(".txt");

        // temporary variables
        char buffer[265];
        TMesh::ScalarType num;

        // use a file stream to read from files
        ifstream input(file.c_str());

        // check stream state
        if (!input.is_open()) {
            if (toLog)
                toLog("Warning. Not able to open file '" + file + "' of vertices curvature. ");
            return false;
        }

        // ignore first two lines (header)
        input.getline(buffer, 256);
        input.getline(buffer, 256);

        // for each vertex
        for (int j = 0; j < triMeshes[i]->VN(); j++) {
            // get k1
            input >> num;
            assert(input.good());
            // sat k1
            triMeshes[i]->vert[j].K1() = num;

            // ignore space
//            input.ignore();

            // get k2
            input >> num;
            assert(input.good());
            // sat k1
            triMeshes[i]->vert[j].K2() = num;

            // ignore space
//            input.ignore();

            // get PD1
            for (int k = 0; k < 3; k++) {
                // get coord
                input >> num;
                assert(input.good());

                // set coord
                triMeshes[i]->vert[j].PD1()[k] = num;

                // ignore space
//                input.ignore();
            }

            // get PD2
            for (int k = 0; k < 3; k++) {
                // get coord
                input >> num;
                assert(input.good());

                // set coord
                triMeshes[i]->vert[j].PD2()[k] = num;

                // ignore space
//                input.ignore();
            }

            // ignore other chars in current line
            input.getline(buffer, 256);
        }

        // close the file
        input.close();

        return true;
    }
};

#endif // ANIMATION_H
