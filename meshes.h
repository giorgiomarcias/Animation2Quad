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

#ifndef MESHES_H
#define MESHES_H

// vcg mesh imports
#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>

using namespace vcg;

/// define a triangle mesh type
class TVertex;
class TFace;
class TUsedTypes : public UsedTypes< Use<TVertex>::AsVertexType, Use<TFace>::AsFaceType > {};
class TVertex: public Vertex<   TUsedTypes,
                                vertex::Coord3f,
                                vertex::Normal3f,
                                vertex::Curvaturef,
                                vertex::CurvatureDirf,
                                vertex::VFAdj,
                                vertex::Qualityf,
                                vertex::BitFlags
                            > {};
class TFace : public Face   <   TUsedTypes,
                                face::VertexRef,
                                face::VFAdj,
                                face::FFAdj,
                                face::Normal3f,
                                face::CurvatureDirf,
                                face::Qualityf,
                                face::Color4b,
                                face::BitFlags
                            > {};
class TMesh : public tri::TriMesh< std::vector<TVertex>, std::vector<TFace> > {};

/// define a polygon mesh type
class PVertex;
class PFace;
class PUsedTypes : public UsedTypes< Use<PVertex>::AsVertexType, Use<PFace>::AsFaceType > {};
class PVertex : public Vertex<  PUsedTypes,
                                vertex::Coord3f,
                                vertex::Normal3f,
                                vertex::Mark,
                                vertex::BitFlags
                             > {};
class PFace : public Face<  PUsedTypes,
                            face::PolyInfo,
                            face::PFVAdj,
                            face::PFFAdj,
                            face::PFHAdj,
                            face::Normal3f,
                            face::BitFlags
                         > {};
class PMesh : public tri::TriMesh< std::vector<PVertex>, std::vector<PFace> > {};

/// define a simple enumeration of types
enum MeshTypes {
    NO_MESHTYPE,
    TRI_MESHTYPE,
    POLY_MESHTYPE
};

#endif // MESHES_H
