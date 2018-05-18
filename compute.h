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

#ifndef COMPUTE_H
#define COMPUTE_H

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>
//#include <vcg/space/triangle2.h>
//#include <vcg/space/triangle3.h>
//#include <QTime>

using namespace std;

/**
 * @brief getTriangle2TriangleMap Get the map from a trinagle t1 to another one t2.
 * @param t1 The first triangle to be mapped into t2.
 * @param t2 The triangle into which t1 will be mapped.
 * @param A The map to be computed.
 */
template < typename ScalarType, template <typename CoordType> class Point, int dim >
void getTriangle2TriangleMap(const Point<ScalarType> * const t1, const Point<ScalarType> * const t2, Eigen::Matrix<ScalarType,dim,dim> &A) {
    // construct the linear system with triangle t1
    Eigen::Matrix<ScalarType,dim*dim,dim*dim> M = Eigen::Matrix<ScalarType,dim*dim,dim*dim>::Zero();
    for (int k = 0; k < dim; k++)
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                M(k*dim+i,i*dim+j) = t1[k][j];
            }

    // construct the rhs of the linear system with triangle t2
    Eigen::Matrix<ScalarType,dim*dim,1> b;
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            b(i*dim+j) = t2[i][j];

    // create a solver and use it to solve for solution x
    Eigen::FullPivLU< Eigen::Matrix<ScalarType,dim*dim,dim*dim> > solver(M);
    Eigen::Matrix<ScalarType,dim*dim,1> x = solver.solve(b);

    // construct the map A with the solution x
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            A(i,j) = x(i*dim+j);
}

/**
 * @brief getSVD Gets the svd of a given square matrix (A = USV).
 * @param A Input matrix.
 * @param U Output U matrix.
 * @param V Output V matrix.
 * @param S Output vector (singular values).
 */
template < typename ScalarType, int dim >
void getSVD(const Eigen::Matrix<ScalarType,dim,dim> &A, Eigen::Matrix<ScalarType,dim,dim> &U, Eigen::Matrix<ScalarType,dim,dim> &V, Eigen::Matrix<ScalarType,dim,1> &S) {
    // create and compute svd
    Eigen::JacobiSVD< Eigen::Matrix<ScalarType,dim,dim> > svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // get the matrices
    U = svd.matrixU();
    V = svd.matrixV();
    S = svd.singularValues();
}

/**
 * @brief averageVectors
 * @param V
 * @param w
 * @param N
 * @return
 */
template < typename ValueType >
Eigen::Matrix<ValueType,2,1> averageVectors(const Eigen::Matrix<ValueType,Eigen::Dynamic,2> &V, const Eigen::Matrix<ValueType,Eigen::Dynamic,1> &w, const int N) {
    // check parameter
    assert(V.rows() == w.rows() && N > 0);

    // create a vector of angles
    Eigen::Matrix<ValueType,Eigen::Dynamic,1> angles(V.rows(), 1);

    // make angle mod 2pi/N by multiplying times N
    for (int i = 0; i < V.rows(); i++)
        angles(i) = std::atan2(V(i,1), V(i,0)) * N;

    // create vector of directions
    Eigen::Matrix<ValueType,Eigen::Dynamic,2> VV(V.rows(), 2);

    // compute directions
    for (int i = 0; i < V.rows(); i++) {
        VV(i,0) = std::cos(angles(i));
        VV(i,1) = std::sin(angles(i));
    }

    // average vector
    Eigen::Matrix<ValueType,2,1> R;
    R(0) = R(1) = ValueType(0);

    // compute average of the unit vectors
    for (int i = 0; i < VV.rows(); i++)
        R += VV.row(i) * w(i);
    R /= VV.rows();

    // scale them back
    ValueType a = std::atan2(R(1), R(0)) / N;
    R(0) = std::cos(a);
    R(1) = std::sin(a);

    return R;
}

//void test() {
//    vcg::Point3d op[3], oq[3], b1[3], b2[3], u;
//    vcg::Point2d p[2], q[2];
//    qsrand(QTime::currentTime().msec());
//    // generate t1 and t2
//    for (int i = 0; i < 3; i++) {
//        op[i] = vcg::Point3d(qrand()%10, qrand()%10, qrand()%10);
//        oq[i] = vcg::Point3d(qrand()%10, qrand()%10, qrand()%10);
//    }

//    // centrate t1 and t2
//    for (int i = 0; i < 3; i++) {
//        op[i] -= op[0];
//        oq[i] -= oq[0];
//    }

//    // print t1
//    std::cout << "t1:" << std::endl;
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 3; j++)
//            std::cout << op[j][i] << "\t";
//        std::cout << std::endl;
//    }
//    //print t2
//    std::cout << "t2:" << std::endl;
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 3; j++)
//            std::cout << oq[j][i] << "\t";
//        std::cout << std::endl;
//    }

//    // new reference base for t1
//    double den1 = (op[1]-op[0]).Norm();
//    b1[0] = (op[1] - op[0]) / den1;
//    b1[2] = ((op[1] - op[0]) ^ (op[2] - op[0])).Normalize();
//    b1[1] = b1[2] ^ b1[0];

//    // express t1 in new reference
//    p[0] = vcg::Point2d(den1,0);
//    p[1] = vcg::Point2d((op[2]-op[0]) * b1[0], (op[2]-op[0]) * b1[1]);

//    // new reference base for t2
//    double den2 = (oq[1]-oq[0]).Norm();
//    b2[0] = (oq[1] - oq[0]) / den2;
//    b2[2] = ((oq[1] - oq[0]) ^ (oq[2] - oq[0])).Normalize();
//    b2[1] = b2[2] ^ b2[0];

//    // express t2 in new reference
//    q[0] = vcg::Point2d(den2,0);
//    q[1] = vcg::Point2d((oq[2]-oq[0]) * b2[0], (oq[2]-oq[0]) * b2[1]);

//    // print t1
//    std::cout << "t1 (local):" << std::endl;
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 2; j++)
//            std::cout << p[j][i] << "\t";
//        std::cout << std::endl;
//    }
//    //print t2
//    std::cout << "t2 (local):" << std::endl;
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 2; j++)
//            std::cout << q[j][i] << "\t";
//        std::cout << std::endl;
//    }

//    // compute A (in 2D)
//    Eigen::Matrix2d A;
//    getTriangle2TriangleMap(p, q, A);

//    // print A
//    std::cout << "A:" << std::endl << A << std::endl;

//    // print A*t1 (should be equal to t2)
//    Eigen::Vector2d s[2], temp;
//    std::cout << "t1 mapped to t2:" << std::endl;
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 2; j++)
//            temp(j) = p[i][j];
//        s[i] = A * temp;
//    }
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 2; j++)
//            std::cout << s[j](i) << "\t";
//        std::cout << std::endl;
//    }

//    // compute P
//    Eigen::Matrix2d P, U;
//    U << A(0,0)-A(1,1), A(0,1)+A(1,0), A(1,0)+A(0,1), A(1,1)-A(0,0);
//    U = U / 2;
//    std::cout << "U:" << std::endl << U << std::endl;
//    double det = sqrt(fabs(U.determinant()));
//    std::cout << "Det: " << det << std::endl;
//    U = U / (det > 0 ? det : 1);
//    std::cout << "U:" << std::endl << U << std::endl;
//    std::cout << "U':" << std::endl << U.transpose() << std::endl;
//    P = U.transpose() * A;

//    // print P
//    std::cout << "P:" << std::endl << P << std::endl;

//    // compute svd (in 2D) of P
//    Eigen::Matrix2d V;
//    Eigen::Vector2d S;
//    getSVD(P, U, V, S);

//    // print S
//    std::cout << "S:" << std::endl << S << std::endl;

//    // print dot product of the columns of U
//    std::cout << "<U(:,0),U(:,1)>: " << U.col(0).dot(U.col(1)) << std::endl;

//    // map A at two cols of U
//    Eigen::Vector2d U0 = A * U.col(0);
//    Eigen::Vector2d U1 = A * U.col(1);

//    // print U0 and U1
//    std::cout << "U0 = A*U(:,0):" << std::endl << U0 << std::endl;
//    std::cout << "U1 = A*U(:,1):" << std::endl << U1 << std::endl;

//    // print dot product of U0 and U1
//    std::cout << "<U0,U1>: " << U0.dot(U1) << std::endl;

//    // trasform in glabal space
//    u = b1[0]*U(0,0) + b1[1]*U(1,0);
//    std::cout << "u:" << std::endl;
//    for (int i = 0; i < 3; i++)
//        std::cout << u[i] << std::endl;

//    // check u is coplanar with t1 or t2
//    std::cout << "u is coplanar with t1 (should be 0): " << u*b1[2] << std::endl;
//}

#endif // COMPUTE_H
