/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <cmath>
#include <algorithm>
#include<iostream>


using namespace easy3d;


/// convert a 3 by 3 matrix of type 'Matrix<double>' to mat3
mat3 to_mat3(Matrix<double> &M) {
    mat3 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}


/// convert M of type 'matN' (N can be any positive integer) to type 'Matrix<double>'
template<typename mat>
Matrix<double> to_Matrix(const mat &M) {
    const int num_rows = M.num_rows();
    const int num_cols = M.num_columns();
    Matrix<double> result(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}

//===========================================================================================
/// my code
//void Normalize(const std::vector<vec3> &points, std::vector<vec3> &nor_points, Matrix<double> &T) {
//    //第一步得到所有特征点的均值,并将所有点的均值为0
//    float meanX = 0;
//    float meanY = 0;
//    for (int i = 0; i < points.size(); i++) {
//        meanX += points[i].x;
//        meanY += points[i].y;
//    }
//    meanX /= points.size();
//    meanY /= points.size();
//
//    //第二步将所有点到原点的距离为根号2
//    float meanDevX = 0;
//    float meanDevY = 0;
//    for (int i = 0; i < points.size(); i++) {
//        nor_points[i].x = points[i].x - meanX;
//        nor_points[i].y = points[i].y - meanY;
//        meanDevX += fabs(nor_points[i].x); //fabs是求一个实数的绝对值  点到原点距离的累加
//        meanDevY += fabs(nor_points[i].y);
//    }
//    meanDevX /= points.size(); //点到原点距离的平均值
//    meanDevY /= points.size();
//    for (int i = 0; i < points.size(); i++) {
//        nor_points[i].x /= meanDevX;
//        nor_points[i].y /= meanDevY;
//    }
//
//    //用于还原特征点到原始的坐标系,获得矩阵   |sX  0  -meanx*sX|  用于取逆     x    快速还原
//    //                                  |0   sY -meany*sY|          *   y
//    //                                  |0   0      1    |              1
//    float sX = 1.0 / meanDevX;
//    float sY = 1.0 / meanDevY;
//    T(0, 0) = sX;
//    T(0, 1) = 0;
//    T(0, 2) = -meanX * sX;
//    T(1, 0) = sY;
//    T(1, 1) = sY;
//    T(1, 2) = -meanY * sY;
//    T(2, 0) = 0;
//    T(2, 1) = 0;
//    T(2, 2) = 1;

//    Matrix<double> T (3, 3, 0.0);
//    T[0][0] = sX;
//    T[1][1] = sY;
//    T[0][2] = -meanX * sX;
//    T[1][2] = -meanY * sY;
//    T[2][2] = 1;
//}

Matrix<double> normalized_matrix(std::vector<vec3> points, double size) {
    double sum_x = 0.0, sum_y = 0.0, ave_x, ave_y;
    for (int i = 0; i < size; i++) {
        sum_x = sum_x + points[i].x;
        sum_y = sum_y + points[i].y;
    }
    ave_x = sum_x / size;
    ave_y = sum_y / size;
//    std::cout << "sum_x\n" << ave_x << '\n';
//    std::cout << "sum_y\n" << ave_y << '\n';
    double temp = 0.0;
    // variance
    for (int i = 0; i < size; i++) {
        temp = temp + pow((points[i].x - ave_x), 2) + pow((points[i].y - ave_y), 2);
    }
    // std
    double std = pow(temp / size, 0.5);
    // √2
    double scale = pow(2, 0.5) / std;
//    std::cout << "scale\n" << scale << '\n';
    Matrix<double> T(3, 3, 0.0);
//    T[0][0] = scale;
//    T[1][0] = 0;
//    T[2][1] = 0;
//    T[0][1] = 0;
//    T[1][1] = scale;
//    T[2][2] = 0;
//    T[0][2] = -scale * ave_x;
//    T[1][2] = -scale * ave_y;
//    T[2][2] = 1;

    T[0][0] = scale;
    T[1][0] = 0;
    T[2][1] = 0;
    T[0][1] = 0;
    T[1][1] = scale;
    T[2][2] = 0;
    T[0][2] = -scale * ave_x;
    T[1][2] = -scale * ave_y;
    T[2][2] = 1;
    std::cout << "T\n" << T << '\n';
    return T;
}

std::vector<vec3> normalization(Matrix<double> T, std::vector<vec3> points, double size) {
    std::vector<vec3> nor_points;
    for (int i = 0; i < size; i++) {
        nor_points.push_back(to_mat3(T) * points[i]);
    }
//    std::cout << "new: \n" << nor_points << "\n";
    return nor_points;
}

std::vector<vec3> get_point3d(std::vector<vec3> points0, std::vector<vec3> points1,
                              Matrix<double> &K, Matrix<double> &M1,
                              Matrix<double> &R, Matrix<double> &t,
                              double size) {
    //extrinsic matrix
    //R1t1
    Matrix<double> Rt(3, 4, 0.0);
    Rt.set_column(R.get_column(0), 0);
    Rt.set_column(R.get_column(1), 1);
    Rt.set_column(R.get_column(2), 2);
    Rt.set_column(t.get_column(0), 3);
    std::cout << "Rt \n" << Rt << '\n';
    // second camera matrix
    Matrix<double> M2(3, 4, 0.0);
    M2 = K * Rt;
    std::cout << "M2 \n" << K * Rt << '\n';
    //create AA
    Matrix<double> AA(4, 4, 0.0);
    std::vector<vec3> point_vec3;
    for (int i = 0; i < size; i++) {
        AA.set_row(points0[i].x * M1.get_row(2) - M1.get_row(0), 0);
        AA.set_row(points0[i].y * M1.get_row(2) - M1.get_row(1), 1);
        AA.set_row(points1[i].x * M2.get_row(2) - M2.get_row(0), 2);
        AA.set_row(points1[i].y * M2.get_row(2) - M2.get_row(1), 3);
        //SVD
//        std::cout << "---- Forth SVD of AA \n" << '\n';
        Matrix<double> Ua(4, 4, 0.0);   // initialized with 0s
        Matrix<double> Sa(4, 4, 0.0);   // initialized with 0s
        Matrix<double> Va(4, 4, 0.0);   // initialized with 0s
        svd_decompose(AA, Ua, Sa, Va);
//        std::cout << "Va \n" << Va << '\n';
        // get last column
        std::vector<double> Va_last_col = Va.get_row(3);
//        std::cout << "Va_last_col \n" << Va_last_col << '\n';
        // change last element to 1
        vec3 one_point;
        one_point[0] = Va_last_col[0] / Va_last_col[3];
        one_point[1] = Va_last_col[1] / Va_last_col[3];
        one_point[2] = Va_last_col[2] / Va_last_col[3];
//        one_point[3] = Va_last_col[3] / Va_last_col[3];
//        std::cout << "one_point \n" << one_point << '\n';
        point_vec3.push_back(one_point);
    }
//    std::cout << "point_vec3 \n" << point_vec3 << '\n';
    return point_vec3;
}


int count_z(std::vector<vec3> points) {
    int cnt = 0;
    for (int i = 0; i < points.size(); i++) {
        if (points[i][2] > 0.0) {
            cnt = cnt + 1;
        }
    }
    return cnt;
}


//===========================================================================================

/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'.
 */
bool Triangulation::triangulation(
        float fx, float fy,     /// input: the focal lengths (same for both cameras)
        float cx, float cy,     /// input: the principal point (same for both cameras)
        const std::vector<vec3> &points_0,    /// input: image points (in homogenous coordinates) in the 1st image.
        const std::vector<vec3> &points_1,    /// input: image points (in homogenous coordinates) in the 2nd image.
        std::vector<vec3> &points_3d,         /// output: reconstructed 3D points
        mat3 &R,   /// output: recovered rotation of 2nd camera (used for updating the viewer and visual inspection)
        vec3 &t    /// output: recovered translation of 2nd camera (used for updating the viewer and visual inspection)
) const {
//  ===== = ======================================================================================
//
//  1.1 Normalization
    double size = points_0.size();
//    std::vector<vec3> points_0_nor(points_0.begin(), points_0.end());
//    std::vector<vec3> points_1_nor(points_1.begin(), points_1.end());
//    Matrix<double> T0(3,3,0.0);
//    Matrix<double> T1(3,3,0.0);
//    Normalize(points_0, points_0_nor, T0);
//    Normalize(points_0, points_1_nor, T1);

    Matrix<double> T0 = normalized_matrix(points_0, size);
    Matrix<double> T1 = normalized_matrix(points_1, size);
    std::vector<vec3> points_0_nor = normalization(T0, points_0, size);
    std::vector<vec3> points_1_nor = normalization(T1, points_1, size);

    for (auto &i : points_1_nor) {
        std::cout << "points_1_nor\n" << i << '\n';
    }

//  1.2 Linear solution
    Matrix<double> A(int(size), 9, 1.0);
    for (int i = 0; i < int(size); i++) {
        A[i][0] = points_0_nor[i].x * points_1_nor[i].x;
        A[i][1] = points_0_nor[i].y * points_1_nor[i].x;
        A[i][2] = points_1_nor[i].x;
        A[i][3] = points_0_nor[i].x * points_1_nor[i].y;
        A[i][4] = points_0_nor[i].y * points_1_nor[i].y;
        A[i][5] = points_1_nor[i].y;
        A[i][6] = points_0_nor[i].x;
        A[i][7] = points_0_nor[i].y;
    }

    std::cout << "---- first SVD \n" << '\n';
    Matrix<double> U(int(size), int(size), 0.0);   // initialized with 0s
    Matrix<double> S(int(size), 9, 0.0);   // initialized with 0s
    Matrix<double> V(9, 9, 0.0);   // initialized with 0s
    svd_decompose(A, U, S, V);
    const std::vector<double> f_hat = V.get_column(9 - 1); //9 rows, 1 col
    std::cout << "f_hat\n" << f_hat << '\n';
//    std::cout << "f_hat.size()\n" << f_hat.size() << '\n';
    Matrix<double> F_hat(3, 3, 0.0);
    F_hat[0][0] = f_hat[0];
    F_hat[0][1] = f_hat[1];
    F_hat[0][2] = f_hat[2];
    F_hat[1][0] = f_hat[3];
    F_hat[1][1] = f_hat[4];
    F_hat[1][2] = f_hat[5];
    F_hat[2][0] = f_hat[6];
    F_hat[2][1] = f_hat[7];
    F_hat[2][2] = f_hat[8];
    std::cout << "F_hat\n" << F_hat << '\n';

//  1.3 Constraint enforcement
    std::cout << "---- Second SVD \n" << '\n';
    Matrix<double> Uf(3, 3, 0.0);   // initialized with 0s
    Matrix<double> Sf(3, 3, 0.0);   // initialized with 0s
    Matrix<double> Vf(3, 3, 0.0);   // initialized with 0s
    svd_decompose(F_hat, Uf, Sf, Vf);
    std::cout << "Sf \n" << Sf << '\n';
    Sf(2, 2) = 0;
    std::cout << "Sf_after \n" << Sf << '\n';

    Matrix<double> F;
    F = Uf * Sf * Vf;
    std::cout << "F\n" << F << '\n';

//  1.4 Denormalization
    Matrix<double> F_de;
    F_de = T1.transpose() * F * T0;
    std::cout << "F_de\n" << F_de << '\n';

//  =================================================================================
//  2.1 get Essential matrix E
    Matrix<double> K(3, 3, 0.0);  //  intrinsic matrix
    K[0][0] = fx;
    K[1][0] = 0;
    K[2][0] = 0;
    K[0][1] = 0;
    K[1][1] = fy;
    K[2][1] = 0;
    K[0][2] = cx;
    K[1][2] = cy;
    K[2][2] = 1;
    std::cout << "intrinsic matrix K\n" << K << '\n';
    Matrix<double> E;
    E = K.transpose() * F_de * K;
    std::cout << "E\n" << E << '\n';

//  2.2 get camera matrix
    std::cout << "---- Third SVD of Essential matrix E \n" << '\n';
    Matrix<double> Ue(3, 3, 0.0);   // initialized with 0s
    Matrix<double> Se(3, 3, 0.0);   // initialized with 0s
    Matrix<double> Ve(3, 3, 0.0);   // initialized with 0s
    svd_decompose(E, Ue, Se, Ve);
    std::cout << "Se \n" << Se << '\n';

    //set W and Z
    Matrix<double> W(3, 3, 0.0);
    W[1][0] = 1.0;
    W[0][1] = -1;
    W[2][2] = 1;
    std::cout << "W \n" << W << '\n';
    Matrix<double> Z(3, 3, 0.0);
    Z[1][0] = -1.0;
    Z[0][1] = 1;
    std::cout << "Z \n" << Z << '\n';

    // R1 and R2
    Matrix<double> R1(3, 3, 0.0);
    R1 = determinant(Ue * W * Ve) * Ue * W * Ve;
    std::cout << "R1 \n" << R1 << '\n';
    Matrix<double> R2(3, 3, 0.0);
    R2 = determinant(Ue * W.transpose() * Ve) * Ue * W.transpose() * Ve;
    std::cout << "R2 \n" << R2 << '\n';

    // t1 and t2
    Matrix<double> temp(3, 1, 0.0);
    temp[0][2] = 1.0;
    std::cout << "temp \n" << temp << '\n';
//    std::cout << "Ue.last_column \n" << Ue.get_column(2) << '\n';
    Matrix<double> t1(3, 1, 0.0);
    t1 = Ue * temp;
    std::cout << "t1 \n" << t1 << '\n';
    Matrix<double> t2(3, 1, 0.0);
    t2 = -Ue * temp;
    std::cout << "t2 \n" << t2 << '\n';

    //extrinsic matrix
    //R1t1

//====================================================================================== = ==
//  3.1
    // first camera matrix
    Matrix<double> M1(3, 4, 0.0);
    M1[0][0] = 1;
    M1[1][1] = 1;
    M1[2][2] = 1;
    std::cout << "M1 \n" << M1 << '\n';
    // calculate number of positive_z
    int cnt11 = count_z(get_point3d(points_0, points_1, K, M1, R1, t1, size));
    int cnt12 = count_z(get_point3d(points_0, points_1, K, M1, R1, t2, size));
    int cnt21 = count_z(get_point3d(points_0, points_1, K, M1, R2, t1, size));
    int cnt22 = count_z(get_point3d(points_0, points_1, K, M1, R2, t2, size));
    std::cout << "cnt11 \n" << cnt11 << '\n';
    std::cout << "cnt12 \n" << cnt12 << '\n';
    std::cout << "cnt21 \n" << cnt21 << '\n';
    std::cout << "cnt22 \n" << cnt22 << '\n';


    int max_cnt = std::max(cnt22, std::max(cnt21, std::max(cnt11, cnt12)));
    if (max_cnt == cnt11) {
        R = to_mat3(R1);
        t = vec3(t1[0][0], t1[1][0], t1[2][0]);
        points_3d = get_point3d(points_0, points_1, K, M1, R1, t1, size);
    }
    if (max_cnt == cnt12) {
        R = to_mat3(R1);
        t = vec3(t2[0][0], t2[1][0], t2[2][0]);
        points_3d = get_point3d(points_0, points_1, K, M1, R1, t2, size);
    }
    if (max_cnt == cnt21) {
        R = to_mat3(R2);
        t = vec3(t1[0][0], t1[1][0], t1[2][0]);
        points_3d = get_point3d(points_0, points_1, K, M1, R2, t1, size);
    }
    if (max_cnt == cnt22) {
        R = to_mat3(R2);
        t = vec3(t2[0][0], t2[1][0], t2[2][0]);
        points_3d = get_point3d(points_0, points_1, K, M1, R2, t2, size);
    }

//  ===== = ======================================================================================

    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for each sub-task. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or feel free to put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any data structure and function offered by Easy3D, in particular the following two\n"
                 "\tfiles for vectors and matrices:\n"
                 "\t    - easy3d/core/mat.h  Fixed-size matrices and related functions.\n"
                 "\t    - easy3d/core/vec.h  Fixed-size vectors and related functions.\n"
                 "\tFor matrices with unknown sizes (e.g., when handling an unknown number of corresponding points\n"
                 "\tstored in a file, where their sizes can only be known at run time), a dynamic-sized matrix data\n"
                 "\tstructure is necessary. In this case, you can use the templated 'Matrix' class defined in\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\tPlease refer to the corresponding header files for more details of these data structures.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please refer to\n"
                 "\t'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations. \n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (original code framework + your implementation).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without any modification.\n\n"
              << std::flush;

    /// Easy3D provides fixed-size matrix types, e.g., mat2 (2x2), mat3 (3x3), mat4 (4x4), mat34 (3x4).
    /// To use these matrices, their sizes should be known to you at the compile-time (i.e., when compiling your code).
    /// Once defined, their sizes can NOT be changed.
    /// In 'Triangulation/matrix.h', another templated 'Matrix' type is also provided. This type can have arbitrary
    /// dimensions and their sizes can be specified at run-time (i.e., when executing your program).
    /// Below are a few examples showing some of these data structures and related APIs.

    /// ----------- fixed-size matrices

    /// define a 3 by 4 matrix M (you can also define 3 by 4 matrix similarly)
    mat34 M(1.0f);  /// entries on the diagonal are initialized to be 1 and others to be 0.

    /// set the first row of M
    M.set_row(0, vec4(1, 1, 1, 1));    /// vec4 is a 4D vector.

    /// set the second column of M
    M.set_col(1, vec4(2, 2, 2, 2));

    /// get the 3 rows of M
//    vec4 M1 = M.row(0);
    vec4 M2 = M.row(1);
    vec4 M3 = M.row(2);


    /// ----------- fixed-size vectors

    /// how to quickly initialize a std::vector
    std::vector<double> rows = {0, 1, 2, 3,
                                4, 5, 6, 7,
                                8, 9, 10, 11};
    /// get the '2'-th row of M
    const vec4 b = M.row(2);    // it assigns the requested row to a new vector b

    /// get the '1'-th column of M
    const vec3 c = M.col(1);    // it assigns the requested column to a new vector c

    /// modify the element value at row 2 and column 1 (Note the 0-based indices)
    M(2, 1) = b.x;

    /// apply transformation M on a 3D point p (p is a 3D vector)
    vec3 p(222, 444, 333);
    vec3 proj = M * vec4(p, 1.0f);  // use the homogenous coordinates. result is a 3D vector

    /// the length of a vector
    float len = p.length();
    /// the squared length of a vector
    float sqr_len = p.length2();

    /// the dot product of two vectors
    float dot_prod = dot(p, proj);

    /// the cross product of two vectors
//    vec3 cross_prod = cross(p, proj);
//    std::cout << "before\n" << cross_prod << "\n";
//
//    /// normalize this vector
//    cross_prod.normalize();
//    std::cout << "after\n" << cross_prod << "\n";

    /// a 3 by 3 matrix (all entries are intentionally NOT initialized for efficiency reasons)
    mat3 eF;
    /// ... here you compute or initialize F.
    /// compute the inverse of K
    mat3 invF = inverse(eF);

    /// ----------- dynamic-size matrices

    /// define a non-fixed size matrix
//    Matrix<double> A(2, 3, 0.0); // all entries initialized to 0.0.

    /// set its first row by a 3D vector (1.1, 2.2, 3.3)
//    A.set_row({1.1, 2.2, 3.3}, 0);   // here "{ 1.1, 2.2, 3.3 }" is of type 'std::vector<double>'
//    std::cout << "A\n" << A << "\n";
//    std::cout << "A[0][1]\n" << A[0][1] << "\n";

    /// get the last column of a matrix
    std::vector<double> last_column = A.get_column(A.cols() - 1);

    // TODO: delete all above demo code in the final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you to check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       However, there are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}
