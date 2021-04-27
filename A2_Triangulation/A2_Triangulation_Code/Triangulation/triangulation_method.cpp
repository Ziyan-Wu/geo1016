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
Matrix<double> normalized_matrix(std::vector<vec3> points, double size) {
//    1. cooordinate / 1
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
    for (int i = 0; i < size; i++) {
        temp += distance(vec2(points[i]), vec2(ave_x, ave_y));
//        temp = temp + pow((points[i].x - ave_x), 2) + pow((points[i].y - ave_y), 2);
    }
    temp /= size;
    // use sqrt not pow(0.5)
    double scale = std::sqrt(2) / temp;
    Matrix<double> T(3, 3, 0.0);
    T[0][0] = scale;
    T[1][0] = 0;
    T[2][0] = 0;
    T[0][1] = 0;
    T[1][1] = scale;
    T[2][1] = 0;
    T[0][2] = -scale * ave_x;
    T[1][2] = -scale * ave_y;
    T[2][2] = 1;
    std::cout << "T\n" << T << '\n';
    return T;
}

std::vector<vec2> normalization(Matrix<double> T, std::vector<vec3> points, double size) {
    std::vector<vec2> nor_points;
    for (int i = 0; i < size; i++) {
        vec3 p = to_mat3(T) * points[i];
        // convert to non-Homogeneous coordinates
        vec2 p2(p.x / p.z, p.y / p.z);
        nor_points.push_back(p2);
    }
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
        std::vector<double> Va_last_col = Va.get_column(3);
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


int count_z(std::vector<vec3> points, Matrix<double> &R, Matrix<double> &t) {
    int cnt = 0;
    for (int i = 0; i < points.size(); i++) {
        mat4 temp_R(to_mat3(R));
//        mat4 temp_t(to_mat3(t));
        mat4 temp_t = mat4::identity();
        temp_t(0, 3) = t(0, 0);
        temp_t(1, 3) = t(1, 0);
        temp_t(2, 3) = t(2, 0);

        vec3 p = temp_t * temp_R * points[i];
        if ((points[i][2] > 0.0) && (p.z > 0.0)) {
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
//  check
//  1.1 Normalization
    double size = points_0.size();
    Matrix<double> T0 = normalized_matrix(points_0, size);
    Matrix<double> T1 = normalized_matrix(points_1, size);
    std::vector<vec2> points_0_nor = normalization(T0, points_0, size);
    std::vector<vec2> points_1_nor = normalization(T1, points_1, size);
    // print points_0_nor
//    for (auto &i : points_1_nor) {
//        std::cout << "points_0_nor\n" << i << '\n';
//    }

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
//    std::cout << "A \n" << A << '\n';

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
    F = Uf * Sf * Vf.transpose();
    std::cout << "F\n" << F << '\n';

//  1.4 Denormalization
    Matrix<double> F_de;
    F_de = T1.transpose() * F * T0;
    std::cout << "F_de\n" << F_de << '\n';
    // F is up to scale, divide all elements by the last one
    F_de(0, 0) = F_de(0, 0) / F_de(2, 2);
    F_de(0, 1) = F_de(0, 1) / F_de(2, 2);
    F_de(0, 2) = F_de(0, 2) / F_de(2, 2);
    F_de(1, 0) = F_de(1, 0) / F_de(2, 2);
    F_de(1, 1) = F_de(1, 1) / F_de(2, 2);
    F_de(1, 2) = F_de(1, 2) / F_de(2, 2);
    F_de(2, 0) = F_de(2, 0) / F_de(2, 2);
    F_de(2, 1) = F_de(2, 1) / F_de(2, 2);
    F_de(2, 2) = F_de(2, 2) / F_de(2, 2);
    std::cout << "F_de after up to scale\n" << F_de << '\n';


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
    E = K.transpose() * F_de * K; // K.transpose() not K
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
    R1 = determinant(Ue * W * Ve) * Ue * W * Ve.transpose(); // Ve.transpose() not Ve
    std::cout << "R1 \n" << R1 << '\n';
    Matrix<double> R2(3, 3, 0.0);
    R2 = determinant(Ue * W.transpose() * Ve) * Ue * W.transpose() * Ve.transpose(); // Ve.transpose() not Ve
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
    M1 = K * M1;
    std::cout << "M1 \n" << M1 << '\n';
    // calculate number of positive_z
    int cnt11 = count_z(get_point3d(points_0, points_1, K, M1, R1, t1, size), R1, t1);
    int cnt12 = count_z(get_point3d(points_0, points_1, K, M1, R1, t2, size), R1, t2);
    int cnt21 = count_z(get_point3d(points_0, points_1, K, M1, R2, t1, size), R2, t1);
    int cnt22 = count_z(get_point3d(points_0, points_1, K, M1, R2, t2, size), R2, t2);

    std::vector<vec3> points_3d_2;

    std::cout << "cnt11 \n" << cnt11 << '\n';
    std::cout << "cnt12 \n" << cnt12 << '\n';
    std::cout << "cnt21 \n" << cnt21 << '\n';
    std::cout << "cnt22 \n" << cnt22 << '\n';

    int max_cnt = std::max(cnt22, std::max(cnt21, std::max(cnt11, cnt12)));
    if (max_cnt == cnt11) {
        R = to_mat3(R1);
        t = vec3(t1[0][0], t1[1][0], t1[2][0]);
        points_3d = get_point3d(points_0, points_1, K, M1, R1, t1, size);
        std::cerr << "cnt11 " << cnt11 << std::endl;
        std::cerr << "R " << R << std::endl;
        std::cerr << "t " << t << std::endl;
    }
    if (max_cnt == cnt12) {
        R = to_mat3(R1);
        t = vec3(t2[0][0], t2[1][0], t2[2][0]);
        points_3d = get_point3d(points_0, points_1, K, M1, R1, t2, size);
        std::cerr << "cnt12 " << cnt12 << std::endl;
        std::cerr << "R " << R << std::endl;
        std::cerr << "t " << t << std::endl;
    }
    if (max_cnt == cnt21) {
        R = to_mat3(R2);
        t = vec3(t1[0][0], t1[1][0], t1[2][0]);
        points_3d = get_point3d(points_0, points_1, K, M1, R2, t1, size);
        std::cerr << "cnt21 " << cnt21 << std::endl;
        std::cerr << "R " << R << std::endl;
        std::cerr << "t " << t << std::endl;
    }
    if (max_cnt == cnt22) {
        R = to_mat3(R2);
        t = vec3(t2[0][0], t2[1][0], t2[2][0]);
        points_3d = get_point3d(points_0, points_1, K, M1, R2, t2, size);
        std::cerr << "cnt22 " << cnt22 << std::endl;
        std::cerr << "R " << R << std::endl;
        std::cerr << "t " << t << std::endl;
    }

//     print 3d points
//    for (auto &i : points_3d) {
//        std::cout << "points_3d\n" << i << '\n';
//    }
    return points_3d.size() > 0;
}
