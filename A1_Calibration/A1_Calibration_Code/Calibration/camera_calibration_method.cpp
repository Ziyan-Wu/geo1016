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

#include "camera_calibration.h"
#include "matrix_algo.h"


using namespace easy3d;


float cotan(float i) { return(1 / std::tan(i)); }
/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 *
 * @param points_3d   An array of 3D points.
 * @param points_2d   An array of 2D points.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 *           - fx and fy: the focal length (in our slides, we use 'alpha' and 'beta'),
 *           - cx and cy: the principal point (in our slides, we use 'u0' and 'v0'),
 *           - skew:      the skew factor ('-alpha * cot_theta')
 *           - R:         the 3x3 rotation matrix encoding camera orientation.
 *           - t:         a 3D vector encoding camera location.
 */
bool CameraCalibration::calibration(
    const std::vector<vec3>& points_3d,
    const std::vector<vec2>& points_2d,
    float& fx, float& fy,
    float& cx, float& cy,
    float& skew,
    mat3& R,
    vec3& t)
{
    std::cout << "TODO: I am going to implement calibration() ..." << std::endl;

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points don't match) done
    // chaeck pts number (both need >= 6) and correspondence number (match)
    if (points_3d.size() != points_2d.size() || points_3d.size() < 6 || points_2d.size() < 6) return false;

    // point_2d should follow the rule of image coordinates (coord x/y should be int and >= 0)
    for (auto pt2d : points_2d){
        auto x = pt2d.x;
        auto y = pt2d.y;
        if (x < 0 ||
            y < 0 ||
            floor(x) - x > 0.001 ||
            floor(y) - y > 0.001){
            return false;
        }
    }

    // construct the P matrix (so P * m = 0).
    Matrix<double> P(2 * points_2d.size(), 12, 0.0);
    for (int i = 0; i < points_2d.size(); i++) {
        std::vector<double> p = { points_3d[i][0],points_3d[i][1],points_3d[i][2],1 };
        std::vector<double> p2 = { points_2d[i][0],points_2d[i][1] };
        for (int j = 0; j < 4; j++) {
            P[2 * i][j] = p[j];
            P[2 * i][j + 4] = 0;
            P[2 * i][j + 8] = -p2[0] * p[j];
            P[2 * i + 1][j] = 0;
            P[2 * i + 1][j + 4] = p[j];
            P[2 * i + 1][j + 8] = -p2[1] * p[j];
        }
    }
    std::cout << "P: \n" << P << std::endl;



    // solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    Matrix<double> U(2 * points_2d.size(), 2 * points_2d.size(), 0.0);   // initialized with 0s
    Matrix<double> S(2 * points_2d.size(), 12, 0.0);   // initialized with 0s
    Matrix<double> V(12, 12, 0.0);   // initialized with 0s
    svd_decompose(P, U, S, V);
    // Check 1: U is orthogonal, so U * U^T must be identity
    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    std::cout << "V: \n" << V << std::endl;

    // Check 3: S must be a diagonal matrix
    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
    std::cout << "P - U * S * V^T: \n" << P - U * S * transpose(V) << std::endl;


    Matrix<double> M(3, 4, 0.0);
    //std::vector<double> v(12);
    for (int i = 0; i < 3; i++) {

        M[i][0] = V[i * 4][11];
        M[i][1] = V[i * 4 + 1][11];
        M[i][2] = V[i * 4 + 2][11];
        M[i][3] = V[i * 4 + 3][11];
    }

    //check if the M is correct
    for (int i = 0; i < points_3d.size(); i++) {
        Matrix<double> pt(4, 1, 0.0);
        pt[0][0] = points_3d[i][0];
        pt[1][0] = points_3d[i][1];
        pt[2][0] = points_3d[i][2];
        pt[3][0] = 1;
        Matrix<double> multi = M * pt;
        std::cout << "\t" << i << ": (" << points_3d_[i] << ") <-> (" << multi[0][0] / multi[2][0] << " " << multi[1][0] / multi[2][0] << ")" << std::endl;
    }

    //extract intrinsic parameters from M.
    vec3 a1(M[0][0], M[0][1], M[0][2]);
    vec3 a2(M[1][0], M[1][1], M[1][2]);
    vec3 a3(M[2][0], M[2][1], M[2][2]);
    vec3 b(M[0][3], M[1][3], M[2][3]);
    std::cout << "a1: " << a1 << '\n';
    std::cout << "a3: " << a3 << '\n';
    std::cout << "a2: " << a2 << '\n';
    std::cout << "b: " << b << '\n';

    float ro = 1 / a3.length();
    cx = pow(ro, 2) * (dot(a1, a3));
    cy = pow(ro, 2) * (dot(a2, a3));
    float cos = -(dot(cross(a1, a3), cross(a2, a3))) / (cross(a1, a3).length() * cross(a2, a3).length());
    float angle = std::acos(cos);
    fx = pow(ro, 2) * cross(a1, a3).length() * std::sin(angle);
    fy = pow(ro, 2) * cross(a2, a3).length() * std::sin(angle);
    skew = -fx * cotan(angle);
    std::cout << "ro: \n" << ro << std::endl;


    //extract extrinsic parameters from M.
    mat3 K;
    K[0] = fx;
    K[1] = 0;
    K[2] = 0;
    K[3] = -fx * double(cotan(angle));
    K[4] = fy / std::sin(angle);
    K[5] = 0;
    K[6] = cx;
    K[7] = cy;
    K[8] = 1;
    std::cout << "K: \n" << K << '\n';

    vec3 r1 = cross(a2, a3) / cross(a2, a3).length();
    vec3 r3 = ro * a3;
    vec3 r2 = cross(r3, r1);


    for (int i = 0; i < 3; i++) {
        R[3 * i] = r1[i];
        R[3 * i + 1] = r2[i];
        R[3 * i + 2] = r3[i];
    }
    std::cout << "R: \n" << R << std::endl;
    mat3 invK = inverse(K);
    std::cout << "K: \n" << K << std::endl;
    std::cout << "invK: \n" << invK << std::endl;
    t = invK * b * ro;
    std::cout << "t: \n" << t << '\n';
    return true;
    // TODO: delete the above code in you final submission (which are just examples).

}

















