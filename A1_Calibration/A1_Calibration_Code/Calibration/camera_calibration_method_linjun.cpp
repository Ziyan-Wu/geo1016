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
#include <cmath> 


using namespace easy3d;

float length(std::vector<float>& v) {
    float sum=0;
    for (int i=0; i < v.size(); i++){
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

float dot(std::vector<float>& a, std::vector<float>& b) {
    float sum = 0;
    for (int i=0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

std::vector<float> cross(std::vector<float>& a, std::vector<float>& b) {
    std::vector<float> result(3);
    result[0] = a[1] * b[2] - b[1] * a[2];
    result[1] = a[2] * b[0] - b[2] * a[0];
    result[2] = a[0] * b[1] - b[0] * a[1];
    return result;
    
}


float cotan(float i) { return(1 / std::tan(i)); }

/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 *
 * @param points_3d   A array of 3D points.
 * @param points_2d   A array of 2D points.
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
    if ((points_2d.size() != points_3d.size())&&(points_2d.size()>=6)) {
        return false;
    }
    else {
        // TODO: construct the P matrix (so P * m = 0).
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

        // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
        //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
        //             should be very close to your input images points.
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
        std::cout << "M - U * S * V^T: \n" << U * S * transpose(V) << std::endl;

        // TODO: extract intrinsic parameters from M.
        Matrix<double> A(3, 3, 0.0);
        Matrix<double> b(3, 1, 0.0);
        Matrix<double> M(3, 4, 0.0);
        //std::vector<double> v(12);
        for (int i = 0; i < 3; i++) {
            A[i][0] = V[i * 4][11];
            A[i][1] = V[i * 4 + 1][11];
            A[i][2] = V[i * 4 + 2][11];
            b[i][0] = V[i * 4 + 3][11];

            M[i][0] = V[i * 4][11];
            M[i][1] = V[i * 4 + 1][11];
            M[i][2] = V[i * 4 + 2][11];
            M[i][3] = V[i * 4 + 3][11];
        }
        
        //check if the M is correct
        for (int i=0; i < points_3d.size(); i++) {
            Matrix<double> pt(4, 1, 0.0);
            pt[0][0] = points_3d[i][0];
            pt[1][0] = points_3d[i][1];
            pt[2][0] = points_3d[i][2];
            pt[3][0] = 1;
            Matrix<double> multi = M * pt;
            std::cout << "\t" << i << ": (" << points_3d_[i] << ") <-> (" << multi[0][0]/multi[2][0]<<" "<< multi[1][0] / multi[2][0]<< ")" << std::endl;
        }

        std::cout << "A: \n" << A << std::endl;
        std::cout << "b: \n" << b << std::endl;
        std::vector<float> a1(3);
        std::vector<float> a2(3);
        std::vector<float> a3(3);
        for (int i = 0; i < 3; i++) {
            a1[i] = A[0][i];
            a2[i] = A[1][i];
            a3[i] = A[2][i];
        }


        float ro = 1 / length(a3);
        cx = pow(ro, 2) * (dot(a1, a3));
        cy = pow(ro, 2) * (dot(a2, a3));
        float cos = (dot(cross(a1, a3), cross(a2, a3))) / (length(cross(a1, a3)) * length(cross(a2, a3)));
        float angle = std::acos(cos);
        fx = pow(ro, 2) * length(cross(a1, a3)) * std::sin(angle);
        fy = pow(ro, 2) * length(cross(a2, a3)) * std::sin(angle);
        skew = -fx * cotan(angle);
        std::cout << "ro: \n" << ro << std::endl;

        // TODO: extract extrinsic parameters from M.
        Matrix<double> K(3, 3, 0.0);
        K[0][0] = fx;
        K[0][1] = fx * double(cotan(angle));
        K[0][2] = cx;
        K[1][1] = fy / std::sin(angle);
        K[1][2] = cy;
        K[2][2] = 1;


        std::vector<float> r1(3);
        std::vector<float> r2(3);
        std::vector<float> r3(3);
        r1 = cross(a2, a3);

        for (int i = 0; i < 3; i++) {
            r1[i] /= length(cross(a2, a3));
        }

        for (int i = 0; i < 3; i++) {
            r3[i] = ro * a3[i];
        }

        r2 = cross(r3, r1);
        //R[2] = 12;

        for (int i = 0; i < 3; i++) {
            R[3 * i] = r1[i];
            R[3 * i + 1] = r2[i];
            R[3 * i + 2] = r3[i];
        }
        std::cout << "R: \n" << R << std::endl;
        Matrix<double> invK(3, 3, 0.0);
        inverse(K, invK);
        std::cout << "K: \n" << K << std::endl;
        std::cout << "invK: \n" << invK << std::endl;
        Matrix<double> T = invK * b;

        std::cout << "T: \n" << T << std::endl;
        for (int i = 0; i < 3; i++) {
            t[i] = T[i][0] * ro;
        }
        return true;
    }
    // TODO: return true when testing your algorithm and in you final submission.
    //return false;


    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can do append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
    // https://en.cppreference.com/w/cpp/container/vector
    //std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 6, 1, 6, 2, 1, 5, 2, 3, 2, 7, 8, 1};
    //array.push_back(5); // append 5 to the array (so the size will increase by 1).
    //array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).



}

















