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

#define PI acos(-1)


using namespace easy3d;


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
        const std::vector<vec3> &points_3d,
        const std::vector<vec2> &points_2d,
        float &fx, float &fy,
        float &cx, float &cy,
        float &skew,
        mat3 &R,
        vec3 &t) {
    std::cout << std::endl;
    std::cout << "TODO: I am going to implement the calibration() function in the following file:" << std::endl
              << "\t" << __FILE__ << std::endl;
    std::cout << "TODO: After implementing the calibration() function, I will disable all unrelated output ...\n\n";

    // check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    if ((points_2d_.size() == points_3d_.size()) || (points_2d_.size() + points_3d_.size()) >= 30) {

        // 1.
        // construct the P matrix (so P * m = 0).
        int n = points_3d_.size();
        const int Pm = 2 * n, Pn = 12;

        std::vector<double> array_P;
        for (int i = 0; i < n; i++) {
            array_P.push_back(points_3d[i][0]);
            array_P.push_back(points_3d[i][1]);
            array_P.push_back(points_3d[i][2]);
            array_P.push_back(1);
            array_P.push_back(0);
            array_P.push_back(0);
            array_P.push_back(0);
            array_P.push_back(0);
            array_P.push_back(-points_2d[i][0] * points_3d[i][0]);
            array_P.push_back(-points_2d[i][0] * points_3d[i][1]);
            array_P.push_back(-points_2d[i][0] * points_3d[i][2]);
            array_P.push_back(-points_2d[i][0] * 1);
            array_P.push_back(0);
            array_P.push_back(0);
            array_P.push_back(0);
            array_P.push_back(0);
            array_P.push_back(points_3d[i][0]);
            array_P.push_back(points_3d[i][1]);
            array_P.push_back(points_3d[i][2]);
            array_P.push_back(1);
            array_P.push_back(-points_2d[i][1] * points_3d[i][0]);
            array_P.push_back(-points_2d[i][1] * points_3d[i][1]);
            array_P.push_back(-points_2d[i][1] * points_3d[i][2]);
            array_P.push_back(-points_2d[i][1] * 1);
        }
        Matrix<double> P(Pm, Pn, array_P.data());
        std::cout << "P matrix: \n" << P << '\n';

        // 2.
        // solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
        //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
        //             should be very close to your input images points.
        Matrix<double> U(Pm, Pm, 0.0);   // initialized with 0s
        Matrix<double> S(Pm, Pn, 0.0);   // initialized with 0s
        Matrix<double> V(Pn, Pn, 0.0);   // initialized with 0s
        // Compute the SVD decomposition of P
        svd_decompose(P, U, S, V);
//        std::cout<<"V matrix: \n"<< V<<'\n';
        const std::vector<double> m = V.get_column(Pn - 1);
        std::cout << "m matrix: \n" << m << '\n';
//        std::cout << "m1 matrix: \n" << m.size() << '\n';
//        std::cout << "m2 matrix: \n" << V.get_column(Pn - 1).size() << '\n';

        // 3.
        // extract intrinsic parameters from M.
        vec3 a1(m[0], m[1], m[2]);
        vec3 a2(m[4], m[5], m[6]);
        vec3 a3(m[8], m[9], m[10]);
        vec3 b(m[3], m[7], m[11]);
        std::cout << "a1: " << a1 << '\n';
        std::cout << "a3: " << a3 << '\n';
        std::cout << "a2: " << a2 << '\n';
        std::cout << "b: " << b << '\n';

//        std::cout<<a3<<'\n';
        double costheta = -dot(cross(a1, a3), cross(a2, a3)) / (cross(a1, a3).length() * cross(a2, a3).length());
        double theta = acos(costheta) * 180 / PI;
        std::cout << "theta: " << theta << '\n';
        double sintheta = sin(theta); //pow(1 - pow(costheta, 2), 0.5);
        std::cout << "sintheta: " << sintheta << '\n';
        double p = 1 / a3.length();
        std::cout << "p: " << p << '\n';

        cx = pow(p, 2) * dot(a1, a3);
        cy = pow(p, 2) * dot(a2, a3);
        fx = pow(p, 2) * cross(a1, a3).length() * sintheta;
        fy = pow(p, 2) * cross(a2, a3).length() * sintheta;
        std::cout << "cx   : " << cx << '\n';
        std::cout << "cy   : " << cy << '\n';
        std::cout << "fx: " << fx << '\n';
        std::cout << "fy : " << fy << '\n';

        // 4.
        // extract extrinsic parameters from M.
        vec3 r1 = cross(a2, a3) / cross(a2, a3).length();
        vec3 r3 = p * a3; // +-???
        vec3 r2 = cross(r3, r1);
        std::cout << "r1: \n" << r1 << '\n';
        std::cout << "r3: \n" << r3 << '\n';
        std::cout << "r2: \n" << r2 << '\n';
        R;
        R[0] = r1[0];
        R[1] = r2[0];
        R[2] = r3[0];
        R[3] = r1[1];
        R[4] = r2[1];
        R[5] = r3[1];
        R[6] = r1[2];
        R[7] = r2[2];
        R[8] = r3[2];
        std::cout << "R: \n" << R << '\n';

        skew = -fx * costheta / sintheta;
        mat3 K;
        K[0] = fx;
        K[1] = 0;
        K[2] = 0;
        K[3] = -fx * costheta / sintheta;
        K[4] = fy / sintheta;
        K[5] = 0;
        K[6] = cx;
        K[7] = cy;
        K[8] = 1;
        std::cout << "K: \n" << K << '\n';
        mat3 invK = inverse(K);
        std::cout << "inverse(K): \n" << invK << '\n';
        std::cout << "K * invK: \n" << K * invK << std::endl;
        t = inverse(K) * b * skew;
        std::cout << "t: \n" << t << '\n';
    }

    // TODO: uncomment the line below to return true when testing your algorithm and in you final submission.
    return true;
}

















