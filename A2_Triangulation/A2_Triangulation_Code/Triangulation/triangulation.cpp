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

#include <easy3d/viewer/drawable_points.h>
#include <easy3d/viewer/camera.h>
#include <easy3d/viewer/texture.h>
#include <easy3d/viewer/primitives.h>
#include <easy3d/fileio/resources.h>
#include <easy3d/fileio/point_cloud_io.h>

#include <3rd_party/glfw/include/GLFW/glfw3.h>	// for the KEYs


using namespace easy3d;


Triangulation::Triangulation(const std::string &title, 
                             const std::string &image_point_file_0,
                             const std::string &image_point_file_1)
        : Viewer(title)
        , texture_(nullptr)
{
    auto load_image_points = [](const std::string image_point_file) -> std::vector<vec3> {
        PointCloud cloud;
        if (easy3d::load_xyz(image_point_file, cloud))
            return cloud.points();
        return std::vector<vec3>();
    };

    image_0_points_ = load_image_points(image_point_file_0);
    image_1_points_ = load_image_points(image_point_file_1);

    if (image_0_points_.empty())
        LOG(ERROR) << "failed loading image points from file: " << image_point_file_0;
    else
        LOG(ERROR) << image_0_points_.size() << " image points loaded from file 0";

    if (image_1_points_.empty())
        LOG(ERROR) << "failed loading image points from file: " << image_point_file_1;
    else
        LOG(ERROR) << image_1_points_.size() << " image points loaded from file 1";

    if (image_0_points_.size() != image_1_points_.size())
        LOG(ERROR) << "numbers of image points do not match (" << image_0_points_.size() << " != " << image_1_points_.size() << ")";

    if (image_0_points_.size() < 8)
        LOG(ERROR) << "no enough point pairs (" << image_0_points_.size() << " loaded, but at least " << 8 << " needed)";
}


std::string Triangulation::usage() const {
    return ("\n============================== Usage ====================================\n"
            "\tCorresponding images points should have been loaded (unless you see    \n"
            "\t\tsome error messages). If no errors, you can press the 'space' key    \n"
            "\t\tto run your 'triangulation()' method.                                \n"
            "\tYou need to implement 'triangulation()' in 'triangulation_method.cpp'. \n"
            "-------------------------------------------------------------------------\n"
            "\tIf reconstruction was successful, the recovered 3D points will be      \n"
            "\tvisualized for you. You can use your mouse to interact with the viewer:\n"
            "\t\t- Left button: rotate                                                \n"
            "\t\t- Right button: move                                                 \n"
            "\t\t- Wheel: zoom in/out                                                 \n"
            "-------------------------------------------------------------------------\n"
            "\tSome other tools provided by the viewer:                               \n"
            "\t\t- 's': take a snapshot of the visualization                          \n"
            "\t\t- 'Ctrl' + 's': save the reconstructed 3D points into an '.xyz' file \n"
            "\t\t   (required in your final submission)                               \n"
            "-------------------------------------------------------------------------\n");
}


void Triangulation::post_draw() {
    Viewer::post_draw();
    if (!texture_) {
        texture_ = Texture::create(resource::directory() + "/data/images.png");
        if (!texture_)
            return;
    }

    int w = width() * dpi_scaling();
    int h = height() * dpi_scaling();

    int tex_w = texture_->width();
    int tex_h = texture_->height();
    float image_as = tex_w / static_cast<float>(tex_h);
    float viewer_as = width() / static_cast<float>(height());
    if (image_as < viewer_as) {// thin
        tex_h = static_cast<int>(height() * 0.4f);
        tex_w = static_cast<int>(tex_h * image_as);
    }
    else {
        tex_w = static_cast<int>(width() * 0.4f);
        tex_h = static_cast<int>(tex_w / image_as);
    }

    const Rect quad(20 * dpi_scaling(), (20 + tex_w) * dpi_scaling(), 20 * dpi_scaling(), (20 + tex_h) * dpi_scaling());
    opengl::draw_quad_filled(quad, texture_->id(), w, h, -0.9f);
}


void Triangulation::cleanup() {
    if (texture_)
        delete texture_;
    Viewer::cleanup();
}


bool Triangulation::key_press_event(int key, int modifiers) {
    if (key == GLFW_KEY_SPACE) {

        if (image_0_points_.size() != image_1_points_.size()) {
            LOG(ERROR) << "numbers of points do not match (" << image_0_points_.size() << " != "
                       << image_1_points_.size() << ")";
            return EXIT_FAILURE;
        }

        if (image_0_points_.size() < 8) {
            LOG(ERROR) << "no enough point pairs (" << image_0_points_.size() << " loaded, but at least " << 8
                       << " needed)";
            return false;
        }

        const float fx = 1000;  /// TODO: tune to see how it affects the reconstruction
        const float fy = 1000;  /// TODO: tune to see how it affects the reconstruction
        const float cx = 320;   /// TODO: tune to see how it affects the reconstruction
        const float cy = 240;   /// TODO: tune to see how it affects the reconstruction
        std::vector<vec3> points_3d;
        mat3 R;
        vec3 t;
        bool success = triangulation(fx, fy, cx, cy, image_0_points_, image_1_points_, points_3d, R, t);
        if (success) {
            if (points_3d.empty()) {
                LOG(ERROR) << "triangulation() has returned 'true', but no 3D points returned";
                return false;
            }

            PointCloud* cloud = dynamic_cast<PointCloud*>(current_model());
            if (cloud)
                cloud->clear();
            else {
                cloud = new PointCloud;
                cloud->set_name("triangulation.xyz");
                add_model(cloud, false);
            }

            for (std::size_t i=0; i<points_3d.size(); ++i)
                cloud->add_vertex(points_3d[i]);
            std::cout << "reconstructed model has " << cloud->n_vertices() << " points" << std::endl;

            PointsDrawable* drawable = cloud->points_drawable("vertices");
            if (!drawable) {
                drawable = cloud->add_points_drawable("vertices");
                drawable->set_default_color(vec4(0.8f, 0.3f, 0.4f, 1.0f));
                drawable->set_point_size(5.0f);
                drawable->set_impostor_type(PointsDrawable::SPHERE);
            }
            drawable->update_vertex_buffer(points_3d);

            // update view using the recovered R and r
            camera()->set_from_calibration(fx, fy, 0.0, cx, cy, R, t, true);

            fit_screen();

            return true;
        } else {
            LOG(ERROR) << "triangulation failed";
            return false;
        }
    }
    else
        return Viewer::key_press_event(key, modifiers);
}