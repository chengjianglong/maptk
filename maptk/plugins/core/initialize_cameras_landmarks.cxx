/*ckwg +29
 * Copyright 2014-2016 by Kitware, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither name of Kitware, Inc. nor the names of any contributors may be used
 *    to endorse or promote products derived from this software without specific
 *    prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * \file
 * \brief Implementation of core camera and landmark initialization algorithm
 */

#include "initialize_cameras_landmarks.h"

#include <deque>
#include <iterator>

#include <vital/vital_foreach.h>

#include <vital/exceptions.h>
#include <vital/io/eigen_io.h>

#include <vital/algo/estimate_essential_matrix.h>
#include <vital/algo/triangulate_landmarks.h>
#include <vital/algo/bundle_adjust.h>
#include <vital/algo/optimize_cameras.h>

#include <maptk/plugins/core/triangulate_landmarks.h>
#include <maptk/metrics.h>
#include <maptk/match_matrix.h>
#include <maptk/triangulate.h>
#include <maptk/transform.h>


using namespace kwiver::vital;

namespace {
inline bool is_power_of_two(unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}


/// detect bad tracks
std::set<track_id_t>
detect_bad_tracks(const camera_map::map_camera_t& cams,
                  const landmark_map::map_landmark_t& lms,
                  const std::vector<track_sptr>& trks,
                  double error_tol = 5.0)
{
  typedef landmark_map::map_landmark_t lm_map_t;
  std::set<track_id_t> to_remove;
  VITAL_FOREACH(const lm_map_t::value_type& p, lms)
  {
    landmark_map::map_landmark_t lm_single;
    lm_single.insert(p);
    double rmse = kwiver::maptk::reprojection_rmse(cams, lm_single, trks);
    if( rmse > error_tol )
    {
      std::cerr << "remove track "<<p.first<<" with rmse "<< rmse <<std::endl;
      to_remove.insert(p.first);
    }
  }
  return to_remove;
}

/// remove landmarks with IDs in the set
void
remove_landmarks(const std::set<track_id_t>& to_remove,
                 landmark_map::map_landmark_t& lms)
{
  VITAL_FOREACH(const track_id_t& tid, to_remove)
  {
    lms.erase(tid);
  }
}


/// remove landmarks with IDs in the set
void
remove_tracks(const std::set<track_id_t>& to_remove,
              std::vector<track_sptr>& trks)
{
  std::vector<track_sptr> kept_tracks;
  VITAL_FOREACH(const track_sptr& t, trks)
  {
    const track_id_t& tid = t->id();
    if( !to_remove.count(tid) )
    {
      kept_tracks.push_back(t);
    }
  }
  trks.swap(kept_tracks);
}

}


namespace kwiver {
namespace maptk {
namespace core {


/// Private implementation class
class initialize_cameras_landmarks::priv
{
public:
  /// Constructor
  priv()
  : verbose(false),
    init_from_last(false),
    retriangulate_all(false),
    next_frame_max_distance(0),
    base_camera(),
    e_estimator(),
    camera_optimizer(),
    // use the core triangulation as the default, users can change it
    lm_triangulator(new maptk::core::triangulate_landmarks()),
    bundle_adjuster()
  {
  }

  priv(const priv& other)
  : verbose(other.verbose),
    init_from_last(other.init_from_last),
    retriangulate_all(other.retriangulate_all),
    next_frame_max_distance(other.next_frame_max_distance),
    base_camera(other.base_camera),
    e_estimator(!other.e_estimator ? algo::estimate_essential_matrix_sptr()
                                   : other.e_estimator->clone()),
    camera_optimizer(!other.camera_optimizer ? algo::optimize_cameras_sptr()
                                             : other.camera_optimizer->clone()),
    lm_triangulator(!other.lm_triangulator ? algo::triangulate_landmarks_sptr()
                                           : other.lm_triangulator->clone()),
    bundle_adjuster(!other.bundle_adjuster ? algo::bundle_adjust_sptr()
                                           : other.bundle_adjuster->clone())
  {
  }

  /// Construct and initialized camera for \a frame
  camera_sptr init_camera(frame_id_t frame, frame_id_t last_frame,
                          const camera_map::map_camera_t& cams,
                          const std::vector<track_sptr>& trks,
                          const landmark_map::map_landmark_t& lms) const;

  /// Re-triangulate all landmarks for provided tracks
  void retriangulate(landmark_map::map_landmark_t& lms,
                     const camera_map::map_camera_t& cams,
                     const std::vector<track_sptr>& trks,
                     const std::set<landmark_id_t>& new_lm_ids) const;

  /// Estimate the translation scale using a 2d-3d correspondence
  double estimate_t_scale(const vector_3d& KRp,
                          const vector_3d& Kt,
                          const vector_2d& pt2d) const;

  /// Compute a valid left camera from an essential matrix
  /**
   *  There for four valid left camera possibilities for any essential
   *  matrix (assuming the right camera is the identity camera).
   *  This function select the left camera such that a corresponding
   *  pair of points triangulates in front of both cameras.
   */
  vital::simple_camera
  extract_valid_left_camera(const essential_matrix_d& e,
                            const vector_2d& left_pt,
                            const vector_2d& right_pt) const;

  bool verbose;
  bool init_from_last;
  bool retriangulate_all;
  unsigned int next_frame_max_distance;
  vital::simple_camera base_camera;
  vital::algo::estimate_essential_matrix_sptr e_estimator;
  vital::algo::optimize_cameras_sptr camera_optimizer;
  vital::algo::triangulate_landmarks_sptr lm_triangulator;
  vital::algo::bundle_adjust_sptr bundle_adjuster;
};


/// Construct and initialized camera for \a frame
camera_sptr
initialize_cameras_landmarks::priv
::init_camera(frame_id_t frame, frame_id_t last_frame,
              const camera_map::map_camera_t& cams,
              const std::vector<track_sptr>& trks,
              const landmark_map::map_landmark_t& lms) const
{
  typedef landmark_map::map_landmark_t lm_map_t;
  // extract coresponding image points and landmarks
  std::vector<vector_2d> pts_right, pts_left;
  std::vector<landmark_sptr> pts_lm;
  for(unsigned int i=0; i<trks.size(); ++i)
  {
    pts_right.push_back(trks[i]->find(last_frame)->feat->loc());
    pts_left.push_back(trks[i]->find(frame)->feat->loc());
    lm_map_t::const_iterator li = lms.find(trks[i]->id());
    if( li != lms.end() )
    {
      pts_lm.push_back(li->second);
    }
    else
    {
      pts_lm.push_back(landmark_sptr());
    }
  }

  // compute the essential matrix from the corresponding points
  camera_map::map_camera_t::const_iterator ci = cams.find(last_frame);
  if( ci == cams.end() )
  {
    throw invalid_value("Camera for last frame not provided.");
  }
  camera_sptr prev_cam = ci->second;
  camera_intrinsics_sptr cal_right = prev_cam->intrinsics();
  const camera_intrinsics_sptr cal_left = base_camera.get_intrinsics();
  std::vector<bool> inliers;
  essential_matrix_sptr E_sptr = e_estimator->estimate(pts_right, pts_left,
                                                       cal_right, cal_left,
                                                       inliers, 2.0);
  const essential_matrix_d E(*E_sptr);

  unsigned num_inliers = static_cast<unsigned>(std::count(inliers.begin(),
                                                          inliers.end(), true));
  if( this->verbose )
  {
    std::cout << "E matrix num inliers = " << num_inliers
              << "/" << inliers.size() << std::endl;
  }

  // get the first inlier index
  unsigned int inlier_idx = 0;
  for(; inlier_idx < inliers.size() && !inliers[inlier_idx]; ++inlier_idx);

  // get the first inlier correspondence to
  // disambiguate essential matrix solutions
  vector_2d left_pt = cal_left->unmap(pts_left[inlier_idx]);
  vector_2d right_pt = cal_right->unmap(pts_right[inlier_idx]);

  // compute the corresponding camera rotation and translation (up to scale)
  vital::simple_camera cam = extract_valid_left_camera(E, left_pt, right_pt);
  cam.set_intrinsics(cal_left);

  // compute the scale from existing landmark locations (if available)
  matrix_3x3d prev_R(prev_cam->rotation());
  vector_3d prev_t = prev_cam->translation();
  matrix_3x3d R = matrix_3x3d(cam.get_rotation());
  vector_3d t = cam.translation();
  std::vector<double> scales;
  scales.reserve(num_inliers);
  for(unsigned int i=0; i<inliers.size(); ++i)
  {
    if( !inliers[i] || !pts_lm[i] )
    {
      continue;
    }
    vector_3d pt3d = prev_R * pts_lm[i]->loc() + prev_t;
    const vector_2d& pt2d = cal_left->unmap(pts_left[i]);
    scales.push_back(estimate_t_scale(R*pt3d, t, pt2d));
  }
  // find the median scale
  double median_scale = 1.0;
  if( !scales.empty() )
  {
    size_t n = scales.size() / 2;
    std::nth_element(scales.begin(), scales.begin()+n, scales.end());
    median_scale = scales[n];
  }
  if( this->verbose )
  {
    std::cout << " median scale = "<< median_scale<<std::endl;
    if( !scales.empty() )
    {
      std::sort(scales.begin(), scales.end());
      std::cout << "    min scale = " << scales.front() << '\n'
                << "    max scale = " << scales.back() << std::endl;
    }
  }

  // adjust pose relative to the previous camera
  vector_3d new_t = cam.get_rotation() * prev_cam->translation()
                  + median_scale * cam.translation();
  cam.set_rotation(cam.get_rotation() * prev_cam->rotation());
  cam.set_translation(new_t);

  return cam.clone();
}


/// Re-triangulate all landmarks for provided tracks
void
initialize_cameras_landmarks::priv
::retriangulate(landmark_map::map_landmark_t& lms,
                const camera_map::map_camera_t& cams,
                const std::vector<track_sptr>& trks,
                const std::set<landmark_id_t>& new_lm_ids) const
{
  typedef landmark_map::map_landmark_t lm_map_t;
  lm_map_t init_lms;
  VITAL_FOREACH(const track_sptr& t, trks)
  {
    const track_id_t& tid = t->id();
    if( !this->retriangulate_all &&
        new_lm_ids.find(tid) == new_lm_ids.end() )
    {
      continue;
    }
    lm_map_t::const_iterator li = lms.find(tid);
    if( li == lms.end() )
    {
      landmark_sptr lm(new landmark_d(vector_3d(0,0,0)));
      init_lms[static_cast<landmark_id_t>(tid)] = lm;
    }
    else
    {
      init_lms.insert(*li);
    }
  }

  landmark_map_sptr lm_map(new simple_landmark_map(init_lms));
  camera_map_sptr cam_map(new simple_camera_map(cams));
  track_set_sptr tracks(new simple_track_set(trks));
  this->lm_triangulator->triangulate(cam_map, tracks, lm_map);

  // detect and remove landmarks with large triangulation error
  std::set<track_id_t> to_remove = detect_bad_tracks(cams,
                                                     lm_map->landmarks(),
                                                     trks, 5.0);
  VITAL_FOREACH(const lm_map_t::value_type& p, lm_map->landmarks())
  {
    lms[p.first] = p.second;
  }
  remove_landmarks(to_remove, lms);
}


/// Estimate the translation scale using a 2d-3d correspondence
double
initialize_cameras_landmarks::priv
::estimate_t_scale(const vector_3d& KRp,
                   const vector_3d& Kt,
                   const vector_2d& pt2d) const
{
  vector_3d a = KRp;
  vector_3d b = Kt;
  a.x() = pt2d.x() * a.z() - a.x();
  b.x() = pt2d.x() * b.z() - b.x();
  a.y() = pt2d.y() * a.z() - a.y();
  b.y() = pt2d.y() * b.z() - b.y();
  double cx = a.x()*b.z() - a.z()*b.x();
  double cy = a.y()*b.z() - a.z()*b.y();
  return (a.x()*cx + a.y()*cy) / -(b.x()*cx + b.y()*cy);
}

/// Compute a valid left camera from an essential matrix
vital::simple_camera
initialize_cameras_landmarks::priv
::extract_valid_left_camera(const essential_matrix_d& e,
                            const vector_2d& left_pt,
                            const vector_2d& right_pt) const
{
  /// construct an identity right camera
  const vector_3d t = e.translation();
  rotation_d R = e.rotation();

  std::vector<vector_2d> pts;
  pts.push_back(right_pt);
  pts.push_back(left_pt);

  std::vector<vital::simple_camera> cams(2);
  const vital::simple_camera& left_camera = cams[1];

  // option 1
  cams[1] = vital::simple_camera(R.inverse()*-t, R);
  vector_3d pt3 = triangulate_inhomog(cams, pts);
  if( pt3.z() > 0.0 && left_camera.depth(pt3) > 0.0 )
  {
    return left_camera;
  }

  // option 2, with negated translation
  cams[1] = vital::simple_camera(R.inverse()*t, R);
  pt3 = triangulate_inhomog(cams, pts);
  if( pt3.z() > 0.0 && left_camera.depth(pt3) > 0.0 )
  {
    return left_camera;
  }

  // option 3, with the twisted pair rotation
  R = e.twisted_rotation();
  cams[1] = vital::simple_camera(R.inverse()*-t, R);
  pt3 = triangulate_inhomog(cams, pts);
  if( pt3.z() > 0.0 && left_camera.depth(pt3) > 0.0 )
  {
    return left_camera;
  }

  // option 4, with negated translation
  cams[1] = vital::simple_camera(R.inverse()*t, R);
  pt3 = triangulate_inhomog(cams, pts);
  if( pt3.z() > 0.0 && left_camera.depth(pt3) > 0.0 )
  {
    return left_camera;
  }
  // should never get here
  return vital::simple_camera();
}


/// Constructor
initialize_cameras_landmarks
::initialize_cameras_landmarks()
: d_(new priv)
{
}


/// Copy Constructor
initialize_cameras_landmarks
::initialize_cameras_landmarks(const initialize_cameras_landmarks& other)
: d_(new priv(*other.d_))
{
}


/// Destructor
initialize_cameras_landmarks
::~initialize_cameras_landmarks()
{
}


/// Get this algorithm's \link vital::config_block configuration block \endlink
vital::config_block_sptr
initialize_cameras_landmarks
::get_configuration() const
{
  // get base config from base class
  vital::config_block_sptr config =
      vital::algo::initialize_cameras_landmarks::get_configuration();

  const camera_intrinsics_sptr K = d_->base_camera.get_intrinsics();

  config->set_value("verbose", d_->verbose,
                    "If true, write status messages to the terminal showing "
                    "debugging information");

  config->set_value("init_from_last", d_->init_from_last,
                    "If true, and a camera optimizer is specified, initialize "
                    "the camera using the closest exiting camera and optimize");

  config->set_value("retriangulate_all", d_->retriangulate_all,
                    "If true, re-triangulate all landmarks observed by a newly "
                    "initialized camera.  Otherwise, only triangulate or "
                    "re-triangulate landmarks that are marked for initialization.");

  config->set_value("next_frame_max_distance", d_->next_frame_max_distance,
                    "Limit the selection of the next frame to initialize to "
                    "within this many frames of an already initialized frame. "
                    "If no valid frames are found, double the search range "
                    "until a valid frame is found. "
                    "A value of zero disables this limit");

  config->set_value("base_camera:focal_length", K->focal_length(),
                    "focal length of the base camera model");

  config->set_value("base_camera:principal_point", K->principal_point().transpose(),
                    "The principal point of the base camera model \"x y\".\n"
                    "It is usually safe to assume this is the center of the "
                    "image.");

  config->set_value("base_camera:aspect_ratio", K->aspect_ratio(),
                    "the pixel aspect ratio of the base camera model");

  config->set_value("base_camera:skew", K->skew(),
                    "The skew factor of the base camera model.\n"
                    "This is almost always zero in any real camera.");

  // nested algorithm configurations
  vital::algo::estimate_essential_matrix
      ::get_nested_algo_configuration("essential_mat_estimator",
                                      config, d_->e_estimator);
  vital::algo::optimize_cameras
      ::get_nested_algo_configuration("camera_optimizer",
                                      config, d_->camera_optimizer);
  vital::algo::triangulate_landmarks
      ::get_nested_algo_configuration("lm_triangulator",
                                      config, d_->lm_triangulator);
  vital::algo::bundle_adjust
      ::get_nested_algo_configuration("bundle_adjuster",
                                      config, d_->bundle_adjuster);
  return config;
}


/// Set this algorithm's properties via a config block
void
initialize_cameras_landmarks
::set_configuration(vital::config_block_sptr config)
{
  const camera_intrinsics_sptr K = d_->base_camera.get_intrinsics();

  // Set nested algorithm configurations
  vital::algo::estimate_essential_matrix
      ::set_nested_algo_configuration("essential_mat_estimator",
                                      config, d_->e_estimator);
  vital::algo::optimize_cameras
      ::set_nested_algo_configuration("camera_optimizer",
                                      config, d_->camera_optimizer);
  vital::algo::triangulate_landmarks
      ::set_nested_algo_configuration("lm_triangulator",
                                      config, d_->lm_triangulator);
  vital::algo::bundle_adjust
      ::set_nested_algo_configuration("bundle_adjuster",
                                      config, d_->bundle_adjuster);

  d_->verbose = config->get_value<bool>("verbose",
                                        d_->verbose);

  d_->init_from_last = config->get_value<bool>("init_from_last",
                                               d_->init_from_last);

  d_->retriangulate_all = config->get_value<bool>("retriangulate_all",
                                                  d_->retriangulate_all);

  d_->next_frame_max_distance =
      config->get_value<unsigned int>("next_frame_max_distance",
                                      d_->next_frame_max_distance);

  vital::config_block_sptr bc = config->subblock("base_camera");
  simple_camera_intrinsics K2(bc->get_value<double>("focal_length",
                                                    K->focal_length()),
                              bc->get_value<vector_2d>("principal_point",
                                                       K->principal_point()),
                              bc->get_value<double>("aspect_ratio",
                                                    K->aspect_ratio()),
                              bc->get_value<double>("skew",
                                                    K->skew()));
  d_->base_camera.set_intrinsics(K2.clone());
}


/// Check that the algorithm's currently configuration is valid
bool
initialize_cameras_landmarks
::check_configuration(vital::config_block_sptr config) const
{
  if (config->get_value<std::string>("camera_optimizer", "") != ""
      && !vital::algo::optimize_cameras
              ::check_nested_algo_configuration("camera_optimizer", config))
  {
    return false;
  }
  if (config->get_value<std::string>("bundle_adjuster", "") != ""
      && !vital::algo::bundle_adjust
              ::check_nested_algo_configuration("bundle_adjuster", config))
  {
    return false;
  }
  return vital::algo::estimate_essential_matrix
             ::check_nested_algo_configuration("essential_mat_estimator",
                                               config)
      && vital::algo::triangulate_landmarks
             ::check_nested_algo_configuration("lm_triangulator", config);
}

namespace
{

/// Extract valid cameras and cameras to initialize.
/**
 * If \a cameras is NULL then return empty cam_map and frame_ids unchanged.
 * If not NULL, return frame_ids containing IDs of all NULL cameras and
 * return cam_map containing all valid cameras.
 * \param [in]     cameras the camera map object to extract from
 * \param [in,out] frame_ids the set of all frames (input),
 *                           the set of frame to initialize (output)
 * \param [out]    cam_map the extract map of valid camera
 */
void extract_cameras(const camera_map_sptr& cameras,
                     std::set<frame_id_t>& frame_ids,
                     camera_map::map_camera_t& cam_map)
{
  cam_map.clear();
  if( !cameras )
  {
    return;
  }

  typedef camera_map::map_camera_t map_cam_t;
  map_cam_t all_cams = cameras->cameras();

  // Find the set of all cameras that need to be initialized
  std::set<frame_id_t> new_frames;
  VITAL_FOREACH(const map_cam_t::value_type& p, all_cams)
  {
    if(p.second)
    {
      cam_map.insert(p);
    }
    else if( frame_ids.count(p.first) )
    {
      new_frames.insert(p.first);
    }
  }
  frame_ids = new_frames;
}


/// Extract valid landmarks and landmarks to initialize.
/**
 * If \a landmarks is NULL then return empty lm_map and track_ids unchanged.
 * If not NULL, return track_ids containing IDs of all NULL landmark and
 * return lm_map containing all valid landmarks.
 * \param [in]     landmarks the landmark map object to extract from
 * \param [in,out] track_ids the set of all tracks (input),
 *                           the set of landmarks to initialize (output)
 * \param [out]    lm_map the extract map of valid landmarks
 */
void extract_landmarks(const landmark_map_sptr& landmarks,
                       std::set<track_id_t>& track_ids,
                       landmark_map::map_landmark_t& lm_map)
{
  lm_map.clear();
  if( !landmarks )
  {
    return;
  }

  typedef landmark_map::map_landmark_t map_landmark_t;
  map_landmark_t all_lms = landmarks->landmarks();

  // Find the set of all landmarks that need to be initialized
  std::set<track_id_t> new_landmarks;
  VITAL_FOREACH(const map_landmark_t::value_type& p, all_lms)
  {
    if(p.second)
    {
      lm_map.insert(p);
    }
    else if( track_ids.count(p.first) )
    {
      new_landmarks.insert(p.first);
    }
  }
  track_ids = new_landmarks;
}


/// Find the closest frame number with an existing camera
frame_id_t
find_closest_camera(const frame_id_t& frame,
                    const camera_map::map_camera_t& cams)
{
  typedef vital::camera_map::map_camera_t map_cam_t;
  frame_id_t other_frame = cams.rbegin()->first;
  map_cam_t::const_iterator ci = cams.lower_bound(frame);
  if( ci == cams.end() )
  {
    other_frame = cams.rbegin()->first;
  }
  else
  {
    other_frame = ci->first;
    if (ci != cams.begin() &&
        (other_frame-frame) >= (frame-(--ci)->first))
    {
      other_frame = ci->first;
    }
  }
  return other_frame;
}


/// Find the subset of new_frames within dist frames of a camera in cams
std::set<frame_id_t>
find_nearby_new_frames(const std::set<frame_id_t>& new_frames,
                       const camera_map::map_camera_t& cams,
                       unsigned int dist)
{
  std::set<frame_id_t> nearby;
  VITAL_FOREACH(camera_map::map_camera_t::value_type p, cams)
  {
    frame_id_t f = p.first < dist ? 0 : p.first - dist;
    for(; f < p.first + dist; ++f)
    {
      nearby.insert(f);
    }
  }
  std::set<frame_id_t> new_nearby;
  std::set_intersection(nearby.begin(), nearby.end(),
                        new_frames.begin(), new_frames.end(),
                        std::inserter(new_nearby, new_nearby.begin()));
  return new_nearby;
}


/// find the best pair of camera indices to start with
void
find_best_initial_pair(const Eigen::SparseMatrix<unsigned int>& mm,
                       int& i, int& j)
{
  typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> vectorXu;
  const int cols = mm.cols();
  const vectorXu d = mm.diagonal();

  // compute the maximum off-diagonal value
  unsigned int global_max_matches = 0;
  for( int k=0; k<cols; ++k)
  {
    for(Eigen::SparseMatrix<unsigned int>::InnerIterator it(mm, k); it; ++it)
    {
      if(it.row() > k && it.value() > global_max_matches)
      {
        global_max_matches = it.value();
      }
    }
  }
  const unsigned int threshold = std::max(global_max_matches / 2, 20u);

  std::cout <<"global max "<<global_max_matches << std::endl;
  std::cout <<"threshold "<<threshold << std::endl;
  for(int x=1; x<cols; ++x)
  {
    unsigned int max_matches = 0;
    int max_i = 0, max_j = 0;
    for(int y=0; y<cols-x; ++y)
    {
      unsigned int matches = mm.coeff(x+y,y);
      if( matches > max_matches )
      {
        max_matches = matches;
        max_i = y;
        max_j = x+y;
      }
    }
    std::cout << "max matches at "<<x<<" is "<< max_matches << " at "<< max_i << ", "<< max_j<<std::endl;
    if( max_matches < threshold )
    {
      break;
    }
    i = max_i;
    j = max_j;
  }
}


// find the frame in the set \p new_frame_ids that sees the most
// landmarks in \p lms in the track set \p tracks.
frame_id_t
next_best_frame(const track_set_sptr tracks,
                const vital::landmark_map::map_landmark_t& lms,
                const std::set<frame_id_t>& new_frame_ids)
{
  const std::vector<track_sptr> trks = tracks->tracks();
  typedef std::map<frame_id_t, unsigned int> frame_map_t;
  frame_map_t vis_count;
  VITAL_FOREACH(const track_sptr& t, trks)
  {
    if( lms.find(t->id()) != lms.end() )
    {
      const std::set<frame_id_t> t_frames = t->all_frame_ids();
      VITAL_FOREACH(const frame_id_t& fid, t_frames)
      {
        if( new_frame_ids.find(fid) == new_frame_ids.end() )
        {
          continue;
        }
        frame_map_t::iterator fmi = vis_count.find(fid);
        if( fmi == vis_count.end() )
        {
          vis_count.insert(std::pair<frame_id_t, unsigned int>(fid,1));
        }
        else
        {
          ++fmi->second;
        }
      }
    }
  }

  // check if remaining new frames see no existing landmarks
  if( vis_count.empty() )
  {
    std::cout << "remaining frames do not see any existing landmarks" << std::endl;
    return *new_frame_ids.begin();
  }

  // find the maximum observation
  unsigned int max_count = 0;
  frame_id_t best_frame = 0;
  VITAL_FOREACH(const frame_map_t::value_type& p, vis_count)
  {
    if(p.second > max_count)
    {
      max_count = p.second;
      best_frame = p.first;
    }
  }
  std::cout << "frame "<< best_frame << " sees "<< max_count << " landmarks"<<std::endl;
  return best_frame;
}


} // end anonymous namespace


/// Initialize the camera and landmark parameters given a set of tracks
void
initialize_cameras_landmarks
::initialize(camera_map_sptr& cameras,
             landmark_map_sptr& landmarks,
             track_set_sptr tracks) const
{
  if( !tracks )
  {
    throw invalid_value("Some required input data is NULL.");
  }
  if( !d_->e_estimator )
  {
    throw invalid_value("Essential matrix estimator not initialized.");
  }
  if( !d_->lm_triangulator )
  {
    throw invalid_value("Landmark triangulator not initialized.");
  }
  typedef vital::camera_map::map_camera_t map_cam_t;
  typedef vital::landmark_map::map_landmark_t map_landmark_t;

  // Extract the existing cameras and camera ids to be initialized
  std::set<frame_id_t> frame_ids = tracks->all_frame_ids();
  map_cam_t cams;
  extract_cameras(cameras, frame_ids, cams);
  std::set<frame_id_t> new_frame_ids(frame_ids);

  // Extract the existing landmarks and landmark ids to be initialized
  std::set<track_id_t> track_ids = tracks->all_track_ids();
  map_landmark_t lms;
  extract_landmarks(landmarks, track_ids, lms);
  std::set<landmark_id_t> new_lm_ids(track_ids.begin(), track_ids.end());

  std::vector<track_sptr> trks = tracks->tracks();

  if(new_frame_ids.empty() && new_lm_ids.empty())
  {
    //nothing to initialize
    return;
  }

  // initialize landmarks if there are already at least two cameras
  if(cams.size() > 2 && !new_lm_ids.empty())
  {
    map_landmark_t init_lms;
    VITAL_FOREACH(const landmark_id_t& lmid, new_lm_ids)
    {
      landmark_sptr lm(new landmark_d(vector_3d(0,0,0)));
      init_lms[static_cast<landmark_id_t>(lmid)] = lm;
    }

    landmark_map_sptr lm_map(new simple_landmark_map(init_lms));
    camera_map_sptr cam_map(new simple_camera_map(cams));
    d_->lm_triangulator->triangulate(cam_map, tracks, lm_map);

    VITAL_FOREACH(const map_landmark_t::value_type& p, lm_map->landmarks())
    {
      lms[p.first] = p.second;
    }
  }

  std::vector<frame_id_t> mm_frames(frame_ids.begin(), frame_ids.end());
  Eigen::SparseMatrix<unsigned int> mm = match_matrix(tracks, mm_frames);
  int init_i=0,init_j=0;
  find_best_initial_pair(mm, init_i, init_j);
  std::cout << "using frames "<< mm_frames[init_i] << " and " << mm_frames[init_j] <<std::endl;

  if(cams.empty())
  {
    // first frame, initialize to base camera
    frame_id_t f = mm_frames[init_i];
    new_frame_ids.erase(f);
    cams[f] = d_->base_camera.clone();
  }

  while( !new_frame_ids.empty() )
  {
    frame_id_t f;
    if( cams.size() == 1 )
    {
      f = mm_frames[init_j];
    }
    else
    {
      unsigned int search_range = d_->next_frame_max_distance;
      if( search_range < 1 )
      {
        f = next_best_frame(tracks, lms, new_frame_ids);
      }
      else
      {
        std::set<frame_id_t> nearby;
        const frame_id_t max_frame = tracks->last_frame();
        while( nearby.empty() && search_range < max_frame )
        {
          nearby = find_nearby_new_frames(new_frame_ids, cams, search_range);
          search_range *= 2;
        }
        f = next_best_frame(tracks, lms, nearby);
      }
    }
    new_frame_ids.erase(f);

    // get the closest frame number with an existing camera
    frame_id_t other_frame = find_closest_camera(f, cams);
    if(d_->verbose)
    {
      std::cout <<"frame "<< f <<" uses reference "<< other_frame <<std::endl;
    }

    // get the subset of tracks that have features on frame f
    track_set_sptr ftracks = tracks->active_tracks(static_cast<int>(f));
    // get the subset of tracks that also  have features on the other frame
    ftracks = ftracks->active_tracks(static_cast<int>(other_frame));

    // find existing landmarks for these tracks
    map_landmark_t flms;
    std::vector<track_sptr> trks = ftracks->tracks();
    VITAL_FOREACH(const track_sptr& t, trks)
    {
      map_landmark_t::const_iterator li = lms.find(t->id());
      if( li != lms.end() )
      {
        flms.insert(*li);
      }
    }

    if( d_->init_from_last && d_->camera_optimizer && flms.size() > 3)
    {
      cams[f] = cams[other_frame]->clone();
    }
    else if( trks.size() > 10 )
    {
      cams[f] = d_->init_camera(f, other_frame, cams, trks, flms);
    }
    else
    {
      break;
    }

    // optionally optimize the new camera
    if( d_->camera_optimizer && flms.size() > 3)
    {
      camera_map::map_camera_t opt_cam_map;
      opt_cam_map[f] = cams[f];
      camera_map_sptr opt_cams(new simple_camera_map(opt_cam_map));
      landmark_map_sptr landmarks(new simple_landmark_map(flms));
      track_set_sptr tracks(new simple_track_set(trks));
      d_->camera_optimizer->optimize(opt_cams, tracks, landmarks);
      cams[f] = opt_cams->cameras()[f];
    }


    // triangulate (or re-triangulate) points seen by the new camera
    d_->retriangulate(lms, cams, trks, new_lm_ids);

    if(d_->verbose)
    {
      camera_map::map_camera_t new_cam_map;
      new_cam_map[f] = cams[f];
      std::vector<double> rpe = maptk::reprojection_errors(new_cam_map, lms, trks);
      if( rpe.empty() )
      {
        std::cerr << "no landmark projections for new camera" << std::endl;
      }
      else
      {
        std::sort(rpe.begin(), rpe.end());
        std::cerr << "new camera reprojections - median: "<<rpe[rpe.size()/2]
                  << " max: " << rpe.back() << std::endl;
      }
    }

    if( d_->bundle_adjuster && cams.size() >= 4 &&
        is_power_of_two(static_cast<unsigned int>(cams.size())) )
    {
      camera_map_sptr ba_cams(new simple_camera_map(cams));
      landmark_map_sptr ba_lms(new simple_landmark_map(lms));
      double init_rmse = maptk::reprojection_rmse(cams, lms, trks);
      std::cerr << "initial reprojection RMSE: " << init_rmse << std::endl;

      d_->bundle_adjuster->optimize(ba_cams, ba_lms, tracks);
      cams = ba_cams->cameras();
      lms = ba_lms->landmarks();

      // detect tracks/landmarks with large error and remove them
      std::set<track_id_t> to_remove = detect_bad_tracks(cams, lms, trks, 5.0);
      remove_landmarks(to_remove, lms);
      std::vector<track_sptr> all_trks = tracks->tracks();
      remove_tracks(to_remove, all_trks);
      tracks = track_set_sptr(new simple_track_set(all_trks));
      double final_rmse = maptk::reprojection_rmse(cams, lms, trks);
      std::cerr << "final reprojection RMSE: " << final_rmse << std::endl;
      std::cout << "updated focal length "<<cams.begin()->second->intrinsics()->focal_length() <<std::endl;
    }

    if(d_->verbose)
    {
      camera_map_sptr ba_cams(new simple_camera_map(cams));
      landmark_map_sptr ba_lms(new simple_landmark_map(lms));
      double curr_rmse = maptk::reprojection_rmse(cams, lms, trks);
      std::cerr << "current reprojection RMSE: " << curr_rmse << std::endl;

      std::cout << "frame "<<f<<" - num landmarks = "<< lms.size() << std::endl;
    }
  }

  // try depth reversal at the end
  if( d_->bundle_adjuster )
  {
    camera_map_sptr ba_cams(new simple_camera_map(cams));
    landmark_map_sptr ba_lms(new simple_landmark_map(lms));
    double init_rmse = maptk::reprojection_rmse(cams, lms, trks);
    std::cerr << "initial reprojection RMSE: " << init_rmse << std::endl;

    d_->bundle_adjuster->optimize(ba_cams, ba_lms, tracks);
    map_cam_t cams1 = ba_cams->cameras();
    map_landmark_t lms1 = ba_lms->landmarks();
    double final_rmse1 = maptk::reprojection_rmse(cams1, lms1, trks);
    std::cerr << "final reprojection RMSE: " << final_rmse1 << std::endl;

    // reverse cameras and optimize again
    camera_map_sptr ba_cams2(new simple_camera_map(cams1));
    landmark_map_sptr ba_lms2(new simple_landmark_map(lms1));
    necker_reverse(ba_cams2, ba_lms2);
    d_->lm_triangulator->triangulate(ba_cams2, tracks, ba_lms2);
    init_rmse = maptk::reprojection_rmse(ba_cams2->cameras(), ba_lms2->landmarks(), trks);
    std::cerr << "flipped initial reprojection RMSE: " << init_rmse << std::endl;
    d_->bundle_adjuster->optimize(ba_cams2, ba_lms2, tracks);
    map_cam_t cams2 = ba_cams2->cameras();
    map_landmark_t lms2 = ba_lms2->landmarks();
    double final_rmse2 = maptk::reprojection_rmse(cams2, lms2, trks);
    std::cerr << "flipped final reprojection RMSE: " << final_rmse2 << std::endl;

    if(final_rmse1 < final_rmse2)
    {
      cams = ba_cams->cameras();
      lms = ba_lms->landmarks();
    }
    else
    {
      cams = ba_cams2->cameras();
      lms = ba_lms2->landmarks();
    }

    // if using bundle adjustment, remove landmarks with large error
    // after optimization
    std::set<track_id_t> to_remove = detect_bad_tracks(cams, lms, trks, 1.0);
    remove_landmarks(to_remove, lms);
  }

  cameras = camera_map_sptr(new simple_camera_map(cams));
  landmarks = landmark_map_sptr(new simple_landmark_map(lms));
}


} // end namespace core
} // end namespace maptk
} // end namespace kwiver
