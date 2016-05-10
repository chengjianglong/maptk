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
 * \The source is implemented by Chengjiang Long (chengjiang.long@ktiware.com)
 * \file
 * \brief Sparse Bungle Adjustment utility
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <exception>
#include <string>
#include <vector>

#include <vital/vital_foreach.h>
#include <vital/config/config_block.h>
#include <vital/config/config_block_io.h>

#include <vital/algorithm_plugin_manager.h>
#include <vital/exceptions.h>
#include <vital/vital_types.h>
#include <vital/util/get_paths.h>

#include <kwiversys/SystemTools.hxx>
#include <kwiversys/CommandLineArguments.hxx>
#include <maptk/version.h>


#include <vcl_iostream.h>
#include <vcl_vector.h>
#include <vnl/vnl_vector.h>

#include <rrel/rrel_linear_regression.h>
#include <rrel/rrel_orthogonal_regression.h>
#include <rrel/rrel_lms_obj.h>
#include <rrel/rrel_ran_sam_search.h>
#include <rrel/rrel_irls.h>
#include <rrel/rrel_ransac_obj.h>
#include <rrel/rrel_trunc_quad_obj.h>
#include <rrel/rrel_m_est_obj.h>
#include <rrel/rrel_tukey_obj.h>
#include <rrel/rrel_muset_obj.h>


typedef kwiversys::SystemTools     ST;

static kwiver::vital::config_block_sptr default_config()
{

  kwiver::vital::config_block_sptr config = kwiver::vital::config_block::empty_config("bundle_adjust_tracks_tool");

  config->set_value("input_track_file", "",
                    "Path an input file containing feature tracks");

  config->set_value("filtered_track_file", "",
                    "Path to write a file containing filtered feature tracks");

  config->set_value("image_list_file", "",
                    "Path to the input image list file used to generated the "
                    "input tracks.");

  config->set_value("input_pos_files", "",
                    "A directory containing the input POS files, or a text file"
                    "containing a newline-separated list of POS files.\n"
                    "\n"
                    "This is optional, leave blank to ignore.\n"
                    "\n"
                    "This is mutually exclusive with the input_krtd_files "
                    "option for system initialization, and shadowed by the "
                    "input_reference_points option when using an "
                    "st_estimator.");

  config->set_value("input_krtd_files", "",
                    "A directory containing input KRTD camera files, or a text "
                    "file containing a newline-separated list of KRTD files.\n"
                    "\n"
                    "This is optional, leave blank to ignore.\n"
                    "\n"
                    "This is mutually exclusive with input_pos_files option "
                    "for system initialization, and shadowed by the "
                    "input_reference_points_file when using an st_estimator.");

  config->set_value("input_reference_points_file", "",
                    "File containing reference points to use for reprojection "
                    "of results into the geographic coordinate system.\n"
                    "\n"
                    "This option is NOT mutually exclusive with input_*_files "
                    "options when using an st_estimator. When both this and "
                    "another input files option are specified, use of the "
                    "reference file is given priority over the input "
                    "cameras.\n"
                    "\n"
                    "Reference points file format (lm=landmark, tNsM=track N state M):\n"
                    "\tlm1.x lm1.y lm1.z t1s1.frame t1s1.x t1s1.y t1s2.frame t1s2.x t1s2.y ...\n"
                    "\tlm2.x lm2.y lm2.z t2s1.frame t2s1.x t2s1.y t2s2.frame t2s2.x t2s2.y ...\n"
                    "\t...\n"
                    "\n"
                    "At least 3 landmarks must be given, with at least 2 "
                    "track states recorded for each landmark, for "
                    "transformation estimation to converge, however more of "
                    "each is recommended.\n"
                    "\n"
                    "Landmark z position, or altitude, should be provided in meters.");

  config->set_value("initialize_unloaded_cameras", "true",
                    "When loading a subset of cameras, should we optimize only the "
                    "loaded cameras or also initialize and optimize the unspecified cameras");

  config->set_value("output_ply_file", "output/landmarks.ply",
                    "Path to the output PLY file in which to write "
                    "resulting 3D landmark points");

  config->set_value("output_pos_dir", "output/pos",
                    "A directory in which to write the output POS files.");

  config->set_value("output_krtd_dir", "output/krtd",
                    "A directory in which to write the output KRTD files.");

  config->set_value("min_track_length", "50",
                    "Filter the input tracks keeping those covering "
                    "at least this many frames.");

  config->set_value("camera_sample_rate", "1",
                    "Sub-sample the cameras for by this rate.\n"
                    "Set to 1 to use all cameras, "
                    "2 to use every other camera, etc.");

  config->set_value("necker_reverse_input", "false",
                    "Apply a Necker reversal to the initial cameras and landmarks");

  config->set_value("base_camera:focal_length", "1.0",
                    "focal length of the base camera model");

  config->set_value("base_camera:principal_point", "640 480",
                    "The principal point of the base camera model \"x y\".\n"
                    "It is usually safe to assume this is the center of the "
                    "image.");

  config->set_value("base_camera:aspect_ratio", "1.0",
                    "the pixel aspect ratio of the base camera model");

  config->set_value("base_camera:skew", "0.0",
                    "The skew factor of the base camera model.\n"
                    "This is almost always zero in any real camera.");

  config->set_value("ins:rotation_offset", "0 0 0 1",
                    "A quaternion used to offset rotation data from POS files "
                    "when updating cameras. This option is only relevent if a "
                    "value is give to the input_pos_files option.");

  return config;
}


// ------------------------------------------------------------------
static bool check_config(kwiver::vital::config_block_sptr config)
{
  bool config_valid = true;

#define MAPTK_CONFIG_FAIL(msg) \
  std::cerr << "Config Check Fail: " << msg << std::endl; \
  config_valid = false

  if (!config->has_value("input_ply_file"))
  {
    MAPTK_CONFIG_FAIL("Not given a landmarks file path");
  }
  else if (! ST::FileExists( config->get_value<std::string>("input_ply_file"), true ) )
  {
    MAPTK_CONFIG_FAIL("Given landmarks file path doesn't point to an existing file.");
  }

  if (!config->has_value("ransac_scale"))
  {
	  MAPTK_CONFIG_FAIL("Not given a ransac prior scale");
  }

  if (!config->has_value("output_plane_file"))
  {
    MAPTK_CONFIG_FAIL("Not given an file to store the estimated parameter for the fitting plane");
  }


#undef MAPTK_CONFIG_FAIL

  return config_valid;
}



#define print_config(config)                                            \
  do                                                                    \
  {                                                                     \
    VITAL_FOREACH( kwiver::vital::config_block_key_t key, config->available_values() ) \
    {                                                                   \
      std::cerr << "\t"                                                 \
           << key << " = " << config->get_value<kwiver::vital::config_block_key_t>(key) \
           << std::endl;                                                \
    }                                                                   \
  } while (false)


static int fit_plane(int argc, char const* argv[])
{
  static bool        opt_help(false);
  static std::string opt_config;
  static std::string opt_out_config;

  kwiversys::CommandLineArguments arg;

  arg.Initialize( argc, argv );
  typedef kwiversys::CommandLineArguments argT;

  arg.AddArgument( "--help",        argT::NO_ARGUMENT, &opt_help, "Display usage information" );
  arg.AddArgument( "--config",      argT::SPACE_ARGUMENT, &opt_config, "Configuration file for tool" );
  arg.AddArgument( "-c",            argT::SPACE_ARGUMENT, &opt_config, "Configuration file for tool" );
  arg.AddArgument( "--output-config", argT::SPACE_ARGUMENT, &opt_out_config,
                   "Output a configuration. This may be seeded with a configuration file from -c/--config." );
  arg.AddArgument( "-o",            argT::SPACE_ARGUMENT, &opt_out_config,
                   "Output a configuration. This may be seeded with a configuration file from -c/--config." );

    if ( ! arg.Parse() )
  {
    std::cerr << "Problem parsing arguments" << std::endl;
    return EXIT_FAILURE;
  }

  if ( opt_help )
  {
    std::cout
      << "USAGE: " << argv[0] << " [OPTS]\n\n"
      << "Options:"
      << arg.GetHelp() << std::endl;
    return EXIT_SUCCESS;
  }

  // register the algorithm implementations
  std::string rel_plugin_path = kwiver::vital::get_executable_path() + "/../lib/maptk";
  kwiver::vital::algorithm_plugin_manager::instance().add_search_path(rel_plugin_path);
  kwiver::vital::algorithm_plugin_manager::instance().register_plugins();

  // Set config to algo chain
  // Get config from algo chain after set
  // Check config validity, store result
  //
  // If -o/--output-config given,
  //   output config result and notify of current (in)validity
  // Else error if provided config not valid.

  // Set up top level configuration w/ defaults where applicable.
  kwiver::vital::config_block_sptr config = kwiver::vital::config_block::empty_config();

  // If -c/--config given, read in confg file, merge in with default just generated
  if( ! opt_config.empty() )
  {
    const std::string prefix = kwiver::vital::get_executable_path() + "/..";
    config->merge_config(kwiver::vital::read_config_file(opt_config, "maptk",
                                                         MAPTK_VERSION, prefix));
  }

  //std::cerr << "[DEBUG] Config AFTER set:" << std::endl;
  //print_config(config);

  bool valid_config = check_config(config);

  if(!valid_config)
  {
    std::cerr << "ERROR: Configuration not valid." << std::endl;
    return EXIT_FAILURE;
  }
  
  std::vector< vnl_vector<double> > mylandmarks;
  std::string input_ply_file = config->get_value<std::string>("input_ply_file");
  std::cout << input_ply_file << std::endl;
  std::ifstream ifs(input_ply_file.c_str(), std::ifstream::in);
  std::string line;
  if(ifs.is_open())
  {
      int lno = 0;
	  int vno;
	  int startlno = 1000000;
	  while(ifs.good())
	  {
		  lno++;
		  std::getline(ifs, line);
		  if(line.find("element vertex") != std::string::npos)
		  {
			  std::string vnostr = line.substr(15);
			  vno = stoi(vnostr);
			  //std::cout << "Vertex # = " << vno << std::endl;
		  }

		  if(line == "end_header")
		  {
			  startlno = lno + 1;
		  }
		  
		  if(lno >= startlno && lno < vno+startlno)
		  {
              //std::cout << lno << ": " << line << std::endl;
			  std::stringstream ss(line);
			  std::string token;
			  std::vector<std::string> valstrs;
			  while(std::getline(ss, token, ' '))
			  {
				  valstrs.push_back(token);
				  //std::cout << token << std::endl;
			  }

			  vnl_vector<double> vertex(3);
			  for(int i=0; i<3; i++)
			  {
				  vertex[i] = stod(valstrs[i]);
			  }

			  mylandmarks.push_back(vertex);
		  }
	  }
  }

  ifs.close();

  for(int i=0; i<mylandmarks.size(); i++)
  {
      std::cout << "i = " << i << ": ";
	  for(int j=0; j<mylandmarks[i].size(); j++)
	  {
	      std::cout << " " << mylandmarks[i][j];
	  }
	  std::cout <<  std::endl;
  }

  // We wish to fit a line to the (x,y) points, assuming x is the
  // independent variable and y is the dependent variable. This is a
  // linear regression problem. We wish to estimate y = ax + b, but
  // the data is of the form (x,y), not (1,x,y), so we set
  // use_intercept=true.
  //
  //bool use_intercept = true;
  rrel_orthogonal_regression *lr = new rrel_orthogonal_regression(mylandmarks);

  // This controls the verbosity of the search techniques.
  int trace_level = 0;

  // These next three parameters are used in the random sampling
  // searches, not in the IRLS searches.

  // The maximum fraction of the data that is expected to be gross outliers.
  double max_outlier_frac = 0.5;

  // The desired probability of finding the correct fit.
  double desired_prob_good = 0.99;

  // The number of different populations in the data set. For most
  // problems, the data is from one source (surface, etc.), so this
  // will be 1.
  int max_pops = 1;

  // Now we try different objective function/search technique
  // combinations to solve this linear regression problem.

  std::string output_plane_file = config->get_value<std::string>("output_plane_file");
  std::cout << "The parameters of the fitting plane are stored in the file:  " <<  output_plane_file << std::endl;
  std::ofstream ofs(output_plane_file.c_str(), std::ifstream::out);
  
  //
  // Pure least-squares.
  //
  // Most problems implement the weighted_least_squares_fit()
  // function, which means that the IRLS search technique can be used
  // on those problems. It also means that we can do a simple LS by
  // calling the function directly without providing a weight vector.
  //
  {
    vnl_vector<double> ls_params;
    vnl_matrix<double> ls_norm_covar;
    if ( !lr->weighted_least_squares_fit( ls_params, ls_norm_covar ) )
      vcl_cout << "Regression failed!!\n";
    else
      vcl_cout << "Regression succeeded.\n"
               << "estimate = " << ls_params[0] << " * x_0 + " << ls_params[1] << " * x_1 + " 
			   << ls_params[2] << " * x_2 + " << ls_params[3] << vcl_endl;

	ofs << "Estimated by least-squares: " << std::endl;
	for(int i=0; i<ls_params.size(); i++)
	{
		vcl_cout << "ls_params[" << i << "] = " << ls_params[i] << vcl_endl;

		if(i<ls_params.size()-1)
		{
			ofs << ls_params[i] << " ";
		}
		else
		{
			ofs << ls_params[i] << std::endl << std::endl;
		}

	}
    vcl_cout << vcl_endl;
  }

  //
  //  Least Median of Squares
  //
  {
    int num_sam_inst = lr->num_samples_to_instantiate();
    rrel_objective* lms = new rrel_lms_obj( num_sam_inst );
    rrel_ran_sam_search* ransam = new rrel_ran_sam_search;
    ransam->set_sampling_params( max_outlier_frac, desired_prob_good, max_pops);
    ransam->set_trace_level(trace_level);

    if ( !ransam->estimate( lr, lms) )
      vcl_cout << "LMS failed!!\n";
    else
      vcl_cout << "LMS succeeded.\n"
               << "estimate = " << ransam->params()[0] << " * x_0 + " << ransam->params()[1] << " * x_1 + "
			   << ransam->params()[2] << " * x_2 + " << ransam->params()[3] << "\n"
               << "scale = " << ransam->scale() << vcl_endl;


	ofs << "Estimated by least median of squares: " << std::endl;
	for(int i=0; i<ransam->params().size(); i++)
	{
		vcl_cout << "ransam->params()[" << i << "] = " << ransam->params()[i]<< vcl_endl;

		if(i<ransam->params().size()-1)
		{
			ofs << ransam->params()[i] << " ";
		}
		else
		{
			ofs << ransam->params()[i] << std::endl << std::endl;
		}
	}
    vcl_cout << vcl_endl;

	delete ransam;
    delete lms;
  }

  //
  //  RANSAC
  //
  {
    rrel_ransac_obj* ransac = new rrel_ransac_obj( );
    rrel_ran_sam_search* ransam = new rrel_ran_sam_search;
    ransam->set_sampling_params( max_outlier_frac, desired_prob_good, max_pops);
    ransam->set_trace_level(trace_level);
  
	double ransac_scale = config->get_value<double>("ransac_scale");
	std::cout << "RANSAC_PRIOR_scale = " << ransac_scale << std::endl;
    lr->set_prior_scale(ransac_scale);

    if ( !ransam->estimate( lr, ransac) )
      vcl_cout << "RANSAC failed!!\n";
    else
      vcl_cout << "RANSAC succeeded.\n"
               << "estimate = " << ransam->params()[0] << " * x_0 + " << ransam->params()[1] << " * x_1 + " 
			   << ransam->params()[2] << " * x_2 + " << ransam->params()[3] << "\n"
               << "scale = " << ransam->scale() << vcl_endl;

	
	ofs << "Estimated by RANSAC: " << std::endl;
	for(int i=0; i<ransam->params().size(); i++)
	{
		vcl_cout << "ransam->params()[" << i << "] = " << ransam->params()[i]<< vcl_endl;

		if(i<ransam->params().size()-1)
		{
			ofs << ransam->params()[i] << " ";
		}
		else
		{
			ofs << ransam->params()[i] << std::endl;
		}
	}
    vcl_cout << vcl_endl;

    delete ransac;
    delete ransam;
  }
  
  ofs.close();


  return EXIT_SUCCESS;
}


int main(int argc, char const* argv[])
{
  try
  {
    return fit_plane(argc, argv);
  }
  catch (std::exception const& e)
  {
    std::cerr << "Exception caught: " << e.what() << std::endl;

    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cerr << "Unknown exception caught" << std::endl;

    return EXIT_FAILURE;
  }
}
