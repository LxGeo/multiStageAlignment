#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"
#include "stitching/vector_on_raster_stitcher.h"
#include "gdal_algs_wrap/gdal_proximity.h"
#include "numcpp/stats.h"
#include <nlohmann/json.hpp>
#include "template_matcher_alignment.h"

namespace LxGeo
{
	using namespace LxSpatialOps;
	namespace MSA
	{

		class MultiStageAlignment : public BaseOperation {

		public:
			std::pair<double, double> xy_constants;
			std::string dhm_completness_column_name = "DHM_COMP";
			std::string dhm_height_column_name = "DHM_H";
			std::string tm_height_column_name = "TM_H";
			std::string correlation_value_column_name = "CORR";
			std::string fitness_column_name = "FITT";
			std::string success_column_name = "SUCCESS";
			std::string method_column_name = "METHOD";
			bool keep_geometries = true;
			
			
		public:
			MultiStageAlignment(
				const std::unordered_map<std::string, std::string>& view1_paths_map,
				const std::unordered_map<std::string, std::string>& view2_paths_map,
				const std::string& couple_file_path,
				const std::string& aligned_vector_1,
				const std::string& aligned_vector_2,
				const std::string proximity_1,
				const std::string proximity_2,
				bool _keep_geometries = true
			) {

				//TODO add check same crs

				add_raster_input_dataset(view1_paths_map.at("ortho"), "ortho_1");
				add_raster_input_dataset(view1_paths_map.at("proba"), "proba_1");
				add_raster_input_dataset(view1_paths_map.at("dsm"), "dsm_1");
				add_vector_input_dataset(view1_paths_map.at("vector"), "vector_1");

				add_raster_input_dataset(view2_paths_map.at("ortho"), "ortho_2");
				add_raster_input_dataset(view2_paths_map.at("proba"), "proba_2");
				add_raster_input_dataset(view2_paths_map.at("dsm"), "dsm_2");
				add_vector_input_dataset(view2_paths_map.at("vector"), "vector_2");

				add_vector_output_dataset(aligned_vector_1, { "vector_1" }, ext_intersection, WriteMode::overwrite, "al_v1");
				add_vector_output_dataset(aligned_vector_2, { "vector_2" }, ext_intersection, WriteMode::overwrite, "al_v2");

				add_raster_output_dataset(proximity_1, { "proba_1" }, ext_intersection, "proximity_1");
				add_raster_output_dataset(proximity_2, { "proba_2" }, ext_intersection, "proximity_2");

				std::ifstream f(couple_file_path);
				nlohmann::json data = nlohmann::json::parse(f);
				std::vector<double> loaded_constants = data["xy_cst"];
				if (loaded_constants.size() != 2)
					throw std::runtime_error("Error reading xy_constants from couple file");
				xy_constants = std::make_pair(loaded_constants[0], loaded_constants[1]);

				init_cpd(OperationDiveStrategy::zoom, 1000.0, 100.0);

				keep_geometries = _keep_geometries;

			}

			GeoImage<cv::Mat> proba2proximity(const GeoImage<cv::Mat>& proba_gimg, size_t band_idx = 2) {
				// functor used to apply argmax on proba map and return a binary map of the respective band
				const std::function<cv::Mat(cv::Mat)> binarizer_fn = [&band_idx](const cv::Mat& in_img) {
					std::vector<cv::Mat> planes;
					cv::split(in_img, planes);
					cv::Mat b_is_larger = cv::Mat::ones(in_img.size(), CV_8UC1);
					for (int b_idx = 0; b_idx < planes.size(); b_idx++) {
						cv::Mat binaryImage;
						cv::compare(planes[band_idx], planes[b_idx], binaryImage, cv::CMP_GE);
						cv::bitwise_and(b_is_larger, binaryImage, b_is_larger);						
					}
					return b_is_larger;
				};
				VirtualGeoImage<cv::Mat> binary_contour = VirtualGeoImage<cv::Mat>(proba_gimg, binarizer_fn);
				GeoImage<cv::Mat> proximity_gimg = proximity_raster(binary_contour, 1);
				return proximity_gimg;
			}
			
			template <std::ranges::input_range Range>
			void assign_heights(Range& in_g_vector, const GeoImage<cv::Mat>& height_raster) {
				float null_value = (height_raster.no_data.has_value()) ? height_raster.no_data.value() : FLT_MAX; //value_or
				RasterPixelsStitcher rps = RasterPixelsStitcher(height_raster);
				for (auto& gwa : in_g_vector) {
					auto stitched_pixels = rps.readPolygonPixels<float>(gwa.get_definition(), RasterPixelsStitcherStartegy::filled_polygon);
					auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
					if (!stats.empty()) {
						gwa.set_double_attribute(dhm_height_column_name, stats.percentile(90));
						double null_percentage =  double(stats.count_null()) / double(stats.count()+ stats.count_null());
						gwa.set_double_attribute(dhm_completness_column_name, 1.0-null_percentage);

					}
					else {
						gwa.set_double_attribute(dhm_height_column_name, 0);
						gwa.set_double_attribute(dhm_completness_column_name, 0);
					}
				}
			}

			template <std::ranges::input_range Range>
			void assign_fitness(Range& in_g_vector, const GeoImage<cv::Mat>& proximity_raster, const std::pair<double,double>& xy_constants, std::string& height_column_name) {

				auto fitness_from_stats_functor = [this](numcpp::DetailedStats<float>& stats)->float {
					if (stats.empty())
						return 1.0;
					auto sigmoid = [](const double& value, double offset = 0.0, double scale = 1.0, double exp_scale = 1.0) {return scale / (1.0 + std::exp(-exp_scale * value)) + offset; };
					double alpha = 0.3, gamma = 0.3;
					double w_objective = alpha * sigmoid(stats.mean(), -1.0, 2.0, 0.3) + gamma * sigmoid(stats.stdev(), -1.0, 2.0, 0.3);
					return w_objective / (alpha + gamma);
				};

				float null_value = (proximity_raster.no_data.has_value()) ? proximity_raster.no_data.value() : FLT_MAX; //value_or
				RasterPixelsStitcher rps = RasterPixelsStitcher(proximity_raster);
				for (auto& gwa : in_g_vector) {
					double h = -gwa.get_double_attribute(height_column_name);
					auto translate_matrix = bg::strategy::transform::translate_transformer<double, 2, 2>(h * xy_constants.first, h * xy_constants.second);
					auto translated_geometry = translate_geometry(gwa.get_definition(), translate_matrix);

					auto stitched_pixels = rps.readPolygonPixels<float>(translated_geometry, RasterPixelsStitcherStartegy::contours);
					auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
					gwa.set_double_attribute(fitness_column_name, fitness_from_stats_functor(stats));
					
				}
			}

			template <std::ranges::input_range Range>
			void move_geometries(Range& in_g_vector, SIGN direction) {
				for (auto& gwa : in_g_vector) {
					std::string height_column_name = (gwa.get_string_attribute(method_column_name) == "1") ? dhm_height_column_name : tm_height_column_name;
					double h = gwa.get_double_attribute(height_column_name) * direction;
					auto translate_matrix = bg::strategy::transform::translate_transformer<double, 2, 2>(h * xy_constants.first, h * xy_constants.second);
					gwa.set_double_attribute("disp_x", h * xy_constants.first);
					gwa.set_double_attribute("disp_y", h * xy_constants.second);
					auto translated_geometry = translate_geometry(gwa.get_definition(), translate_matrix);
					gwa.set_definition(translated_geometry);
				}
			}

			/*
			Expecting a polygon geovector & respective dsm & target proximity map to use for confidence assignment && xy_constants for heihgt-disparity conversion
			Updates respective geovector by adding height_column_name, fitness_column_name, success_column_name to each geometry
			*/
			template <std::ranges::input_range Range>
			void height_alignment_op(Range& in_geovector, const GeoImage<cv::Mat>& respective_dsm, const GeoImage<cv::Mat>& target_proximity_map,
				const std::pair<double, double>& xy_constants) {
				
				assign_heights(in_geovector, respective_dsm);
				for (auto& gwa : in_geovector ) {
					gwa.set_string_attribute(method_column_name, "1");
				}
				
				assign_fitness(in_geovector, target_proximity_map, xy_constants, dhm_height_column_name);

				auto success_filter_predicate = [this](const Geometries_with_attributes<Boost_Polygon_2>& feature) {
					double success_value = (1.0 - feature.get_double_attribute(fitness_column_name)) * feature.get_double_attribute(dhm_completness_column_name);
					return success_value > 0.5;
				};

				for (auto& gwa : std::views::filter(in_geovector, success_filter_predicate)) {
					gwa.set_int_attribute(success_column_name, 1);
				}

			}

			template <std::ranges::input_range Range>
			void template_matching_alignment_op(Range& in_g_vector, const GeoImage<cv::Mat>& template_ortho, const GeoImage<cv::Mat>& search_ortho,
				const GeoImage<cv::Mat>& target_proximity_map, const std::pair<double, double>& xy_constants) {
				using namespace templateMatchingAlignment;
				double disparity_to_height_scale_factor = 1 * std::abs(template_ortho.geotransform[1]) / numcpp::norm(std::vector<double>({ xy_constants.first, xy_constants.second }));
				double rotation_angle = -DEGS(std::atan2(xy_constants.second, xy_constants.first));
				
				int left_search_pixels = 200, right_search_pixels = 0;

				// Create epipolar images
				auto epi_template_geoimage = rotate(template_ortho, rotation_angle);
				auto epi_search_geoimage = rotate(search_ortho, rotation_angle);
				std::string disp_column_name = "DISP";
				polygonsMatcher::estimate_disp_1d(in_g_vector, epi_template_geoimage, epi_search_geoimage,
					TemplateMatchingMethod::ccoef_N, left_search_pixels, right_search_pixels,
					disp_column_name, correlation_value_column_name
				);
				
				
				for (auto& gwa : in_g_vector) {
					gwa.set_string_attribute(method_column_name, "2");
					gwa.set_double_attribute(tm_height_column_name, gwa.get_double_attribute(disp_column_name)*disparity_to_height_scale_factor);
				}

				auto confident_filter_predicate = [this](const Geometries_with_attributes<Boost_Polygon_2>& feature) { // template matching cofidence
					return feature.get_double_attribute(correlation_value_column_name) > 0.3; // depends on TM_* method
				};

				auto polygon_successful_height_assigned = std::views::filter(in_g_vector, confident_filter_predicate);
				assign_fitness(polygon_successful_height_assigned, target_proximity_map, xy_constants, tm_height_column_name); // alignment confidence requieres height attribute set

				auto good_confidence_filter_predicate = [this](const Geometries_with_attributes<Boost_Polygon_2>& feature) {
					
					double success_value = (1.0 - feature.get_double_attribute(fitness_column_name)) * feature.get_double_attribute(correlation_value_column_name);
					return success_value > 0.5;
				};

				for (auto& gwa : std::views::filter(in_g_vector, good_confidence_filter_predicate)) {
					gwa.set_int_attribute(success_column_name, 1);
				}
			}

			template <std::ranges::input_range Range>
			void preset_geovector(Range& in_g_vector) {
				for (auto& gwa : in_g_vector) {
					gwa.set_string_attribute(method_column_name, "");
					gwa.set_double_attribute(tm_height_column_name, 0.0);
					gwa.set_double_attribute(dhm_height_column_name, 0.0);
					gwa.set_double_attribute(fitness_column_name, 1.0);
					gwa.set_double_attribute(correlation_value_column_name, 0.0);
					gwa.set_double_attribute(dhm_completness_column_name, 0.0);
					gwa.set_int_attribute(success_column_name, 0);
				}
			}

			/**
			The opperation consists of estimating disparity (height) to every polygon using different methods in order:
			1) SGBM disparity 2) Template matching disparity
			Where in each method:
				a) Height is assigned is assigned for every polygon alongside a confidence value of the assignment
				-In the case of SGBM: a low confidence value will be assigned to geometries with unreliable statistcs for example (high variance of values) or (High null value percentage)
				-In the case of template matching a low confidence value will be assigned to geometries with low correlation of patches.
				b) Geometries will be filtered using a the height assignemnt confidence, then another confidence attribute will be assigned based on geometry superpostion after alignment on the probability map.
				-This will be the same for both.
				c) Geometries will be filtered based on the last computed confidence where good confidence geometries will be ignored in the next steps.

			We end up with initial no moved geometries with the following fields:
			METHOD: index of the method used to succefully align a geometry
			HEIGHT: The value of height using the respective method
			CONF: Confidence value of alignment using the respective method

			Finally a translation will be applied using the HEIGHT attribute.
			
			*/
			ViewPair op(ViewPair& in_view_pair) override {
				ViewPair out_view;

				std::future<GeoImage<cv::Mat>> f_contour_proximity_1 = std::async(std::launch::async, [this, &in_view_pair] {return proba2proximity(in_view_pair.raster_views["proba_1"]); });
				std::future<GeoImage<cv::Mat>> f_contour_proximity_2 = std::async(std::launch::async, [this, &in_view_pair] {return proba2proximity(in_view_pair.raster_views["proba_2"]); });

				GeoVector<Boost_Polygon_2>& polygons1 = boost::get<GeoVector<Boost_Polygon_2>>(in_view_pair.vector_views["vector_1"]);
				//auto polygons1 = polygons1g | std::views::filter([](auto& gwa) {return gwa.get_int_attribute("BUILD_ID") == 12; });
				GeoVector<Boost_Polygon_2>& polygons2 = boost::get<GeoVector<Boost_Polygon_2>>(in_view_pair.vector_views["vector_2"]);
				preset_geovector(polygons1); preset_geovector(polygons2);

				GeoImage<cv::Mat> contour_proximity_1 = f_contour_proximity_1.get(), contour_proximity_2 = f_contour_proximity_2.get();
				out_view.raster_views["proximity_1"] = GeoImage<cv::Mat>(contour_proximity_1);
				out_view.raster_views["proximity_2"] = GeoImage<cv::Mat>(contour_proximity_2);
				
				//height_alignment_op(polygons1, in_view_pair.raster_views["dsm_1"], contour_proximity_2, xy_constants);
				//height_alignment_op(polygons2, in_view_pair.raster_views["dsm_2"], contour_proximity_1, { -xy_constants.first, -xy_constants.second });
				
				auto failure_filter_predicate = [this](const Geometries_with_attributes<Boost_Polygon_2>& feature) {
					return feature.get_int_attribute(success_column_name)!=1;
				};

				auto polygon1_remaining = polygons1 | std::views::filter(failure_filter_predicate);
				auto polygon2_remaining = polygons2 | std::views::filter(failure_filter_predicate);

				template_matching_alignment_op(polygon1_remaining, in_view_pair.raster_views["ortho_1"], in_view_pair.raster_views["ortho_2"], contour_proximity_2, xy_constants);
				template_matching_alignment_op(polygon2_remaining, in_view_pair.raster_views["ortho_2"], in_view_pair.raster_views["ortho_1"], contour_proximity_1, {-xy_constants.first, -xy_constants.second});
				
				if (!keep_geometries) {
					std::cout << "Moving second layer geometries!" << std::endl;
					move_geometries(polygons1, SIGN::NEGATIVE);
					std::cout << "Moving first layer geometries!" << std::endl;
					move_geometries(polygons2, SIGN::POSITIVE);
				}

				out_view.valid_geometries_indices["al_v1"].idx = 0;
				out_view.valid_geometries_indices["al_v2"].idx = 0;
				out_view.vector_views.emplace("al_v1" , std::move(polygons1));
				out_view.vector_views.emplace("al_v2" , std::move(polygons2));
				return out_view;
			}


		};
	}
}