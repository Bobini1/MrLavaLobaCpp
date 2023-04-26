//
// Created by bobini on 26.04.23.
//

#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <iostream>
#include <fstream>
#include "Params.hpp"
auto getDesc() {
    auto desc = boost::program_options::options_description("Generic options");
    desc.add_options()
      ("help", "produce help message")
        ("version", "print version string")
          ("input_file", boost::program_options::value<std::vector<std::string>>(), "set input file");

    auto inputDataProperties = boost::program_options::options_description("Input data options");
    inputDataProperties.add_options()
      ("run_name", boost::program_options::value<std::string>(), "name of the run (used to save the parameters and the output)")
        ("source", boost::program_options::value<std::string>(), "File name of ASCII digital elevation model")
          ("vent_flag", boost::program_options::value<unsigned int>(), "This flag select how multiple initial coordinates are treated:\n"
                                                                       "vent_flag = 0  => the initial lobes are on the vents coordinates\n"
                                                                       "and the flows start initially from the first vent, "
                                                                       "then from the second and so on.\n"
                                                                       "vent_flag = 1  => the initial lobes are chosen randomly from the vents "
                                                                       "coordinates and each vent has the same probability\n"
                                                                       "vent_flag = 2  => the initial lobes are on the polyline connecting "
                                                                       "the vents and all the point of the polyline "
                                                                       "have the same probability\n"
                                                                       "vent_flag = 3  => the initial lobes are on the polyline connecting "
                                                                       "the vents and all the segments of the polyline "
                                                                       "have the same probability\n"
                                                                       "vent_flag = 4  => the initial lobes are on multiple "
                                                                       "fissures and all the point of the fissures "
                                                                       "have the same probability\n"
                                                                       "vent_flag = 5  => the initial lobes are on multiple "
                                                                       "fissures and all the fissures "
                                                                       "have the same probability\n"
                                                                       "vent_flag = 6  => the initial lobes are on the polyline connecting "
                                                                       "the vents and the probability of "
                                                                       "each segment is fixed by \"fissure probabilities\"\n"
                                                                       "vent_flag = 7  => the initial lobes are on multiple "
                                                                       "fissures and the probability of "
                                                                       "each fissure is fixed by \"fissure_probabilities\"\n"
                                                                       "vent_flag = 8  => the initial lobes are chosen randomly from the vents "
                                                                       "coordinates and the probability of each vent ")
            ("Fissure probabilities", boost::program_options::value<std::vector<int>>(), "this values defines the probabilities of the different segments of the polyline "
                                                                                         "or of the different fissures.")
              ("hazard_flag", boost::program_options::value<bool>(), "Whether a raster map is saved where the values "
                                                                     "represent the probability of a cell to be covered. ")
                ("masking_treshold", boost::program_options::value<double>(), "Fraction of the volume emplaced or the area invaded (according to the flag "
                                                                              "flag_threshold) used to save the run_name_*_masked.asc files."
                                                                              "In this way we cut the \"tails\" with a low thickness (*_thickness_masked.asc) "
                                                                              "and a low probability (*_hazard_masked.asc). The file is used also "
                                                                              "to check the convergence of the solutions increasing the number of flows.\n"
                                                                              "The full outputs are saved in the files run_name_*_full.asc\n"
                                                                              "The masked files are saved only when masking_thresold < 1.")
                  ("m_flows", boost::program_options::value<int>(), "Number of flows")
                    ("min_n_lobes", boost::program_options::value<int>(), "Minimum number of lobes generated for each flow")
                      ("max_n_lobes", boost::program_options::value<int>(), "Maximum number of lobes generated for each flow")
                        ("volume_flag", boost::program_options::value<bool>(), "total volume is read in input, and the "
                                                                               "thickness or the area of the lobes are evaluated according to the "
                                                                               "flag fixed_dimension_flag and the relationship V = n*area*T. "
                                                                               "Otherwise the thickness and the area are read in input and the total "
                                                                               "volume is evaluated (V = n*area*T).")
                          ("total_volume", boost::program_options::value<double>(), "Total volume (this value is used only when volume_flag = true) set \"true\" to be confirmed.")
                            ("fixed_dimension_flag", boost::program_options::value<bool>(), "This flag select which dimension of the lobe is fixed when volume_flag=true. "
                                                                                            "If true, the thickness of the lobes is assigned. "
                                                                                            "Otherwise, the area of the lobes is assigned")
                              ("lobe_area", boost::program_options::value<double>(), "Area of each lobe ( only effective when volume_flag = false or fixed_dimension_flag = true )")
                                ("lobe_area", boost::program_options::value<double>(), "Thickness of each lobe ( only effective when volume_flag = 0 or fixed_dimension_flag  2 )")
                                  ("thickness_ratio", boost::program_options::value<double>(), "Ratio between the thickness of the first lobe of the flow and the thickness of the "
                                                                                               "last lobe.\n"
                                                                                               "thickness_ratio < 1   => the thickness increases with lobe \"age\"\n"
                                                                                               "thickness_ratio = 1   => all the lobes have the same thickness\n"
                                                                                               "thickness_ratio > 1   => the thickness decreases with lobe \"age\"")
                                    ("thickening_parameter", boost::program_options::value<double>(), "This parameter (between 0 and 1) allows for a thickening of the flow giving "
                                                                                                      "controlling the modification of the slope due to the presence of the flow.\n"
                                                                                                      "thickening_parameter = 0  => minimum thickening (maximum spreading)\n"
                                                                                                      "thickening_parameter = 1  => maximum thickening produced in the output\n"
                                                                                                      "default thickening_parameter = 0.2\n"
                                                                                                      "if you reduce this, the impact of the lava flow is lessened in the computation of the slope, "
                                                                                                      "but the thickness is still correct. this allows for \"channel\" flow, if = 1, "
                                                                                                      "then sublava flow would not happen. ")
                                      ("lobe_exponent", boost::program_options::value<double>(), "This parameter (between 0 and 1) allows for a thickening of the flow giving "
                                                                                                 "controlling the modification of the slope due to the presence of the flow.\n"
                                                                                                 "lobe_exponent = 0  => minimum thickening (maximum spreading)\n"
                                                                                                 "lobe_exponent = 1  => maximum thickening produced in the output\n"
                                                                                                 "default lobe_exponent = 0.2\n"
                                                                                                 "if you reduce this, the impact of the lava flow is lessened in the computation of the slope, "
                                                                                                 "but the thickness is still correct. this allows for \"channel\" flow, if = 1, "
                                                                                                 "then sublava flow would not happen. ")
                                        ("lobe_exponent", boost::program_options::value<double>(), "This parameter (between 0 and 1) allows for a thickening of the flow giving "
                                                                                                   "controlling the modification of the slope due to the presence of the flow.\n"
                                                                                                   "lobe_exponent = 0  => minimum thickening (maximum spreading)\n"
                                                                                                   "lobe_exponent = 1  => maximum thickening produced in the output\n"
                                                                                                   "default lobe_exponent = 0.2\n"
                                                                                                   "if you reduce this, the impact of the lava flow is lessened in the computation of the slope, "
                                                                                                   "but the thickness is still correct. this allows for \"channel\" flow, if = 1, "
                                                                                                   "then sublava flow would not happen.")
                                          ("lobe_exponent", boost::program_options::value<double>(), "Lobe_exponent is associated to the probability that a new lobe will "
                                                                                                     "be generated by a young or old (far or close to the vent when the "
                                                                                                     "flag start_from_dist_flag=true) lobe. The closer is lobe_exponent to 0 the "
                                                                                                     "larger is the probability that the new lobe will be generated from a "
                                                                                                     "younger lobe.\n"
                                                                                                     "lobe_exponent = 1  => there is a uniform probability distribution "
                                                                                                     "assigned to all the existing lobes for the choice "
                                                                                                     "of the next lobe from which a new lobe will be "
                                                                                                     "generated. \n"
                                                                                                     "lobe_exponent = 0  => the new lobe is generated from the last one.")
                                            ("max_slope_prob", boost::program_options::value<double>(), "max_slope_prob is related to the probability that the direction of "
                                                                                                        "the new lobe is close to the maximum slope direction:\n"
                                                                                                        "max_slope_prob = 0 => all the directions have the same probability;\n"
                                                                                                        "max_slope_prob > 0 => the maximum slope direction has a larger "
                                                                                                        "probaiblity, and it increases with increasing "
                                                                                                        "value of the parameter;\n"
                                                                                                        "max_slope_prob = 1 => the direction of the new lobe is the maximum "
                                                                                                        "slope direction.")
                                              ("inertial_exponent", boost::program_options::value<double>(), "Inertial exponent: \n"
                                                                                                             "inertial_exponent = 0 => the max probability direction for the new lobe is the "
                                                                                                             "max slope direction;\n"
                                                                                                             "inertial_exponent > 0 => the max probability direction for the new lobe takes "
                                                                                                             "into account also the direction of the parent lobe and "
                                                                                                             "the inertia increaes with increasing exponent")
      ;
    auto demCrop = boost::program_options::options_description("Parameters to crop the DEM file. This reduces the computational "
                                                               "time when the DEM is large and the area covered by the flow "
                                                               "is a lot smaller");
    demCrop.add_options()
      ("east_to_vent", boost::program_options::value<double>(), "Distance from the vent to the eastern border of the DEM")
        ("west_to_vent", boost::program_options::value<double>(), "Distance from the vent to the western border of the DEM")
          ("north_to_vent", boost::program_options::value<double>(), "Distance from the vent to the northern border of the DEM")
            ("south_to_vent", boost::program_options::value<double>(), "Distance from the vent to the southern border of the DEM");

    auto ventCoords = boost::program_options::options_description("Coordinates of the vents");
    ventCoords.add_options()
      ("x_vent", boost::program_options::value<std::vector<int>>(), "Vent coordinates. Two values define a fissure instead.")
        ("y_vent", boost::program_options::value<std::vector<int>>(), "Vent coordinates. Two values define a fissure instead.");

    auto ventCoordsExtra = boost::program_options::options_description("Extra coordinates of the vents. "
                                                                       "these coordinates are used when multiple fissures are defined:\n"
                                                                       "the first one goes from (x_vent[0],y_vent[0]) to (x_vent_end[0],y_vent_end[0])\n"
                                                                       "the second one goes from (x_vent[1],y_vent[1]) to (x_vent_end[1],y_vent_end[1])");
    ventCoordsExtra.add_options()
      ("x_vent_end", boost::program_options::value<std::vector<int>>(), "Vent end x coordinates.")
        ("y_vent_end", boost::program_options::value<std::vector<int>>(), "Vent end y coordinates.");

    ventCoords.add(ventCoordsExtra);

    auto advancedOptions = boost::program_options::options_description("Advanced options");
    advancedOptions.add_options()
      ("restart_files", boost::program_options::value<std::vector<std::string>>(), "this is a list of names of exiting *_thickness_*.asc files to be used if "
                                                                                   "you want to start a simulation considering the existence of previous "
                                                                                   "flows. Leave the field empy [] if you don't want to use it. "
                                                                                   "If you want to use more than one file separate them with a comma. "
                                                                                   "Example: restart_file = ['flow_001_thickness_full.asc','flow_002_thickness_full.asc]")
        ("saveshape_flag", boost::program_options::value<bool>(), "If saveshape_flag = true then the output with the lobes is saved on a shapefile")
          ("saveraster_flag", boost::program_options::value<bool>(), "If saveraster_flag = true then the raster output is saved as a *.asc file")
            ("flag_threshold", boost::program_options::value<int>(), "Flag to select if it is cutted the volume of the area\n"
                                                                   "flag_threshold = 1  => volume\n"
                                                                   "flag_threshold = 2  => area")
              ("force_max_length", boost::program_options::value<bool>(), "Flag for maximum distances (number of chained lobes) from the vent")
                ("max_length", boost::program_options::value<int>(), "Maximum distances (number of chained lobes) from the vent. "
                                                                   "This parameter is used only when force_max_length = 1")
                  ("n_init", boost::program_options::value<int>(), "Number of repetitions of the first lobe (useful for initial spreading)")
                    ("n_check_loop", boost::program_options::value<bool>(), "This parameter is to avoid that a chain of loop get stuck in a hole. It is "
                                                                          "active only when n_check_loop>0. When it is greater than zero the code "
                                                                          "check if the last n_check_loop lobes are all in a small box. If this is the "
                                                                          "case then the slope modified by flow is evaluated, and the hole is filled.")
                      ("start_from_dist_flag", boost::program_options::value<bool>(), "This flag controls which lobes have larger probability:\n"
                                                                                      "start_from_dist_flag = true  => the lobes with a larger distance from "
                                                                                      "the vent have a higher probability\n"
                                                                                      "start_form_dist_flag = false  => the younger lobes have a higher "
                                                                                      "probability")
                        ("dist_fact", boost::program_options::value<bool>(), "This factor is to choose where the center of the new lobe will be:\n"
                                                                             "dist_fact = 0  => the center of the new lobe is on the border of the "
                                                                             "previous one;\n"
                                                                             "dist fact > 0  => increase the distance of the center of the new lobe "
                                                                             "from the border of the previous one;\n"
                                                                             "dist_fact = 1  => the two lobes touch in one point only.")
                          ("npoints", boost::program_options::value<int>(), "Number of points for the lobe profile")
                            ("aspect_ratio_coeff", boost::program_options::value<double>(), "This parameter affect the shape of the lobes. The larger is this parameter, "
                                                                                            "the larger is the effect of a small slope on the eccentricity of the lobes:\n"
                                                                                            "aspect_ratio_coeff = 0  => the lobe is always a circle\n"
                                                                                            "aspect_ratio_coeff > 0  => the lobe is an ellipse, with the aspect ratio "
                                                                                            "increasing with the slope ")
                              ("max_aspect_ratio", boost::program_options::value<double>(), "Maximum aspect ratio of the lobes")
                                ("shape_name", boost::program_options::value<std::string>(), "Shapefile name (use '' if no shapefile is present)");

    auto betaOptions = boost::program_options::options_description("Beta options: The number of lobes of the flow is defined accordingly to a random uniform "
                                                                   "distribution or to a beta law, as a function of the flow number.\n"
                                                                   "a_beta, b_beta = 0  => n_lobes is sampled randomly in [min_n_lobes,max_n_lobes]\n"
                                                                   "a_beta, b_beta > 0  => n_lobes = min_n_lobes + 0.5 * ( max_n_lobes - min_n_lobes ) "
                                                                   "* beta(flow/n_flows,a_beta,b_beta)");
    betaOptions.add_options()
      ("a_beta", boost::program_options::value<double>(), "Parameter a of the beta distribution")
        ("b_beta", boost::program_options::value<double>(), "Parameter b of the beta distribution");

    advancedOptions.add(betaOptions);

    inputDataProperties.add(demCrop);
    inputDataProperties.add(ventCoords);

    auto allOptions = boost::program_options::options_description();
    allOptions.add(desc).add(inputDataProperties).add(advancedOptions);
    return allOptions;
}

Params::Params(int argc, char** argv)
: desc(getDesc())
{
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    // read config file if present
    if (vm.contains("input_file")) {
        auto ifs = std::ifstream(vm["input_file"].as<std::string>().c_str());
        if (!ifs) {
            throw std::runtime_error("can not open config file " + vm["input_file"].as<std::string>());
        }
        boost::program_options::store(boost::program_options::parse_config_file(ifs, desc), vm);
        boost::program_options::notify(vm);
    }
}
auto
Params::contains(const std::string& key) const -> bool
{
    return vm.contains(key);
}
auto
Params::getVm() const -> const boost::program_options::variables_map&
{
    return vm;
}
auto
Params::operator[](const std::string& key) const
  -> const boost::program_options::variable_value&
{
    return vm[key];
}