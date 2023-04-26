#include <iostream>
#include <string>

#include <fmt/core.h>
#include <boost/program_options.hpp>
#include <fstream>

auto main(int argc, char** argv) -> int
{
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

    auto ventCoordsExtra = boost::program_options::options_description("Extra coordinates of the vents "
                                                                       "this coordinates are used when multiple fissures are defined:\n"
                                                                       "the first one goes from (x_vent[0],y_vent[0]) to (x_vent_end[0],y_vent_end[0])\n"
                                                                       "the second one goes from (x_vent[1],y_vent[1]) to (x_vent_end[1],y_vent_end[1])");
    ventCoordsExtra.add_options()
      ("x_vent_end", boost::program_options::value<std::vector<int>>(), "Vent end coordinates.")
      ("y_vent_end", boost::program_options::value<std::vector<int>>(), "Vent end coordinates.");

    ventCoords.add(ventCoordsExtra);

    inputDataProperties.add(demCrop);
    inputDataProperties.add(ventCoords);

    auto allOptions = boost::program_options::options_description();
    allOptions.add(desc).add(inputDataProperties);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, allOptions), vm);
    boost::program_options::notify(vm);
    // read config file if present
    if (vm.contains("config")) {
      auto ifs = std::ifstream(vm["config"].as<std::string>().c_str());
      if (!ifs) {
        std::cout << "can not open config file: " << vm["config"].as<std::string>() << "\n";
        return 0;
      }
      boost::program_options::store(boost::program_options::parse_config_file(ifs, allOptions), vm);
      boost::program_options::notify(vm);
    }

    if (vm.contains("help")) {
      std::cout << allOptions << "\n";
      return 0;
    }

    if (vm.contains("version")) {
      std::cout << "Version: " << "0.1.0" << "\n";
      return 0;
    }

    return 0;
}
