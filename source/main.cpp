#include <iostream>
#include <string>

#include <fmt/core.h>
#include <boost/program_options.hpp>
#include <fstream>
#include <Eigen/Dense>
#include <filesystem>
#include "Params.hpp"

auto
calculateCumFissLength(
  const std::vector<int>& xVent,
  const std::vector<int>& yVent,
  const std::optional<std::vector<int>>& xVentEnd,
  const std::optional<std::vector<int>>& yVentEnd,
  int ventFlag,
  int nVents,
  const std::optional<std::vector<double>>& fissureProbabilities)
  -> Eigen::Vector2d
{
    auto xVentEndMode =
      xVentEnd.has_value() && !xVentEnd->empty() && ventFlag > 3;
    auto [firstJ, cumFissLength] = [&]() {
        if (xVentEndMode) {
            return std::make_pair(0, Eigen::Vector2d::Zero(nVents + 1).eval());
        }
        return std::make_pair(1, Eigen::Vector2d::Zero(nVents).eval());
    }();

    for (auto j = firstJ; j < nVents; ++j) {
        if (xVentEndMode) {
            auto deltaXVent = xVentEnd->at(j) - xVent[j];
            auto deltaYVent = yVentEnd->at(j) - yVent[j];
            cumFissLength[j + 1] =
              cumFissLength[j] +
              std::sqrt(deltaXVent * deltaXVent + deltaYVent * deltaYVent);
        } else {
            auto deltaXVent = xVent[j] - xVent[j - 1];
            auto deltaYVent = yVent[j] - yVent[j - 1];
            cumFissLength[j] =
              cumFissLength[j - 1] +
              std::sqrt(deltaXVent * deltaXVent + deltaYVent * deltaYVent);
        }
        return cumFissLength;
    }

    if (fissureProbabilities.has_value() && ventFlag > 5) {
        auto fissureProbabilitiesCum =
          Eigen::Vector2d::Zero(fissureProbabilities->size()).eval();
        fissureProbabilitiesCum[0] = (*fissureProbabilities)[0];
        for (auto j = 1; j < fissureProbabilities->size(); ++j) {
            fissureProbabilitiesCum[j] =
              fissureProbabilitiesCum[j - 1] + fissureProbabilities->at(j);
        }
        if (ventFlag == 8) {
            // cumulative sum
            cumFissLength = fissureProbabilitiesCum;
        } else if (ventFlag > 5) {
            // keep first element, copy the rest from fissureProbabilitiesCum
            for (auto j = 1; j < fissureProbabilities->size(); ++j) {
                cumFissLength[j] = fissureProbabilitiesCum[j];
            }
        }
    }

    if (nVents > 1) {
        cumFissLength = cumFissLength / cumFissLength[nVents - 1];
    }

    return cumFissLength;
}

auto
createBackup(const std::string& fileName) -> void
{
    auto i = 0;
    auto condition = true;
    auto baseName = fileName;
    while (condition) {
        auto runName = fmt::format("{}_{:03}", baseName, i);
        auto backupFile = fmt::format("{}.bak", runName);
        condition = std::filesystem::exists(backupFile);
        ++i;
    }
    std::filesystem::copy_file(fileName, fmt::format("{}.bak", baseName));
}

auto
main(int argc, char** argv) -> int
{

    auto params = Params{ argc, argv };

    if (params.contains("help")) {
        std::cout << params << "\n";
        return 0;
    }

    if (params.contains("version")) {
        std::cout << "Version: "
                  << "0.1.0"
                  << "\n";
        return 0;
    }

    auto fillingParameter = 1.0 - params["filling_parameter"].as<double>();
    auto nVents = params["x_vent"].as<std::vector<int>>().size();

    auto ventFlag = params["vent_flag"].as<int>();

    auto xVent = params["x_vent"].as<std::vector<int>>();
    auto yVent = params["y_vent"].as<std::vector<int>>();
    auto xVentEnd = params.getOpt<std::vector<int>>("x_vent_end");
    auto yVentEnd = params.getOpt<std::vector<int>>("y_vent_end");
    auto fissureProbabilities =
      params.getOpt<std::vector<double>>("fissure_probabilities");

    auto cumFissLength = calculateCumFissLength(
      xVent, yVent, xVentEnd, yVentEnd, ventFlag, nVents, fissureProbabilities);

    if (params.contains("input_file")) {
        createBackup(params["input_file"].as<std::string>());
    }

    auto aBeta = params["a_beta"].as<double>();
    auto bBeta = params["b_beta"].as<double>();

    auto maxNLobes = params["max_n_lobes"].as<int>();
    auto minNLobes = params["min_n_lobes"].as<int>();
    auto allocNLobes = [&]() {
        if (aBeta == 0.0 && bBeta == 0.0) {
            return params["max_n_lobes"].as<int>();
        }
        auto nFlows = params["n_flows"].as<int>();
        auto xBeta =
          Eigen::VectorXd::LinSpaced(nFlows, 0, nFlows - 1) / (nFlows - 1);

        // TODO: check math!
        auto betaPdf =
          (xBeta.array().pow(aBeta) * (1.0 - xBeta.array()).pow(bBeta))
            .matrix();

        return static_cast<int>(
          std::round(minNLobes +
                     0.5 * (maxNLobes -
                            minNLobes) *
                       betaPdf.maxCoeff()));
    }();

    auto angle = Eigen::VectorXd::Zero(allocNLobes).eval();
    auto x = Eigen::VectorXd::Zero(allocNLobes).eval();
    auto y = Eigen::VectorXd::Zero(allocNLobes).eval();
    auto x1 = Eigen::VectorXd::Zero(allocNLobes).eval();
    auto x2 = Eigen::VectorXd::Zero(allocNLobes).eval();
    auto h = Eigen::VectorXd::Zero(allocNLobes).eval();

    auto distInt = Eigen::VectorXi::Constant(allocNLobes, -1).eval();
    auto descendents = Eigen::VectorXi::Zero(allocNLobes).eval();
    auto parent = Eigen::VectorXi::Zero(allocNLobes).eval();
    auto alfaInertial = Eigen::VectorXd::Zero(allocNLobes).eval();


    if (params["volume_flag"].as<bool>()) {
        auto lobeArea = params["lobe_area"].as<double>();
        auto totalVolume = params["total_volume"].as<double>();
        auto nFlows = params["n_flows"].as<int>();

        if (!params["fixed_dimension_flag"].as<bool>()) {
            // avg_lobe_thickness = total_volume / ( n_flows * lobe_area * 0.5 * ( min_n_lobes + max_n_lobes ) )
            auto avgLobeThickness =
              totalVolume /
              (nFlows * lobeArea *
               0.5 * (minNLobes +
                      maxNLobes));
        }
        lobeArea = totalVolume /
                   (nFlows *
                    params["avg_lobe_thickness"].as<double>() *
                    0.5 * (minNLobes +
                           maxNLobes));
    }

    auto npoints = params["npoints"].as<int>();
    auto t = Eigen::VectorXd::LinSpaced(npoints, 0, 2.0 * M_PI);
    auto xCircle = t.array().cos();
    auto yCircle = t.array().sin();

    // # Parse the header using a loop and
    // # the built-in linecache module
    // hdr = [getline(source, i) for i in range(1,7)]
    // values = [float(h.split(" ")[-1].strip()) for h in hdr]
    // del hdr

    auto source = std::ifstream{ params["input_file"].as<std::string>() };
    auto line = std::string{};
    std::getline(source, line);
    auto values = std::array<double, 6>{};
    for (auto i = 0; i < 6; ++i) {
        std::getline(source, line);
        auto lastSpace = line.find_last_of(" ");
        values[i] = std::stod(line.substr(lastSpace + 1));
    }


    auto [cols, rows, lx, ly, cell, nd] = values;

    cols = std::round(cols);
    rows = std::round(rows);

    auto crop_flag = params.contains("west_to_vent") && params.contains("east_to_vent") &&
                     params.contains("south_to_vent") && params.contains("north_to_vent");

    auto start = std::chrono::high_resolution_clock::now();






    return 0;
}
