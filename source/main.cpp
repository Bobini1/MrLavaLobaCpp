#include <iostream>
#include <string>

#include <fmt/core.h>
#include <boost/program_options.hpp>
#include <fstream>
#include <Eigen/Dense>
#include <filesystem>
#include <random>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/KroneckerProduct>

#include "Params.hpp"
#include "rtnorm.hpp"

auto
calculateCumFissLength(const std::vector<int>& xVent,
                       const std::vector<int>& yVent,
                       int nVents) -> Eigen::VectorXd
{
    auto cumFissLength = Eigen::VectorXd::Zero(nVents).eval();
    for (auto j = 1; j < nVents; ++j) {
        auto deltaXVent = xVent[j] - xVent[j - 1];
        auto deltaYVent = yVent[j] - yVent[j - 1];
        cumFissLength[j] =
          cumFissLength[j - 1] +
          std::sqrt(deltaXVent * deltaXVent + deltaYVent * deltaYVent);
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

template<typename Floating>
auto signedFmod(Floating x, Floating y) -> Floating
{
    auto result = std::fmod(x, y);
    if (result < 0) {
        result += y;
    }
    return result;
}

auto
loadTxt(const std::string& filename,
        size_t rows,
        size_t cols,
        size_t skipRows = 0) -> Eigen::MatrixXd
{
    auto file = std::ifstream(filename);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    auto line = std::string();
    auto matrix = Eigen::MatrixXd::Zero(rows, cols).eval();
    auto i = 0;
    while (std::getline(file, line)) {
        if (i >= skipRows) {
            auto j = 0;
            auto lineStream = std::stringstream(line);
            auto value = std::string();
            while (std::getline(lineStream, value, ' ')) {
                matrix(i - skipRows, j) = std::stod(value);
                ++j;
            }
        }
        ++i;
    }
    return matrix;
}

auto
ellipse(double xc,
        double yc,
        double ax1,
        double ax2,
        double angle,
        Eigen::MatrixXd X_circle,
        Eigen::MatrixXd Y_circle) -> std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
{
    auto cos_angle = std::cos(angle * M_PI / 180);
    auto sin_angle = std::sin(angle * M_PI / 180);

    auto X = (ax1 * X_circle).eval();
    auto Y = (ax2 * Y_circle).eval();

    auto xe = (xc + X.array() * cos_angle - Y.array() * sin_angle).eval();
    auto ye = (yc + X.array() * sin_angle + Y.array() * cos_angle).eval();

    return { xe, ye };
}

auto
localIntersection(const Eigen::MatrixXd& XsLocal,
                  const Eigen::MatrixXd& YsLocal,
                  double xcE,
                  double ycE,
                  double ax1,
                  double ax2,
                  double angle,
                  const Eigen::VectorXd& xv,
                  const Eigen::VectorXd& yv,
                  double nv2) -> Eigen::MatrixXd
{
    auto nxCell = XsLocal.rows();
    auto nyCell = XsLocal.cols();

    auto c = std::cos(angle * M_PI / 180.0);
    auto s = std::sin(angle * M_PI / 180.0);

    auto c1 = c / ax1;
    auto s1 = s / ax1;

    auto c2 = c / ax2;
    auto s2 = s / ax2;

    Eigen::VectorXd xvLocal = xv.array() - xcE;
    Eigen::VectorXd yvLocal = yv.array() - ycE;

    Eigen::MatrixXd XsLocalReshaped = XsLocal.transpose();
    Eigen::MatrixXd YsLocalReshaped = YsLocal.transpose();
    Eigen::VectorXd XsLocal1d = Eigen::Map<const Eigen::VectorXd>(
      XsLocalReshaped.data(), XsLocalReshaped.size());
    Eigen::VectorXd YsLocal1d = Eigen::Map<const Eigen::VectorXd>(
      YsLocalReshaped.data(), YsLocalReshaped.size());


    Eigen::VectorXd c1xvPS1yv = c1 * xvLocal.array() + s1 * yvLocal.array();
    Eigen::VectorXd c2yvMS2yv = c2 * yvLocal.array() - s2 * xvLocal.array();

    Eigen::VectorXd term1 = (c1 * c1 + s2 * s2) * XsLocal1d.array().square();
    Eigen::VectorXd term2 =
      (2.0 * c1 * s1 - 2.0 * c2 * s2) * XsLocal1d.array() * YsLocal1d.array();
    Eigen::MatrixXd term3 = Eigen::kroneckerProduct(
      XsLocal1d, 2.0 * c1 * c1xvPS1yv - 2.0 * s2 * c2yvMS2yv);
    term3.resize(XsLocal.size(), term3.size() / XsLocal.size());
    Eigen::VectorXd term4 = (c2 * c2 + s1 * s1) * YsLocal1d.array().square();
    Eigen::MatrixXd term5 = Eigen::kroneckerProduct(
      YsLocal1d, 2.0 * c2 * c2yvMS2yv + 2.0 * s1 * c1xvPS1yv);
    term5.resize(XsLocal.size(), term5.size() / XsLocal.size());
    Eigen::VectorXd term6 =
      c1xvPS1yv.array().square() + c2yvMS2yv.array().square();

    Eigen::VectorXd term124 = term1 + term2 + term4;
    Eigen::MatrixXd term356 = (term3 + term5).reshaped(term3.cols(), term3.rows()).colwise() + term6;

    Eigen::MatrixXd termTot = term356.rowwise() + term124.transpose();

    Eigen::MatrixXd inside = (termTot.array() <= 1.0).cast<double>();

    Eigen::VectorXd areaFract1d = inside.colwise().sum();

    areaFract1d /= nv2;

    Eigen::MatrixXd areaFract =
      Eigen::Map<const Eigen::MatrixXd>(areaFract1d.data(), nxCell, nyCell);

    return areaFract;
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

    auto fillingParameter = 1.0 - params["thickening_parameter"].as<double>();
    auto nVents = params["x_vent"].as<std::vector<int>>().size();

    auto xVent = params["x_vent"].as<std::vector<int>>();
    auto yVent = params["y_vent"].as<std::vector<int>>();

    auto cumFissLength = calculateCumFissLength(xVent, yVent, nVents);

    // skip shapefile stuff

    // skip backup stuff

    auto aBeta = params["a_beta"].as<double>();
    auto bBeta = params["b_beta"].as<double>();

    auto maxNLobes = params["max_n_lobes"].as<int>();
    auto minNLobes = params["min_n_lobes"].as<int>();
    auto allocNLobes = params["max_n_lobes"].as<int>(); // skip if

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

    auto lobeArea = params["lobe_area"].as<double>();
    auto totalVolume = params["total_volume"].as<double>();
    auto nFlows = params["n_flows"].as<int>();
    auto avgLobeThickness =
      totalVolume / (nFlows * lobeArea * 0.5 * (minNLobes + maxNLobes));

    auto npoints = params["npoints"].as<int>();
    auto t = Eigen::VectorXd::LinSpaced(npoints, 0, 2.0 * M_PI).eval();
    auto xCircle = t.array().cos().eval();
    auto yCircle = t.array().sin().eval();

    auto source = params["source"].as<std::string>();
    auto values = std::array<double, 6>{};

    {
        auto sourceStream = std::ifstream{ source };
        if (!sourceStream) {
            std::cerr << "Could not open source file: " << source << "\n";
            std::cerr << "reason: " << std::strerror(errno) << "\n";
            std::cout << "File path " << source << " at absolute location "
                      << std::filesystem::absolute(source)
                      << " does not exist\n";
            return 1;
        }
        auto line = std::string{};
        for (auto i = 0; i < 6; ++i) {
            std::getline(sourceStream, line);
            auto lastSpace = line.find_last_of(' ');
            values[i] = std::stod(line.substr(lastSpace + 1));
        }
    }

    auto [cols, rows, lx, ly, cell, nd] = values;

    cols = std::round(cols);
    rows = std::round(rows);

    auto arr = loadTxt(params["source"].as<std::string>(), rows, cols, 6);

    auto nx = arr.cols();
    Eigen::VectorXd xs =
      lx + cell * (0.5 + Eigen::VectorXd::LinSpaced(nx, 0, nx - 1).array());
    auto xmin = xs.minCoeff();
    auto xmax = xs.maxCoeff();

    auto ny = arr.rows();
    Eigen::VectorXd ys =
      ly + cell * (0.5 + Eigen::VectorXd::LinSpaced(ny, 0, ny - 1).array());
    auto ymin = ys.minCoeff();
    auto ymax = ys.maxCoeff();

    std::cout << xmin << " " << ymin << "\n";

    auto Xs = Eigen::MatrixXd::Zero(ny, nx).eval();
    auto Ys = Eigen::MatrixXd::Zero(ny, nx).eval();
    for (auto i = 0; i < ny; ++i) {
        Xs.row(i) = xs;
    }
    for (auto i = 0; i < nx; ++i) {
        Ys.col(i) = ys;
    }

    auto Zs = Eigen::MatrixXd::Zero(ny, nx).eval();
    for (auto i = 0; i < ny; ++i) {
        Zs.row(i) = arr.row(ny - i - 1);
    }

    // skip restart files

    auto nv = 20;
    auto xvYvInit = Eigen::VectorXd::LinSpaced(nv, -0.5 * cell, 0.5 * cell).eval();
    auto xv = Eigen::MatrixXd::Zero(nv, nv).eval();
    auto yv = Eigen::MatrixXd::Zero(nv, nv).eval();
    for (auto i = 0; i < nv; ++i) {
        xv.col(i) = xvYvInit;
    }
    for (auto i = 0; i < nv; ++i) {
        yv.row(i) = xvYvInit;
    }
    xv = xv.reshaped(nv * nv, 1);
    yv = yv.reshaped(nv * nv, 1);

    auto nv2 = nv * nv;

    auto Ztot = Zs;
    auto ZtotTemp = Zs;

    auto nTest = 0;

    // random device
    auto rd = std::random_device{};
    // random number generator
    auto gen = std::mt19937(rd());

    auto gslGen = gsl_rng_alloc(gsl_rng_mt19937);

    auto Xs1d = Eigen::VectorXd::Map(Xs.data(), Xs.size());
    auto Ys1d = Eigen::VectorXd::Map(Ys.data(), Ys.size());
    auto nxy = Xs1d.size();

    auto points = Eigen::MatrixXd::Zero(nxy, 2).eval();
    auto Zflow = Eigen::MatrixXd::Zero(ny, nx).eval();

    auto maxAspectRatio = params["max_aspect_ratio"].as<double>();
    auto maxSemiaxis = std::sqrt(lobeArea * maxAspectRatio / M_PI);
    auto maxCells = std::ceil(2.0 * maxSemiaxis / cell) + 4;
    auto maxCellsInt = static_cast<int>(maxCells);

    std::cout << "max_semiaxis " << maxSemiaxis << "\n";

    auto jtopArray = Eigen::VectorXi::Zero(allocNLobes).eval();
    auto jbottomArray = Eigen::VectorXi::Zero(allocNLobes).eval();

    auto irightArray = Eigen::VectorXi::Zero(allocNLobes).eval();
    auto ileftArray = Eigen::VectorXi::Zero(allocNLobes).eval();

    auto zhazard = Eigen::MatrixXd::Zero(ny, nx).eval();
    auto zhazardTemp = Eigen::MatrixXd::Zero(ny, nx).eval();

    Eigen::MatrixXd zdist = Zflow.array() + 9999.0;

    std::cout << "End pre-processing\n\n";

    auto flowsCounter = 0;
    auto nLobesTot = 0;

    auto thicknessRatio = params["thickness_ratio"].as<double>();

    auto nFlowsCounter = params["n_flows_counter"].as<int>();

    auto nInit = params["n_init"].as<int>();

    auto maxSlopeProb = params["max_slope_prob"].as<double>();

    auto aspectRatioCoeff = params["aspect_ratio_coeff"].as<double>();

    auto lobeExponent = params["lobe_exponent"].as<double>();

    auto inertialExponent = params["inertial_exponent"].as<double>();

    auto distFact = params["dist_fact"].as<double>();

    auto nLobesCounter = params["n_lobes_counter"].as<int>();

    auto estRemTime = std::string{};

    auto start = std::chrono::steady_clock::now();

    for (size_t flow = 0; flow < nFlows; flow++) {
        auto ZflowLocalArray =
          Eigen::Tensor<int, 3>(allocNLobes, maxCellsInt, maxCellsInt).eval();
        auto iFirstCheck = 0;

        flowsCounter++;

        auto nLobes =
          std::uniform_real_distribution<double>(minNLobes, maxNLobes)(gen);

        nLobesTot += nLobes;

        auto thicknessMin =
          2.0 * thicknessRatio / (thicknessRatio + 1.0) * avgLobeThickness;
        auto deltaLobeThickness =
          2.0 * (avgLobeThickness - thicknessMin) / (nLobes - 1.0);

        auto lastPercentage5 =
          static_cast<int>(std::round(flow * 20.0 / nFlows));
        auto lastPercentage =
          static_cast<int>(std::round(flow * 100.0 / nFlows));

        std::cout << "\r[" << std::string(lastPercentage5, '=')
                  << std::string(20 - lastPercentage5, ' ') << "] "
                  << lastPercentage << "% " << estRemTime << std::flush;

        if (flowsCounter == nFlowsCounter) {
            flowsCounter = 0;

            ZtotTemp = Ztot;
        }

        auto lobesCounter = 0;

        for (size_t i = 0; i < nInit; i++) {
            if (nFlows == 1) {
                auto lastPercentage =
                  static_cast<int>(std::round(i * 20.0 / (nInit - 1)) * 5);

                std::cout << "\r[" << std::string(lastPercentage / 5, '=')
                          << std::string(20 - lastPercentage / 5, ' ') << "] "
                          << lastPercentage << "%" << std::flush;
            }

            auto iVent = static_cast<int>(std::floor(flow * nVents / nFlows));

            x(i) = xVent[iVent];
            y(i) = yVent[iVent];

            distInt(i) = 0;
            descendents(i) = 0;

            auto xi = (x(i) - xmin) / cell;
            auto yi = (y(i) - ymin) / cell;

            auto ix = static_cast<int>(std::floor(xi));
            auto iy = static_cast<int>(std::floor(yi));

            auto xiFract = xi - ix;
            auto yiFract = yi - iy;

            auto FxTest =
              (xiFract * (Ztot(iy + 1, ix + 1) - Ztot(iy + 1, ix)) +
               (1.0 - xiFract) * (Ztot(iy, ix + 1) - Ztot(iy, ix))) /
              cell;
            auto FyTest =
              (yiFract * (Ztot(iy + 1, ix + 1) - Ztot(iy, ix + 1)) +
               (1.0 - yiFract) * (Ztot(iy + 1, ix) - Ztot(iy, ix))) /
              cell;

            auto maxSlopeAngle = signedFmod(
              180.0 + (180.0 * std::atan2(FyTest, FxTest) / M_PI), 360.0);
            auto slope = std::sqrt(std::pow(FxTest, 2) + std::pow(FyTest, 2));

            auto slopeDeg = 180.0 * std::atan(slope) / M_PI;

            auto randAngleNew = [&]() {
                if (slopeDeg > 0) {
                    auto sigma = (1.0 - maxSlopeProb) / maxSlopeProb *
                                 (90.0 - slopeDeg) / slopeDeg;
                    return rtnorm(gslGen, -180.0, 180.0, 0.0, sigma).first;
                }
                auto rand =
                  std::uniform_real_distribution<double>(0.0, 1.0)(gen);
                return 360.0 * std::abs(rand - 0.5);
            }();
            angle(i) = maxSlopeAngle + randAngleNew;

            auto aspectRatio =
              std::min(maxAspectRatio, 1.0 + aspectRatioCoeff * slope);

            x1(i) = std::sqrt(lobeArea / M_PI) * std::sqrt(aspectRatio);
            x2(i) = std::sqrt(lobeArea / M_PI) / std::sqrt(aspectRatio);

            auto [xe, ye] =
              ellipse(x(i), y(i), x1(i), x2(i), angle(i), xCircle, yCircle);

            auto minXe = xe.minCoeff();
            auto maxXe = xe.maxCoeff();

            auto minYe = ye.minCoeff();
            auto maxYe = ye.maxCoeff();

            auto iLeft = static_cast<int>(std::round(
              std::distance(xs.begin(),
                            std::find_if(xs.begin(),
                                         xs.end(),
                                         [&](auto x) { return x > minXe; })) -
              2));
            auto iRight = static_cast<int>(std::round(
              std::distance(xs.begin(),
                            std::find_if(xs.begin(),
                                         xs.end(),
                                         [&](auto x) { return x > maxXe; })) +
              2));

            auto jBottom = static_cast<int>(std::round(
              std::distance(ys.begin(),
                            std::find_if(ys.begin(),
                                         ys.end(),
                                         [&](auto y) { return y > minYe; })) -
              2));
            auto jTop = static_cast<int>(std::round(
              std::distance(ys.begin(),
                            std::find_if(ys.begin(),
                                         ys.end(),
                                         [&](auto y) { return y > maxYe; })) +
              2));

            Eigen::MatrixXd xsLocal =
              Xs.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft);
            Eigen::MatrixXd ysLocal =
              Ys.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft);

            Eigen::MatrixXd areaFract = localIntersection(xsLocal,
                                                          ysLocal,
                                                          x(i),
                                                          y(i),
                                                          x1(i),
                                                          x2(i),
                                                          angle(i),
                                                          xv,
                                                          yv,
                                                          nv2);

            auto zflowLocal = areaFract;
            auto zflowLocalInt =
              areaFract.unaryExpr([](auto x) { return std::ceil(x); })
                .cast<int>().eval();

            auto lobeThickness = thicknessMin + (i - 1) * deltaLobeThickness;

            Zflow.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft) +=
              lobeThickness * zflowLocal;
            ZtotTemp.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft) =
              Zs.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft) +
              fillingParameter *
                Zflow.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft);

            auto zdistLocal = zflowLocalInt
                                .unaryExpr([&](auto x) {
                                    return x * distInt(i) + 9999 * (x == 0);
                                })
                                .eval();
            zdist.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft) =
              zdistLocal.cast<double>().cwiseMin(
                zdist.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft));

            jtopArray(i) = jTop;
            jbottomArray(i) = jBottom;
            irightArray(i) = iRight;
            ileftArray(i) = iLeft;

            lobesCounter++;
        }
        for (int i = nInit; i < nLobes; i++) {
            auto idx0 = std::uniform_real_distribution<double>(0.0, 1.0)(gen);
            auto idx1 = std::pow(idx0, lobeExponent);

            auto idx2 = i * idx1;
            auto idx3 = std::floor(idx2);
            auto idx = static_cast<int>(idx3);

            parent(i) = idx;
            distInt(i) = distInt(idx) + 1;

            auto last = i;

            for (int j = 0; j < distInt(idx) + 1; j++) {
                auto previous = parent(last);
                descendents(previous) = descendents(previous) + 1;
                last = previous;
            }

            auto xi = (x(idx) - xmin) / cell;
            auto yi = (y(idx) - ymin) / cell;

            auto ix = static_cast<int>(std::floor(xi));
            auto iy = static_cast<int>(std::floor(yi));

            auto ix1 = ix + 1;
            auto iy1 = iy + 1;

            if ((ix <= 1) || (ix1 >= nx - 1) || (iy <= 1) || (iy1 >= ny - 1)) {
                break;
            }

            auto xiFract = xi - ix;
            auto yiFract = yi - iy;

            auto fxLobe = (xiFract * (Ztot(iy1, ix1) - Ztot(iy1, ix)) +
                           (1.0 - xiFract) * (Ztot(iy, ix1) - Ztot(iy, ix))) /
                          cell;
            auto fyLobe = (yiFract * (Ztot(iy1, ix1) - Ztot(iy, ix1)) +
                           (1.0 - yiFract) * (Ztot(iy1, ix) - Ztot(iy, ix))) /
                          cell;

            auto slope = std::sqrt(std::pow(fxLobe, 2) + std::pow(fyLobe, 2));
            auto maxSlopeAngle =
              signedFmod(180 + (180 * std::atan2(fyLobe, fxLobe) / M_PI), 360.0);

            auto slopeDeg = 180.0 * std::atan(slope) / M_PI;

            auto randAngleNew = [&]() {
                if (slopeDeg > 0.0) {
                    auto sigma = (1.0 - maxSlopeProb) / maxSlopeProb *
                                 (90.0 - slopeDeg) / slopeDeg;
                    return rtnorm(gslGen, -180, 180, 0, sigma).first;
                }
                return std::uniform_real_distribution<double>(-180.0,
                                                              180.0)(gen);
            }();

            auto newAngle = maxSlopeAngle + randAngleNew;

            constexpr auto deg2rad = M_PI / 180.0;

            auto cosAngle1 = std::cos(angle(idx) * deg2rad);
            auto sinAngle1 = std::sin(angle(idx) * deg2rad);

            auto cosAngle2 = std::cos(newAngle * deg2rad);
            auto sinAngle2 = std::sin(newAngle * deg2rad);

            alfaInertial(i) = std::pow(
              1.0 - std::pow(2.0 * std::atan(slope) / M_PI, inertialExponent),
              1.0 / inertialExponent);

            auto cosAngleAvg =
              (1.0 - alfaInertial(i)) * cosAngle2 + alfaInertial(i) * cosAngle1;
            auto sinAngleAvg =
              (1.0 - alfaInertial(i)) * sinAngle2 + alfaInertial(i) * sinAngle1;

            auto angleSigned = 180 * std::atan2(sinAngleAvg, cosAngleAvg) / M_PI;
            auto angleAvg =
              signedFmod(angleSigned, 360.0);

            newAngle = angleAvg;

            auto a = std::tan(deg2rad * (newAngle - angle(idx)));

            auto xt = [&]() {
                auto xt = std::sqrt(std::pow(x1(idx), 2) * std::pow(x2(idx), 2) /
                                           (std::pow(x2(idx), 2) + std::pow(x1(idx), 2) * std::pow(a, 2)));
                if (std::cos(deg2rad * (newAngle - angle(idx))) <= 0) {
                    xt = -xt;
                }
                return xt;
            }();

            auto yt = a * xt;

            auto deltaX = xt * cosAngle1 - yt * sinAngle1;
            auto deltaY = xt * sinAngle1 + yt * cosAngle1;

            xi = (x(idx) + deltaX - xmin) / cell;
            yi = (y(idx) + deltaY - ymin) / cell;

            ix = static_cast<int>(std::floor(xi));
            iy = static_cast<int>(std::floor(yi));

            ix1 = ix + 1;
            iy1 = iy + 1;

            if ((ix <= 1) || (ix1 >= nx - 1) || (iy <= 1) || (iy1 >= ny - 1)) {
                break;
            }

            xiFract = xi - ix;
            yiFract = yi - iy;

            fxLobe = (xiFract * (Ztot(iy1, ix1) - Ztot(iy1, ix)) +
                      (1.0 - xiFract) * (Ztot(iy, ix1) - Ztot(iy, ix))) /
                     cell;
            fyLobe = (yiFract * (Ztot(iy1, ix1) - Ztot(iy, ix1)) +
                      (1.0 - yiFract) * (Ztot(iy1, ix) - Ztot(iy, ix))) /
                     cell;

            //         slope = np.sqrt(np.square(Fx_lobe)+np.square(Fy_lobe))
            //        aspect_ratio = min(max_aspect_ratio,1.0 +
            //        aspect_ratio_coeff * slope)

            slope = std::sqrt(std::pow(fxLobe, 2) + std::pow(fyLobe, 2));
            auto aspectRatio =
              std::min(maxAspectRatio, 1.0 + aspectRatioCoeff * slope);

            auto newx1 = std::sqrt(lobeArea / M_PI) * std::sqrt(aspectRatio);
            auto newx2 = std::sqrt(lobeArea / M_PI) / std::sqrt(aspectRatio);

            auto v1 = std::sqrt(std::pow(deltaX, 2) + std::pow(deltaY, 2));

            auto v2 = v1 + newx1;

            auto v = (v1 * (1.0 - distFact) + v2 * distFact) / v1;

            auto xNew = x(idx) + v * deltaX;
            auto yNew = y(idx) + v * deltaY;

            angle(idx) = newAngle;
            x1(idx) = newx1;
            x2(idx) = newx2;
            x(idx) = xNew;
            y(idx) = yNew;

            auto [xe, ye] = ellipse(
              x(idx), y(idx), x1(idx), x2(idx), angle(idx), xCircle, yCircle);

            auto minXe = xe.minCoeff();
            auto maxXe = xe.maxCoeff();

            auto minYe = ye.minCoeff();
            auto maxYe = ye.maxCoeff();

            auto iParent = parent(idx);

            auto iLeft = [&]() {
                if (minXe < xs(0)) {
                    return 0l;
                } else if (minXe >= xs(nx - 1)) {
                    return nx - 1;
                } else {
                    return static_cast<long>(std::distance(
                             xs.data(),
                             std::upper_bound(
                               xs.data(), xs.data() + nx, minXe))) -
                           2;
                }
            }();

            auto iRight = [&]() {
                if (maxXe < xs(0)) {
                    return 0l;
                } else if (maxXe >= xs(nx - 1)) {
                    return nx - 1;
                } else {
                    return static_cast<long>(std::distance(
                             xs.data(),
                             std::upper_bound(
                               xs.data(), xs.data() + nx, maxXe))) +
                           2;
                }
            }();

            auto jBottom = [&]() {
                if (minYe < ys(0)) {
                    return 0l;
                } else if (minYe >= ys(ny - 1)) {
                    return ny - 1;
                } else {
                    return static_cast<long>(std::distance(
                             ys.data(),
                             std::upper_bound(
                               ys.data(), ys.data() + ny, minYe))) -
                           2;
                }
            }();

            auto jTop = [&]() {
                if (maxYe < ys(0)) {
                    return 0l;
                } else if (maxYe >= ys(ny - 1)) {
                    return ny - 1;
                } else {
                    return static_cast<long>(std::distance(
                             ys.data(),
                             std::upper_bound(
                               ys.data(), ys.data() + ny, maxYe))) +
                           2;
                }
            }();

            Eigen::MatrixXd xsLocal = Xs.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft);
            Eigen::MatrixXd ysLocal = Ys.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft);

            auto areaFract = localIntersection(xsLocal,
                                               ysLocal,
                                               x(idx),
                                               y(idx),
                                               x1(idx),
                                               x2(idx),
                                               angle(idx),
                                               xv,
                                               yv,
                                               nv2);

            Eigen::MatrixXd zFlowLocal = areaFract;
            Eigen::MatrixXd zDistLocal =
              zFlowLocal * distInt(idx) + (9999 * zFlowLocal.array() == 0).cast<double>().matrix();

            zdist.block(
              jBottom, iLeft, jTop - jBottom, iRight - iLeft) =
              zdist
                .block(jBottom, iLeft, jTop - jBottom, iRight - iLeft)
                .cwiseMin(zDistLocal);

            auto lobeThickness = thicknessMin + (i - 1) * deltaLobeThickness;

            Zflow.block(
              jBottom, iLeft, jTop - jBottom, iRight - iLeft) +=
              lobeThickness * zFlowLocal;
            ZtotTemp.block(
              jBottom, iLeft, jTop - jBottom, iRight - iLeft) =
              Zs.block(jBottom, iLeft, jTop - jBottom, iRight - iLeft) +
              fillingParameter *
                Zflow.block(
                  jBottom, iLeft, jTop - jBottom, iRight - iLeft);

            jtopArray(idx) = jTop;
            jbottomArray(idx) = jBottom;
            irightArray(idx) = iRight;
            ileftArray(idx) = iLeft;

            lobesCounter++;

            if (lobesCounter == nLobesCounter) {
                lobesCounter = 0;
                Ztot = ZtotTemp;
            }
        }

        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
          std::chrono::steady_clock::now() - start);
        auto estimated =
          std::ceil(elapsed.count() * nFlows / (flow + 1) - elapsed.count());
        estRemTime = std::to_string(estimated);
    }

    std::cout << "\r[" << std::string(20, '=') << "] 100%" << std::endl;

    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - start);

    std::cout << "Total number of lobes " << nLobesTot
              << " Average number of lobes "
              << static_cast<long>(std::round(1.0 * nLobesTot / nFlows))
              << std::endl;
    std::cout << "Time elapsed " << elapsed.count() << " sec." << std::endl;
    std::cout << "Saving files" << std::endl;

    auto runName = params["run_name"].as<std::string>();
    auto outputFull = runName + "_thickness_full.asc";
    {
        std::ofstream zflowFile(outputFull);
        zflowFile << "ncols " << Zflow.cols() << std::endl;
        zflowFile << "nrows " << Zflow.rows() << std::endl;
        zflowFile << "xllcorner " << lx << std::endl;
        zflowFile << "yllcorner " << ly << std::endl;
        zflowFile << "cellsize " << cell << std::endl;
        zflowFile << "NODATA_value 0" << std::endl;

        zflowFile << Zflow.colwise().reverse().format(
          Eigen::IOFormat(5,
                          Eigen::DontAlignCols,
                          " ",
                          "\n",
                          "",
                          "",
                          "",
                          ""));
    }
    std::cout << outputFull << " saved" << std::endl;

    auto maxLobes =
      static_cast<long>(std::floor(Zflow.maxCoeff() / avgLobeThickness));

    auto maskingThreshold = params["masking_threshold"].as<double>();

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> maskedZflow;
    for (int i = 0; i < 10 * maxLobes; i++) {
        maskedZflow = (Zflow.array() < i * 0.1 * avgLobeThickness).eval();

        auto totalZflow = Zflow.sum();

        auto volumeFraction = maskedZflow.sum() / totalZflow;

        auto coverageFraction = volumeFraction;

        if (coverageFraction < maskingThreshold) {
            std::cout << "Total volume " << cell * cell * totalZflow
                      << " Masked volume " << cell * cell * maskedZflow.sum()
                      << " Volume fraction " << coverageFraction << std::endl;



            auto outputMasked = runName + "_thickness_masked.asc";
            {
                std::ofstream zflowMaskedFile(outputMasked);
                zflowMaskedFile << "ncols " << Zflow.cols() << std::endl;
                zflowMaskedFile << "nrows " << Zflow.rows() << std::endl;
                zflowMaskedFile << "xllcorner " << lx << std::endl;
                zflowMaskedFile << "yllcorner " << ly << std::endl;
                zflowMaskedFile << "cellsize " << cell << std::endl;
                zflowMaskedFile << "NODATA_value 0" << std::endl;

                zflowMaskedFile
                  << ((1 - maskedZflow.cast<double>().array()) * Zflow.array())
                       .colwise()
                       .reverse()
                       .format(Eigen::IOFormat(5,
                                               Eigen::DontAlignCols,
                                               " ",
                                               "\n",
                                               "",
                                               "",
                                               "",
                                               ""));
            }

            std::cout << outputMasked << " saved" << std::endl;

            break;
        }
    }

    auto outputDist = runName + "_dist_full.asc";
    {
        std::ofstream zdistFile(outputDist);
        zdistFile << "ncols " << zdist.cols() << std::endl;
        zdistFile << "nrows " << zdist.rows() << std::endl;
        zdistFile << "xllcorner " << lx << std::endl;
        zdistFile << "yllcorner " << ly << std::endl;
        zdistFile << "cellsize " << cell << std::endl;
        zdistFile << "NODATA_value 0" << std::endl;

        zdistFile << zdist.colwise().reverse().format(
          Eigen::IOFormat(Eigen::FullPrecision,
                          Eigen::DontAlignCols,
                          " ",
                          "\n",
                          "",
                          "",
                          "",
                          ""));
    }

    std::cout << outputDist << " saved" << std::endl;

    zdist = ((1 - maskedZflow.array()).cast<double>() * zdist.array()).eval();

    auto outputDistMasked = runName + "_dist_masked.asc";
    {
        std::ofstream zdistMaskedFile(outputDistMasked);
        zdistMaskedFile << "ncols " << zdist.cols() << std::endl;
        zdistMaskedFile << "nrows " << zdist.rows() << std::endl;
        zdistMaskedFile << "xllcorner " << lx << std::endl;
        zdistMaskedFile << "yllcorner " << ly << std::endl;
        zdistMaskedFile << "cellsize " << cell << std::endl;
        zdistMaskedFile << "NODATA_value 0" << std::endl;

        zdistMaskedFile << zdist.matrix().colwise().reverse().format(
          Eigen::IOFormat(Eigen::FullPrecision,
                          Eigen::DontAlignCols,
                          " ",
                          "\n",
                          "",
                          "",
                          "",
                          ""));
    }

    gsl_rng_free(gslGen);

    return 0;
}
