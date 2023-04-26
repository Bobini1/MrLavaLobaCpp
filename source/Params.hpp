//
// Created by bobini on 26.04.23.
//

#ifndef MRLAVALOBA_PARAMS_HPP
#define MRLAVALOBA_PARAMS_HPP

#include <boost/program_options/variables_map.hpp>
#include <optional>
class Params
{
    boost::program_options::variables_map vm;
    boost::program_options::options_description desc;
  public:
    Params(int argc, char** argv);
    auto getVm() const -> const boost::program_options::variables_map&;
    template<typename T>
    auto get(const std::string& key) const -> T
    {
        return vm[key].as<T>();
    }
    auto operator[](const std::string& key) const -> const boost::program_options::variable_value&;
    template<typename T>
    auto getOpt(const std::string& key) const -> std::optional<T> {
        if (vm.contains(key)) {
            return vm[key].as<T>();
        }
        return std::nullopt;

    }
    [[nodiscard]] auto contains(const std::string& key) const -> bool;

    friend auto operator<<(std::ostream& out, const Params& params) -> std::ostream& {
        out << params.desc;
        return out;
    }
};

#endif // MRLAVALOBA_PARAMS_HPP
