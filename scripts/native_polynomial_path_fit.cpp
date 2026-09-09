// Native driver for OpenSim's PolynomialPathFitter, plus an independent sampler.
//
// Why this file exists: the engine does not run the model's GeometryPaths.
// `native_mechanical_stream.cpp` replaces all 80 with fitted multivariate
// polynomials read from a separate FunctionBasedPathSet file.  An isotropic
// stature change can be pushed through those polynomials exactly (path length
// is homogeneous of degree one in the geometry), but an ANISOTROPIC change --
// a wider pelvis, a different shoulder:hip ratio -- cannot: it moves 56 of 98
// muscle paths by between 0.01% and 6.23% and leaves 42 alone, and a single
// coefficient factor is wrong for every muscle but at most one.  Those paths
// have to be REFITTED.  PolynomialPathFitter does that and is built in
// libosimActuators; only a driver was missing.
//
// Two subcommands:
//
//   fit     run PolynomialPathFitter on a model and a coordinate trajectory,
//           writing <output-dir>/<modelname>_FunctionBasedPathSet.xml.
//
//   sample  evaluate path lengths and moment arms at every row of a coordinate
//           table and write them as CSV.  Runs against the model's own
//           GeometryPaths, or -- with --pathset -- against a FunctionBasedPathSet
//           substituted in by ModOpReplacePathsWithFunctionBasedPaths.  This is
//           the independent read-back: the comparison between the two is
//           computed in Python from CSV, not taken from the fitter's own log.
//
// Every numeric result this program produces is written to stdout or to a file;
// nothing is judged here.  The gates live in scripts/verify_path_refitting.py.
#include <OpenSim/OpenSim.h>
#include <OpenSim/Actuators/PolynomialPathFitter.h>
#include <OpenSim/Actuators/ModelOperators.h>
#include <OpenSim/Simulation/TableProcessor.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenSim;

namespace {

struct Args {
    std::map<std::string, std::string> single;
    std::multimap<std::string, std::string> repeated;

    bool has(const std::string& k) const { return single.count(k) > 0; }
    std::string get(const std::string& k, const std::string& fallback = "") const {
        auto it = single.find(k);
        return it == single.end() ? fallback : it->second;
    }
    double num(const std::string& k, double fallback) const {
        auto it = single.find(k);
        return it == single.end() ? fallback : std::stod(it->second);
    }
    int integer(const std::string& k, int fallback) const {
        auto it = single.find(k);
        return it == single.end() ? fallback : std::stoi(it->second);
    }
    std::vector<std::string> all(const std::string& k) const {
        std::vector<std::string> out;
        auto range = repeated.equal_range(k);
        for (auto it = range.first; it != range.second; ++it) out.push_back(it->second);
        return out;
    }
};

Args parseArgs(int argc, char** argv, int from) {
    Args args;
    for (int i = from; i < argc; ++i) {
        std::string token(argv[i]);
        if (token.rfind("--", 0) != 0) {
            throw std::runtime_error("Expected --flag, got: " + token);
        }
        std::string key = token.substr(2);
        std::string value = "1";
        auto eq = key.find('=');
        if (eq != std::string::npos) {
            value = key.substr(eq + 1);
            key = key.substr(0, eq);
        } else if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
            value = argv[++i];
        }
        args.repeated.emplace(key, value);
        args.single[key] = value;
    }
    return args;
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string item;
    std::istringstream stream(s);
    while (std::getline(stream, item, sep)) out.push_back(item);
    return out;
}

// Keep every `stride`-th row of a table, counting from row 0.  The upstream
// example does exactly this to `coordinates.sto` before fitting; replicating it
// is what makes the reproduction gate a like-for-like comparison.
void decimate(TimeSeriesTable& table, int stride) {
    if (stride <= 1) return;
    auto times = table.getIndependentColumn();
    for (int i = static_cast<int>(times.size()) - 1; i >= 0; --i) {
        if (i % stride != 0) table.removeRow(times[i]);
    }
}

// A coordinate the trajectory does not mention still has to be given values;
// the fitter requires a column for every Coordinate in the model.  `--fill-sine`
// reproduces the upstream example's treatment of the toe joints.
void appendSineColumn(TimeSeriesTable& table, const std::string& column,
        double amplitude, double frequency, double phase) {
    Sine sine(amplitude, frequency, phase);
    auto times = table.getIndependentColumn();
    SimTK::Vector values(static_cast<int>(table.getNumRows()), 0.0);
    SimTK::Vector t(1, 0.0);
    for (int i = 0; i < static_cast<int>(table.getNumRows()); ++i) {
        t.set(0, times[i]);
        values.set(i, sine.calcValue(t));
    }
    table.appendColumn(column, values);
}

void appendConstantColumn(TimeSeriesTable& table, const std::string& column,
        double value) {
    SimTK::Vector values(static_cast<int>(table.getNumRows()), value);
    table.appendColumn(column, values);
}

TimeSeriesTable loadValues(const Args& args, const Model& model) {
    TimeSeriesTable table(args.get("coordinates"));
    decimate(table, args.integer("row-stride", 1));

    // --fill-sine <coordinate-path>,<amplitude>,<frequency>,<phase>
    for (const auto& spec : args.all("fill-sine")) {
        auto parts = split(spec, ',');
        if (parts.size() != 4) throw std::runtime_error("--fill-sine needs path,amp,freq,phase");
        appendSineColumn(table, parts[0], std::stod(parts[1]), std::stod(parts[2]),
                std::stod(parts[3]));
    }
    // --fill-constant <coordinate-path>,<value>
    for (const auto& spec : args.all("fill-constant")) {
        auto parts = split(spec, ',');
        if (parts.size() != 2) throw std::runtime_error("--fill-constant needs path,value");
        appendConstantColumn(table, parts[0], std::stod(parts[1]));
    }
    // Any coordinate still missing is filled at its default value, and the fact
    // is printed rather than absorbed -- a silently defaulted coordinate is a
    // coordinate the fit never explored.
    if (args.has("fill-missing-with-default")) {
        const auto& labels = table.getColumnLabels();
        for (const auto& coordinate : model.getComponentList<Coordinate>()) {
            std::string column = coordinate.getAbsolutePathString() + "/value";
            if (std::find(labels.begin(), labels.end(), column) == labels.end()) {
                std::cout << "filled-with-default " << column << " "
                          << coordinate.getDefaultValue() << "\n";
                appendConstantColumn(table, column, coordinate.getDefaultValue());
            }
        }
    }
    return table;
}

int runFit(const Args& args) {
    Model model(args.get("model"));
    model.initSystem();

    PolynomialPathFitter fitter;
    fitter.setModel(ModelProcessor(model));
    fitter.setCoordinateValues(TableProcessor(loadValues(args, model)));
    fitter.setOutputDirectory(args.get("output-dir", "."));
    fitter.setMaximumPolynomialOrder(args.integer("max-order", 5));
    fitter.setMinimumPolynomialOrder(args.integer("min-order", 2));
    fitter.setNumSamplesPerFrame(args.integer("samples-per-frame", 10));
    fitter.setUseStepwiseRegression(args.integer("stepwise", 1) != 0);
    fitter.setPathLengthTolerance(args.num("path-length-tolerance", 1e-3));
    fitter.setMomentArmTolerance(args.num("moment-arm-tolerance", 1e-3));
    if (args.has("threads")) fitter.setNumParallelThreads(args.integer("threads", 4));
    if (args.has("latin-hypercube")) fitter.setLatinHypercubeAlgorithm(args.get("latin-hypercube"));

    auto bounds = split(args.get("global-bounds", "-30,30"), ',');
    fitter.setGlobalCoordinateSamplingBounds(
            SimTK::Vec2(std::stod(bounds[0]), std::stod(bounds[1])));
    // --bounds <coordinate-path>,<lo>,<hi>   (degrees; repeatable)
    for (const auto& spec : args.all("bounds")) {
        auto parts = split(spec, ',');
        if (parts.size() != 3) throw std::runtime_error("--bounds needs path,lo,hi");
        fitter.appendCoordinateSamplingBounds(
                parts[0], SimTK::Vec2(std::stod(parts[1]), std::stod(parts[2])));
    }

    fitter.run();
    std::cout << "fit-complete " << args.get("output-dir", ".") << "/"
              << model.getName() << "_FunctionBasedPathSet.xml\n";
    return 0;
}

int runSample(const Args& args) {
    Model model = args.has("pathset")
            ? [&] {
                  ModelProcessor processor(args.get("model"));
                  processor.append(ModOpReplacePathsWithFunctionBasedPaths(
                          args.get("pathset")));
                  return processor.process();
              }()
            : Model(args.get("model"));
    SimTK::State state = model.initSystem();

    TimeSeriesTable table(args.get("coordinates"));
    decimate(table, args.integer("row-stride", 1));
    const auto& labels = table.getColumnLabels();

    // Which coordinates get a moment arm.  Default: every coordinate the model
    // holds, so nothing is quietly excluded; --moment-arm-coordinates narrows it.
    std::vector<std::string> armCoordinates;
    if (args.has("moment-arm-coordinates")) {
        armCoordinates = split(args.get("moment-arm-coordinates"), ',');
    } else {
        for (const auto& coordinate : model.getComponentList<Coordinate>()) {
            armCoordinates.push_back(coordinate.getAbsolutePathString());
        }
    }

    std::vector<const PathActuator*> actuators;
    for (const auto& actuator : model.getComponentList<PathActuator>()) {
        actuators.push_back(&actuator);
    }
    std::sort(actuators.begin(), actuators.end(),
            [](const PathActuator* a, const PathActuator* b) {
                return a->getAbsolutePathString() < b->getAbsolutePathString();
            });

    std::ofstream out(args.get("output"));
    out.precision(17);
    out << "row,time,actuator,quantity,coordinate,value\n";

    for (size_t row = 0; row < table.getNumRows(); ++row) {
        double time = table.getIndependentColumn()[row];
        state.setTime(time);
        for (size_t c = 0; c < labels.size(); ++c) {
            std::string label = labels[c];
            const std::string suffix = "/value";
            if (label.size() <= suffix.size() ||
                    label.compare(label.size() - suffix.size(), suffix.size(), suffix) != 0) {
                continue;
            }
            std::string coordinatePath = label.substr(0, label.size() - suffix.size());
            if (!model.hasComponent<Coordinate>(coordinatePath)) continue;
            model.getComponent<Coordinate>(coordinatePath)
                    .setValue(state, table.getRowAtIndex(row)[static_cast<int>(c)], false);
        }
        model.assemble(state);
        model.realizePosition(state);

        for (const auto* actuator : actuators) {
            std::string name = actuator->getAbsolutePathString();
            out << row << "," << time << "," << name << ",length,,"
                << actuator->getLength(state) << "\n";
            for (const auto& coordinatePath : armCoordinates) {
                if (!model.hasComponent<Coordinate>(coordinatePath)) continue;
                auto& coordinate = const_cast<Coordinate&>(
                        model.getComponent<Coordinate>(coordinatePath));
                double arm = actuator->computeMomentArm(state, coordinate);
                out << row << "," << time << "," << name << ",moment_arm,"
                    << coordinatePath << "," << arm << "\n";
            }
        }
    }
    out.close();
    std::cout << "sample-complete " << args.get("output") << " rows=" << table.getNumRows()
              << " actuators=" << actuators.size() << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: native_polynomial_path_fit {fit|sample} --flag value ...\n";
        return 2;
    }
    try {
        Logger::setLevel(Logger::Level::Info);
        std::string command(argv[1]);
        Args args = parseArgs(argc, argv, 2);
        if (command == "fit") return runFit(args);
        if (command == "sample") return runSample(args);
        std::cerr << "Unknown subcommand: " << command << "\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "ERROR " << e.what() << "\n";
        return 1;
    }
}
