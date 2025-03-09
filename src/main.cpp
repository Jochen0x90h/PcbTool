#include <fstream>
#include <filesystem>
#include <iostream>
#include <regex>
#include "kicad.hpp"


namespace fs = std::filesystem;


// remove quotes from string
std::string clean(std::string_view s) {
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"')
        return std::string(s.substr(1, s.size() - 2));
    return std::string(s);
}

bool match(const std::string &str, const std::list<std::regex> &list) {
    for (auto regex : list) {
        if (std::regex_match(str, regex))
            return true;
    }
    return false;
}

int main(int argc, const char **argv) {
    if (argc < 2)
        return 1;

    std::list<std::regex> hideReference;
    std::list<std::regex> removeNet;
    fs::path path;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "-hr") {
            // hide reference
            ++i;
            hideReference.emplace_back(argv[i]);
        } else if (arg == "-rn") {
            // remove net
            ++i;
            removeNet.emplace_back(argv[i]);
        } else {
            path = argv[i];
        }
    }

    if (path.empty() || (hideReference.empty() && removeNet.empty())) {
        std::cout << "Usage: pcb-tool [options] pcb.kicad_pcb" << std::endl;
        std::cout << "    -hr <hide reference>, supports regular expressions, e.g. R*" << std::endl;
        std::cout << "    -rn <remove net>, supports regular expressions, e.g. GND" << std::endl;
    }

    // read .kicad_pcb file
    std::ifstream in(path.string());
    if (!in) {
        // error
        return 1;
    }
    kicad::Container file;
    kicad::readFile(in, file);
    in.close();

    for (auto element1 : file.elements) {
        auto container1 = dynamic_cast<kicad::Container *>(element1);
        if (container1) {
            // check if it is a footprint
            if (container1->id == "footprint") {
                auto footprint = container1;

                // get footprint name
                auto footprintName = clean(footprint->getValue(0));
                //std::cout << "Footprint: " << footprintName << std::endl;

                // get reference and hide if requested
                std::string reference;
                for (auto element2 : container1->elements) {
                    auto container2 = dynamic_cast<kicad::Container *>(element2);
                    if (container2) {
                        if (container2->id == "property") {
                            auto property = container2;
                            if (property->getValue(0) == "\"Reference\"") {
                                reference = clean(property->getValue(1));

                                // check if reference matches a reference to hide
                                if (match(reference, hideReference)) {
                                    std::cout << "Hide reference: " << reference << std::endl;

                                    auto hide = property->find("hide");
                                    property->erase(hide);

                                    // add hide attribute
                                    property->add("hide", "yes");
                                }
                            }
                        }
                    }
                }

                // remove nets from pads if requested
                for (auto element2 : container1->elements) {
                    auto container2 = dynamic_cast<kicad::Container *>(element2);
                    if (container2) {
                        if (container2->id == "pad") {
                            auto pad = container2;

                            auto net = pad->find("net");
                            if (net != nullptr) {
                                auto netName = clean(net->getValue(1));

                                // check if net matches a net to remove
                                if (match(netName, removeNet)) {
                                    // erase net
                                    std::cout << "Erase net " << netName << " of " << reference << std::endl;
                                    pad->erase(net);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // write .kicad_pcb file
    std::ofstream out(path.string());
    if (!out) {
        // error
        return 1;
    }
    kicad::writeFile(out, file);
    out.close();

    return 0;
}