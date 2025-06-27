#include <fstream>
#include <filesystem>
#include <iostream>
#include <regex>
#include "kicad.hpp"


namespace fs = std::filesystem;

constexpr float pif = 3.1416f;


template <typename T>
struct Vector2 {
    T x;
    T y;
};

template <typename T1, typename T2>
inline auto operator +(const Vector2<T1> &a, const Vector2<T2> &b) {
    return Vector2<decltype(a.x + b.x)>(a.x + b.x, a.y + b.y);
}

template <typename T1, typename T2>
inline auto operator -(const Vector2<T1> &a, const Vector2<T2> &b) {
    return Vector2<decltype(a.x - b.x)>(a.x - b.x, a.y - b.y);
}

using float2 = Vector2<float>;

float distance(float2 a, float2 b) {
    float2 d = a - b;
    return sqrt(d.x * d.x + d.y * d.y);
}


struct Orientation {
    float2 position;
    float rotation;

    /// @brief Transform a point from local to global coordinates
    /// @param p Point
    /// @return Point in global coordinates
    float2 toGlobal(float2 p) {
        float r = this->rotation * pif / 180.0f;
        float s = sin(r);
        float c = cos(r);

        return {this->position.x + c * p.x + s * p.y, this->position.y + c * p.y - s * p.x};
    }

    /// @brief Transform a point from global to local coordinates
    /// @param p Point
    /// @return Point in local coordinates
    float2 toLocal(float2 p) {
        float r = -this->rotation * pif / 180.0f;
        float s = sin(r);
        float c = cos(r);

        return {(p.x - this->position.x) * c + (p.y - this->position.y) * s,
            (p.y - this->position.y) * c - (p.x - this->position.x) * s};
    }

    /// @brief Transform an orientation from local to global coordinates
    /// @param o Orientation
    /// @return Orientation in global coordinates
    Orientation toGlobal(Orientation o) {
        float r = this->rotation * pif / 180.0f;
        float s = sin(r);
        float c = cos(r);

        return {{this->position.x + c * o.position.x + s * o.position.y,
            this->position.y + c * o.position.y - s * o.position.x},
            this->rotation + o.rotation};
    }

    /// @brief Transform an orientation from global to local coordinates
    /// @param o Orientation
    /// @return Orientation in local coordinates
    Orientation toLocal(Orientation o) {
        float r = -this->rotation * pif / 180.0f;
        float s = sin(r);
        float c = cos(r);

        return {{(o.position.x - this->position.x) * c + (o.position.y - this->position.y) * s,
            (o.position.y - this->position.y) * c - (o.position.x - this->position.x) * s},
            o.rotation - this->rotation};
    }
};

struct Pad {
    std::string name;
    kicad::Container *at;
    float2 position;
    int netNumber;
};

struct Footprint {
    std::string name;
    kicad::Container *at;
    Orientation orientation;
    std::string reference;
    std::string value;
    std::vector<Pad> pads;

    void place(const Orientation &o) {
        this->at->setFloat(0, o.position.x);
        this->at->setFloat(1, o.position.y);
        this->at->setFloat(2, o.rotation);
        for (auto &pad : pads) {
            pad.at->setFloat(2, o.rotation);
        }
    }
};

struct Template {
    std::list<Footprint> refFootprints;
    std::list<Footprint> placementFootprints;
};


// remove quotes from string
std::string clean(std::string_view s) {
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"')
        return std::string(s.substr(1, s.size() - 2));
    return std::string(s);
}

// check if a string matches at least one regluar expression in the list
bool match(const std::string &str, const std::list<std::regex> &list) {
    for (auto regex : list) {
        if (std::regex_match(str, regex))
            return true;
    }
    return false;
}


void getFootprints(kicad::Container &file, std::list<Footprint>& refFootprints,
    std::list<Footprint> &placementFootprints)
{
    for (auto element1 : file.elements) {
        auto container1 = dynamic_cast<kicad::Container *>(element1);
        if (container1) {
            // check if it is a footprint
            if (container1->id == "footprint") {
                auto footprint = container1;

                // get footprint name
                auto footprintName = footprint->getString(0);
                //std::cout << "Footprint: " << footprintName << std::endl;

                // iterate over properties of footprint
                kicad::Container *at = nullptr;
                float2 position = {0, 0};
                float rotation = 0;
                std::string reference;
                std::string value;
                std::vector<Pad> pads;
                for (auto element2 : footprint->elements) {
                    auto container2 = dynamic_cast<kicad::Container *>(element2);
                    if (container2) {
                        if (container2->id == "at") {
                            at = container2;
                            position.x = at->getFloat(0);
                            position.y = at->getFloat(1);
                            rotation = at->getFloat(2);
                        } else if (container2->id == "property") {
                            auto property = container2;
                            std::string propertyName = property->getString(0);
                            if (propertyName == "Reference") {
                                // get reference and hide if requested
                                reference = property->getString(1);
                            } else if (propertyName == "Value") {
                                // get reference and hide if requested
                                value = property->getString(1);
                            }
                        }
                    }
                }

                bool addToRefFootprints = at != nullptr && (value.starts_with("TLV3544IPW") || value.starts_with("TSV912AIDSG"));
                bool addToPlacementFootprints = at != nullptr && footprintName == "local:0402";

                if (addToRefFootprints || addToPlacementFootprints) {
                    // iterate over pads of footprint
                    for (auto element2 : footprint->elements) {
                        auto container2 = dynamic_cast<kicad::Container *>(element2);
                        if (container2) {
                            if (container2->id == "pad") {
                                auto pad = container2;

                                // get pad name (typically a number)
                                std::string padName = pad->getString(0);

                                // get pad position
                                auto at = pad->find("at");
                                //auto at = pad->findFloat2("at");
                                auto x = at->getFloat(0);
                                auto y = at->getFloat(1);

                                // get net number
                                int netNumber = -1;
                                auto net = pad->find("net");
                                if (net != nullptr) {
                                    netNumber = net->getInt(0);
                                    //auto netName = net->getString(1);
                                }
                                pads.push_back({padName, at, {x, y}, netNumber});
                            }
                        }
                    }
                }

                if (addToRefFootprints) {
                    refFootprints.push_back({footprintName, at, {position, rotation}, reference, value, pads});
                } else if (addToPlacementFootprints) {
                    placementFootprints.push_back({footprintName, at, {position, rotation}, reference, value, pads});
                }
            }
        }
    }

}



int main(int argc, const char **argv) {
    if (argc < 2)
        return 1;

/*    Orientation o{{10.0f, 10.0f}, 90.0f};
    float2 p{0.0f, 1.0f};
    float2 q = o.toGlobal(p);
    float2 r = o.toLocal(q);*/

    // schematic options
    std::map<std::string, std::string> exchangeFootprint;

    // pcb options
    std::map<int, double> exchangeSegmentWidth;
    std::map<int, double> exchangeViaSize;
    std::list<std::regex> hideReference;
    std::list<std::regex> removeNet;
    std::list<std::regex> moveFootprint;

    std::list<Template> templates;

    // paths of files to process
    std::list<fs::path> paths;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "-xfp") {
            // exchange footprint
            i += 2;
            exchangeFootprint[argv[i - 1]] = argv[i];
        } else if (arg == "-rsw") {
            // exchange segment width
            i += 2;
            try {
                double oldWidth = std::stod(argv[i - 1]);
                double newWidth = std::stod(argv[i]);
                if (newWidth > 0)
                    exchangeSegmentWidth[int(std::round(oldWidth * 1000.0))] = newWidth;
            } catch (std::exception &) {
                // invalid width
            }
        } else if (arg == "-rvs") {
            // exchange via size
            i += 2;
            try {
                double oldSize = std::stod(argv[i - 1]);
                double newSize = std::stod(argv[i]);
                if (newSize > 0)
                    exchangeViaSize[int(std::round(oldSize * 1000.0))] = newSize;
            } catch (std::exception &) {
                // invalid width
            }
        } else if (arg == "-hr") {
            // hide reference
            ++i;
            hideReference.emplace_back(argv[i]);
        } else if (arg == "-rn") {
            // remove net
            ++i;
            removeNet.emplace_back(argv[i]);
        } else if (arg == "-mfp") {
            ++i;
            moveFootprint.emplace_back(argv[i]);
        } else if (arg == "-t") {
            // template
            ++i;
            fs::path path = argv[i];

            // read template
            std::ifstream in(path.string());
            if (in) {
                kicad::Container file;
                kicad::readFile(in, file);
                in.close();

                // get footprints from template
                auto &t = templates.emplace_back();
                getFootprints(file, t.refFootprints, t.placementFootprints);
                if (t.refFootprints.empty()) {
                    templates.pop_back();
                    std::cout << "Warning: Template " << path.string() << " has no reference footprint" << std::endl;
                } else if (t.refFootprints.size() > 1) {
                    std::cout << "Warning: Template " << path.string() << " has multiple reference footprints" << std::endl;
                }
            }
        } else {
            paths.emplace_back(argv[i]);
        }
    }

    bool hasExchangeSegmentWidth = !exchangeSegmentWidth.empty();
    bool hasExchangeViaSize = !exchangeViaSize.empty();
    if (paths.empty() || argc < 4) {
        std::cout << "Usage: pcb-tool [options] schematic.kicad_sch layout.kicad_pcb ..." << std::endl;
        std::cout << "    -xfp <old footprint> <new footprint>           Exchange footprint (only schematic)" << std::endl;
        std::cout << "    -xsw <old segment width> <new segment width>   Exchange segment width (only pcb)" << std::endl;
        std::cout << "    -xvs <old via size> <new via size>             Exchange segment width (only pcb)" << std::endl;
        std::cout << "    -hr <reference>                                Hide reference, supports regular expressions, e.g. R.* (only pcb)" << std::endl;
        std::cout << "    -rn <net>                                      Remove net, supports regular expressions, e.g. GND (only pcb)" << std::endl;
        std::cout << "    -mfp <footprint>                               Move footprint to border (only pcb)" << std::endl;
        return 1;
    }

    for (auto &path : paths) {
        bool isSchematic = path.extension() == ".kicad_sch";
        bool isPcb = path.extension() == ".kicad_pcb";

        if (!isSchematic && !isPcb) {
            std::cout << "Error: File type not supported" << std::endl;
            return 1;
        }

        // read .kicad_pcb file
        std::cout << "Read " << path.string() << std::endl;
        std::ifstream in(path.string());
        if (!in) {
            // error
            return 1;
        }
        kicad::Container file;
        kicad::readFile(in, file);
        in.close();

        if (isSchematic) {
            for (auto element1 : file.elements) {
                auto container1 = dynamic_cast<kicad::Container *>(element1);
                if (container1) {
                    // check if it is a symbol
                    if (container1->id == "symbol") {
                        // iterate over properties
                        for (auto element2 : container1->elements) {
                            auto container2 = dynamic_cast<kicad::Container *>(element2);
                            if (container2) {
                                if (container2->id == "property") {
                                    auto property = container2;
                                    if (property->getString(0) == "Footprint") {
                                        // replace footprint
                                        std::string footprint = property->getString(1);
                                        auto it = exchangeFootprint.find(footprint);
                                        if (it != exchangeFootprint.end()) {
                                            std::cout << "Replace footprint " << it->first << " by " << it->second << std::endl;
                                            property->elements[1] = new kicad::Value('"' + it->second + '"');
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (isPcb) {
            float2 moveFootprintPosition = {};
            for (auto element1 : file.elements) {
                auto container1 = dynamic_cast<kicad::Container *>(element1);
                if (container1) {
                    // check if it is a footprint
                    if (container1->id == "footprint") {
                        auto footprint = container1;

                        // get footprint name
                        auto footprintName = footprint->getString(0);
                        //std::cout << "Footprint: " << footprintName << std::endl;

                        // iterate over properties of footprint
                        float2 position = {0, 0};
                        float rotation = 0;
                        std::string reference;
                        std::string value;
                        for (auto element2 : footprint->elements) {
                            auto container2 = dynamic_cast<kicad::Container *>(element2);
                            if (container2) {
                                if (container2->id == "at") {
                                    auto at = container2;

                                    // move footprint away
                                    if (match(footprintName, moveFootprint)) {
                                        at->setFloat(0, moveFootprintPosition.x);
                                        at->setFloat(1, moveFootprintPosition.y);
                                        //at->setFloat(2, 0);
                                        moveFootprintPosition.y += 2;
                                    }

                                    position.x = at->getFloat(0);
                                    position.y = at->getFloat(1);
                                    rotation = at->getFloat(2);
                                } else if (container2->id == "property") {
                                    auto property = container2;
                                    std::string propertyName = property->getString(0);
                                    if (propertyName == "Reference") {
                                        // get reference and hide if requested
                                        reference = property->getString(1);

                                        // check if reference matches a reference to hide
                                        if (match(reference, hideReference)) {
                                            std::cout << "Hide reference: " << reference << std::endl;

                                            auto hide = property->find("hide");
                                            property->erase(hide);

                                            // add hide attribute
                                            property->add("hide", "yes");
                                        }
                                    } else if (propertyName == "Value") {
                                        // get reference and hide if requested
                                        value = property->getString(1);
                                    }
                                }
                            }
                        }

                        // iterate over pads of footprint
                        for (auto element2 : footprint->elements) {
                            auto container2 = dynamic_cast<kicad::Container *>(element2);
                            if (container2) {
                                if (container2->id == "pad") {
                                    auto pad = container2;

                                    // get pad name (typcally a number)
                                    //std::string padName = clean(pad->getValue(0));

                                    auto net = pad->find("net");
                                    if (net != nullptr) {
                                        auto netName = net->getString(1);

                                        // remove net from pad if requested
                                        if (match(netName, removeNet)) {
                                            // erase net
                                            std::cout << "Erase net " << netName << " of " << reference << std::endl;
                                            pad->erase(net);
                                        }
                                    }
                                }
                            }
                        }
                    } else if (hasExchangeSegmentWidth && container1->id == "segment") {
                        auto segment = container1;
                        auto width = segment->find("width");
                        if (width != nullptr) {
                            try {
                                double oldWidth = width->getDouble(0);
                                auto it = exchangeSegmentWidth.find(int(std::round(oldWidth * 1000.0)));
                                if (it != exchangeSegmentWidth.end()) {
                                    width->setDouble(0, it->second);
                                    //width->elements[0] = new kicad::Value(std::to_string(it->second));
                                }
                            } catch (std::exception &) {
                            }
                        }
                    } else if (hasExchangeViaSize && container1->id == "via") {
                        auto via = container1;
                        auto size = via->find("size");
                        if (size != nullptr) {
                            try {
                                double oldSize = size->getDouble(0);
                                auto it = exchangeViaSize.find(int(std::round(oldSize * 1000.0)));
                                if (it != exchangeViaSize.end()) {
                                    //size->elements[0] = new kicad::Value(std::to_string(it->second));
                                    size->setDouble(0, it->second);
                                }
                            } catch (std::exception &) {
                            }
                        }
                    }
                }
            }

            // automatic placement of resistors and capacitors
            if (!templates.empty()) {
                // collection of footprints for processing
                std::list<Footprint> refFootprints;
                std::list<Footprint> placementFootprints;
                getFootprints(file, refFootprints, placementFootprints);

                for (auto &refFootprint : refFootprints) {

                    // generate net positions and places
                    std::map<int, std::vector<float2>> netPositions;
                    std::list<Orientation> places;
                    /*for (auto &refFootprint : refFootprints)*/ {
                        // add net positions
                        for (auto &pad : refFootprint.pads) {
                            float2 p = refFootprint.orientation.toGlobal(pad.position);
                            //std::cout << "pad " << pad.name << " " << p.x << " " << p.y << std::endl;
                            netPositions[pad.netNumber].push_back(p);
                        }

                        // select template
                        for (auto &t : templates) {
                            auto &trf = t.refFootprints.front();
                            if (trf.name == refFootprint.name && refFootprint.value.starts_with(trf.value)) {
                                // copy places from template
                                for (auto &tpf : t.placementFootprints) {
                                    places.push_back(refFootprint.orientation.toGlobal(trf.orientation.toLocal(tpf.orientation)));
                                }
                                break;
                            }
                        }
                    }

                    while (!placementFootprints.empty() && !places.empty()) {
                        // find best placement (try all placement footprints on all available places)
                        float bestDistance = std::numeric_limits<float>::max();
                        std::list<Footprint>::iterator bestFootprint;
                        std::list<Orientation>::iterator bestPlace;
                        int bestFlip;
                        for (auto footrpintIt = placementFootprints.begin(); footrpintIt != placementFootprints.end(); ++footrpintIt) {
                            auto &footprint = *footrpintIt;

                            // test all places
                            for (auto placeIt = places.begin(); placeIt != places.end(); ++placeIt) {
                                auto &place = *placeIt;

                                // also test flipped by 180 degrees
                                for (int flip = 0; flip <= 1; ++flip) {
                                    // give capacitors a "head start" over resistors
                                    float d = footprint.reference.starts_with("C") ? -1.0f : 0.0f;

                                    // calculate distance sum
                                    for (auto &pad : footprint.pads) {
                                        auto orientation = place;
                                        orientation.rotation += flip ? 180.0f : 0.0f;
                                        auto padPosition = orientation.toGlobal(pad.position);

                                        float padDistance = 100.0f;
                                        if (pad.netNumber >= 0) {
                                            for (auto &netPosition : netPositions[pad.netNumber]) {
                                                padDistance = std::min(padDistance, distance(padPosition, netPosition));
                                            }
                                        }
                                        d += padDistance;
                                    }

                                    // check if smallest distance sum so far
                                    if (d < bestDistance) {
                                        bestDistance = d;
                                        bestFootprint = footrpintIt;
                                        bestPlace = placeIt;
                                        bestFlip = flip;
                                    }
                                }
                            }
                        }

                        if (bestDistance > 190.0f)
                            break;

                        // place
                        auto orientation = *bestPlace;
                        orientation.rotation += bestFlip ? 180.0f : 0.0f;
                        bestFootprint->place(orientation);
                        //bestFootprint->at->setFloat(0, orientation.position.x);
                        //bestFootprint->at->setFloat(1, orientation.position.y);
                        //bestFootprint->at->setFloat(2, orientation.rotation);

                        // remove
                        places.erase(bestPlace);
                        placementFootprints.erase(bestFootprint);
                    }
                }

                // move unplaced footprints away
                /*Orientation o = {{-5, 0}, 0};
                for (auto &footprint : placementFootprints) {
                    footprint.place(o);
                    o.position.y += 2;
                }*/
            }
        }

        // write .kicad_sch or .kicad_pcb file
        std::cout << "Write " << path.string() << std::endl;
        std::ofstream out(path.string());
        if (!out) {
            // error
            return 1;
        }
        kicad::writeFile(out, file);
        out.close();
    }
    return 0;
}