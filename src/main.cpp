#include <xlnt/xlnt.hpp> // https://github.com/tfussell/xlnt
#include <nlohmann/json.hpp>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <regex>
#include <set>
#include "kicad.hpp"


namespace fs = std::filesystem;
using json = nlohmann::json;
constexpr double pi = 3.1416f;


// check if a string matches at least one regluar expression in the list
bool match(const std::string &str, const std::list<std::regex> &list) {
    for (auto regex : list) {
        if (std::regex_match(str, regex))
            return true;
    }
    return false;
}

// get type from reference, e.g. "C20" -> "C"
std::string getType(std::string_view str) {
    for (int i = 0; i < str.length(); ++i) {
        char ch = str[i];
        if (ch >= '0' && ch <= '9')
            return std::string(str.substr(0, i));
    }
    return std::string(str);
}

// schematic

struct NetVoltage {
    std::regex net;
    double voltage;
};

std::optional<double> match(const std::string &str, const std::list<NetVoltage> &netVoltages) {
    for (auto &netVoltage : netVoltages) {
        if (std::regex_match(str, netVoltage.net)) {
            return netVoltage.voltage;
        }
    }
    return {};
}

// get the value of a cell in an excel, also when it is part of a merged range
std::string getCellValue(const xlnt::worksheet &sheet, int row, int col) {
    xlnt::cell_reference ref(col, row);

    // check merged ranges
    auto mergedRanges = sheet.merged_ranges();
    for (const auto& range : mergedRanges) {
        auto tl = range.top_left();
        auto br = range.bottom_right();
        if (tl.row() <= row && tl.column() <= col && br.row() >= row && br.column() >= col) {
            ref = tl;
            break;
        }
    }

    return sheet[ref].to_string();
}

struct Sheet {
    fs::path path;
    kicad::Container file;
    std::string uuidPath;
};

struct Properties {
    std::string manufacturer;
    std::string mpn;
    std::string description;
    std::string datasheet;
    std::string lcscPn;
};

void addProperty(kicad::Container *symbol, kicad::Container *property, std::string_view name, std::string_view value,
    double x, double y)
{
    if (property == nullptr) {
        property = symbol->add("property");
        property->setString(0, name);
        property->setString(1, value);
        property->add("at")->addNumber(x).addNumber(y).addNumber(0);
        property->add("effects")->add("hide")->addTag("yes");
    } else {
        property->setString(1, value);
    }
}



struct Symbol {
    Symbol(std::string_view reference) : reference(reference) {}

    std::string reference;

    double minVoltage = std::numeric_limits<double>::max();
    double maxVoltage = -std::numeric_limits<double>::max();
    bool fixedVoltage = false;

    std::set<int> nets;

    bool setVoltage(double minVoltage, double maxVoltage) {
        if (this->fixedVoltage)
            return false;

        bool changed = false;
        if (minVoltage < this->minVoltage) {
            this->minVoltage = minVoltage;
            changed = true;
        }
        if (maxVoltage > this->maxVoltage) {
            this->maxVoltage = maxVoltage;
            changed = true;
        }
        return changed;
    }

    bool setVoltage(double voltage) {
        return setVoltage(voltage, voltage);
    }

    bool voltageValid() {
        return this->minVoltage != std::numeric_limits<double>::max()
            && this->maxVoltage != -std::numeric_limits<double>::max();
    }

    double voltage() {
        return this->maxVoltage - this->minVoltage;
    }
};


// pcb

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

using double2 = Vector2<double>;

double distance(double2 a, double2 b) {
    double2 d = a - b;
    return sqrt(d.x * d.x + d.y * d.y);
}


struct Orientation {
    double2 position;
    double rotation;

    /// @brief Transform a point from local to global coordinates
    /// @param p Point
    /// @return Point in global coordinates
    double2 toGlobal(double2 p) {
        double r = this->rotation * pi / 180.0;
        double s = sin(r);
        double c = cos(r);

        return {this->position.x + c * p.x + s * p.y, this->position.y + c * p.y - s * p.x};
    }

    /// @brief Transform a point from global to local coordinates
    /// @param p Point
    /// @return Point in local coordinates
    double2 toLocal(double2 p) {
        double r = -this->rotation * pi / 180.0;
        double s = sin(r);
        double c = cos(r);

        return {(p.x - this->position.x) * c + (p.y - this->position.y) * s,
            (p.y - this->position.y) * c - (p.x - this->position.x) * s};
    }

    /// @brief Transform an orientation from local to global coordinates
    /// @param o Orientation
    /// @return Orientation in global coordinates
    Orientation toGlobal(Orientation o) {
        double r = this->rotation * pi / 180.0;
        double s = sin(r);
        double c = cos(r);

        return {{this->position.x + c * o.position.x + s * o.position.y,
            this->position.y + c * o.position.y - s * o.position.x},
            this->rotation + o.rotation};
    }

    /// @brief Transform an orientation from global to local coordinates
    /// @param o Orientation
    /// @return Orientation in local coordinates
    Orientation toLocal(Orientation o) {
        double r = -this->rotation * pi / 180.0;
        double s = sin(r);
        double c = cos(r);

        return {{(o.position.x - this->position.x) * c + (o.position.y - this->position.y) * s,
            (o.position.y - this->position.y) * c - (o.position.x - this->position.x) * s},
            o.rotation - this->rotation};
    }
};

struct Pad {
    std::string name;
    kicad::Container *at;
    double2 position;
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
        this->at->setNumber(0, o.position.x);
        this->at->setNumber(1, o.position.y);
        this->at->setNumber(2, o.rotation);
        for (auto &pad : pads) {
            pad.at->setNumber(2, o.rotation);
        }
    }
};

struct Template {
    std::list<Footprint> refFootprints;
    std::list<Footprint> placementFootprints;
};



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
                double2 position = {0, 0};
                double rotation = 0;
                std::string reference;
                std::string value;
                std::vector<Pad> pads;
                for (auto element2 : footprint->elements) {
                    auto container2 = dynamic_cast<kicad::Container *>(element2);
                    if (container2) {
                        if (container2->id == "at") {
                            at = container2;
                            position.x = at->getNumber(0);
                            position.y = at->getNumber(1);
                            rotation = at->getNumber(2);
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
                                auto x = at->getNumber(0);
                                auto y = at->getNumber(1);

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
    // schematic options
    std::map<std::string, std::string> exchangeFootprint;
    std::list<NetVoltage> netVoltages;
    std::map<std::string, std::pair<double, double>> symbolVoltages;
    std::map<std::string, Properties> importedBom;
    bool annotateSchematic = false;

    // pcb options
    std::map<int, double> exchangeSegmentWidth;
    std::map<int, double> exchangeViaSize;
    std::list<std::regex> hideReference;
    std::list<std::regex> removeNet;
    std::list<std::regex> moveFootprint;

    // schematic and pcb options
    std::list<std::regex> removeProperty;

    std::list<Template> templates;

    // paths of files to process
    std::list<fs::path> paths;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "-xfp") {
            // exchange footprint
            i += 2;
            exchangeFootprint[argv[i - 1]] = argv[i];
        } else if (arg == "-nv") {
            // define net voltage
            i += 2;
            //netVoltages[argv[i - 1]] = std::stod(argv[i]);
            netVoltages.emplace_back(std::regex(argv[i - 1]), std::stod(argv[i]));
        } else if (arg == "-sv") {
            // define symbol voltage
            i += 3;
            symbolVoltages[argv[i - 2]] = {std::stod(argv[i - 1]), std::stod(argv[i])};
        } else if (arg == "-ijb") {
            // import JLCPCB BOM
            ++i;
            fs::path xlsPath = argv[i];
            xlnt::workbook wb;
            wb.load(xlsPath);
            auto sheet = wb.active_sheet();
            for (uint32_t row = 6; row <= sheet.highest_row(); ++row) {
                auto references = getCellValue(sheet, row, 1);
                auto manufacturer = getCellValue(sheet, row, 8);
                auto mpn = getCellValue(sheet, row, 7);
                auto description = getCellValue(sheet, row, 10);
                auto datasheet = getCellValue(sheet, row, 14);
                auto lcscPn = getCellValue(sheet, row, 13);

                // split references (designators) by comma
                std::stringstream ss(references);
                while (ss.good()) {
                    std::string reference;
                    std::getline(ss, reference, ',');
                    importedBom[reference] = {manufacturer, mpn, description, datasheet, lcscPn};
                }
            }
        } else if (arg == "-as") {
            // annotate schematic
            annotateSchematic = true;
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
            // move footprint to top-left corner
            ++i;
            moveFootprint.emplace_back(argv[i]);
        } else if (arg == "-rp") {
            // remove property
            ++i;
            removeProperty.emplace_back(argv[i]);
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

    // project
    std::map<std::string, std::string> sheetUuidMap;

    // voltage propagation
    bool defineVoltages = !netVoltages.empty();
    std::list<Symbol> symbols;
    std::map<int, std::set<Symbol *>> netToSymbols;
    std::map<std::string, Symbol *> referenceToSymbol;

    bool hasExchangeSegmentWidth = !exchangeSegmentWidth.empty();
    bool hasExchangeViaSize = !exchangeViaSize.empty();

    if (paths.empty()) {
        std::cout << "Usage: pcb-tool [options] schematic.kicad_sch layout.kicad_pcb ..." << std::endl;
        std::cout << "    -xfp <old footprint> <new footprint>           Exchange footprint (only schematic)" << std::endl;
        std::cout << "    -ijb <bom file>                                Import downloaded JLCPCB BOM Excel file (only schematic)" << std::endl;
        std::cout << "    -as                                            Annotate schematic using sheet number as first number (only schematic)" << std::endl;
        std::cout << "    -xsw <old segment width> <new segment width>   Exchange segment width (only pcb)" << std::endl;
        std::cout << "    -xvs <old via size> <new via size>             Exchange segment width (only pcb)" << std::endl;
        std::cout << "    -hr <reference>                                Hide reference, supports regular expressions, e.g. R.* (only pcb)" << std::endl;
        std::cout << "    -rn <net>                                      Remove net, supports regular expressions, e.g. GND (only pcb)" << std::endl;
        std::cout << "    -mfp <footprint>                               Move footprint to border (only pcb)" << std::endl;
        return 1;
    }

    for (auto &path : paths) {
        auto ext = path.extension();
        bool isProject = ext == ".kicad_pro";
        bool isSchematic = ext == ".kicad_sch";
        bool isPcb = ext == ".kicad_pcb";

        if (!isProject && !isSchematic && !isPcb) {
            std::cout << "Error: File type of " << path.string() << " not supported" << std::endl;
            return 1;
        }

        std::cout << "Read " << path.string() << std::endl;
        if (isProject) {
            // read project
            std::ifstream is(path.string());
            if (is.is_open()) {
                try {
                    json j = json::parse(is);
                    json sheets = j.at("sheets");
                    for (int i = 0; i < sheets.size(); ++i) {
                        json sheet = sheets.at(i);
                        std::string uuid = sheet.at(0).get<std::string>();
                        std::string name = sheet.at(1).get<std::string>();
                        sheetUuidMap[name] = uuid;
                    }
                } catch (std::exception &e) {
                }
            }
        } else if (isSchematic) {
            std::string projectName = path.stem().string();
            std::map<int, Sheet> sheets;

            // read root sheet
            std::ifstream in(path.string());
            if (!in) {
                // error
                std::cout << "Can't read file " << path.string() << std::endl;
                return 1;
            }
            auto &rootSheet = sheets.emplace(-1, path).first->second;
            kicad::readFile(in, rootSheet.file);
            in.close();

            // set UUID path of root sheet
            rootSheet.uuidPath = '/' + sheetUuidMap["Root"];

            // load other sheets
            int rootPage =  -1;
            for (auto container1 : rootSheet.file) {
                if (container1->id == "sheet_instances") {
                    auto path = container1->find("path");
                    if (path != nullptr) {
                        auto page = path->find("page");
                        if (page != nullptr) {
                            rootPage = page->getInt(0, -1);
                        }
                    }
                } else if (container1->id == "sheet") {
                    // sheet
                    std::string sheetName;
                    fs::path sheetPath;
                    int sheetPage = -1;
                    auto properties = container1;
                    for (auto property : *properties) {
                        if (property->id == "property") {
                            auto propertyName = property->getString(0);
                            if (propertyName == "Sheetname") {
                                sheetName = property->getString(1);
                            } else if (propertyName == "Sheetfile") {
                                sheetPath = path.parent_path() / property->getString(1);
                            }
                        } else if (property->id == "instances") {
                            auto project = property->find("project");
                            if (project != nullptr) {
                                auto path = project->find("path");
                                if (path != nullptr) {
                                    auto page = path->find("page");
                                    if (page != nullptr) {
                                        sheetPage = page->getInt(0, -1);
                                    }
                                }
                            }
                        }
                    }

                    // read sheet
                    if (!sheetName.empty() && !sheetPath.empty() && sheetPage >= 0) {
                        std::cout << "Read sheet " << sheetPath.string() << " page " << sheetPage << std::endl;
                        std::ifstream in(sheetPath.string());
                        if (!in) {
                            // error
                            std::cout << "Can't read sheet " << sheetPath.string() << std::endl;
                            return 1;
                        }
                        auto &sheet = sheets.emplace(sheetPage, sheetPath).first->second;
                        kicad::readFile(in, sheet.file);
                        in.close();

                        // set UUID path of sheet
                        sheet.uuidPath = rootSheet.uuidPath + '/' + sheetUuidMap[sheetName];
                    }
                }
            }

            // assign root page number
            auto node = sheets.extract(-1);
            node.key() = rootPage;
            sheets.insert(std::move(node));

            // varialbes for schematic annotation
            std::map<std::string, std::string> oldToNewReference;
            std::map<std::string, int> nextReference;

            // process each sheet
            for (auto &p : sheets) {
                int sheetPage = p.first;
                auto &sheet = p.second;

                // symbol references sorted by y for schematic annotation
                std::multimap<std::pair<double, double>, std::pair<kicad::Value *, kicad::Value *>> sortedReferences;

                for (auto container1 : sheet.file) {
                    // check if it is a symbol
                    if (container1->id == "symbol") {
                        auto symbol = container1;

                        // iterate over properties
                        std::string reference;
                        int unit = -1;
                        double x = 0;
                        double y = 0;
                        kicad::Value *ref1 = nullptr;
                        kicad::Value *ref2 = nullptr;
                        kicad::Container *manufacturer = nullptr;
                        kicad::Container *mpn = nullptr;
                        kicad::Container *description = nullptr;
                        kicad::Container *datasheet = nullptr;
                        kicad::Container *lcscPn = nullptr;
                        kicad::Container *voltage = nullptr;
                        for (auto property : *symbol) {
                            if (property->id == "at") {
                                // position
                                x = property->getNumber(0);
                                y = property->getNumber(1);
                            } else if (property->id == "unit") {
                                // unit index (parts with multiple units such as op amps or gates)
                                if (unit == -1)
                                    unit = property->getInt(0);
                            } else if (property->id == "property") {
                                // property
                                auto propertyName = property->getString(0);

                                if (propertyName == "Reference") {
                                    // reference (e.g. "U1")
                                    if (reference.empty())
                                        reference = property->getString(1);
                                    ref1 = dynamic_cast<kicad::Value *>(property->elements[1]);
                                } else if (propertyName == "Footprint") {
                                    // replace footprint
                                    auto footprint = property;
                                    auto footprintName = footprint->getString(1);
                                    auto it = exchangeFootprint.find(footprintName);
                                    if (it != exchangeFootprint.end()) {
                                        std::cout << "Replace footprint " << it->first << " by " << it->second << std::endl;
                                        footprint->setString(1, it->second);
                                    }
                                } else if (propertyName == "Description") {
                                    description = property;
                                } else if (propertyName == "Datasheet") {
                                    datasheet = property;
                                } else {
                                    // check if property should be removed (only allowed for the following custom properties)
                                    if (match(propertyName, removeProperty)) {
                                        //next = symbol->elements.erase(it);
                                        property->action = kicad::Element::Action::DELETE;
                                    } else {
                                        if (propertyName == "Manufacturer") {
                                            manufacturer = property;
                                        } else if (propertyName == "MPN") {
                                            mpn = property;
                                        } else if (propertyName == "LCSC PN") {
                                            lcscPn = property;
                                        } else if (propertyName == "Voltage") {
                                            voltage = property;
                                        }
                                    }
                                }
                            } else if (property->id == "instances") {
                                for (auto project : *property) {
                                    // remove unused instances
                                    //if (project->action == kicad::Element::Action::NONE)
                                    //    project->action = kicad::Element::Action::DELETE;

                                    // check if UUID path of instance matches
                                    auto path = project->find("path");
                                    if (path != nullptr && path->getString(0) == sheet.uuidPath) {
                                        auto r = path->find("reference");
                                        if (r != nullptr) {
                                            // set project name
                                            project->setString(0, projectName);

                                            // get reference
                                            reference = r->getString(0);
                                            std::cout << "Instance " << reference << std::endl;
                                            ref2 = dynamic_cast<kicad::Value *>(r->elements[0]);

                                            // keep instance
                                            project->action = kicad::Element::Action::KEEP;
                                        }

                                        // get unit index
                                        auto u = path->find("unit");
                                        if (u != nullptr) {
                                            unit = u->getInt(0);
                                        }
                                    }
                                }
                            }
                        }

                        if (!reference.empty() && reference[0] != '#') {
                            // add properties from imported bom
                            {
                                auto it = importedBom.find(reference);
                                if (it != importedBom.end()) {
                                    addProperty(symbol, manufacturer, "Manufacturer", it->second.manufacturer, x, y);
                                    addProperty(symbol, mpn, "MPN", it->second.mpn, x, y);
                                    addProperty(symbol, description, "Description", it->second.description, x, y);
                                    addProperty(symbol, datasheet, "Datasheet", it->second.datasheet, x, y);
                                    addProperty(symbol, lcscPn, "LCSC PN", it->second.lcscPn, x, y);
                                }
                            }

                            // add propagated voltage (propagation is done using pcb's)
                            {
                                auto it = referenceToSymbol.find(reference);
                                if (it != referenceToSymbol.end()) {
                                    auto sym = it->second;
                                    if (sym->voltageValid()) {
                                        std::stringstream ss;
                                        ss << sym->voltage();
                                        addProperty(symbol, voltage, "Voltage", ss.str(), x, y);
                                    }
                                }
                            }

                            if (annotateSchematic) {
                                sortedReferences.insert({{y, x}, {ref1, ref2}});
                            }
                        }
                    }
                }

                // delete all elements marked as DELETE
                sheet.file.sweep();

                // annotate schematic (assign references)
                if (annotateSchematic) {
                    // update next reference indices
                    for (auto &p : nextReference) {
                        if (p.second < sheetPage) {
                            p.second = sheetPage;
                        }
                    }

                    for (auto &p : sortedReferences) {
                        auto ref1 = p.second.first;
                        auto ref2 = p.second.second;
                        auto reference = (ref2 != nullptr ? ref2 : ref1)->getString();

                        // check if new reference already exists
                        if (!oldToNewReference.contains(reference)) {
                            // create new reference
                            auto type = getType(reference);
                            int index = sheetPage;
                            if (nextReference.contains(type)) {
                                index = nextReference[type];
                            }

                            std::string newReference = type + std::to_string(index);
                            oldToNewReference[reference] = newReference;

                            nextReference[type] = index + 1;
                        }

                        if (ref1 != nullptr)
                            ref1->setString(oldToNewReference[reference]);
                        if (ref2 != nullptr)
                            ref2->setString(oldToNewReference[reference]);
                    }
                }
            }

            // print mapping from old to new reference
            for (auto &p : oldToNewReference) {
                std::cout << p.first << " -> " << p.second << std::endl;
            }

            // write sheets
            for (auto &p : sheets) {
                auto &sheet = p.second;
                std::cout << "Write " << sheet.path.string() << std::endl;
                std::ofstream out(sheet.path.string());
                if (!out) {
                    // error
                    return 1;
                }
                kicad::writeFile(out, sheet.file);
                out.close();
            }
        }

        if (isPcb) {
            // read pcb file
            std::ifstream in(path.string());
            if (!in) {
                // error
                std::cout << "Can't read file " << path.string() << std::endl;
                return 1;
            }
            kicad::Container file;
            kicad::readFile(in, file);
            in.close();

            symbols.clear();
            netToSymbols.clear();

            double2 moveFootprintPosition = {};
            for (auto container1 : file) {
                // check if it is a footprint
                if (container1->id == "footprint") {
                    auto footprint = container1;

                    // get footprint name
                    auto footprintName = footprint->getString(0);
                    //std::cout << "Footprint: " << footprintName << std::endl;

                    // iterate over properties of footprint
                    double2 position = {0, 0};
                    double rotation = 0;
                    std::string reference;
                    std::string value;
                    bool excludeFromBom = false;
                    for (auto element2 : footprint->elements) {
                        auto property = dynamic_cast<kicad::Container *>(element2);
                        if (property) {
                            if (property->id == "at") {
                                auto at = property;

                                // move footprint away
                                if (match(footprintName, moveFootprint)) {
                                    at->setNumber(0, moveFootprintPosition.x);
                                    at->setNumber(1, moveFootprintPosition.y);
                                    //at->setNumber(2, 0);
                                    moveFootprintPosition.y += 2;
                                }

                                position.x = at->getNumber(0);
                                position.y = at->getNumber(1);
                                rotation = at->getNumber(2);
                            } else if (property->id == "property") {
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
                                        property->add("hide")->addTag("yes");
                                    }
                                } else if (propertyName == "Value") {
                                    // get reference and hide if requested
                                    value = property->getString(1);
                                }
                            } else if (property->id == "attr") {
                                //doNotPopulate = property->contains("dnp");
                                excludeFromBom = property->contains("exclude_from_bom");
                                //throughHole = property->contains("through_hole");
                            }

                        }
                    }

                    // define voltage
                    Symbol *symbol = nullptr;
                    if (defineVoltages && !excludeFromBom) {
                        symbol = &symbols.emplace_back(reference);
                        if (symbolVoltages.contains(reference)) {
                            // set voltage if it is defined for the symbol
                            auto &voltage = symbolVoltages[reference];
                            symbol->setVoltage(voltage.first, voltage.second);
                            symbol->fixedVoltage = true;
                        }
                        referenceToSymbol[reference] = symbol;
                    }

                    // iterate over pads of footprint
                    for (auto container2 : *footprint) {
                        if (container2->id == "pad") {
                            auto pad = container2;

                            // get pad name (typcally a number)
                            //std::string padName = clean(pad->getValue(0));

                            auto net = pad->find("net");
                            if (net != nullptr) {
                                auto netNumber = net->getInt(0);
                                auto netName = net->getString(1);

                                // define voltage
                                if (symbol != nullptr) {
                                    auto v = match(netName, netVoltages);
                                    if (v) {
                                        // set voltage if it is defined for the net
                                        symbol->setVoltage(*v);
                                    } else {
                                        // add net unless it is an input
                                        auto pinType = pad->findString("pintype");
                                        if (pinType != "input") {
                                            symbol->nets.insert(netNumber);
                                        }
                                        netToSymbols[netNumber].insert(symbol);
                                    }
                                }

                                // remove net from pad if requested
                                if (match(netName, removeNet)) {
                                    // erase net
                                    std::cout << "Erase net " << netName << " of " << reference << std::endl;
                                    pad->erase(net);
                                }
                            }
                        }
                    }
                } else if (hasExchangeSegmentWidth && container1->id == "segment") {
                    auto segment = container1;
                    auto width = segment->find("width");
                    if (width != nullptr) {
                        try {
                            double oldWidth = width->getNumber(0);
                            auto it = exchangeSegmentWidth.find(int(std::round(oldWidth * 1000.0)));
                            if (it != exchangeSegmentWidth.end()) {
                                width->setNumber(0, it->second);
                            }
                        } catch (std::exception &) {
                        }
                    }
                } else if (hasExchangeViaSize && container1->id == "via") {
                    auto via = container1;
                    auto size = via->find("size");
                    if (size != nullptr) {
                        try {
                            double oldSize = size->getNumber(0);
                            auto it = exchangeViaSize.find(int(std::round(oldSize * 1000.0)));
                            if (it != exchangeViaSize.end()) {
                                size->setNumber(0, it->second);
                            }
                        } catch (std::exception &) {
                        }
                    }
                }
            }

            // define voltages
            bool changed;
            do {
                changed = false;
                for (auto &symbol : symbols) {
                    for (auto net : symbol.nets) {
                        auto &symbols2 = netToSymbols[net];
                        for (auto symbol2 : symbols2) {
                            double oldMin = symbol2->minVoltage;
                            double oldMax = symbol2->maxVoltage;
                            bool valid = symbol2->voltageValid();
                            bool c = symbol2->setVoltage(symbol.minVoltage, symbol.maxVoltage);
                            if (c) {
                                std::cout << symbol.reference << " changes voltage of " << symbol2->reference;
                                if (valid)
                                    std::cout << " from [" << oldMin << ' ' << oldMax << "]";
                                std::cout << " to [" << symbol2->minVoltage << ' ' << symbol2->maxVoltage << ']' << std::endl;
                            }

                            changed |= c;
                        }
                    }
                }
            } while (changed);

            // automatic placement of resistors and capacitors
            if (!templates.empty()) {
                // collection of footprints for processing
                std::list<Footprint> refFootprints;
                std::list<Footprint> placementFootprints;
                getFootprints(file, refFootprints, placementFootprints);

                for (auto &refFootprint : refFootprints) {

                    // generate net positions and places
                    std::map<int, std::vector<double2>> netPositions;
                    std::list<Orientation> places;
                    /*for (auto &refFootprint : refFootprints)*/ {
                        // add net positions
                        for (auto &pad : refFootprint.pads) {
                            double2 p = refFootprint.orientation.toGlobal(pad.position);
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
                        double bestDistance = std::numeric_limits<double>::max();
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
                                    double d = footprint.reference.starts_with("C") ? -1.0f : 0.0f;

                                    // calculate distance sum
                                    for (auto &pad : footprint.pads) {
                                        auto orientation = place;
                                        orientation.rotation += flip ? 180.0f : 0.0f;
                                        auto padPosition = orientation.toGlobal(pad.position);

                                        double padDistance = 100.0f;
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
                        //bestFootprint->at->setNumber(0, orientation.position.x);
                        //bestFootprint->at->setNumber(1, orientation.position.y);
                        //bestFootprint->at->setNumber(2, orientation.rotation);

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

            // write pcb file
            std::cout << "Write " << path.string() << std::endl;
            std::ofstream out(path.string());
            if (!out) {
                // error
                return 1;
            }
            kicad::writeFile(out, file);
            out.close();
        }
    }
    return 0;
}