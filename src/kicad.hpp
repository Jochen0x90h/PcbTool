#pragma once

#include <istream>
#include <ostream>
#include <list>
#include <map>
#include <string>
#include <vector>


/// @brief classes and functions for reading and writing kicad files
///
namespace kicad {

std::string toString(std::string_view value);


class Element {
public:

    virtual ~Element() {}

    virtual int count() = 0;

    virtual void write(std::ostream &s, int indent) = 0;
};


class Value : public Element {
public:
    Value() = default;
    Value(std::string_view value) : value(value) {}

    virtual ~Value();

    int count() override;

    void write(std::ostream &s, int indent) override;

    std::string value;
};


class Container : public Element {
public:
    template <typename T>
    struct Value2 {
        T x;
        T y;
    };


    Container() = default;
    Container(std::string_view id) : id(id) {}

    virtual ~Container();

    void clear();

    int count() override;

    void add(Element *element) {this->elements.push_back(element);}


    Container &addValue(std::string_view value);

    Container &addValue(int value) {
        return addValue(std::to_string(value));
    }

    Container &addValue(double value) {
        return addValue(std::to_string(value));
    }

    template <typename T>
    Container &addValues(const T &value) {
        return addValue(value);
    }

    template <typename T, typename ...Args>
    Container &addValues(const T &value, Args ...args) {
        addValue(value);
        return addValues(args...);
    }

    /// @brief Add a new container
    /// @param id id of container
    /// @return the new container
    Container *add(std::string_view id);

    template <typename T>
    Container *add(std::string_view id, const T &value) {
        auto container = add(id);
        container->addValue(value);
        return container;
    }

    template <typename T, typename ...Args>
    Container *add(std::string_view id, const T &value, Args ...args) {
        auto container = add(id);
        container->addValue(value).addValues(args...);
        return container;
    }

    Value *getValue(int index) {
        // check index
        if (unsigned(index) >= this->elements.size())
            return nullptr;

        // check if element is of type Value
        return dynamic_cast<Value *>(this->elements[index]);
    }

    /// @brief Get a tag at given index (string without quotes ).
    /// @param index Index of the element to get
    /// @param defaultValue Value to return if index is out of bounds or element is not of type Value
    /// @return tag at given index
    std::string getTag(int index, const std::string &defaultValue = {});

    /// @brief Get a string at given index, removing quotes if there are any.
    /// @param index Index of the element to get
    /// @param defaultValue Value to return if index is out of bounds or element is not of type Value
    /// @return string value at given index
    std::string getString(int index, const std::string &defaultValue = {});

    /// @brief Get an integer at given index, returning defaultValue if index is out of bounds or element is not of type Value.
    /// @param index Index of the element to get
    /// @param defaultValue Value to return if index is out of bounds or element is not of type Value
    /// @return integer value at given index
    int getInt(int index, int defaultValue = 0);

    /// @brief Get a double at given index, returning defaultValue if index is out of bounds or element is not of type Value.
    /// @param index Index of the element to get
    /// @param defaultValue Value to return if index is out of bounds or element is not of type Value
    /// @return double value at given index
    double getDouble(int index, double defaultValue = 0.0);
    float getFloat(int index, float defaultValue = 0.0f) {return float(getDouble(index, defaultValue));}


    /// @brief Set a double at a given index
    /// @param index Index of the element to set
    /// @param value Value to set
    void setDouble(int index, double value);
    void setFloat(int index, float value) {setDouble(index, value);}


    /// @brief Check if the container contains a value
    /// @param value to find
    /// @return true if value is an element of the container
    bool contains(const std::string &value);

    /// @brief Find element container with given id
    /// @param id id of sub-container to find
    /// @return container or nullptr if not found or not of type Container
    Container *find(const std::string &id);

    /// @brief Find element container with given id and return its first value (of type kicad::Value)
    /// @param id id of sub-container to find
    /// @return value or empty string if not found
    std::string findString(const std::string &id);

    /// @brief Find element container with given id and return its first and second value (of type kicad::Value)
    /// @param id id of sub-container to find
    /// @return values or empty strings if not found
    Value2<std::string> findString2(const std::string &id);

    Value2<float> findFloat2(const std::string &id);

    /// @brief Erase element
    /// @param element Element to erase
    void erase(Element *element);


    void write(std::ostream &s, int indent) override;
    static void newLine(std::ostream &s, int indent);


    std::string id;
    std::vector<Element *> elements;
};


/// @brief Read a kicad file
/// @param buffer buffer of an open file or network socket that is in ready state
void readFile(std::istream &s, Container &kicad);

/// @brief Write a kicad file
/// @param buffer buffer of an open file or network socket that is in ready state
inline void writeFile(std::ostream &s, Container &kicad) {
    kicad.write(s, 0);
}

} // namespace kicad
