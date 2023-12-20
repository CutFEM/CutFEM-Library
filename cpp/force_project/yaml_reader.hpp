/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/

#ifndef CUTFEM_FORCE_PROJECT_YAML_READER_HPP
#define CUTFEM_FORCE_PROJECT_YAML_READER_HPP

#include <ryml_all.hpp>

namespace ryaml {
using namespace std;

template <class CharContainer> size_t file_get_contents(const char *filename, CharContainer *v) {
    ::FILE *fp = ::fopen(filename, "rb");

    C4_CHECK_MSG(fp != nullptr, "could not open file");
    ::fseek(fp, 0, SEEK_END);
    long sz = ::ftell(fp);
    v->resize(static_cast<typename CharContainer::size_type>(sz));
    if (sz) {
        ::rewind(fp);
        size_t ret = ::fread(&(*v)[0], 1, v->size(), fp);
        C4_CHECK(ret == (size_t)sz);
    }
    ::fclose(fp);
    return v->size();
}

/// @brief Load a file from disk into an existing CharContainer
template <class CharContainer> CharContainer file_get_contents(const char *filename) {
    CharContainer cc;
    file_get_contents(filename, &cc);
    return cc;
}
} // namespace ryaml

/// @brief Check if a node of a yaml tree has a specific child
/// @details If the child does not exist, the program will stop and write the details in
/// the log file
/// @param t the node to check
/// @param str the name of the child
void checkChildExist(ryml::ConstNodeRef t, std::string str);

void checkHasVal(ryml::ConstNodeRef t);

/// @brief Class that encapsulate a node of a yaml tree.
/// @details This class is used to access the data of a yaml tree and to check if the
/// required path or names exists.
class YamlReaderNode {

  public:
    /// @brief Build a YamlReaderNode from a yaml constRefNode and a string
    /// @param t the constant reference to the yaml node
    /// @param s the name of the node to access
    YamlReaderNode(ryml::ConstNodeRef t, std::string s);

    /// @brief Return a YamlReaderNode of the child s of the current node, if it exist.
    /// Otherwise it will stop the program
    YamlReaderNode operator[](std::string s);
    ryml::ConstNodeRef operator[](int i) const { return yaml_node[i]; };

    /// @brief Deserialize a value from the yaml node
    /// @tparam T
    /// @param v
    template <typename T> void operator>>(T &v) { yaml_node >> v; }

    ryml::ConstNodeRef get() const { return yaml_node; }

    bool has_val() const { return yaml_node.has_val(); }

    template <typename T> std::string val() const;

    std::string key() const { return std::string(yaml_node.key().str, yaml_node.key().len); }

    size_t num_children() const { return yaml_node.num_children(); }

    bool has_child(std::string s) const;

  private:
    /// @brief main tree
    ryml::ConstNodeRef yaml_node;
};

/// @brief Class handling Yaml format
class YamlReader {
  protected:
    std::string filename_;
    ryml::Tree tree_;

  public:
    YamlReader(const std::string &f) : filename_(f) {

        std::string contents = ryaml::file_get_contents<std::string>(filename_.c_str());
        tree_                = ryml::parse_in_arena(ryml::to_csubstr(contents));

        LOG_INFO << " Open file " << filename_ << " for reading" << logger::endl;
    }
    void showData() const { std::cout << tree_ << std::endl; }

    const ryml::Tree &tree() const { return tree_; }

    YamlReaderNode operator[](std::string s) const;

    friend class FileReader;
};

#endif
