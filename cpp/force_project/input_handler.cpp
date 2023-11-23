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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <mutex>
#include <chrono>
#include <map>

#include "../common/logger.hpp"
#include "input_handler.hpp"

void checkChildExist(ryml::ConstNodeRef t, std::string str) {
    bool has_child = t.has_child(c4::to_csubstr(str));
    if (!has_child) {
        if (t.has_key()) {
            LOG_CRITICAL << "The key [ " << std::string(t.key().str, t.key().len) << " ] doesnt have a child [ " << str
                         << " ]" << logger::endl;
        } else {
            LOG_CRITICAL << "The root node doesnt have a child [ " << str << " ]" << logger::endl;
        }
        exit(EXIT_FAILURE);
    }
}

YamlReaderNode::YamlReaderNode(ryml::ConstNodeRef t, std::string str) {
    checkChildExist(t, str);
    yaml_node = t[c4::to_csubstr(str)];
}

YamlReaderNode YamlReaderNode::operator[](std::string str) { return YamlReaderNode(yaml_node, str); }

YamlReaderNode YamlReader::operator[](std::string s) const { return YamlReaderNode(tree_, s); }
