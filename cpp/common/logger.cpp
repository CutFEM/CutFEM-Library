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
#include <thread>
#include <chrono>
#include <map>
#include "logger.hpp"

CutFEMLogger &operator<<(CutFEMLogger &logger, const std::string &message) {
    if (message == logger::endl) {
        logger.flushMessage();
    } else {
        logger.addMessage(message);
    }
    return logger;
}

CutFEMLogger &LOG(const Severity severity) { return CutFEMLogger::get()(severity); }

CutFEMLogger &LOG(int line, const char *source_file, const Severity severity) {
    return CutFEMLogger::get()(line, source_file, severity);
}
