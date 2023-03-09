#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <sstream>
#include "logger.hpp"

Logger &operator<<(Logger &logger, const std::string &message) {
    if (message == logger::endl) {
        logger.flushMessage();
    } else {
        logger.addMessage(message);
    }
    return logger;
}

Logger &LOG(const Severity severity) { return Logger::get()(severity); }
Logger &LOG(int line, const char *source_file, const Severity severity) {
    return Logger::get()(line, source_file, severity);
}
