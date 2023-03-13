
#ifndef CUTFEM_COMMON_LOGGER_HPP
#define CUTFEM_COMMON_LOGGER_HPP

enum class Severity {
    Trace,
    Debug,
    Info,
    Warning,
    Error,
    Critical,
};

enum class LogOutput {
    Console,
    File,
    Everywhere,
};
constexpr const char *capitalize(Severity severity) noexcept {
    switch (severity) {
    case Severity::Trace:
        return "TRACE";
    case Severity::Debug:
        return "DEBUG";
    case Severity::Info:
        return "INFO";
    case Severity::Warning:
        return "WARNING";
    case Severity::Error:
        return "ERROR";
    case Severity::Critical:
        return "CRITICAL";
    default:
        return "UNKNOWN";
    }
}
namespace logger {
const std::string endl = "\n";
}

class Logger {
  private:
    std::mutex m_logMutex;
    std::mutex m_threadMsgMutex;
    const char *filename_separator = "/\\";

  public:
    Logger() {}

    ~Logger() { close(); }

  public:
    static Logger &get() {
        static Logger self;
        return self;
    }

    Logger(const Logger &)            = delete;
    Logger &operator=(const Logger &) = delete;

    void init(std::string &filename, const LogOutput output = LogOutput::Everywhere) {
        m_filename = filename;
        m_output   = output;
        open();
    }

    static void initialize(std::string filename, const LogOutput output = LogOutput::Everywhere) {

        Logger &logger_instance    = get();
        logger_instance.m_filename = filename;
        logger_instance.m_output   = output;
        if (output == LogOutput::Everywhere || output == LogOutput::File) {
            try {
                logger_instance.open();
            } catch (const std::exception &exc) {
                std::cout << "[logger.hpp] - " << exc.what() << " Logging to console only." << std::endl;
                logger_instance.m_output = LogOutput::Console;
            }
        }
    }

    void open() {

        if (m_streamOpen) {
            return;
        }

        if (m_output == LogOutput::Everywhere || m_output == LogOutput::File) {
            m_file.open(m_filename, std::ios::out);
            m_streamOpen = m_file.is_open();

            if (!m_streamOpen) {
                throw std::runtime_error("Couldn't open log file.");
            }
        }
    }

    void close() {
        if (!m_streamOpen) {
            return;
        }

        if (m_output == LogOutput::Everywhere || m_output == LogOutput::File) {
            m_file.close();
        }
    }

    void flush() {
        if (!m_streamOpen) {
            return;
        } else {
            std::unique_lock<std::mutex> locker(m_threadMsgMutex);
            for (const auto &elem : m_threadsMessageMap) {
                log(elem.second.rdbuf(), m_threadsSeverityMap[elem.first]);
            }
            m_threadsMessageMap.clear();
            m_threadsSeverityMap.clear();
        }

        if (m_output == LogOutput::Everywhere || m_output == LogOutput::File) {
            m_file.flush();
        }
    }

    void log(const std::string &msg, const Severity severity = Severity::Debug) {

        std::string output_message =
            wrapString(timestamp()) + wrapString("#" + threadId()) + " " + wrapString(capitalize(severity)) + " " + msg;

        std::unique_lock<std::mutex> locker(m_logMutex);

        if (m_output == LogOutput::Console) {
            std::cout << msg << std::endl;
        }

        else if (m_output == LogOutput::File) {
            m_file << output_message << std::endl;
        }

        else {
            std::cout << msg << std::endl;
            m_file << output_message << std::endl;
        }
    }

    void log(int line, std::string source_file, const std::string &msg, const Severity severity) {

        std::size_t found    = source_file.find_last_of(filename_separator);
        std::string location = source_file.substr(found + 1) + "@" + std::to_string(line);

        std::string output_message = wrapString(timestamp()) + wrapString("#" + threadId()) + wrapString(location) +
                                     " " + wrapString(capitalize(severity)) + " " + msg;

        std::unique_lock<std::mutex> locker(m_logMutex);

        if (m_output == LogOutput::Console) {
            std::cout << msg << std::endl;
        }

        else if (m_output == LogOutput::File) {
            m_file << output_message << std::endl;
        }

        else {
            std::cout << msg << std::endl;
            m_file << output_message << std::endl;
        }
    }

    template <class T> void log(const T &t, const Severity severity = Severity::Debug) {
        std::stringstream ss;
        ss << t;
        log(ss.str(), severity);
    }

    template <class T>
    void log(const int line, const std::string source_file, const T &t, const Severity severity = Severity::Debug) {
        std::stringstream ss;
        ss << t;
        log(line, source_file, ss.str(), severity);
    }

    void addMessage(const std::string &msg) {
        std::unique_lock<std::mutex> locker(m_threadMsgMutex);
        m_threadsMessageMap[threadId()] << msg;
    }

    void flushMessage() {
        const std::string threadID = threadId();
        std::unique_lock<std::mutex> locker(m_threadMsgMutex);

        if (m_threadsFilesMap[threadID].empty()) {
            log(m_threadsMessageMap[threadID].rdbuf(), m_threadsSeverityMap[threadID]);
        } else {
            log(m_threadsLineMap[threadID], m_threadsFilesMap[threadID], m_threadsMessageMap[threadID].rdbuf(),
                m_threadsSeverityMap[threadID]);
            m_threadsLineMap.erase(threadID);
            m_threadsFilesMap.erase(threadID);
        };
        m_threadsMessageMap.erase(threadID);
        m_threadsSeverityMap.erase(threadID);
    }

    Logger &operator()(const Severity severity) {
        std::unique_lock<std::mutex> locker(m_threadMsgMutex);
        m_threadsSeverityMap[threadId()] = severity;
        m_severity                       = severity;
        return get();
    }

    Logger &operator()(int line, const char *source_file, const Severity severity) {
        std::unique_lock<std::mutex> locker(m_threadMsgMutex);
        m_threadsSeverityMap[threadId()] = severity;
        m_severity                       = severity;
        m_threadsLineMap[threadId()]     = line;
        m_threadsFilesMap[threadId()]    = source_file;
        return get();
    }

  private:
    std::string m_filename;
    std::fstream m_file;
    LogOutput m_output  = LogOutput::Everywhere;
    Severity m_severity = Severity::Info;

    std::map<std::string, Severity> m_threadsSeverityMap;
    std::map<std::string, std::stringstream> m_threadsMessageMap;
    std::map<std::string, int> m_threadsLineMap;
    std::map<std::string, std::string> m_threadsFilesMap;

    bool m_streamOpen = false;

    static std::string timestamp() {
        auto now  = std::chrono::system_clock::now();
        time_t tt = std::chrono::system_clock::to_time_t(now);
        auto ms   = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

        std::stringstream ss;
        ss << std::put_time(std::localtime(&tt), "%F %T");
        ss << "." << std::setfill('0') << std::setw(3) << ms.count();

        return ss.str();
    }

    static std::string threadId() {
        std::stringstream ss;
        std::thread::id id = std::this_thread::get_id();
        ss << std::this_thread::get_id();
        return ss.str();
    }

    static std::string wrapString(const std::string &wrap_str, const std::string &prefix = "[",
                                  const std::string &sufix = "]") {
        return prefix + wrap_str + sufix;
    }
};

Logger &operator<<(Logger &logger, const std::string &message);

template <class T> Logger &operator<<(Logger &logger, const T &t) {
    std::stringstream ss;
    ss << t;
    logger.addMessage(ss.str());
    return logger;
}

Logger &LOG(const Severity severity = Severity::Info);
Logger &LOG(int line, const char *source_file, const Severity severity);

#define LOG_TRACE (LOG(__LINE__, __FILE__, Severity::Trace))
#define LOG_DEBUG (LOG(__LINE__, __FILE__, Severity::Debug))
#define LOG_INFO (LOG(__LINE__, __FILE__, Severity::Info))
#define LOG_WARNING (LOG(__LINE__, __FILE__, Severity::Warning))
#define LOG_ERROR (LOG(__LINE__, __FILE__, Severity::Error))
#define LOG_CRITICAL (LOG(__LINE__, __FILE__, Severity::Critical))

// Below is to test tread saftey.
// Cannot be there otherwise one get duplicated symbol. Has to be in a cpp file
// void threadFunc() {
//   LOG().log("TreadFunc");
//   LOG(Severity::Error) << "Interesting log" << 11 << logger::endl;
// }

// void threadFunc2() {
//   static int thread_num = 0;
//   for (int i = 0; i < 10; i++) {
//     LOG() << "tread #" << ++thread_num;
//     LOG() << "iteration #" << i + 1 << logger::endl;
//     LOG_DEBUG << "Macro test with thread#" << ++thread_num << logger::endl;
//   }
// }

// void testLoggerTreadSaftey() {

//   Logger::initialize("testLoggerTreadSaftey.txt");
//   std::vector<std::thread> threads;

//   for (int i = 0; i < 10; i++) {
//     threads.emplace_back(std::thread(i % 2 ? threadFunc : threadFunc2));
//   }
//   threadFunc();

//   for (auto &th : threads) {
//     if (th.joinable()) {
//       th.join();
//     }
//   }
// }

#endif