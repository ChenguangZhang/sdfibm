#include "logger.h"
#include <fstream>

namespace sdfibm {

Logger*  Logger::m_this = nullptr;
const std::string Logger::m_logfilename = "cloud.log";
std::ofstream Logger::m_logfile;

}
