#ifndef UTILS_H
#define UTILS_H

#include <ctime>
#include <iomanip>
#include "types.h"

namespace sdfibm{

namespace AsciiColors{
    static const std::string COLOR_NORMAL ="\033[0m";
    static const std::string COLOR_RED    ="\033[0;31;49m";
    static const std::string COLOR_GREEN  ="\033[0;32;49m";
    static const std::string COLOR_YELLOW ="\033[0;33;49m";
    static const std::string COLOR_BLUE   ="\033[0;34;49m";
    static const std::string COLOR_MAGENTA="\033[0;35;49m";
    static const std::string COLOR_CYAN   ="\033[0;36;49m";
    static const std::string COLOR_WHITE  ="\033[0;37;49m";
    static const std::string COLOR_INFO   ="\033[0;32;49m";
    static const std::string COLOR_WARNING="\033[4;33;49m";
    static const std::string COLOR_ERROR  ="\033[1;31;43m";
}

inline void PrintInfo(const std::string& message)
{
    std::cout << AsciiColors::COLOR_INFO << message << AsciiColors::COLOR_NORMAL << std::endl;
}
inline void PrintWarning(const std::string& message)
{
    std::cout << AsciiColors::COLOR_WARNING << message << AsciiColors::COLOR_NORMAL << std::endl;
}
inline void PrintError(const std::string& message)
{
    std::cout << AsciiColors::COLOR_ERROR << message << AsciiColors::COLOR_NORMAL << std::endl;
}

inline void Quit(const std::string& msg)
{
    std::cout << msg << std::endl;
    std::exit(1);
}

// below requires rather new c++ version
inline std::string GetTimeStringNew()
{
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::stringstream ss;
    ss << std::put_time(&tm, "%h %d %Y %H:%M:%S");
    return ss.str();
}

inline std::string GetTimeString()
{
    char buf[80];
    time_t t = time(NULL);
    struct tm *tstruct = localtime(&t);
    strftime(buf, sizeof(buf), "%y/%m/%d %X", tstruct);
    return buf;
}

}
#endif
