#ifndef LOGGER_H
#define LOGGER_H
// adapted from: https://cppcodetips.wordpress.com/2014/01/02/a-simple-logger-class-in-c/

#include <fstream>
#include <iostream>
#include <cstdarg>
#include <string>
#include <time.h>
#include "utils.h"

#define LOG(m) Logger::GetLogger()->log(m)
#define LOGF Logger::m_logfile

namespace sdfibm{


class Logger
{
public:
    void log(const std::string& msg)
    {
        m_logfile << '[' << GetTimeString() << "]\t" << msg << '\n';
    }
    static Logger* GetLogger()
    {
        if (m_this == nullptr)
        {
            m_this = new Logger();
            m_logfile.open(m_logfilename, std::ios::out);
        }
        return m_this;
    }
    ~Logger()
    {
        m_logfile.close();
    }
    static const std::string m_logfilename;
    static std::ofstream m_logfile;

private:
    Logger(){}
    Logger(const Logger&){}
    Logger& operator=(const Logger&){return *this;}

    static Logger* m_this;
};

}
#endif