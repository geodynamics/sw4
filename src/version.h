#ifndef EW_VERSION_H
#define EW_VERSION_H

#include <string>
#include <iostream>

namespace ewversion {

extern const char* madeby;
extern const char* when;
extern const char* hostname;
extern const char* basedir;
extern const char* optimization;
extern const char* compiler;
extern const char* version;

std::string getVersionInfo();
};

#endif
