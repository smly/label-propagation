#include "utils.h"

namespace utils {
int stoi (const std::string str)
{
    int ret;
    std::stringstream ss;
    ss << str;
    ss >> ret;
    return ret;
}
double stod (const std::string str)
{
    double ret;
    std::stringstream ss;
    ss << str;
    ss >> ret;
    return ret;
}
}
