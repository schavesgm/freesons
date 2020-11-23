#pragma once

#include <algorithm>
#include <sstream>

namespace Argparser {

    // Check if a flag is present in the command line
    bool Check_Flag(char** beg, char** end, const std::string& flag) {
        return std::find(beg, end, flag) != end;
    }

    // Get the option after some given flag as a string
    char* Get_Option(char** beg, char** end, const std::string& flag)
    {
        // Look for the flag inside begin and end pointers
        char** flag_pos = std::find(beg, end, flag);

        // Dereference the next memory address if its not the end
        if (flag_pos != end && ++flag_pos != end) {
            return *flag_pos;
        }
        
        // Return nullptr if nothing was found
        return nullptr;
    }


    // Cast a char pointer to a given type
    template <typename T>
    T Cast_To(const char* flag_option)
    {
        // Variable to hold the converted value
        T casted_value;
        
        // Transform the input to the given cast value
        std::istringstream ss(flag_option);

        // Move the data into the variable
        ss >> casted_value;

        return casted_value;
    }

};
