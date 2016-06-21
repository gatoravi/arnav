/*  common.h -- misc utilities

    Copyright (c) 2015, The Griffith Lab

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef COMMON_H_
#define COMMON_H_

#include <iostream>
#include <sstream>

using namespace std;

namespace common {
    //Convert a number to a string
    template <typename T>
        string num_to_str(T num) {
            stringstream ss;
            ss << num;
            return ss.str();
    }

    //Convert a number-string to a uint32
    inline uint32_t str_to_uint32(string num) {
            stringstream ss;
            uint32_t num_uint;
            ss << num;
            ss >> num_uint;
            return num_uint;
    }

    //Convert a number-string to a double
    inline double str_to_double(string num) {
            stringstream ss;
            double d1;
            ss << num;
            ss >> d1;
            return d1;
    }
}

#endif
