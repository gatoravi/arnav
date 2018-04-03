/*  common.h -- misc utilities

    Copyright (c) 2015, The Griffith Lab

    Author: Avinash Ramu <aramu@genome.wustl.edu>

    This file is part of bcall.

    bcall is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bcall is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bcall.  If not, see <http://www.gnu.org/licenses/>. */

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
