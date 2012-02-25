///////////////////////////////////////////////////////////////////////////////
//
// POSSUM Population Genetic Summary Statistic(s) Calculator.
//
// Copyright 2010 Jeet Sukumaran.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POSSUM_CLI_H
#define POSSUM_CLI_H

#include <iostream>
#include <sstream>
#include <vector>
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"

namespace possum {

///////////////////////////////////////////////////////////////////////////////
// Support Functions

// returns identification string to be displayed on splash/title/logs
std::string get_program_identification(const char * name);

// Process string of numbers specifying sub-population sizes.
std::vector<unsigned int> parse_subpopulation_sizes(const std::string& pop_sizes_str);

// Reports error.
void report_invalid_subpopulation_sizes(const std::vector<unsigned int>& subpop_sizes,
    unsigned int sum_subpop_sizes, unsigned int expected_total);

// Converts from one simple streamable type to another.
template <typename T, typename U>
T to_scalar(U from) {
    std::ostringstream o;
    o << from;
    std::istringstream i(o.str());
    T target;
    i >> target;
    if (i.fail() || !i.eof()) {
        throw std::runtime_error("Failed to convert '" + o.str() + "' to integer");
    }
    return target;
}

///////////////////////////////////////////////////////////////////////////////
// DataFormats

class DataFormatSpecification {

    public:
        typedef std::map<std::string, MultiFormatReader::DataFormatType> StringFormatMapType;

        // constructor
        DataFormatSpecification();

        // display names of available formats
        void show_help_formats(ostream& out=std::cout) const;

        // returns NCL format id for given format string
        MultiFormatReader::DataFormatType get_format_from_name(const std::string& data_format_name) const;

    private:
        DataFormatSpecification::StringFormatMapType     name_format_map_;

};

///////////////////////////////////////////////////////////////////////////////
// ResultsFormatter

class ResultsFormatter {

    public:
        ResultsFormatter(bool write_field_name=true,
                bool extended_rows=false,
                const char * field_separator="\t",
                bool skip_missing=false);

        ~ResultsFormatter();

        void set_write_field_names(bool wfn) {
            this->write_field_names_ = wfn;
        }

        void set_extended_rows(bool ext) {
            this->extended_rows_ = ext;
        }

        void set_field_separator(const char * sep) {
            this->field_separator_ = sep;
        }

        void set_skip_missing(bool skip) {
            this->skip_missing_ = skip;
        }

        // adds a string to the list of field/column names if not already in list
        void add_field_name(const char* name);

        // as above, but registers a default value to be written.
        template <typename T>
        void add_field_name(const char* name, const T& default_value) {
            this->add_field_name(name);
            this->default_values_.insert(std::make_pair(std::string(name), to_scalar<std::string>(default_value)));
        }

        // adds `name` to the list of field/column names if not already in list,
        // converts `value` to string and add associates this with the field `name`
        template <typename T>
        void add_field_value(const char* name, const T& value) {
            this->add_field_name(name);
            this->field_values_.insert(std::make_pair(std::string(name), to_scalar<std::string>(value)));
        }

        // adds `name` to the list of field/column names if not already in list,
        // converts `value` to string and add associates this with the field `name`
        template <typename T>
        void add_field_value(const std::string& name, const T& value) {
            this->add_field_value(name.c_str(), value);
        }

        void write_results(ostream& out) const;

    private:
        std::vector<std::string>                field_names_;
        std::map<std::string, std::string>      field_values_;
        std::map<std::string, std::string>      default_values_;
        bool                                    write_field_names_;
        bool                                    extended_rows_;
        std::string                             field_separator_;
        bool                                    skip_missing_;

}; // ResultsFormatter

} // possum namespace

#endif
