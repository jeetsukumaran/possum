///////////////////////////////////////////////////////////////////////////////
//
// MOPOGEN Population Genetic Summary Statistic(s) Calculator.
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

#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include "config.h"
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"
#include "possum_defs.hpp"
#include "textutil.hpp"
#include "cli.hpp"

#ifdef HAVE_CONFIG_H
#	include <config.h>
#   include "possum_info.h"
#else
#	include "win_config.h"
#endif

namespace possum {

///////////////////////////////////////////////////////////////////////////////
// Support Functions

// returns identification string to be displayed on splash/title/logs
std::string get_program_identification(const char * name) {
    std::ostringstream s;
    if (name == NULL) {
        name = PACKAGE_NAME;
    }
    s << name << " v" << PACKAGE_VERSION;
#if defined(BUILDDESC)
    s << " (" << BUILDDESC << ")";
#else
    s << " (" << __DATE__ << ")";
#endif
    return  s.str();
}

// Process string of numbers specifying sub-population sizes.
std::vector<unsigned int> parse_subpopulation_sizes(const std::string& pop_sizes_str) {
    std::vector<unsigned int> pop_sizes;
    if (pop_sizes_str.size() == 0) {
        return pop_sizes;
    }
    unsigned int total_size = 0;
    std::vector<std::string> tokens = textutil::split(pop_sizes_str, ",");
    pop_sizes.reserve(tokens.size());
    for (unsigned long i = 0; i < tokens.size(); ++i) {
        unsigned int psz = to_scalar<unsigned int>(textutil::strip(tokens[i]));
        pop_sizes.push_back(psz);
        total_size += psz;
    }
    return pop_sizes;
}

void report_invalid_subpopulation_sizes(const std::vector<unsigned int>& subpop_sizes,
    unsigned int sum_subpop_sizes, unsigned int expected_total) {
    std::cerr << "Expecting sum of samples from sub-populations to equal ";
    std::cerr << expected_total << ", but instead found ";
    std::cerr << sum_subpop_sizes << ": ";
    for (std::vector<unsigned int>::const_iterator x = subpop_sizes.begin();
            x != subpop_sizes.end();
            ++x)  {
        std::cerr << *x << " ";
    }
    std::cerr << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// DataFormats

DataFormatSpecification::DataFormatSpecification() {
    this->name_format_map_.insert(std::make_pair("nexus", MultiFormatReader::NEXUS_FORMAT));
    this->name_format_map_.insert(std::make_pair("dnafasta", MultiFormatReader::FASTA_DNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("rnafasta", MultiFormatReader::FASTA_RNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("aafasta", MultiFormatReader::FASTA_AA_FORMAT));
    this->name_format_map_.insert(std::make_pair("dnaphylip", MultiFormatReader::PHYLIP_DNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("rnaphylip", MultiFormatReader::PHYLIP_RNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("aaphylip", MultiFormatReader::PHYLIP_AA_FORMAT));
    this->name_format_map_.insert(std::make_pair("dnaphylipinterleaved", MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("rnaphylipinterleaved", MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("aaphylipinterleaved", MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT));
    this->name_format_map_.insert(std::make_pair("dnarelaxedphylip", MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("rnarelaxedphylip", MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("aarelaxedphylip", MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT));
    this->name_format_map_.insert(std::make_pair("dnarelaxedphylipinterleaved", MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("rnarelaxedphylipinterleaved", MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("discreterelaxedphylipinterleaved", MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT));
    this->name_format_map_.insert(std::make_pair("dnaaln", MultiFormatReader::ALN_DNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("rnaaln", MultiFormatReader::ALN_RNA_FORMAT));
    this->name_format_map_.insert(std::make_pair("aaaln", MultiFormatReader::ALN_AA_FORMAT));
}

void DataFormatSpecification::show_help_formats(ostream& out) const {
    out << "Available format specifications:" << std::endl;
    for (DataFormatSpecification::StringFormatMapType::const_iterator nfmi = this->name_format_map_.begin();
            nfmi != this->name_format_map_.end();
            ++nfmi) {
        out << "    " << nfmi->first << std::endl;
    }
}

MultiFormatReader::DataFormatType DataFormatSpecification::get_format_from_name(const std::string& data_format_name) const {
    std::string normalized_name = possum::textutil::lower(data_format_name);
    DataFormatSpecification::StringFormatMapType::const_iterator format_iter = this->name_format_map_.find(normalized_name);
    if (format_iter == this->name_format_map_.end()) {
        std::cerr << "'" << normalized_name << "' is not a recognized data format specification." << std::endl;
        this->show_help_formats(std::cerr);
        exit(1);
    }
    MultiFormatReader::DataFormatType ncl_format(format_iter->second);
    return ncl_format;
}

///////////////////////////////////////////////////////////////////////////////
// ResultsFormatter

ResultsFormatter::ResultsFormatter(bool write_field_name,
        bool extended_rows,
        const char * field_separator,
        bool skip_missing)
            : write_field_names_(write_field_name),
              extended_rows_(extended_rows),
              field_separator_(field_separator),
              skip_missing_(skip_missing) {

}

ResultsFormatter::~ResultsFormatter() {
}

void ResultsFormatter::add_field_name(const char* name) {
    std::vector<std::string>::iterator i = std::find(this->field_names_.begin(), this->field_names_.end(), name);
    if ( i == this->field_names_.end()) {
        this->field_names_.push_back(name);
    }
}

void ResultsFormatter::write_results(ostream& out) const {
    if (this->extended_rows_) {
        for (std::vector<std::string>::const_iterator fname = this->field_names_.begin();
                fname != this->field_names_.end();
                ++fname) {
            std::map<std::string, std::string>::const_iterator field = this->field_values_.find(*fname);
            std::map<std::string, std::string>::const_iterator default_field = this->default_values_.find(*fname);
            if (field != this->field_values_.end()
                    || default_field != this->default_values_.end()
                    || (!this->skip_missing_)) {
                if (this->write_field_names_) {
                    out << *fname;
                    out << this->field_separator_;
                }
                if (field != this->field_values_.end()) {
                    out << field->second;
                } else if (default_field != this->default_values_.end()) {
                    out << default_field->second;
                }
                out << std::endl;
            } // if field has value
        } // for each field
    } else {
        unsigned int field_count = 0;
        if (this->write_field_names_) {
            for (std::vector<std::string>::const_iterator fname = this->field_names_.begin();
                    fname != this->field_names_.end();
                    ++fname) {
                ++field_count;
                std::map<std::string, std::string>::const_iterator field = this->field_values_.find(*fname);
                std::map<std::string, std::string>::const_iterator default_field = this->default_values_.find(*fname);
                if (field != this->field_values_.end()
                        || default_field != this->default_values_.end()
                        || (!this->skip_missing_)) {
                    if (field_count > 1) {
                        out << this->field_separator_;
                    }
                    out << *fname;
                } else {
                    --field_count;
                }
            }
            out << std::endl;
        }
        field_count = 0;
        for (std::vector<std::string>::const_iterator fname = this->field_names_.begin();
                fname != this->field_names_.end();
                ++fname) {
            ++field_count;
            std::map<std::string, std::string>::const_iterator field = this->field_values_.find(*fname);
            std::map<std::string, std::string>::const_iterator default_field = this->default_values_.find(*fname);
            if (field != this->field_values_.end()) {
                if (field_count > 1) {
                    out << this->field_separator_;
                }
                out << field->second;
            } else if (default_field != this->default_values_.end()) {
                if (field_count > 1) {
                    out << this->field_separator_;
                }
                out << default_field->second;
            } else if (!this->skip_missing_) {
                if (field_count > 1) {
                    out << this->field_separator_;
                }
            } else {
                --field_count;
            }
        }
        out << std::endl;
    } // extended_rows/not extended_rows
}

} // namespace possum
