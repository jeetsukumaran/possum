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

#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include "config.h"
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"
#include "possum_defs.hpp"
#include "character.hpp"
#include "cmdopt.hpp"
#include "cli.hpp"

void write_headers(bool extended_rows=false,
        bool report_filename=false,
        unsigned int num_subpops=0,
        bool by_population=false,
        bool by_subpopulation_pair=false,
        bool per_sequence=false,
        bool report_seq_len=true,
        bool report_label=false) {

    std::vector<std::string> field_names;

    if (report_label) {
        field_names.push_back("id");
    }
    if (report_filename) {
        field_names.push_back("file");
    }
    if (report_seq_len) {
        field_names.push_back("seq.len.original");
        field_names.push_back("seq.len.working");
    }
    field_names.push_back("pi");
    if (per_sequence) {
        field_names.push_back("pi.per.seq");
    }
    field_names.push_back("num.seg.sites");
    if (per_sequence) {
        field_names.push_back("num.seg.sites.per.seq");
    }
    field_names.push_back("num.mutations");
    if (per_sequence) {
        field_names.push_back("num.mutations.per.seq");
    }
    field_names.push_back("wattersons.theta");
    if (per_sequence) {
        field_names.push_back("wattersons.theta.per.seq");
    }
    field_names.push_back("tajimas.d");

    // populate report (multi-pop mode)
    if (num_subpops > 1) {

        field_names.push_back("pi.within");
        if (per_sequence) {
            field_names.push_back("pi.within.per.seq");
        }
        field_names.push_back("pi.between");
        if (per_sequence) {
            field_names.push_back("pi.between.per.seq");
        }
        field_names.push_back("pi.total");
        if (per_sequence) {
            field_names.push_back("pi.total.per.seq");
        }
        field_names.push_back("pi.net");
        if (per_sequence) {
            field_names.push_back("pi.net.per.seq");
        }
        field_names.push_back("fst");
        if (by_population) {
            for (unsigned int i = 0; i < num_subpops; ++i) {
                std::string label = possum::compose_subpopulation_label("pi.within", i);
                field_names.push_back(label);
                if (per_sequence) {
                    field_names.push_back(label+".per.seq");
                }
            }
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    std::string label = possum::compose_subpopulation_label("pi.between", i, j);
                    field_names.push_back(label);
                    if (per_sequence) {
                        field_names.push_back(label+".per.seq");
                    }
                }
            }
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    std::string label = possum::compose_subpopulation_label("pi.net", i, j);
                    field_names.push_back(label);
                    if (per_sequence) {
                        field_names.push_back(label+".per.seq");
                    }
                }
            }
        }
        if (by_subpopulation_pair) {
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    std::string label = possum::compose_subpopulation_label("pi", i, j);
                    field_names.push_back(label);
                    if (per_sequence) {
                        field_names.push_back(label+".per.seq");
                    }
                }
            }
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    std::string label = possum::compose_subpopulation_label("wattersons.theta", i, j);
                    field_names.push_back(label);
                    if (per_sequence) {
                        field_names.push_back(label+".per.seq");
                    }
                }
            }
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    std::string label = possum::compose_subpopulation_label("tajimas.d", i, j);
                    field_names.push_back(label);
                }
            }
        }
    }
    std::string sep;
    if (extended_rows) {
        sep = "\n";
    } else {
        sep = "\t";
    }
    for (unsigned int i = 0; i < field_names.size(); ++i) {
        if (i > 0) {
            std::cout << sep;
        }
        std::cout << field_names[i];
    }
    std::cout << std::endl;
}

void process_filestream(std::istream& input,
        MultiFormatReader::DataFormatType ncl_format,
        const std::vector<unsigned int>& pop_sizes,
        bool by_population=false,
        bool by_subpopulation_pair=false,
        bool extended_rows=false,
        bool no_field_names=false,
        bool per_sequence=false,
        const char * filename=NULL,
        bool report_seq_len=true,
        unsigned int sum_subpop_sizes=0,
        const char * id_label=NULL) {

    // load data
    possum::SequenceData seq_data;
    seq_data.set_quiet(true);
    seq_data.set_calc_subpopulation_pair_stats(by_subpopulation_pair);
    seq_data.read(input, ncl_format);

    // for results
    possum::ResultsFormatter rfmt(!no_field_names, extended_rows, "\t", false);

    // check/set population mode
    bool multi_pop_mode = false;
    if (pop_sizes.size() > 0) {
        // multi population mode
        if ( sum_subpop_sizes > 0 && sum_subpop_sizes != seq_data.size() ) {
            possum::report_invalid_subpopulation_sizes(pop_sizes, sum_subpop_sizes, seq_data.size());
        }
        multi_pop_mode = true;
        seq_data.set_subpopulations(pop_sizes);
    }

    // calc stats
    seq_data.calc_stats();

    // populate report (single + multi-pop mode)
    if (id_label != NULL) {
        rfmt.add_field_value("id", id_label);
    }
    if (filename != NULL) {
        rfmt.add_field_value("file", filename);
    }
    if (report_seq_len) {
        rfmt.add_field_value("seq.len.original", seq_data.get_original_seq_len());
        rfmt.add_field_value("seq.len.working", seq_data.get_working_seq_len());
    }
    rfmt.add_field_value("pi", seq_data.get_mean_pdist());
    if (per_sequence) {
        rfmt.add_field_value("pi.per.seq", seq_data.get_avg_pairwise_diffs());
    }
    rfmt.add_field_value("num.seg.sites", static_cast<float>(seq_data.get_num_segregating_sites())/seq_data.get_seq_len());
    if (per_sequence) {
    rfmt.add_field_value("num.seg.sites.per.seq", seq_data.get_num_segregating_sites());
    }
    rfmt.add_field_value("num.mutations", static_cast<float>(seq_data.get_num_mutations())/seq_data.get_seq_len());
    if (per_sequence) {
        rfmt.add_field_value("num.mutations.per.seq", seq_data.get_num_mutations());
    }
    rfmt.add_field_value("wattersons.theta", static_cast<float>(seq_data.get_wattersons_theta())/seq_data.get_seq_len());
    if (per_sequence) {
        rfmt.add_field_value("wattersons.theta.per.seq", seq_data.get_wattersons_theta());
    }
    rfmt.add_field_value("tajimas.d", seq_data.get_tajimas_d());

    // populate report (multi-pop mode)
    if (multi_pop_mode) {
        unsigned int num_subpops = seq_data.get_num_subpops();

        rfmt.add_field_value("pi.within", static_cast<float>(seq_data.get_avg_pairwise_diffs_within())/seq_data.get_seq_len());
        if (per_sequence) {
            rfmt.add_field_value("pi.within.per.seq", seq_data.get_avg_pairwise_diffs_within());
        }
        rfmt.add_field_value("pi.between", static_cast<float>(seq_data.get_avg_pairwise_diffs_between())/seq_data.get_seq_len());
        if (per_sequence) {
            rfmt.add_field_value("pi.between.per.seq", seq_data.get_avg_pairwise_diffs_between());
        }
        rfmt.add_field_value("pi.total", static_cast<float>(seq_data.get_avg_pairwise_diffs_total())/seq_data.get_seq_len());
        if (per_sequence) {
            rfmt.add_field_value("pi.total.per.seq", seq_data.get_avg_pairwise_diffs_total());
        }
        rfmt.add_field_value("pi.net", static_cast<float>(seq_data.get_avg_pairwise_diffs_net())/seq_data.get_seq_len());
        if (per_sequence) {
            rfmt.add_field_value("pi.net.per.seq", seq_data.get_avg_pairwise_diffs_net());
        }
        rfmt.add_field_value("fst", seq_data.get_fst());
        if (by_population) {
            for (unsigned int i = 0; i < num_subpops; ++i) {
                float value = seq_data.get_pairwise_diffs_within(i);
                std::string label = possum::compose_subpopulation_label("pi.within", i);
                rfmt.add_field_value(label, static_cast<float>(value)/seq_data.get_seq_len());
                if (per_sequence) {
                    rfmt.add_field_value(label+".per.seq", value);
                }
            }
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    float value = seq_data.get_pairwise_diffs_between(i, j);
                    std::string label = possum::compose_subpopulation_label("pi.between", i, j);
                    rfmt.add_field_value(label, static_cast<float>(value)/seq_data.get_seq_len());
                    if (per_sequence) {
                        rfmt.add_field_value(label+".per.seq", value);
                    }
                }
            }
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    float value = seq_data.get_pairwise_diffs_net(i, j);
                    std::string label = possum::compose_subpopulation_label("pi.net", i, j);
                    rfmt.add_field_value(label, static_cast<float>(value)/seq_data.get_seq_len());
                    if (per_sequence) {
                        rfmt.add_field_value(label+".per.seq", value);
                    }
                }
            }
        }
        if (by_subpopulation_pair) {
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    float value = seq_data.get_subpopulation_pair_avg_pairwise_diffs(i, j);
                    std::string label = possum::compose_subpopulation_label("pi", i, j);
                    rfmt.add_field_value(label, static_cast<float>(value)/seq_data.get_seq_len());
                    if (per_sequence) {
                        rfmt.add_field_value(label+".per.seq", value);
                    }
                }
            }
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    float value = seq_data.get_subpopulation_pair_wattersons_theta(i, j);
                    std::string label = possum::compose_subpopulation_label("wattersons.theta", i, j);
                    rfmt.add_field_value(label, static_cast<float>(value)/seq_data.get_seq_len());
                    if (per_sequence) {
                        rfmt.add_field_value(label+".per.seq", value);
                    }
                }
            }
            for (unsigned int i = 0; i < num_subpops - 1; ++i) {
                for (unsigned int j = i + 1; j < num_subpops ; ++j) {
                    std::string label = possum::compose_subpopulation_label("tajimas.d", i, j);
                    float value = seq_data.get_subpopulation_pair_tajimas_d(i, j);
                    rfmt.add_field_value(label, value);
                }
            }
        }
    }
    rfmt.write_results(std::cout);
}

void process_filepath(const std::string& filepath,
        MultiFormatReader::DataFormatType ncl_format,
        const std::vector<unsigned int>& pop_sizes,
        bool by_population=false,
        bool by_subpopulation_pair=false,
        bool extended_rows=false,
        bool no_field_names=false,
        bool per_sequence=false,
        bool report_filename=false,
        bool report_seq_len=true,
        unsigned int sum_subpop_sizes=0,
        const char * id_label = NULL) {
    std::ifstream input(filepath.c_str());
    if (!input.good()) {
        std::cerr << "Error opening data file: '" << filepath << "'" << std::endl;
        exit(1);
    }
    const char * filename = NULL;
    if (report_filename) {
        filename = filepath.c_str();
    }
    process_filestream(input,
            ncl_format,
            pop_sizes,
            by_population,
            by_subpopulation_pair,
            extended_rows,
            no_field_names,
            per_sequence,
            filename,
            report_seq_len,
            sum_subpop_sizes,
            id_label);
}

int main(int argc, char* argv[]) {
    possum::DataFormatSpecification data_formats;
    bool help_formats = false;
    std::string data_format_name = "nexus";
    std::string pop_sizes_str;
    std::string files_from;
    bool field_names_only = false;
    bool by_population = false;
    bool by_subpopulation_pair=false;
    bool quiet = false;
    bool no_field_names = false;
    bool extended_rows = false;
    bool per_sequence = false;
    bool from_std_input = false;
    bool report_filenames = false;
    bool report_seq_len=false;
    bool report_label=false;
    const char * id_label=NULL;
    std::string label_str;

    std::string program_ident = possum::get_program_identification("possum-sum");
    cmdopt::OptionParser parser(
            program_ident.c_str(),
            "Calculates summary statistics for molecular population genetic sequences.",
            "possum-sum [options] [FILE [FILE [FILE [...]]]]",
            "General Options:");

    parser.add_switch(&field_names_only, NULL, "--show-fields",
            "write field (column) names and exit, without processing any files or data");
    parser.add_switch(&quiet, "-q", "--quiet",
            "suppress progress messages; only output results");

    parser.add_single_hyphen_switch(&from_std_input, "process data from standard input", "Source Options:");
    parser.add_option<std::string>(&files_from, "-@", "--files-from",
            "process data from filenames given in LISTFILE", "LISTFILE", "Source Options:");

    parser.add_option<std::string>(&data_format_name, "-f", "--format",
            "data format of input files (default = 'nexus'); see 'possum --help-formats' for available choices",
            NULL, "Source Format Options:");
    parser.add_switch(&help_formats, NULL, "--help-formats",
            "show supported input data formats and exit",
            NULL, "Source Format Options:");

    parser.add_option<std::string>(&pop_sizes_str, "-N", "--subpopulation-sizes",
            "comma-separated list of (sub)population sizes if samples are from multiple populations"
            " (samples in the data file are assumed to be grouped by subpopulation and"
            " listed in the same order as specified here)", "'N1,N2,...,Nk'", "Subpopulation Options:");

    parser.add_switch(&by_population, "-p", "--report-subpopulations",
            "report detailed statistics on sub-populations", NULL, "Statistics Options:");
    parser.add_switch(&by_subpopulation_pair, "-P", "--report-subpopulation-pairs",
            "report detailed statistics on pairs of sub-populations", NULL, "Statistics Options:");
    parser.add_switch(&per_sequence, "-s", "--report-per-sequence",
            "report per-sequence results (in addition to per-site results)", NULL, "Statistics Options:");
    parser.add_switch(&report_filenames, "-n", "--report-filenames",
            "include names of data files in output", NULL, "Statistics Options:");
    parser.add_option<std::string>(&label_str, "-i", "--report-id",
            "add identifier column with label", NULL, "Statistics Options:");
    parser.add_switch(&report_seq_len, "-l", "--report-sequence-length",
            "include length of sequences (before and after sites with uncertain/polymorphic characters are dropped) in output", NULL, "Statistics Options:");

    parser.add_switch(&no_field_names, "-v", "--values-only",
            "raw results only: do not write column/field names", NULL, "Formatting Options:");
    parser.add_switch(&extended_rows, "-x", "--extended-format",
            "write results reporting individual statistics by rows instead of by columns", NULL, "Formatting Options:");

    parser.parse(argc, argv);
    std::vector<std::string> args = parser.get_args();

    if (!quiet) {
        std::cerr << program_ident << "\n" << std::endl;
    }

    if (help_formats) {
        data_formats.show_help_formats(std::cout);
        exit(0);
    }

    std::vector<unsigned int> pop_sizes = possum::parse_subpopulation_sizes(pop_sizes_str);
    unsigned int sum_subpop_sizes = 0;
    if (pop_sizes.size() > 0) {
        if (!quiet) {
            std::cerr << "Subpopulation sizes:";
        }
        for (std::vector<unsigned int>::const_iterator ci = pop_sizes.begin(); ci != pop_sizes.end(); ++ci) {
            if (!quiet) {
                std::cerr << " " << *ci;
            }
            sum_subpop_sizes += *ci;
        }
        if (!quiet) {
            std::cerr << "." << std::endl << std::endl;
        }

    }

    if (label_str.size() > 0) {
        report_label = true;
        id_label = label_str.c_str();
    }

    if (field_names_only) {
        write_headers(extended_rows,
                report_filenames,
                pop_sizes.size(),
                by_population,
                by_subpopulation_pair,
                per_sequence,
                report_seq_len,
                report_label);
        exit(0);
    }

    // parse input data source
    if (files_from.size() == 0 && args.size() == 0 && !from_std_input && !field_names_only) {
        parser.write_usage(std::cerr);
        exit(1);
    } else if (from_std_input && args.size() != 0) {
        std::cerr << "Cannot read data from both standard input and data files given as arguments" << std::endl;
        exit(1);
    } else if (from_std_input && files_from.size() != 0) {
        std::cerr << "Cannot read data from both standard input and data files specified by LISTFILE ('-T') argument" << std::endl;
        exit(1);
    } else if (files_from.size() != 0 && args.size() != 0) {
        std::cerr << "Cannot read data from both LISTFILE ('-T') and files given as arguments" << std::endl;
        exit(1);
    }

    MultiFormatReader::DataFormatType ncl_format = data_formats.get_format_from_name(data_format_name);
    if (!quiet) {
        std::cerr << "Reading '" << data_format_name << "' format data." << std::endl;
    }

    if (from_std_input) {
        process_filestream(std::cin,
            ncl_format,
            pop_sizes,
            by_population,
            by_subpopulation_pair,
            extended_rows,
            no_field_names,
            per_sequence,
            NULL,
            report_seq_len,
            sum_subpop_sizes,
            id_label);
    } else if (args.size() > 0) {
        for (unsigned int argi = 0; argi < args.size(); ++argi) {
            bool suppress_headers = true;
            if (!no_field_names && (argi == 0 || extended_rows)) {
                suppress_headers = false;
            }
            process_filepath(args[argi],
                ncl_format,
                pop_sizes,
                by_population,
                by_subpopulation_pair,
                extended_rows,
                suppress_headers,
                per_sequence,
                report_filenames,
                report_seq_len,
                sum_subpop_sizes,
                id_label);
        }
    } else if (files_from.size() > 0) {
        std::ifstream input(files_from.c_str());
        if (!input.good()) {
            std::cerr << "Error opening list file: '" << files_from << "'" << std::endl;
            exit(1);
        }
        if (!quiet) {
            std::cerr << std::endl;
        }
        unsigned int count = 0;
        while (!input.eof()) {
            std::string datafile;
            getline(input, datafile);
            if (datafile.size() > 0) {
                bool suppress_headers = true;
                if (!no_field_names && (count == 0 || extended_rows)) {
                    suppress_headers = false;
                }
                process_filepath(datafile.c_str(),
                    ncl_format,
                    pop_sizes,
                    by_population,
                    by_subpopulation_pair,
                    extended_rows,
                    suppress_headers,
                    per_sequence,
                    report_filenames,
                    report_seq_len,
                    sum_subpop_sizes,
                    id_label);
                count += 1;
            }
        }
    }

}
