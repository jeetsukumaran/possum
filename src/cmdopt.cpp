///////////////////////////////////////////////////////////////////////////////
//
// Command-line options parsing library.
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

#include "cmdopt.hpp"

namespace cmdopt {

///////////////////////////////////////////////////////////////////////////////
// helper function

std::string wrap_text(const std::string& source,
        unsigned line_width,
        unsigned first_line_indent,
        unsigned subsequent_line_indent)  {
    std::string wrapped;
    unsigned col_count = 1;
    unsigned line_count = 1;
    std::string subsequent_line_indent_str(subsequent_line_indent, ' ');
    for (std::string::const_iterator s = source.begin();
            s != source.end();
            ++s, ++col_count) {

        if (*s == '\n') {
            wrapped += "\n";
            col_count = 0;
            line_count += 1;
            continue;
        }

        if (col_count > line_width) {
            std::string::size_type last_break = wrapped.find_last_of("\n");
            std::string::size_type wrap_pos = wrapped.find_last_of(" ");
            if (wrap_pos == std::string::npos or ((last_break != std::string::npos) and (last_break > wrap_pos))) {
                wrapped += "\n";
                col_count = 1;
            } else {
                wrapped.replace(wrap_pos, 1, "\n" + subsequent_line_indent_str);
                col_count = wrapped.size() - wrap_pos;
            }
        }

        if (col_count == 1 and line_count == 1 and first_line_indent > 0) {
            for (unsigned i = 0; i < first_line_indent; ++i) {
                wrapped += ' ';
            }
            col_count += first_line_indent;
        } else if (col_count == 1 and line_count > 1) {
            wrapped += subsequent_line_indent_str;
            col_count += subsequent_line_indent;
        }
        wrapped += *s;

    }

    return wrapped;
}

///////////////////////////////////////////////////////////////////////////////
// OptionArg

OptionArg::OptionArg(const char * help, const char * meta_var)
    : is_switch_(false),
      is_set_(false) {
    if (help != NULL) {
        this->help_ = help;
    }
    if (meta_var != NULL) {
        this->set_meta_var(meta_var);
    }
}

OptionArg::~OptionArg() {}

std::ostream& OptionArg::write_help(std::ostream& out) const {
    std::string help_str;
    help_str += "  ";
    if (this->short_flag_.size() > 0) {
        help_str += this->short_flag_;
        if (not this->is_switch_) {
            help_str += " ";
            if (this->meta_var_.size() == 0) {
                help_str += "VALUE";
            } else {
                help_str += this->meta_var_;
            }
        }
        if (this->long_flag_.size() > 0) {
            help_str += ", ";
        }
    }
    if (this->long_flag_.size() > 0) {
        help_str += this->long_flag_;
        if (not this->is_switch_) {
            help_str += "=";
            if (this->meta_var_.size() == 0) {
                help_str += "VALUE";
            } else {
                help_str += this->meta_var_;
            }
        }
    }
    if (this->help_.size() > 0) {
        if (help_str.size() > CMDOPTS_OPTION_COL_WIDTH-2) {
            help_str += "\n";
        } else {
            while (help_str.size() < CMDOPTS_OPTION_COL_WIDTH) {
                help_str += " ";
            }
        }
        std::string help_msg = this->help_;
        std::string::size_type defval = help_msg.find("%default");
        std::string replace_val = this->current_value_as_string();
        while (defval != std::string::npos) {
            help_msg.replace(defval, 8, replace_val.c_str());
            defval = help_msg.find("%default");
        }
        help_str += help_msg;
        std::string help_desc = wrap_text(help_str, CMDOPTS_LINE_WIDTH, 0, CMDOPTS_OPTION_COL_WIDTH);
        help_str = help_desc;
    }
    out << help_str;
    return out;
}

///////////////////////////////////////////////////////////////////////////////
// Specializations of TypedOptionArg

template <>
void TypedOptionArg<std::string>::process_value_string(const std::string& val_str) {
    *this->store_ = val_str;
    this->set_is_set(true);
}

///////////////////////////////////////////////////////////////////////////////
// OptionParser

OptionParser::OptionParser(
        const char * version,
        const char * description,
        const char * usage,
        const char * default_option_group)
    : show_help_(false),
      show_version_(false),
      single_hyphen_switch_(NULL),
      double_hyphen_switch_(false) {
    if (usage != NULL) {
        this->usage_.assign(usage);
    } else {
        this->usage_ = "%prog [options] [args]";
    }
    if (description != NULL) {
        this->description_.assign(description);
    }
    if (version != NULL) {
        this->version_.assign(version);
    }
    if (default_option_group != NULL) {
        this->default_option_group_ = default_option_group;
    } else {
        this->default_option_group_ = "Options:";
    }
    this->version_option_ = this->add_switch(&this->show_version_, NULL, "--version", "show program's version number and exit");
    this->help_option_ = this->add_switch(&this->show_help_, "-h", "--help",  "show this help message and exit");
}

OptionParser::~OptionParser() {
    for (std::vector<OptionArg *>::iterator oap = this->option_args_.begin();
            oap != this->option_args_.end();
            ++oap) {
        delete *oap;
    }
}

std::ostream& OptionParser::write_usage(std::ostream& out) const {
    if (this->usage_.size() != 0) {
        std::string usage = "Usage: " + this->usage_;
        out << usage << std::endl;
    }
    return out;
}

std::ostream& OptionParser::write_description(std::ostream& out) const {
    if (this->description_.size() != 0) {
        out << wrap_text(this->description_, CMDOPTS_LINE_WIDTH) << std::endl << std::endl;
    }
    return out;
}

std::ostream& OptionParser::write_version(std::ostream& out) const {
    out << this->version_ << std::endl;
    return out;
}

std::ostream& OptionParser::write_help(std::ostream& out) const {
    this->write_usage(out);
    out << std::endl;
    this->write_description(out);

    unsigned idx = 0;
    for (std::vector< const char *>::const_iterator ogni = this->option_group_names_.begin();
            ogni != this->option_group_names_.end();
            ++ogni) {
        std::map< const char *, std::vector<OptionArg *> >::const_iterator ogi = this->option_groups_.find(*ogni);
        const std::vector<OptionArg *>& opts = ogi->second;
        if (opts.size() > 0) {
            ++idx;
            if (idx > 1) {
                out << std::endl;
            }
            out << *ogni << std::endl;
            for (std::vector<OptionArg *>::const_iterator oi = opts.begin();
                    oi != opts.end();
                    ++oi) {
                (*oi)->write_help(out);
                out << std::endl;
            }
        }
    }
    return out;
}

void OptionParser::parse(int argc, char * argv[]) {

    for (int i = 1; i < argc; ++i) {

        if ( (!this->double_hyphen_switch_) && (argv[i][0] == '-') ) {

            if (strlen(argv[i]) == 2 && (strncmp(argv[i], "--", 2) == 0) ) {
                this->double_hyphen_switch_ =  true;
            } else if ((this->single_hyphen_switch_ != NULL)
                  && (strlen(argv[i]) == 1)
                  && (strncmp(argv[i], "-", 1) == 0)) {
                *this->single_hyphen_switch_ = true;
            } else {

                // actual flag processing

                if (strncmp(argv[i], "--", 2) == 0) {

                    // long-flag processing
                    std::string arg_name;
                    std::string arg_value;

                    bool parsing_name = true;
                    for (char *a = argv[i]; *a; ++a) {
                        if (parsing_name) {
                            if (*a == '=') {
                                parsing_name = false;
                            } else {
                                arg_name += *a;
                            }
                        } else {
                            arg_value += *a;
                        }
                    }

                    std::vector< std::string > matches;
                    for (std::map< std::string, OptionArg * >::iterator oai = this->key_opt_map_.begin();
                         oai != this->key_opt_map_.end();
                         ++oai) {
                        const std::string& a = oai->first;
                        if (a == arg_name) {
                            matches.clear();
                            matches.push_back(a);
                            break;
                        }
                        if (a.compare(0, arg_name.size(), arg_name) == 0 ) {
                            matches.push_back(a);
                        }
                    }

                    if (matches.size() == 0) {
                        std::cerr << "Unrecognized option '" << arg_name << "'" << std::endl;
                        exit(1);
                    } else if (matches.size() > 1) {
                        std::cerr << "Multiple matches found for option beginning with '" << arg_name << "':" << std::endl;
                        for (std::vector<std::string>::iterator mi = matches.begin(); mi != matches.end(); ++mi) {
                            std::cerr << *mi << std::endl;
                        }
                        exit(1);
                    }

                    OptionArg& oa = *(this->key_opt_map_[matches[0]]);
                    i = this->handle_value_assignment(oa, arg_name, arg_value, argc, argv, i);

                } else if (argv[i][0] == '-') {
                    // short flag processing
                    std::string arg(argv[i]);
                    std::string arg_name;
                    std::string arg_value;
                    if (arg.size() < 2) {
                        std::cerr << "Incomplete option '" << arg << "'" << std::endl;
                        exit(1);
                    }
                    if (arg.size() == 2) {
                        arg_name = arg;
                    } else {
                        arg_name = arg.substr(0, 2);
                        arg_value = arg.substr(2, arg.size());
                    }
                    std::map< std::string, OptionArg * >::iterator oai = this->key_opt_map_.find(arg_name);
                    if (oai == this->key_opt_map_.end()) {
                        std::cerr << "Unrecognized option '" << arg_name << "'" << std::endl;
                        exit(1);
                    }
                    OptionArg& oa = *(oai->second);
                    if (oa.is_switch()) {
                        this->handle_value_assignment(oa, arg_name, arg_value, argc, argv, i);
                        if (arg_value.size() > 0) {
                            std::string dummy_value;
                            for (std::string::iterator si = arg_value.begin(); si != arg_value.end(); ++si) {
                                std::string inferred_arg_name = "-";
                                inferred_arg_name.push_back(*si);
                                std::map< std::string, OptionArg * >::iterator sopti = this->key_opt_map_.find(inferred_arg_name);
                                if (sopti == this->key_opt_map_.end()) {
                                    std::cerr << "Unrecognized option '" << inferred_arg_name << "' (given in switch group '" << argv[i] << "')" << std::endl;
                                    exit(1);
                                }
                                if (!sopti->second->is_switch()) {
                                    std::cerr << "Option '" << inferred_arg_name << "' (given in switch group '" << argv[i] << "') is not a switch flag" << std::endl;
                                    exit(1);
                                }
                                this->handle_value_assignment(*(sopti->second), inferred_arg_name, dummy_value, argc, argv, i);
                            }
                        }
                    } else {
                        i = this->handle_value_assignment(oa, arg_name, arg_value, argc, argv, i);
                    }
                } // short flag processing
            } // actual flag processing
        } else {
            this->pos_args_.push_back(argv[i]);
        }

        // help option specified
        if (this->show_help_) {
            this->write_help(std::cout);
            exit(0);
        }

        // show version
        if (this->show_version_) {
            this->write_version(std::cout);
            exit(0);
        }
    }
}

OptionArg* OptionParser::get_option_ptr(const char * flag) {
    std::map< std::string, OptionArg * >::iterator oai = this->key_opt_map_.find(flag);
    assert (oai != this->key_opt_map_.end() );
    return oai->second;
}

OptionArg& OptionParser::get_option(const char * flag) {
    return *(get_option_ptr(flag));
}

bool OptionParser::is_set(const char * flag) {
    OptionArg& oa = this->get_option(flag);
    return oa.is_set();
}

void OptionParser::add_option_to_option_group(OptionArg * oa, const char * option_group) {
    if (option_group == NULL) {
        option_group = this->default_option_group_;
    }
    if (this->option_groups_.find(option_group) == this->option_groups_.end()) {
        this->option_group_names_.push_back(option_group);
    }
    this->option_groups_[option_group].push_back(oa);
}

} // namespace cmdopt
