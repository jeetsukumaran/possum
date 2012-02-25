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
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cassert>
#include <map>
#include <set>
#include "character.hpp"

namespace possum {

// Support Functions {{{1
// =============================================================================

std::string compose_subpopulation_label(const std::string& stat_name, unsigned int i, unsigned int j) {
    std::ostringstream label;
    label << stat_name << "." << i + 1 << "." << j + 1;
    return label.str();
}

std::string compose_subpopulation_label(const std::string& stat_name, unsigned int i) {
    std::ostringstream label;
    label << stat_name << "." << i + 1;
    return label.str();
}

unsigned long binomial_coefficient(unsigned long n, unsigned long k) {
    if (k > n)
        return 0;
    if (k > n/2)
        k = n-k; // Take advantage of symmetry
    double accum = 1;
    unsigned i;
    for (i = 1; i <= k; i++)
         accum = accum * (n-k+i) / i;
    return static_cast<unsigned long>(accum + 0.5); // avoid rounding error
}

// 1}}}

// PairwiseDifferences {{{1
// =============================================================================

void PairwiseDifferences::set_num_diffs(unsigned int i, unsigned int j, unsigned int num_diffs) {
    std::pair<unsigned int, unsigned int> taxon_pair = std::make_pair(i, j);
    this->diff_counts_.insert(std::make_pair(taxon_pair, num_diffs));
}

unsigned int PairwiseDifferences::operator()(unsigned int i, unsigned int j) const {
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>::const_iterator idx = this->diff_counts_.find(std::make_pair(i, j));
    if (idx == this->diff_counts_.end()) {
        idx = this->diff_counts_.find(std::make_pair(j, i));
    }
    if (idx == this->diff_counts_.end()) {
        std::ostringstream o;
        o << "Taxon pair (" << i << ", " << j << ") not found";
        throw std::runtime_error(o.str());
    }
    return idx->second;
}

// 1}}}

// CharacterStateMatrix {{{1
// =============================================================================

// lifecycle {{{2
// -----------------------------------------------------------------------------

CharacterStateMatrix::CharacterStateMatrix() {
}

CharacterStateMatrix::~CharacterStateMatrix() {
}

// 2}}}

// metrics {{{2
// -----------------------------------------------------------------------------

bool CharacterStateMatrix::empty() const {
    return this->state_vectors_.empty();
}

CharacterStateMatrix::size_type CharacterStateMatrix::size() const {
    return this->state_vectors_.size();
}

CharacterStateMatrix::size_type CharacterStateMatrix::capacity() const {
    return this->state_vectors_.capacity();
}

CharacterStateMatrix::size_type CharacterStateMatrix::max_size() const {
    return this->state_vectors_.max_size();
}

CharacterStateVectorType& CharacterStateMatrix::at(CharacterStateMatrix::size_type n) {
    return this->state_vectors_.at(n);
}

// 2}}}

// accessors {{{2
// -----------------------------------------------------------------------------

const CharacterStateVectorType& CharacterStateMatrix::at(CharacterStateMatrix::size_type n) const {
    return this->state_vectors_.at(n);
}

CharacterStateVectorType& CharacterStateMatrix::operator[](CharacterStateMatrix::size_type n) {
    return this->state_vectors_[n];
}

const CharacterStateVectorType& CharacterStateMatrix::operator[](CharacterStateMatrix::size_type n) const {
    return this->state_vectors_[n];
}

CharacterStateVectorType& CharacterStateMatrix::back() {
    return this->state_vectors_.back();
}

const CharacterStateVectorType& CharacterStateMatrix::back() const {
    return this->state_vectors_.back();
}

CharacterStateVectorType& CharacterStateMatrix::front() {
    return this->state_vectors_.front();
}

const CharacterStateVectorType& CharacterStateMatrix::front() const {
    return this->state_vectors_.front();
}

// 2}}}

// iteration {{{2
// -----------------------------------------------------------------------------

CharacterStateMatrix::iterator CharacterStateMatrix::begin() {
    return this->state_vectors_.begin();
}

CharacterStateMatrix::const_iterator CharacterStateMatrix::begin() const {
    return this->state_vectors_.begin();
}

CharacterStateMatrix::iterator CharacterStateMatrix::end() {
    return this->state_vectors_.end();
}

CharacterStateMatrix::const_iterator CharacterStateMatrix::end() const {
    return this->state_vectors_.end();
}

CharacterStateMatrix::reverse_iterator CharacterStateMatrix::rbegin() {
    return this->state_vectors_.rbegin();
}

CharacterStateMatrix::const_reverse_iterator CharacterStateMatrix::rbegin() const {
    return this->state_vectors_.rbegin();
}

CharacterStateMatrix::reverse_iterator CharacterStateMatrix::rend() {
    return this->state_vectors_.rend();
}

CharacterStateMatrix::const_reverse_iterator CharacterStateMatrix::rend() const {
    return this->state_vectors_.rend();
}

// 2}}}

// mutators {{{2
// -----------------------------------------------------------------------------

void CharacterStateMatrix::assign(CharacterStateMatrix::size_type n, const CharacterStateVectorType& u) {
    this->state_vectors_.assign(n, u);
}

void CharacterStateMatrix::clear() {
    this->state_vectors_.clear();
}

CharacterStateMatrix::iterator CharacterStateMatrix::erase(CharacterStateMatrix::iterator position) {
    return this->state_vectors_.erase(position);
}

CharacterStateMatrix::iterator CharacterStateMatrix::erase(CharacterStateMatrix::iterator first, CharacterStateMatrix::iterator last) {
    return this->state_vectors_.erase(first, last);
}

void CharacterStateMatrix::push_back(const CharacterStateVectorType& v) {
    this->state_vectors_.push_back(v);
}

CharacterStateVectorType& CharacterStateMatrix::new_vector() {
    this->state_vectors_.push_back(CharacterStateVectorType());
    return this->state_vectors_.back();
}

void CharacterStateMatrix::pop_back() {
    this->state_vectors_.pop_back();
}

CharacterStateMatrix::iterator CharacterStateMatrix::insert(CharacterStateMatrix::iterator position, const CharacterStateVectorType& x) {
    return this->state_vectors_.insert(position, x);
}

void CharacterStateMatrix::insert(CharacterStateMatrix::iterator position, CharacterStateMatrix::size_type n, const CharacterStateVectorType& x ) {
    this->state_vectors_.insert(position, n, x);
}

void CharacterStateMatrix::reserve(CharacterStateMatrix::size_type n) {
    this->state_vectors_.reserve(n);
}

void CharacterStateMatrix::resize(CharacterStateMatrix::size_type n, CharacterStateVectorType c) {
    this->state_vectors_.resize(n, c);
}

void CharacterStateMatrix::swap(CharacterStateMatrix m) {
    this->state_vectors_.swap(m.state_vectors_);
}

CharacterStateMatrix& CharacterStateMatrix::operator=(const CharacterStateMatrix& m) {
    this->state_vectors_ = m.state_vectors_;
    return *this;
}

void CharacterStateMatrix::drop_columns(const std::vector<bool>& to_drop) {
    for (unsigned int ri = 0;
            ri < this->state_vectors_.size();
            ++ri) {
        CharacterStateVectorType& old_vec = this->state_vectors_[ri];
        CharacterStateVectorType new_vec;
        new_vec.reserve(old_vec.size());
        assert(to_drop.size() == old_vec.size());
        for (unsigned int ci = 0; ci < old_vec.size(); ++ci) {
            if (!to_drop[ci]) {
                new_vec.push_back(old_vec[ci]);
            }
        }
        this->state_vectors_[ri].swap(new_vec);
    }
}

// 2}}}

// stats/calcs {{{2
// -----------------------------------------------------------------------------

PairwiseDifferences& CharacterStateMatrix::get_pairwise_diff_counts(PairwiseDifferences& pdiffs) const {
    // TODO: strip off columns with any ambiguous/uncertain or gap characters
    // from this->state_ vectors.
    assert(this->state_vectors_.size() > 0);
    unsigned int seq_len = this->seq_len();
    pdiffs.clear();
    for (unsigned int s1 = 0; s1 < (this->state_vectors_.size()-1); ++s1) {
        for (unsigned int s2 = s1 + 1; s2 <  this->state_vectors_.size(); ++s2) {
            unsigned int num_diffs = 0;
            CharacterStateVectorType seq1 = this->state_vectors_[s1];
            CharacterStateVectorType seq2 = this->state_vectors_[s2];
            if (seq1.size() != seq_len || seq2.size() != seq_len) {
                unsigned int err_seq_idx = 0;
                if (seq1.size() != seq_len) {
                    err_seq_idx = s1+1;
                } else if (seq2.size() != seq_len) {
                    err_seq_idx = s2+1;
                }
                std::cerr << "ERROR: Expecting " << seq_len << " characters in all sequences";
                std::cerr << " (based on first sequence length), but sequence #" << err_seq_idx << " has ";
                std::cerr << seq1.size() << " characters" << std::endl;
                exit(1);
            }
            for (unsigned int cidx = 0; cidx < seq1.size(); ++cidx) {
                if (seq1[cidx] != seq2[cidx]) {
                    ++num_diffs;
                }
            }
            pdiffs.set_num_diffs(s1, s2, num_diffs);
        }
    }
    return pdiffs;
}

// 2}}}

// 1}}}

// SequenceData {{{1
// =============================================================================

// lifecycle {{{2
// -----------------------------------------------------------------------------

SequenceData::SequenceData()
    : original_seq_len_(0),
      working_seq_len_(0),
      avg_pairwise_diffs_(0),
      avg_pairwise_diffs_total_(0),
      avg_pairwise_diffs_within_(0),
      avg_pairwise_diffs_between_(0),
      avg_pairwise_diffs_net_(0),
      num_mutations_(0),
      wattersons_theta_(0),
      tajimas_d_(0),
      use_num_mutations_instead_of_seg_sites_(0),
      calc_subpopulation_pair_stats_(true),
      quiet_(false) {
}

SequenceData::~SequenceData() {

}

// 2}}}

// serialization/deserialization {{{2
// -----------------------------------------------------------------------------

void SequenceData::read(const char * filepath, MultiFormatReader::DataFormatType fmt) {
	assert(filepath);
	if (not this->quiet_) {
	    std::cerr << "Reading '" << filepath << "'" << std::endl;
	}
    std::ifstream input(filepath);
    this->read(input, fmt);
}

void SequenceData::read(std::istream& input, MultiFormatReader::DataFormatType fmt) {
	try {
	    MultiFormatReader nexus_reader;
        nexus_reader.SetWarningOutputLevel(NxsReader::SKIPPING_CONTENT_WARNING);
		try {
			nexus_reader.ReadStream(input, fmt);
			unsigned total_seqs_read = 0;
	        const unsigned ntaxblocks = nexus_reader.GetNumTaxaBlocks();
	        for (unsigned tbi = 0; tbi < ntaxblocks; ++tbi) {
                const NxsTaxaBlock * tb = nexus_reader.GetTaxaBlock(tbi);
                const unsigned ntax = tb->GetNTax();
                const unsigned ncharblocks = nexus_reader.GetNumCharactersBlocks(tb);
                if (ncharblocks == 0) {
                    continue;
                }
                for (unsigned cbi = 0; cbi < ncharblocks; ++cbi) {
                    const NxsCharactersBlock * cb = nexus_reader.GetCharactersBlock(tb, cbi);
                    const NxsDiscreteDatatypeMapper * dm = cb->GetDatatypeMapperForChar(0);
                    if (dm == NULL) {
                        throw NxsNCLAPIException("No DatatypeMapper in SequenceData::read()");
                    }
                    if (cb->IsMixedType()) {
                        throw NxsNCLAPIException("Mixed datatypes are not supported");
                    }
                    const unsigned nchar = cb->GetNCharTotal();
                    this->col_has_gap_or_missing_.assign(nchar, false);
                    for (unsigned tax_idx = 0; tax_idx < ntax; tax_idx++) {
                        const std::string tax_label = NxsString::GetEscaped(tb->GetTaxonLabel(tax_idx));
                        this->taxon_labels_.push_back(tax_label);
                        CharacterStateVectorType& sv = this->state_matrix_.new_vector();
                        sv.reserve(nchar);
                        const NxsDiscreteStateRow & row = cb->GetDiscreteMatrixRow(tax_idx);
                        for (unsigned char_idx = 0; char_idx < nchar; ++char_idx) {
                            CharacterStateType c = row[char_idx];
                            sv.push_back(c);
                            if (c < 0 || c > 3) {
                                this->col_has_gap_or_missing_[char_idx] = true;
                            }
                        } // iterate over chars
                        ++total_seqs_read;
                    } // iterate over taxa
                } // iterate over char blocks
            } // iterate over taxa blocks
            if (not this->quiet_) {
                std::cerr << total_seqs_read << " sequences read.\n";
            }
        } catch (...) {
			nexus_reader.DeleteBlocksFromFactories();
			throw;
        }
        nexus_reader.DeleteBlocksFromFactories();
    } catch (const NxsException &x) {
		std::cerr << "Error:\n " << x.msg << endl;
		if (x.line >= 0)
			std::cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
	}
}

// 2}}}

// subpopulations {{{2
// -----------------------------------------------------------------------------


void SequenceData::set_subpopulations(const std::vector<unsigned int>& pop_sizes) {
    unsigned int total = 0;
    this->subpopulation_taxon_idxs_.clear();
    std::vector<unsigned int> x;
    this->subpopulation_taxon_idxs_.assign(pop_sizes.size(), x);
    for (unsigned int subpop_idx = 0; subpop_idx < pop_sizes.size(); ++subpop_idx) {
        unsigned int pop_size = pop_sizes.at(subpop_idx);
        this->subpopulation_taxon_idxs_.reserve(pop_size);
        for (unsigned int j = 0; j < pop_size; ++j) {
            this->subpopulation_taxon_idxs_.at(subpop_idx).push_back(total + j);
        }
        total += pop_size;
    }
    if (total != this->size()) {
        std::cerr << "Subpopulation numbers do not sum to total number of sequences:";
        std::cerr << " expecting " << this->size() << " but instead found " << total << "." << std::endl;
        exit(1);
    }
}

void SequenceData::report_subpopulations(std::ostream& out) const {
    out << this->subpopulation_taxon_idxs_.size() << " subpopulations." << std::endl;
    for (unsigned int i = 0; i < this->subpopulation_taxon_idxs_.size(); ++i ) {
        out << "[" << i << ": " << this->subpopulation_taxon_idxs_.at(i).size() << " samples]: ";
        for (std::vector<unsigned int>::const_iterator j = this->subpopulation_taxon_idxs_.at(i).begin();
                j != this->subpopulation_taxon_idxs_.at(i).end();
                ++j) {
            out << *j << " ";
        }
        out << std::endl;
    }
}

// 2}}}

// operations {{{2
// -----------------------------------------------------------------------------

void SequenceData::drop_gap_or_uncertain_columns(void) {
    this->state_matrix_.drop_columns(this->col_has_gap_or_missing_);
    this->col_has_gap_or_missing_.assign(this->state_matrix_.seq_len(), false);
}

// 2}}}

// counting/stats (public API) {{{2
// -----------------------------------------------------------------------------

void SequenceData::calc_stats() {
    this->original_seq_len_ = this->state_matrix_.seq_len();
    this->drop_gap_or_uncertain_columns();
    this->working_seq_len_ = this->state_matrix_.seq_len();
    this->state_matrix_.get_pairwise_diff_counts(this->pairwise_diffs_);
    this->calc_diversity_general_();
    this->calc_summary_stats_general_();
    if (!this->subpopulation_taxon_idxs_.empty()) {
        this->calc_diversity_subpopulations_();
        if (this->calc_subpopulation_pair_stats_) {
            this->calc_summary_stats_subpopulation_pairs_();
        }
    }
}

float SequenceData::get_avg_pairwise_diffs() const {
    return this->avg_pairwise_diffs_;
}

float SequenceData::get_avg_pairwise_diffs_within() const {
    return this->avg_pairwise_diffs_within_;
}

float SequenceData::get_avg_pairwise_diffs_between() const {
    return this->avg_pairwise_diffs_between_;
}

float SequenceData::get_avg_pairwise_diffs_total() const {
    return this->avg_pairwise_diffs_total_;
}

float SequenceData::get_avg_pairwise_diffs_net() const {
    return this->avg_pairwise_diffs_net_;
}

unsigned int SequenceData::get_num_segregating_sites() const {
    return this->segregating_sites_.size();
}

float SequenceData::get_mean_pdist() const {
    return this->avg_pairwise_diffs_ / this->working_seq_len_;
}

float SequenceData::get_mean_pdist_within() const {
    return this->avg_pairwise_diffs_within_ / this->working_seq_len_;
}

float SequenceData::get_mean_pdist_between() const {
    return this->avg_pairwise_diffs_between_ / this->working_seq_len_;
}

float SequenceData::get_mean_pdist_total() const {
    return this->avg_pairwise_diffs_total_ / this->working_seq_len_;
}

float SequenceData::get_mean_pdist_net() const {
    return this->avg_pairwise_diffs_net_ / this->working_seq_len_;
}

unsigned int SequenceData::get_num_singleton_sites() const {
    return this->singleton_sites_.size();
}

unsigned int SequenceData::get_num_mutations() const {
    return this->num_mutations_;
}

float SequenceData::get_wattersons_theta() const {
    return this->wattersons_theta_;
}

float SequenceData::get_tajimas_d() const {
    return this->tajimas_d_;
}

float SequenceData::get_fst(unsigned int method) {
    float fst = 0;
    if (method == 0) {
        fst = this->avg_pairwise_diffs_net_ / this->avg_pairwise_diffs_between_;
    } else if (method == 1) {
        fst = this->avg_pairwise_diffs_net_ / (2 * this->avg_pairwise_diffs_within_ + this->avg_pairwise_diffs_net_);
    } else if (method == 2) {
        fst = 1 - this->avg_pairwise_diffs_within_ / this->avg_pairwise_diffs_total_;
    } else {
        std::ostringstream o;
        o << "Unrecognized Fst method: " << method;
        throw std::runtime_error(o.str());
    }
    return fst;
}

float SequenceData::get_subpopulation_fst(unsigned int i, unsigned int j) {
    float fst = 0;
    fst = this->get_pairwise_diffs_net(i, j) / this->get_pairwise_diffs_between(i, j);
    return fst;
}

float SequenceData::get_pairwise_diffs_within(unsigned int i) const {
    std::map< unsigned int, float >::const_iterator xi = this->subpop_diffs_within_.find(i);
    if (xi != this->subpop_diffs_within_.end()) {
        return xi->second;
    } else {
        std::ostringstream o;
        o << "Statistics for subpopulation " << i << " not available";
        throw std::runtime_error(o.str());
    }
}

float SequenceData::get_pairwise_diffs_between(unsigned int i, unsigned int j) const {
    std::map< std::pair<unsigned int, unsigned int>, float>::const_iterator xij = this->subpop_diffs_between_.find(std::make_pair(i, j));
    if (xij == this->subpop_diffs_between_.end()) {
        xij = this->subpop_diffs_between_.find(std::make_pair(j, i));
    }
    if (xij != this->subpop_diffs_between_.end()) {
        return xij->second;
    } else {
        std::ostringstream o;
        o << "Statistics for subpopulation <" << i << "," << j << "> not available";
        throw std::runtime_error(o.str());
    }
}

float SequenceData::get_pairwise_diffs_net(unsigned int i, unsigned int j) const {
    return this->get_pairwise_diffs_between(i, j) - (this->get_pairwise_diffs_within(i) + this->get_pairwise_diffs_within(j))/2;
}

float SequenceData::get_subpopulation_pair_avg_pairwise_diffs(unsigned int i, unsigned int j) {
    std::map< std::pair<unsigned int, unsigned int>, float>::const_iterator xij = this->subpop_pairs_avg_pairwise_diffs_.find(std::make_pair(i, j));
    if (xij == this->subpop_pairs_avg_pairwise_diffs_.end()) {
        xij = this->subpop_pairs_avg_pairwise_diffs_.find(std::make_pair(j, i));
    }
    if (xij != this->subpop_pairs_avg_pairwise_diffs_.end()) {
        return xij->second;
    } else {
        std::ostringstream o;
        o << "Pairwise differences for subpopulation pair <" << i << "," << j << "> not available";
        throw std::runtime_error(o.str());
    }
}

float SequenceData::get_subpopulation_pair_wattersons_theta(unsigned int i, unsigned int j) {
    std::map< std::pair<unsigned int, unsigned int>, float>::const_iterator xij = this->subpop_pairs_theta_.find(std::make_pair(i, j));
    if (xij == this->subpop_pairs_theta_.end()) {
        xij = this->subpop_pairs_theta_.find(std::make_pair(j, i));
    }
    if (xij != this->subpop_pairs_theta_.end()) {
        return xij->second;
    } else {
        std::ostringstream o;
        o << "Watterson's theta for subpopulation pair <" << i << "," << j << "> not available";
        throw std::runtime_error(o.str());
    }
}

float SequenceData::get_subpopulation_pair_tajimas_d(unsigned int i, unsigned int j) {
    std::map< std::pair<unsigned int, unsigned int>, float>::const_iterator xij = this->subpop_pairs_tajimas_d_.find(std::make_pair(i, j));
    if (xij == this->subpop_pairs_tajimas_d_.end()) {
        xij = this->subpop_pairs_tajimas_d_.find(std::make_pair(j, i));
    }
    if (xij != this->subpop_pairs_tajimas_d_.end()) {
        return xij->second;
    } else {
        std::ostringstream o;
        o << "Tajima's D for subpopulation pair <" << i << "," << j << "> not available";
        throw std::runtime_error(o.str());
    }
}

// 2}}}

// counting/stats (support/private API) {{{2
// -----------------------------------------------------------------------------

float SequenceData::sum_diffs_between_(const std::vector<unsigned int>& pop1, const std::vector<unsigned int>& pop2) const {
    unsigned int total_diffs = 0;
    unsigned int total_comparisons = 0;
    for (std::vector<unsigned int>::const_iterator i = pop1.begin(); i != pop1.end(); ++i) {
        for (std::vector<unsigned int>::const_iterator j = pop2.begin(); j != pop2.end(); ++j) {
            if (*i != *j) {
                ++total_comparisons;
                total_diffs += this->pairwise_diffs_(*i, *j);
            }
        }
    }
    return static_cast<float>(total_diffs)/total_comparisons;
}

void SequenceData::calc_diversity_general_() {
    this->avg_pairwise_diffs_ = 0.0;
    unsigned int num_pairwise_diffs_ = 0;
    unsigned int compares = 0;
    for (unsigned int i = 0; i < this->taxon_labels_.size() - 1; ++i) {
        for (unsigned int j = i + 1; j < this->taxon_labels_.size(); ++j) {
            num_pairwise_diffs_ += this->pairwise_diffs_(i, j);
            ++compares;
        }
    }
    this->avg_pairwise_diffs_ = static_cast<float>(num_pairwise_diffs_)/compares;
}

void SequenceData::calc_diversity_subpopulations_() {
    float num_pairwise_diffs_within = 0.0;
    this->avg_pairwise_diffs_within_ = 0.0;
    unsigned int total_comparisons_within = 0;
    for (unsigned int i = 0; i < this->subpopulation_taxon_idxs_.size(); ++i) {
        float nd = this->sum_diffs_between_(this->subpopulation_taxon_idxs_[i], this->subpopulation_taxon_idxs_[i]);
        this->subpop_diffs_within_[i] = nd;
        num_pairwise_diffs_within += nd;
        ++total_comparisons_within;
    }
    this->avg_pairwise_diffs_within_ = num_pairwise_diffs_within / total_comparisons_within;

    float num_pairwise_diffs_between = 0.0;
    this->avg_pairwise_diffs_between_ = 0.0;
    unsigned int total_comparisons_between = 0;
    for (unsigned int i = 0; i < this->subpopulation_taxon_idxs_.size() - 1; ++i) {
        for (unsigned int j = i + 1; j < this->subpopulation_taxon_idxs_.size() ; ++j) {
            float nd = this->sum_diffs_between_(this->subpopulation_taxon_idxs_[i], this->subpopulation_taxon_idxs_[j]);
            num_pairwise_diffs_between += nd;
            this->subpop_diffs_between_[std::make_pair(i, j)] = nd;
            this->subpop_diffs_between_[std::make_pair(j, i)] = nd;
            ++total_comparisons_between;
        }
    }
    this->avg_pairwise_diffs_between_ = num_pairwise_diffs_between / total_comparisons_between;
    this->avg_pairwise_diffs_total_ = (this->avg_pairwise_diffs_within_*total_comparisons_within + 2*this->avg_pairwise_diffs_between_ * total_comparisons_between) / (total_comparisons_within + 2 * total_comparisons_between);
    this->avg_pairwise_diffs_net_ = this->avg_pairwise_diffs_between_ - this->avg_pairwise_diffs_within_;
}

void SequenceData::calc_summary_stats_general_() {

    this->segregating_sites_.clear();
    this->singleton_sites_.clear();
    std::map<int, unsigned int> site_states;
    for (unsigned int site = 0; site < this->working_seq_len_; ++site) {
        site_states.clear();
        for (unsigned int seq = 0; seq < this->state_matrix_.size(); ++seq) {
            site_states[this->state_matrix_[seq][site]] += 1;
        }
        if (site_states.size() > 1) {
            this->segregating_sites_.push_back(site);
        }
        this->num_mutations_ += (site_states.size() - 1);
        if (site_states.size() > 1) {
            for (std::map<int, unsigned int>::const_iterator sti = site_states.begin();
                    sti != site_states.end();
                    ++sti) {
                if (sti->second == 1) {
                    this->singleton_sites_.push_back(site);
                    break;
                }
            }
        }
    }

    float a1 = 0.0;
    float a2 = 0.0;
    this->calc_a1_a2_(this->state_matrix_.size(), a1, a2);

    unsigned int S = 0;
    if (this->use_num_mutations_instead_of_seg_sites_) {
        S = this->num_mutations_;
    } else {
        S = this->segregating_sites_.size();
    }
    this->wattersons_theta_ = static_cast<float>(S) / a1;
    this->tajimas_d_ = this->calc_tajimas_d_(this->state_matrix_.size(), this->avg_pairwise_diffs_, S, a1, a2);

}

void SequenceData::calc_summary_stats_subpopulations_() {

}

void SequenceData::calc_summary_stats_subpopulation_pairs_() {
    this->subpop_pairs_avg_pairwise_diffs_.clear();
    this->subpop_pairs_num_segregating_sites_.clear();
    this->subpop_pairs_num_mutations_.clear();
    for (unsigned int pop_i = 0; pop_i < this->subpopulation_taxon_idxs_.size() - 1; ++pop_i) {
        for (unsigned int pop_j = pop_i + 1; pop_j < this->subpopulation_taxon_idxs_.size() ; ++pop_j) {
            std::pair<unsigned int, unsigned int> subpop_pair_idx = std::make_pair(pop_i, pop_j);
            std::pair<unsigned int, unsigned int> subpop_pair_idx_rev = std::make_pair(pop_j, pop_i);
            this->subpop_pairs_num_segregating_sites_[subpop_pair_idx] = 0;
            this->subpop_pairs_num_mutations_[subpop_pair_idx] = 0;
            std::vector<unsigned int> joined_taxon_idxs(this->subpopulation_taxon_idxs_[pop_i]);
            joined_taxon_idxs.reserve(joined_taxon_idxs.size() + this->subpopulation_taxon_idxs_[pop_j].size());
            std::copy(this->subpopulation_taxon_idxs_[pop_j].begin(), this->subpopulation_taxon_idxs_[pop_j].end(), std::back_inserter(joined_taxon_idxs));

            // calculate num segregating sites etc.
            std::set<int> site_states;
            for (unsigned int site = 0; site < this->working_seq_len_; ++site) {
                site_states.clear();
                for (unsigned int seq = 0; seq < joined_taxon_idxs.size(); ++seq) {
                    site_states.insert(this->state_matrix_[joined_taxon_idxs[seq]][site]);
                    // if (site_states.size() > 1) {
                    //     ++this->subpop_pairs_num_segregating_sites_;
                    //     break;
                    // }
                }
                if (site_states.size() > 1) {
                    ++this->subpop_pairs_num_segregating_sites_[subpop_pair_idx];
                    this->subpop_pairs_num_mutations_[subpop_pair_idx] += site_states.size() - 1;
                }
            }

            float avg_diffs_ij = static_cast<float>((static_cast<float>(this->get_pairwise_diffs_within(pop_i) * this->subpopulation_taxon_idxs_[pop_i].size() * (this->subpopulation_taxon_idxs_[pop_i].size() - 1)) / 2)
                                    + (static_cast<float>(this->get_pairwise_diffs_within(pop_j) * this->subpopulation_taxon_idxs_[pop_j].size() * (this->subpopulation_taxon_idxs_[pop_j].size() - 1)) / 2)
                                    + (this->get_pairwise_diffs_between(pop_i,pop_j) * this->subpopulation_taxon_idxs_[pop_i].size() * this->subpopulation_taxon_idxs_[pop_j].size()) ) / binomial_coefficient(joined_taxon_idxs.size(), 2);
            this->subpop_pairs_avg_pairwise_diffs_[subpop_pair_idx] = avg_diffs_ij;
            this->subpop_pairs_avg_pairwise_diffs_[subpop_pair_idx_rev] = avg_diffs_ij;
            float a1 = 0.0;
            float a2 = 0.0;
            this->calc_a1_a2_(joined_taxon_idxs.size(), a1, a2);
            unsigned int S = 0;
            this->subpop_pairs_num_segregating_sites_[subpop_pair_idx_rev] = this->subpop_pairs_num_segregating_sites_[subpop_pair_idx];
            this->subpop_pairs_num_mutations_[subpop_pair_idx_rev] = this->subpop_pairs_num_mutations_[subpop_pair_idx];
            if (this->use_num_mutations_instead_of_seg_sites_) {
                S = this->subpop_pairs_num_mutations_[subpop_pair_idx];
            } else {
                S = this->subpop_pairs_num_segregating_sites_[subpop_pair_idx];
            }
            this->subpop_pairs_theta_[subpop_pair_idx] = static_cast<float>(S) / a1;
            this->subpop_pairs_theta_[subpop_pair_idx_rev] = this->subpop_pairs_theta_[subpop_pair_idx];
            this->subpop_pairs_tajimas_d_[subpop_pair_idx] = this->calc_tajimas_d_(
                    joined_taxon_idxs.size(),
                    avg_diffs_ij,
                    S,
                    a1,
                    a2);
            this->subpop_pairs_tajimas_d_[subpop_pair_idx_rev] = this->subpop_pairs_tajimas_d_[subpop_pair_idx];
        }
    }
}

float SequenceData::calc_tajimas_d_(unsigned int num_sequences, float num_diffs, unsigned int S, float a1, float a2) {
    // ### VERIFICATION ###
    // ###
    // ### Given: num_sequences = 10, num_pairwise_differences = 3.888889, S = 16
    // ###  i.e.: tajimas_d(10, 3.888889, 16)  == -1.44617198561
    // ###  Then:    a1 == 2.82896825397
    // ###           a2 == 1.53976773117
    // ###           b1 == 0.407407407407
    // ###           b2 == 0.279012345679
    // ###           c1 == 0.0539216450284
    // ###           c2 == 0.0472267720013
    // ###           e1 == 0.0190605338016
    // ###           e2 == 0.0049489277699
    // ###           D ==  -1.44617198561

    float b1 = static_cast<float>(num_sequences+1)/(3*(num_sequences-1));
    float b2 = static_cast<float>(2 * ( (num_sequences * num_sequences) + num_sequences + 3)) / (9 * num_sequences * (num_sequences-1));
    float c1 = b1 - 1.0/a1;
    float c2 = b2 - static_cast<float>(num_sequences+2)/(a1 * num_sequences) + static_cast<float>(a2)/(a1 * a1);
    float e1 = c1/a1;
    float e2 = c2 / ( (a1*a1) + a2 );
    float denominator2 = (e1 * S ) + ((e2 * S) * (S - 1) );
    float denominator = sqrt(denominator2);
    float numerator = static_cast<float>(num_diffs - (static_cast<float>(S)/a1));
    float D = numerator / denominator;
//
//     if (true) {
//         std::cout << "=== Tajima's D Calculations ===" << std::endl;
//         std::cout << "S  = " << S << std::endl;
//         std::cout << "a1 = " << a1 << std::endl;
//         std::cout << "a2 = " << a2 << std::endl;
//         std::cout << "b1 = " << b1 << std::endl;
//         std::cout << "b2 = " << b2 << std::endl;
//         std::cout << "c1 = " << c1 << std::endl;
//         std::cout << "c2 = " << c2 << std::endl;
//         std::cout << "e1 = " << e1 << std::endl;
//         std::cout << "e2 = " << e2 << std::endl;
//         std::cout << "-------------------------------" << std::endl;
//         std::cout << "denominator^2 = " << denominator2 << std::endl;
//         std::cout << "denominator   = " << denominator << std::endl;
//         std::cout << "numerator     = " << numerator << std::endl;
//         std::cout << "===============================" << std::endl;
//     }

    return D;

}

void SequenceData::calc_a1_a2_(unsigned int num_sequences, float& a1, float& a2) {
    a1 = 0;
    a2 = 0;
    for (unsigned int i = 1; i < num_sequences; ++i) {
        a1 += 1.0/i;
        a2 += 1.0/(i*i);
    }
}

// 2}}}

// 1}}}

} // namespace possum
