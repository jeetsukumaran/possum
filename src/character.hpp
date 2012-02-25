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

#ifndef MOPOGEN_CHARACTER_H
#define MOPOGEN_CHARACTER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include "config.h"
#include "ncl/ncl.h"

namespace possum {

// low-level data types
typedef int CharacterStateType;
typedef std::vector<CharacterStateType> CharacterStateVectorType;
typedef std::map<int, CharacterStateVectorType> TaxonStateVectorMapType;

// Support Functions {{{1
// =============================================================================

std::string compose_subpopulation_label(const std::string& stat_name, unsigned int i, unsigned int j);
std::string compose_subpopulation_label(const std::string& stat_name, unsigned int i);

// 1}}}

// PairwiseDifferences {{{1
// =============================================================================

class PairwiseDifferences {

    public:
        PairwiseDifferences() { }
        ~PairwiseDifferences() { }
        bool empty() {
            return this->diff_counts_.empty();
        }
        void clear() {
            this->diff_counts_.clear();
        }
        void set_num_diffs(unsigned int i, unsigned int j, unsigned int num_diffs);
        unsigned int operator()(unsigned int i, unsigned int j) const;

    private:
        std::map<std::pair<unsigned int, unsigned int>, unsigned int>   diff_counts_;

}; // PairwiseDifferences

// 1}}}

// CharacterStateMatrix {{{1
// =============================================================================

class CharacterStateMatrix {

    public:

        // types {{{2
        // -----------------------------------------------------------------------------

        typedef std::vector<CharacterStateVectorType>::size_type                    size_type;
        typedef std::vector<CharacterStateVectorType>::iterator                     iterator;
        typedef std::vector<CharacterStateVectorType>::const_iterator               const_iterator;
        typedef std::vector<CharacterStateVectorType>::reverse_iterator             reverse_iterator;
        typedef std::vector<CharacterStateVectorType>::const_reverse_iterator       const_reverse_iterator;
        typedef std::vector<CharacterStateVectorType>::allocator_type               allocator_type;
        typedef CharacterStateVectorType::size_type                                 subsize_type;

        // 2}}}

        // lifecycle {{{2
        // -----------------------------------------------------------------------------

        CharacterStateMatrix();
        ~CharacterStateMatrix();

        // 2}}}

        // metrics {{{2
        // -----------------------------------------------------------------------------

        bool empty() const;
        size_type size() const;
        size_type capacity() const;
        size_type max_size() const;
        unsigned int seq_len() const {
            if (!this->state_vectors_.empty()) {
                return this->state_vectors_[0].size();
            } else {
                return 0;
            }
        }

        // 2}}}

        // accessors {{{2
        // -----------------------------------------------------------------------------

        CharacterStateVectorType& at(size_type n);
        const CharacterStateVectorType& at(size_type n) const;
        CharacterStateVectorType& operator[](size_type n);
        const CharacterStateVectorType& operator[](size_type n) const;
        CharacterStateVectorType& back();
        const CharacterStateVectorType& back() const;
        CharacterStateVectorType& front();
        const CharacterStateVectorType& front() const;

        // 2}}}

        // iterators {{{2
        // -----------------------------------------------------------------------------

        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;
        reverse_iterator rbegin();
        const_reverse_iterator rbegin() const;
        reverse_iterator rend();
        const_reverse_iterator rend() const;

        // 2}}}

        // mutators {{{2
        // -----------------------------------------------------------------------------

        void assign(size_type n, const CharacterStateVectorType& u);
        void clear();
        iterator erase(iterator position);
        iterator erase(iterator first, iterator last);
        void push_back(const CharacterStateVectorType& v);
        CharacterStateVectorType& new_vector();
        void pop_back();
        iterator insert(iterator position, const CharacterStateVectorType& x);
        void insert (iterator position, size_type n, const CharacterStateVectorType& x);
        template <class InputIterator> void insert(iterator position, InputIterator first, InputIterator last) {
            return this->state_vectors_.insert(position, first, last);
        }
        void reserve(size_type n);
        void resize(size_type n, CharacterStateVectorType c = CharacterStateVectorType());
        void swap(CharacterStateMatrix m);
        CharacterStateMatrix& operator=(const CharacterStateMatrix& m);
        void drop_columns(const std::vector<bool>& to_drop);

        // 2}}}

        // stats/calcs {{{2
        // -----------------------------------------------------------------------------

        PairwiseDifferences get_pairwise_diff_counts() const {
            PairwiseDifferences p;
            return this->get_pairwise_diff_counts(p);
        }
        PairwiseDifferences& get_pairwise_diff_counts(PairwiseDifferences& pdiffs) const;

        // 2}}}

    private:
        std::vector<CharacterStateVectorType>                           state_vectors_;

}; // CharacterStateMatrix

// 1}}}

// SequenceData {{{1
// =============================================================================

class SequenceData {

public:

    // lifecycle {{{2
    // -----------------------------------------------------------------------------

    /**
     * Default constructor.
     */
    SequenceData();

    /**
     * Constructor, parses file of specified format.
     * @param filepath      path to data file
     * @param fmt           format of data
     * @param quiet         whether (quiet=False) or not (quiet=True) to report progress messages
     */
    SequenceData(const char * filepath, MultiFormatReader::DataFormatType fmt, bool quiet=true);

    /**
     * Destructor.
     */
    ~SequenceData();

    // 2}}}

    // deserialization {{{2
    // -----------------------------------------------------------------------------

    /**
     * Parses file of specified format.
     * @param filepath      path to data file
     * @param fmt           format of data
     * @param quiet         whether (quiet=False) or not (quiet=True) to report progress messages
     */
    void read(const char * filepath, MultiFormatReader::DataFormatType fmt);

    /**
     * Parses file of specified format.
     * @param input         stream opened for input
     * @param fmt           format of data
     * @param quiet         whether (quiet=False) or not (quiet=True) to report progress messages
     */
    void read(std::istream& input, MultiFormatReader::DataFormatType fmt);

    // 2}}}

    // metrics {{{2
    // -----------------------------------------------------------------------------

    /**
     * Returns number of sequences.
     */
    unsigned int size(void) const {
        return this->state_matrix_.size();
    }

    // 2}}}

    // getters/setters {{{2
    // -----------------------------------------------------------------------------

    /**
     * Sets the verbosity.
     */
    void set_quiet(bool q) {
        this->quiet_ = q;
    }

    /**
     * Sets whether or not sub-population pair stats are calculated.
     */
    void set_calc_subpopulation_pair_stats(bool p) {
        this->calc_subpopulation_pair_stats_ = p;
    }

    /**
     * Returns length of sequence.
     */
    unsigned int get_seq_len() {
        return this->state_matrix_.seq_len();
    }

    /**
     * Returns length of sequence.
     */
    unsigned int get_original_seq_len() {
        return this->original_seq_len_;
    }

    /**
     * Returns length of sequence.
     */
    unsigned int get_working_seq_len() {
        return this->working_seq_len_;
    }

    /**
     * Returns number of subpopulations.
     */
    unsigned int get_num_subpops() {
        return this->subpopulation_taxon_idxs_.size();
    }

    // 2}}}

    // subpopulation management {{{2
    // -----------------------------------------------------------------------------

    /**
     * Assign subpopulations from a vector of population sizes.
     */
    void set_subpopulations(const std::vector<unsigned int>& pop_sizes);

    /**
     * For debugging.
     */
    void report_subpopulations(std::ostream& out) const;

    // 2}}}

    // operations {{{2
    // -----------------------------------------------------------------------------

    /**
     * Returns number of sequences.
     */
    void drop_gap_or_uncertain_columns(void);
    // 2}}}

    // stats (public API) {{{2
    // -----------------------------------------------------------------------------

    /**
     * Process.
     */
    void calc_stats();

    // nucleotide diversity {{{3
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    /**
     * \pi of entire sample, i.e. the number of differences
     * across all pairwise comparisons in the dataset divided by the number of
     * comparisons.
     */
    float get_avg_pairwise_diffs() const;

    /**
     * Mean within-population \pi of entire sample, i.e. the
     * mean of (the number of differences across all pairwise comparisons
     * within each subpopulation divided by the number of comparisons).
     */
    float get_avg_pairwise_diffs_within() const;

    /**
     * Mean between-population \pi of entire sample, i.e. the
     * mean of (the number of differences across all pairwise comparisons
     * between each subpopulation divided by the number of comparisons).
     */
    float get_avg_pairwise_diffs_between() const;

    /**
     * Mean total-population \pi.
     */
    float get_avg_pairwise_diffs_total() const;

    /**
     * \pi_{between} - \pi_{within}.
     */
    float get_avg_pairwise_diffs_net() const;

    // 3}}}

    // p-distances {{{3
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    /**
     * Average number of differences across all pairwise comparisons in the dataset, divided by sequence length.
     */
    float get_mean_pdist() const;

    /**
     * Average number of differences of all pairwise comparisons within each population, divided by sequence length.
     */
    float get_mean_pdist_within() const;

    /**
     * Average number of differences of all pairwise comparisons between each population, divided by sequence length.
     */
    float get_mean_pdist_between() const;

    /**
     * Mean total-population \pi, divided by sequence length.
     */
    float get_mean_pdist_total() const;

    /**
     * \pi_{between} - \pi_{within} , divided by sequence length.
     */
    float get_mean_pdist_net() const;

    // 3}}}

    // extended subpopulation {{{3
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    /**
     * \pi of subpopulation i, i.e. the number of differences across all
     * pairwise comparisons in the subpopulation divided by the number of
     * comparisons.
     */
    float get_pairwise_diffs(unsigned int i) const;

    /**
     * Mean within-population \pi of subpopulation i, i.e. the
     * mean of (the number of differences across all pairwise comparisons
     * within the subpopulation).
     */
    float get_pairwise_diffs_within(unsigned int i) const;

    /**
     * Mean between-population \pi of subpopulation i and j, i.e. the
     * the number of differences across all pairwise comparisons
     * of individuals between each subpopulation.
     */
    float get_pairwise_diffs_between(unsigned int i, unsigned int j) const;

    /**
     * \pi_{between} - \pi_{within}, for subpopulation i and j.
     */
    float get_pairwise_diffs_net(unsigned int i, unsigned int j) const;

    // 3}}}


    // other stats {{{3
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    /**
     * Number of polymorphisms/segregating sites.
     */
    unsigned int get_num_segregating_sites() const;

    /**
     * Number of sites with at least two types of nucleotides (or amino acids)
     * with, at most, one occurring multiple times
     */
    unsigned int get_num_singleton_sites() const;

    /**
     * Number of distinct states in each site - 1, summed over all sites.
     */
    unsigned int get_num_mutations() const;

    /**
     * Watterson's theta.
     */
    float get_wattersons_theta() const;

    /**
     * Tajima's D.
     */
    float get_tajimas_d() const;

    /**
     * Returns the Fixation index (Fst) of the population, as given by:
     *
     *  Method Id       Formulation                             Authority                           Source
     *  -----------     -----------                             ---------                           ------
     *  0               (pi.between - pi.within) / pi.between   Hudson (1992)                       Wikipedia [given by libsequence doc as: pi.net / (pi.within + pi.net)]
     *  1               pi.net / (2 * pi.within + pi.net)       Slatkin (1993)                      libsequence documentation
     *  2               1 - (pi.within/pi.total)                Hudson, Boos, Kaplan (1992)         libsequence documentation
     */
    float get_fst(unsigned int method=0);

    /**
     * Returns the Fst for a pairwise comparison between two subpopulations.
     */
    float get_subpopulation_fst(unsigned int i, unsigned int j);

    /**
     * Returns \pi for two populations considered together.
     */
    float get_subpopulation_pair_avg_pairwise_diffs(unsigned int i, unsigned int j);

    /**
     * Returns Watterson's theta for two populations considered together.
     */
    float get_subpopulation_pair_wattersons_theta(unsigned int i, unsigned int j);

    /**
     * Returns Tajima's D for two populations considered together.
     */
    float get_subpopulation_pair_tajimas_d(unsigned int i, unsigned int j);


    // 3}}}


    // 2}}}

    // stats (private/support) {{{2
    // -----------------------------------------------------------------------------

private:

    float sum_diffs_between_(const std::vector<unsigned int>& pop1, const std::vector<unsigned int>& pop2) const;
    void calc_diversity_general_();
    void calc_summary_stats_general_();
    void calc_diversity_subpopulations_();
    void calc_summary_stats_subpopulations_();
    void calc_summary_stats_subpopulation_pairs_();
    float calc_tajimas_d_(unsigned int num_sequences, float num_diffs, unsigned int S, float a1=0.0, float a2=0.0);
    void calc_a1_a2_(unsigned int num_sequences, float& a1, float& a2);

    // 2}}}

    // data {{{2
    // -----------------------------------------------------------------------------

private:

    std::vector<std::string>                    taxon_labels_;
    CharacterStateMatrix                        state_matrix_;
    std::vector<bool>                           col_has_gap_or_missing_;
    std::vector< std::vector<unsigned int> >    subpopulation_taxon_idxs_;
    PairwiseDifferences                         pairwise_diffs_;
    unsigned int                                original_seq_len_;
    unsigned int                                working_seq_len_;

    std::map< unsigned int, float >                                 subpop_diffs_within_;
    std::map< std::pair<unsigned int, unsigned int>, float>         subpop_diffs_between_;
    std::map< std::pair<unsigned int, unsigned int>, float>         subpop_pairs_avg_pairwise_diffs_;
    std::map< std::pair<unsigned int, unsigned int>, unsigned int>  subpop_pairs_num_segregating_sites_;
    std::map< std::pair<unsigned int, unsigned int>, unsigned int>  subpop_pairs_num_mutations_;
    std::map< std::pair<unsigned int, unsigned int>, float>         subpop_pairs_theta_;
    std::map< std::pair<unsigned int, unsigned int>, float>         subpop_pairs_tajimas_d_;


    float                                       avg_pairwise_diffs_;
    float                                       avg_pairwise_diffs_total_;
    float                                       avg_pairwise_diffs_within_;
    float                                       avg_pairwise_diffs_between_;
    float                                       avg_pairwise_diffs_net_;
    std::vector<unsigned int>                   segregating_sites_;
    std::vector<unsigned int>                   singleton_sites_;
    unsigned int                                num_mutations_;
    float                                       wattersons_theta_;
    float                                       tajimas_d_;

    bool                                        use_num_mutations_instead_of_seg_sites_;
    bool                                        calc_subpopulation_pair_stats_;
    bool                                        quiet_;


    // 2}}}

};
// 1}}}

} // namespace possum


#endif
