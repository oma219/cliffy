/*
 * File: pfp_lcp_doc.hpp
 * Description: Based heavily on Massmiliano Rossi's pfp_lcp.hpp code, this 
 *              code add the ability to build the document array profiles
 *              based on the prefix-free parse of a text.
 * Date: September 1st, 2022
 */

#ifndef _LCP_PFP_DOC_HH
#define _LCP_PFP_DOC_HH

#include <common.hpp>
#include <iostream>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
#include <gsacak.h>
}

#include <pfp.hpp>
#include <ref_builder.hpp>
#include <deque>
#include <vector>
#include <bits/stdc++.h>
#include <omp.h>
#include <immintrin.h>
#include <x86intrin.h>
#include <emmintrin.h>

// struct for LCP queue, along with method for printing (debugging)
typedef struct
{
    size_t run_num = 0;
    uint8_t bwt_ch = 0;
    size_t doc_num = 0;
    bool is_start = false;
    bool is_end = false;
    size_t lcp_with_prev_suffix = 0;
} queue_entry_t;

std::ostream& operator<<(std::ostream& os, const queue_entry_t& lcp_queue_entry) {
    return os << "[run #: " << lcp_queue_entry.run_num <<
                 ", bwt ch: " << lcp_queue_entry.bwt_ch <<
                 ", doc num: " << lcp_queue_entry.doc_num << 
                 ", is_start: " << lcp_queue_entry.is_start <<
                 ", is_end: " << lcp_queue_entry.is_end <<
                 ", lcp_with_prev: " << lcp_queue_entry.lcp_with_prev_suffix << "]";
}

class pfp_lcp{
public:

    pf_parsing& pf;
    std::vector<size_t> min_s; // Value of the minimum lcp_T in each run of BWT_T
    std::vector<size_t> pos_s; // Position of the minimum lcp_T in each run of BWT_T

    uint8_t head; // character of current BWT run
    size_t length = 0; // length of the current BWT run

    bool rle; // run-length encode the BWT
    size_t total_num_runs = 0;
    size_t NUMCOLSFORTABLE = 0; 

    /*
     * Second alternate approach: avoids usage of predecessor lcp table
     * since it is time-consuming to continue updating a table that has
     * a size of num_docs x num_characters
     */
    pfp_lcp(size_t random_param, size_t random_param2, size_t random_param3, pf_parsing &pfp_, std::string filename, RefBuilder* ref_build, bool rle_ = true) : 
                pf(pfp_),
                min_s(1, pf.n),
                pos_s(1,0),
                head(0),
                num_docs(ref_build->num_docs),
                ch_doc_counters(256, std::vector<size_t>(ref_build->num_docs, 0)),
                ch_doc_encountered(256, std::vector<bool>(ref_build->num_docs, false)),
                predecessor_max_lcp(256, std::vector<size_t>(ref_build->num_docs, ref_build->total_length)),
                queue_pos_per_tuple(256, std::vector<std::deque<size_t>>(ref_build->num_docs, std::deque<size_t>(0))),
                rle(rle_)
                // heads(1, 0)
    {        
        // Opening output files
        std::string outfile = filename + std::string(".lcp");
        if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".ssa");
        if ((ssa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".esa");
        if ((esa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        if (rle) {
            outfile = filename + std::string(".bwt.heads");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

            outfile = filename + std::string(".bwt.len");
            if ((bwt_file_len = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
        } else {
            outfile = filename + std::string(".bwt");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
        }

        outfile = filename + std::string(".sdap");
        if ((sdap_file = fopen(outfile.c_str(), "w")) == nullptr) {
            error("open() file " + outfile + " failed");
        }
        if (fwrite(&ref_build->num_docs, sizeof(size_t), 1, sdap_file) != 1)
            error("SDAP write error: number of documents");

        outfile = filename + std::string(".edap");
        if ((edap_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");
        if (fwrite(&ref_build->num_docs, sizeof(size_t), 1, edap_file) != 1)
            error("SDAP write error: number of documents");

        assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict);

        // variables for bwt/lcp/sa construction
        phrase_suffix_t curr;
        phrase_suffix_t prev;

        // variables for doc profile construction 
        uint8_t prev_bwt_ch = 0;
        size_t curr_run_num = 0;
        size_t pos = 0;
        double avg_queue_length = 0.0;
        size_t total_pos_traversed = 0;
        std::vector<size_t> curr_da_profile (ref_build->num_docs, 0);
        std::vector<bool> docs_to_collect (ref_build->num_docs, false);

        // DEBUGGING
        std::vector<double> times = {0, 0, 0, 0, 0, 0};
        std::vector<size_t> method_wins = {0, 0, 0};
        auto start = std::chrono::system_clock::now();
        auto sec = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();

        inc(curr);
        while (curr.i < pf.dict.saD.size())
        {
            // Make sure current suffix is a valid proper phrase suffix (at least w characters but not whole phrase)
            if(is_valid(curr)){

                // Compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);
                phrase_suffix_t next = curr;

                // Go through suffix array of dictionary and store all phrase ids with same suffix
                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {
                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_suffix.push_back(next);
                    }
                }

                // Hard case: phrases with different BWT characters precediing them
                int_t lcp_suffix = compute_lcp_suffix(curr, prev);

                typedef std::pair<int_t *, std::pair<int_t *, uint8_t>> pq_t;

                // using lambda to compare elements.
                auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                    return *lhs.first > *rhs.first;
                };
                
                // Merge a list of occurrences of each phrase in the BWT of the parse
                std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
                for (auto s: same_suffix)
                {
                    size_t begin = pf.pars.select_ilist_s(s.phrase + 1);
                    size_t end = pf.pars.select_ilist_s(s.phrase + 2);
                    pq.push({&pf.pars.ilist[begin], {&pf.pars.ilist[end], s.bwt_char}});
                }

                size_t prev_occ;
                bool first = true;
                while (!pq.empty())
                {
                    auto curr_occ = pq.top();
                    pq.pop();

                    if (!first)
                    {
                        // Compute the minimum s_lcpP of the the current and previous occurrence of the phrase in BWT_P
                        lcp_suffix = curr.suffix_length + min_s_lcp_T(*curr_occ.first, prev_occ);
                    }
                    first = false;

                    // Update min_s
                    print_lcp(lcp_suffix, j);
                    update_ssa(curr, *curr_occ.first);
                    update_bwt(curr_occ.second.second, 1);
                    update_esa(curr, *curr_occ.first);

                    ssa = (pf.pos_T[*curr_occ.first] - curr.suffix_length) % (pf.n - pf.w + 1ULL);
                    esa = (pf.pos_T[*curr_occ.first] - curr.suffix_length) % (pf.n - pf.w + 1ULL);


                    /* Start of the DA Profiles code */
                    
                    uint8_t curr_bwt_ch = curr_occ.second.second;
                    size_t lcp_i = lcp_suffix;
                    size_t sa_i = ssa;
                    size_t doc_i = ref_build->doc_ends_rank(ssa);

                    // Determine whether current suffix is a run boundary
                    bool is_start = (pos == 0 || curr_bwt_ch != prev_bwt_ch) ? 1 : 0;
                    bool is_end = (pos == ref_build->total_length-1); // only special case, common case is below

                    // Handle scenario where the previous suffix was a end of a run (and now we know 
                    // because we see the next character). So we need to reach into queue.
                    if (pos != 0 && prev_bwt_ch != curr_bwt_ch) 
                        lcp_queue.back().is_end = 1; 
                    
                    if (is_start) {curr_run_num++;}
                    size_t pos_of_LF_i = (sa_i > 0) ? (sa_i - 1) : (ref_build->total_length-1);
                    size_t doc_of_LF_i = ref_build->doc_ends_rank(pos_of_LF_i);

                    // Add the current suffix data to LCP queue 
                    queue_entry_t curr_entry = {curr_run_num-1, curr_bwt_ch, doc_of_LF_i, is_start, is_end, lcp_i};
                    lcp_queue.push_back(curr_entry);
                    ch_doc_counters[curr_bwt_ch][doc_of_LF_i] += 1;
                    ch_doc_encountered[curr_bwt_ch][doc_of_LF_i] = true;

                    // Prepare variables to use during the lcp queue traversal
                    size_t min_lcp = lcp_i; // this is lcp with previous suffix
                    std::fill(docs_to_collect.begin(), docs_to_collect.end(), false);
                    docs_to_collect[doc_of_LF_i] = true;

                    // Re-initialize doc profiles and max lcp with current document (itself)
                    std::fill(curr_da_profile.begin(), curr_da_profile.end(), 0);
                    curr_da_profile[doc_of_LF_i] = ref_build->total_length - pos_of_LF_i;

                    // lambda to check if we haven't found a certain document
                    auto all = [&](std::vector<bool> docs_collected) {
                            bool found_all = true;
                            for (auto elem: docs_collected)
                                found_all &= elem;
                            return found_all;
                    };

                    // Update the vector containing all the lcp values in queue
                    lcp_vals_in_queue.push_back(lcp_i);

                    // Update the queue_pos lists for the current <ch, doc> pair
                    queue_pos_per_tuple[curr_bwt_ch][doc_of_LF_i].push_back(pos);

                    // Put together a list of queue positions to traverse. For
                    // example for <A, 1> with 3 documents, we want to look at all
                    // <A, 0> and <A, 2> records only until the last previous occurrence
                    // of <A, 1>.
                    auto merge_lists = [&](uint8_t curr_row_ch, size_t curr_row_doc) {
                            std::vector<size_t> queue_pos_for_traversal;
                            size_t curr_pos_list_size = queue_pos_per_tuple[curr_bwt_ch][curr_row_doc].size();

                            // lower bound position is the previous occurrence of the same bwt ch and document
                            // pair since we will not have max lcp with any suffixes above that suffix.
                            size_t lower_bound = (curr_pos_list_size >= 2) ?
                                                 queue_pos_per_tuple[curr_bwt_ch][curr_row_doc][curr_pos_list_size-2] :
                                                 0;

                            // Grab the last index for each list, and the documents # with non-zero
                            // number of entries
                            std::vector<size_t> docs_with_values;
                            std::vector<int> iter_pos;
                            std::vector<bool> docs_seen(num_docs, true);
                            for (size_t k = 0; k < num_docs; k++) {
                                iter_pos.push_back(queue_pos_per_tuple[curr_row_ch][k].size()-1);
                                if (queue_pos_per_tuple[curr_row_ch][k].size()) {
                                    docs_with_values.push_back(k);
                                    docs_seen[k] = false;
                                }
                            }
                            
                            // Merge position lists until we no longer find any positions 
                            // above our lower bound position
                            bool still_values_left = true;
                            while (!all(docs_seen)) {
                                //still_values_left = false;
                                for (auto curr_doc: docs_with_values) {
                                    size_t queue_pos = 0;
                                    if (iter_pos[curr_doc] >= 0) 
                                    {
                                        queue_pos = queue_pos_per_tuple[curr_row_ch][curr_doc][iter_pos[curr_doc]];
                                        queue_pos_for_traversal.push_back(queue_pos);
                                        iter_pos[curr_doc]--;
                                        //still_values_left = true;
                                        docs_seen[curr_doc] = true;
                                    }
                                }
                            }    
                            std::sort(queue_pos_for_traversal.begin(), queue_pos_for_traversal.end(), std::greater<size_t>());
                            return queue_pos_for_traversal;
                    };
                   
                    // Determine which queue positions we need to check and 
                    // compute the lcp with.
                    std::vector<size_t> queue_pos_for_traversal = merge_lists(curr_bwt_ch, doc_of_LF_i);
                    ASSERT(queue_pos_for_traversal.front() == pos,
                           "Issue with merge_lists occurred");
                    queue_pos_for_traversal.erase(queue_pos_for_traversal.begin());


                    // Go through all the necessary predecessor suffixes to 
                    // update the profiles
                    bool passed_same_document = false;
                    size_t traversal_count = 0;
                    for (auto queue_pos: queue_pos_for_traversal) {
                        size_t true_pos = queue_pos - num_records_ejected;

                        ASSERT((lcp_queue[true_pos].bwt_ch == curr_bwt_ch), 
                               "Issue occurred in doc profile queue traversal");
                        
                        // Case 1: Cross same character and document pair
                        if (lcp_queue[true_pos].doc_num == doc_of_LF_i)
                            passed_same_document = true;
                        // Case 2: Cross same character, but a different document
                        else 
                        {
                            size_t min_lcp = *std::min_element(lcp_vals_in_queue.begin()+true_pos+1, 
                                                               lcp_vals_in_queue.end());

                            // Case 2a: Update the current DA profile since we haven't seen this document yet
                            if (!docs_to_collect[lcp_queue[true_pos].doc_num]) 
                            {
                                curr_da_profile[lcp_queue[true_pos].doc_num] = min_lcp + 1;
                                docs_to_collect[lcp_queue[true_pos].doc_num] = true;
                            }

                            // Case 2b: Update the predecessor DA profile if its run boundary
                            // and we haven't passed the <ch, doc> pair as current suffix
                            if (!passed_same_document && (lcp_queue[true_pos].is_start or lcp_queue[true_pos].is_end))
                            {
                                size_t start_pos = ref_build->num_docs * true_pos;
                                size_t curr_max_lcp = lcp_queue_profiles[start_pos + doc_of_LF_i];
                                lcp_queue_profiles[start_pos + doc_of_LF_i] = std::max(curr_max_lcp, min_lcp+1);
                            }
                        }
                        traversal_count++;
                        // We can safely exit if we have seen all other documents
                        // and we have passed an occurrence of the same <ch, document> pair
                        if (all(docs_to_collect) && passed_same_document)
                            break;
                    }
                    avg_queue_length += lcp_queue.size();
                    total_pos_traversed += traversal_count;

                    //std::cout << lcp_queue.size() << "  " << queue_pos_for_traversal.size() << std::endl;
                                        
                    // Add the current profile to the vector (should always be multiple of # of docs)
                    for (auto elem: curr_da_profile)
                        lcp_queue_profiles.push_back(elem);

                    assert(lcp_queue_profiles.size() % ref_build->num_docs == 0);
                    assert(lcp_queue_profiles.size() == (lcp_queue.size() * ref_build->num_docs));

                    // Try to trim the LCP queue and adjust the count matrix ...
                    // Method #1: non-heuristic by removing entries with multiples 
                    size_t curr_pos = 0;
                    size_t records_to_remove_method1 = 0;
                    while (curr_pos < lcp_queue.size()) {
                        uint8_t curr_ch = lcp_queue[curr_pos].bwt_ch;
                        size_t curr_doc = lcp_queue[curr_pos].doc_num;
                        assert(ch_doc_counters[curr_ch][curr_doc] >= 1);

                        size_t count = 0;
                        for (size_t i = 0; i < num_docs; i++) {
                            if (i != curr_doc && ch_doc_counters[curr_ch][i] >= 1)
                                count++;
                        }

                        if (ch_doc_counters[curr_ch][curr_doc] >= 1 && count == (num_docs-1))
                            records_to_remove_method1++;
                        // TODO: generalize this to take into account characters that only occur once
                        else if (curr_ch == EndOfDict || (curr_ch != 'A' && curr_ch != 'C'
                                && curr_ch != 'G' && curr_ch != 'T'))
                            records_to_remove_method1++;
                        else
                            break;
                        curr_pos++;
                    }

                    // Makes sure the number of records to remove is not the 
                    // full queue
                    size_t records_to_remove = std::min(records_to_remove_method1, lcp_queue.size()-1);

                    // Remove those top n elements, and update counts
                    update_lcp_queue(records_to_remove, true);
                    
                    /* End of the DA Profiles code (except for some update statements below) */

                    // Update prevs
                    prev_occ = *curr_occ.first;
                    prev_bwt_ch = curr_bwt_ch;

                    // Update pq
                    curr_occ.first++;
                    if (curr_occ.first != curr_occ.second.first)
                        pq.push(curr_occ);

                    j += 1;
                    pos += 1;
                }

                prev = same_suffix.back();
                curr = next;
            }
            else {
                inc(curr);
            }
        }

        // print last BWT char and SA sample
        print_sa();
        print_bwt();
        print_doc_profiles();

        // Close output files
        fclose(ssa_file);
        fclose(esa_file);
        fclose(bwt_file);
        fclose(lcp_file);

        if (rle)
            fclose(bwt_file_len);
    }


    /*
     * First alternate approach: Uses a predecessor lcp table to 
     * reduce the number of positions that need to be traversed 
     * in the lcp queue
     */
    pfp_lcp(bool taxcomp, bool topk, size_t num_cols, pf_parsing &pfp_, std::string filename, RefBuilder* ref_build, bool rle_ = true) : 
                pf(pfp_),
                min_s(1, pf.n),
                pos_s(1,0),
                head(0),
                NUMCOLSFORTABLE(num_cols),
                num_docs(ref_build->num_docs),
                ch_doc_counters(256, std::vector<size_t>(ref_build->num_docs, 0)),
                ch_doc_encountered(256, std::vector<bool>(ref_build->num_docs, false)),
                predecessor_max_lcp(256, std::vector<size_t>(ref_build->num_docs, ref_build->total_length)),
                queue_pos_per_tuple(256, std::vector<std::deque<size_t>>(ref_build->num_docs, std::deque<size_t>(0))),
                rle(rle_),
                use_taxcomp(taxcomp),
                use_topk(topk)
                // heads(1, 0)
    {       
        // Opening output files
        std::string outfile = filename + std::string(".lcp");
        if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".ssa");
        if ((ssa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".esa");
        if ((esa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        if (rle) {
            outfile = filename + std::string(".bwt.heads");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

            outfile = filename + std::string(".bwt.len");
            if ((bwt_file_len = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
        } else {
            outfile = filename + std::string(".bwt");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
        }


        if (taxcomp) {
            outfile = filename + std::string(".taxcomp.sdap");
            if ((sdap_tax = fopen(outfile.c_str(), "w")) == nullptr) 
                error("open() file " + outfile + " failed");
            outfile = filename + std::string(".taxcomp.edap");
            if ((edap_tax = fopen(outfile.c_str(), "w")) == nullptr) 
                error("open() file " + outfile + " failed");

            outfile = filename + std::string(".taxcomp.of.sdap");
            if ((sdap_overtax = fopen(outfile.c_str(), "w")) == nullptr) 
                error("open() file " + outfile + " failed");
            if (fwrite(&ref_build->num_docs, sizeof(size_t), 1, sdap_overtax) != 1)
                error("SDAP write error: number of documents");
            outfile = filename + std::string(".taxcomp.of.edap");
            if ((edap_overtax = fopen(outfile.c_str(), "w")) == nullptr) 
                error("open() file " + outfile + " failed");
            if (fwrite(&ref_build->num_docs, sizeof(size_t), 1, edap_overtax) != 1)
                error("EDAP write error: number of documents");

            tax_sdap_overflow_ptr += sizeof(size_t);
            tax_edap_overflow_ptr += sizeof(size_t);
        } else if (topk) {
            outfile = filename + std::string(".topk.sdap");
            if ((sdap_topk = fopen(outfile.c_str(), "w")) == nullptr) 
                error("open() file " + outfile + " failed");
            outfile = filename + std::string(".topk.edap");
            if ((edap_topk = fopen(outfile.c_str(), "w")) == nullptr) 
                error("open() file " + outfile + " failed");
        } else {
            outfile = filename + std::string(".sdap");
            if ((sdap_file = fopen(outfile.c_str(), "w")) == nullptr) 
                error("open() file " + outfile + " failed");
            if (fwrite(&ref_build->num_docs, sizeof(size_t), 1, sdap_file) != 1)
                error("SDAP write error: number of documents");

            outfile = filename + std::string(".edap");
            if ((edap_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
            if (fwrite(&ref_build->num_docs, sizeof(size_t), 1, edap_file) != 1)
                error("SDAP write error: number of documents");
        }

        assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict);

        // variables for bwt/lcp/sa construction
        phrase_suffix_t curr;
        phrase_suffix_t prev;

        // variables for doc profile construction 
        uint8_t prev_bwt_ch = 0;
        size_t curr_run_num = 0;
        size_t pos = 0;
        std::vector<size_t> curr_da_profile (ref_build->num_docs, 0);
        std::vector<bool> docs_to_collect (ref_build->num_docs, false);

        // DEBUG: variables defined for debugging
        std::vector<double> times = {0, 0, 0, 0, 0, 0};
        std::vector<size_t> method_wins = {0, 0, 0};
        double avg_queue_length = 0.0;
        size_t total_pos_traversed = 0;
        auto start = std::chrono::system_clock::now();
        auto sec = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();

        // DEBUG: output file for storing variables
        // std::ofstream debugging_fd (filename + ".output_data.txt");
        // leaveout_fd.open(filename + ".leave_one_data.txt");

        // Create a backup predecessor max lcp table, and re-initialize with max_lcp
        size_t num_blocks_of_32 = num_docs/32;
        num_blocks_of_32++;

        // Since I am using a int array, I want to make sure I don't have overflow, I am using 
        // a signed int8 since the 512 uint8t load/store were not working.
        uint16_t max_lcp_init = (ref_build->total_length > MAXLCPVALUE) ? MAXLCPVALUE : ref_build->total_length;

        uint16_t predecessor_max_lcp_2[256][num_blocks_of_32 * 32]; 
        for (size_t i = 0; i < 256; i++)
            for (size_t j = 0; j < num_blocks_of_32 * 32; j++)
                predecessor_max_lcp_2[i][j] = max_lcp_init;
        // for (size_t i = 0; i < 256; i++)
        //     for (size_t j = 0; j < num_docs; j++)
        //         predecessor_max_lcp[i][j] = max_lcp_init;

        inc(curr);
        while (curr.i < pf.dict.saD.size())
        {
            // Make sure current suffix is a valid proper phrase suffix (at least w characters but not whole phrase)
            if(is_valid(curr)){

                // Compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);
                phrase_suffix_t next = curr;

                // Go through suffix array of dictionary and store all phrase ids with same suffix
                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {
                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_suffix.push_back(next);
                    }
                }

                // Hard case: phrases with different BWT characters precediing them
                int_t lcp_suffix = compute_lcp_suffix(curr, prev);

                typedef std::pair<int_t *, std::pair<int_t *, uint8_t>> pq_t;

                // using lambda to compare elements.
                auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                    return *lhs.first > *rhs.first;
                };
                
                // Merge a list of occurrences of each phrase in the BWT of the parse
                std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
                for (auto s: same_suffix)
                {
                    size_t begin = pf.pars.select_ilist_s(s.phrase + 1);
                    size_t end = pf.pars.select_ilist_s(s.phrase + 2);
                    pq.push({&pf.pars.ilist[begin], {&pf.pars.ilist[end], s.bwt_char}});
                }

                size_t prev_occ;
                bool first = true;
                while (!pq.empty())
                {
                    auto curr_occ = pq.top();
                    pq.pop();

                    if (!first)
                    {
                        // Compute the minimum s_lcpP of the the current and previous occurrence of the phrase in BWT_P
                        lcp_suffix = curr.suffix_length + min_s_lcp_T(*curr_occ.first, prev_occ);
                    }
                    first = false;

                    // Update min_s
                    print_lcp(lcp_suffix, j);
                    update_ssa(curr, *curr_occ.first);
                    update_bwt(curr_occ.second.second, 1);
                    update_esa(curr, *curr_occ.first);

                    ssa = (pf.pos_T[*curr_occ.first] - curr.suffix_length) % (pf.n - pf.w + 1ULL);
                    esa = (pf.pos_T[*curr_occ.first] - curr.suffix_length) % (pf.n - pf.w + 1ULL);

                    /* Start of the DA Profiles code */
                    uint8_t curr_bwt_ch = curr_occ.second.second;
                    size_t lcp_i = lcp_suffix;
                    size_t sa_i = ssa;
                    size_t doc_i = ref_build->doc_ends_rank(ssa);

                    // Determine whether current suffix is a run boundary
                    bool is_start = (pos == 0 || curr_bwt_ch != prev_bwt_ch) ? 1 : 0;
                    bool is_end = (pos == ref_build->total_length-1); // only special case, common case is below

                    // Handle scenario where the previous suffix was a end of a run (and now we know 
                    // because we see the next character). So we need to reach into queue.
                    if (pos != 0 && prev_bwt_ch != curr_bwt_ch) 
                        lcp_queue.back().is_end = 1; 
                    
                    if (is_start) {curr_run_num++;}
                    size_t pos_of_LF_i = (sa_i > 0) ? (sa_i - 1) : (ref_build->total_length-1);
                    size_t doc_of_LF_i = ref_build->doc_ends_rank(pos_of_LF_i);

                    // DEBUG: creates the debugging file
                    //debugging_fd << pos << std::setw(20) << (curr_run_num-1) << std::setw(20) << curr_bwt_ch << std::setw(20) << sa_i << std::setw(20) << doc_of_LF_i << std::setw(20) << lcp_i << std::endl;

                    // Add the current suffix data to LCP queue 
                    queue_entry_t curr_entry = {curr_run_num-1, curr_bwt_ch, doc_of_LF_i, is_start, is_end, lcp_i};
                    lcp_queue.push_back(curr_entry);
                    ch_doc_counters[curr_bwt_ch][doc_of_LF_i] += 1;
                    ch_doc_encountered[curr_bwt_ch][doc_of_LF_i] = true;

                    // Prepare variables to use during the lcp queue traversal
                    size_t min_lcp = lcp_i; // this is lcp with previous suffix
                    bool passed_same_document = false;
                    std::fill(docs_to_collect.begin(), docs_to_collect.end(), false);
                    docs_to_collect[doc_of_LF_i] = true;

                    // Re-initialize doc profiles and max lcp with current document (itself)
                    std::fill(curr_da_profile.begin(), curr_da_profile.end(), 0);
                    curr_da_profile[doc_of_LF_i] = ref_build->total_length - pos_of_LF_i;

                    // lambda to check if we haven't found a certain document
                    auto all = [&](std::vector<bool> docs_collected) {
                            bool found_all = true;
                            for (auto elem: docs_collected)
                                found_all &= elem;
                            return found_all;
                    };

                    // Update the vector containing all the lcp values in queue
                    lcp_vals_in_queue.push_back(lcp_i);

                    // Update the queue_pos lists for the current <ch, doc> pair
                    queue_pos_per_tuple[curr_bwt_ch][doc_of_LF_i].push_back(pos);

                    // Update the predecessor max lcp structure with the current lcp
                    // so basiscally iterate through all values and take the min
                    // IMPORTANT: this is the old code before SIMDifying.
                    
                    // for (size_t ch_num = 0; ch_num < 256; ch_num++) {
                    //     for (size_t doc_num = 0; doc_num < num_docs; doc_num++) {
                    //         predecessor_max_lcp[ch_num][doc_num] = std::min(predecessor_max_lcp[ch_num][doc_num], std::min(lcp_i, (size_t) 127));
                    //     }
                    // }
                    // predecessor_max_lcp[curr_bwt_ch][doc_of_LF_i] = std::min((size_t) 127, ref_build->total_length - pos_of_LF_i);
                    
                    /* Start of SIMD code */ 

                    // Update the predecessor lcp table, this allows us to compute the maximum
                    // lcp with respect to all the predecessor occurrences of other documents.
                    // For example: if we are at <A, 2>, we will find the maximum lcp with the 
                    // the previous occurrences of <A, 0> and <A, 1> for 3 documents.
                    #if AVX512BW_PRESENT       
                        /*
                         * NOTE: I decided to use the mask store/load because the non-mask
                         * version of the same function was not present. I looked online and 
                         * found this thread, that describes how you can replicate their function
                         * by using masks which is what I did: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95483
                         * This link was helpful in understanding conditional
                         * load masks: https://www.cs.virginia.edu/~cr4bd/3330/F2018/simdref.html
                         */

                        // initialize an lcp_i vector
                        uint16_t lcp_i_vector[32];
                        for (size_t i = 0; i < 32; i++)
                            lcp_i_vector[i] = std::min((size_t) MAXLCPVALUE, lcp_i);

                        // load the array of constants (lcp_i)
                        __m512i arr1, arr2, arr3; 
                        __mmask32 k = ~0; // all 32 bits on, means all 32 values will be written
                        arr2 = _mm512_maskz_loadu_epi16(~0, (const __m512i*) &lcp_i_vector[0]);

                        std::vector<size_t> dna_chars = {65, 67, 71, 84, 85}; // A, C, G, T, U
                        for (size_t ch_num: dna_chars) { // Optimization for DNA
                        //for (size_t ch_num = 0; ch_num < 256; ch_num++) {
                            // use SIMD for all groups of 32
                            for (size_t i = 0; i < (num_blocks_of_32 * 32); i+=32) {
                                // zero-mask, all the set bit positions are loaded
                                arr1 = _mm512_maskz_loadu_epi16(~0, (const __m512i*) &predecessor_max_lcp_2[ch_num][i]); 

                                arr3 = _mm512_min_epu16(arr1, arr2);
                                _mm512_mask_storeu_epi16((__m512i*) &predecessor_max_lcp_2[ch_num][i], k, arr3); 
                            }
                        }
                        // Reset the LCP with respect to the current <ch, doc> pair
                        predecessor_max_lcp_2[curr_bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, ref_build->total_length - pos_of_LF_i);
                    
                    #else                        
                        for (size_t ch_num = 0; ch_num < 256; ch_num++) {
                            for (size_t doc_num = 0; doc_num < num_docs; doc_num++) {
                                predecessor_max_lcp[ch_num][doc_num] = std::min(predecessor_max_lcp[ch_num][doc_num], std::min(lcp_i, (size_t) MAXLCPVALUE));
                            }
                        }
                        // Reset the LCP with respect to the current <ch, doc> pair
                        predecessor_max_lcp[curr_bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, ref_build->total_length - pos_of_LF_i);
                    #endif

                    /* End of SIMD changes */

                    // CHECK CODE: Makes sure the predecessor table is correct
                    // for (size_t i = 0; i < 256; i++) {
                    //     for (size_t j = 0; j < num_docs; j++) {
                            
                    //         if (predecessor_max_lcp[i][j] != (size_t) predecessor_max_lcp_2[i][j]) {
                    //             std::cout << "\nwe found a discrepancy!!!" << std::endl;

                    //             std::cout << "predecessor_max_lcp[" << i << "][" << j << "]\n";
                    //             std::cout << predecessor_max_lcp[i][j] << std:: endl;
                    //             std::cout << unsigned(predecessor_max_lcp_2[i][j]) << std:: endl;
                    //         }
                    //     }
                    // }

                    // Initialize the curr_da_profile with max lcp for predecessor 
                    // occurrences of the same BWT character from another document, and
                    // we check and make sure they occurred to avoid initializing it
                    // with 1 (0 + 1 = 1)
                    for (size_t i = 0; i < num_docs; i++) {
                        if (i != doc_of_LF_i && ch_doc_encountered[curr_bwt_ch][i])
                            curr_da_profile[i] = predecessor_max_lcp_2[curr_bwt_ch][i] + 1; //predecessor_max_lcp[curr_bwt_ch][i] + 1;
                    }

                    // Put together a list of queue positions to traverse. For
                    // example for <A, 1> with 3 documents, we want to look at all
                    // <A, 0> and <A, 2> records only until the last previous occurrence
                    // of <A, 1>.
                    auto merge_lists = [&](uint8_t curr_row_ch, size_t curr_row_doc) {
                            std::vector<size_t> queue_pos_for_traversal;
                            size_t curr_pos_list_size = queue_pos_per_tuple[curr_bwt_ch][curr_row_doc].size();

                            // lower bound position is the previous occurrence of the same bwt ch and document
                            // pair since we will not have max lcp with any suffixes above that suffix.
                            size_t lower_bound = (curr_pos_list_size >= 2) ?
                                                 queue_pos_per_tuple[curr_bwt_ch][curr_row_doc][curr_pos_list_size-2] :
                                                 0;

                            // Grab the last index for each list, and the documents # with non-zero
                            // number of entries
                            std::vector<size_t> docs_with_values;
                            std::vector<int> iter_pos;
                            for (size_t k = 0; k < num_docs; k++) {
                                iter_pos.push_back(queue_pos_per_tuple[curr_row_ch][k].size()-1);
                                if (queue_pos_per_tuple[curr_row_ch][k].size() && k != curr_row_doc)
                                    docs_with_values.push_back(k);
                            }
                            
                            // Merge position lists until we no longer find any positions 
                            // above our lower bound position
                            bool still_values_left = true;
                            while (still_values_left) {
                                still_values_left = false;
                                for (auto curr_doc: docs_with_values) {
                                    size_t queue_pos = 0;
                                    if (iter_pos[curr_doc] >= 0 && 
                                        (queue_pos = queue_pos_per_tuple[curr_row_ch][curr_doc][iter_pos[curr_doc]])>= lower_bound) 
                                    {
                                        queue_pos_for_traversal.push_back(queue_pos);
                                        iter_pos[curr_doc]--;
                                        still_values_left = true;
                                    }
                                }
                            }    
                            return queue_pos_for_traversal;
                    };

                    // Determine which queue positions we need to check and 
                    // compute the lcp with.
                    std::vector<size_t> queue_pos_for_traversal = merge_lists(curr_bwt_ch, doc_of_LF_i);

                    // Go through all the necessary predecessor suffixes to 
                    // update the profiles
                    for (auto queue_pos: queue_pos_for_traversal) {
                        size_t true_pos = queue_pos - num_records_ejected;

                        ASSERT((lcp_queue[true_pos].bwt_ch == curr_bwt_ch), 
                        "Issue occurred in doc profile construction");
                        ASSERT((lcp_queue[true_pos].doc_num != doc_of_LF_i), 
                        "Issue occurred in doc profile construction");

                        size_t start_pos = ref_build->num_docs * true_pos;
                        size_t curr_max_lcp = lcp_queue_profiles[start_pos + doc_of_LF_i];
                        size_t min_lcp = *std::min_element(lcp_vals_in_queue.begin()+true_pos+1, 
                                                           lcp_vals_in_queue.end());
                        
                        lcp_queue_profiles[start_pos + doc_of_LF_i] = std::max(curr_max_lcp, min_lcp+1);
                    }
                    avg_queue_length += lcp_queue.size();
                    total_pos_traversed += queue_pos_for_traversal.size();

                    if (lcp_queue.size() >= MAXQUEUELENGTH)
                        FATAL_ERROR("queue length during construction grew too large");
                    
                    
                    // Add the current profile to the vector (should always be multiple of # of docs)
                    for (auto elem: curr_da_profile)
                        lcp_queue_profiles.push_back(elem);

                    assert(lcp_queue_profiles.size() % ref_build->num_docs == 0);
                    assert(lcp_queue_profiles.size() == (lcp_queue.size() * ref_build->num_docs));


                    /* Try to trim the LCP queue and adjust the count matrix. Below are a 
                     * list of methods to trim the queue ordered based on how optimal they
                     * in terms of reducing the size of the queue.
                     *
                     *     Heuristic Methods:
                     *         #1 - Removes all entries above a small lcp value (reduces queue to 1)
                     *         #2 - Removes all entries if queue is filled with non-relevant characters
                     *              like Ns in DNA database (reduces queue to 1)
                     *         #3 - Fail safe -> keep queue under a certain length at all times (reduces
                     *                           the queue to that certain length)
                     *
                     *     Non-heuristic Methods:  
                     *         #4 - Removes entries in the queue whose document array profile 
                     *              has been finalized (reduces the queue by 0 <= x <= size-1) entries
                     *
                     *     IMPORTANT: We need to know if we used a heuristic or non-heuristic method when
                     *                we go to trim the queue in order to know if we should adjust the 
                     *                count matrix.
                     * 
                     */
                    size_t records_to_remove_heuristic = 0, records_to_remove_non_heuristic = 0;
                    bool heuristic_used = false, non_heuristic_used = false;

                    // Method #1 - Heuristic
                    if (lcp_i <= 8) {
                        heuristic_used = true;
                        records_to_remove_heuristic = lcp_queue.size()-1;
                    // Method #2 - Heuristic 
                    } else if (lcp_queue.size() % 1000 == 0 && lcp_queue.size() >= 1000) {
                        size_t non_usual_chars = 0;
                        auto is_usual_char = [] (auto ch) {return (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'U');};

                        for (auto entry: lcp_queue) {
                            if (!is_usual_char(entry.bwt_ch))
                                non_usual_chars++;
                        }

                        if ((non_usual_chars+0.0)/lcp_queue.size() > 0.90) {
                            records_to_remove_heuristic = lcp_queue.size()-1;
                            heuristic_used = true;
                        }
                    // Method #3 - Heuristic
                    } else if (lcp_queue.size() > 2500) {
                        records_to_remove_heuristic = lcp_queue.size() - 2500;
                        heuristic_used = true;
                    } 

                    /* 
                     * If using heuristics did not yield flushing the queue, let's
                     * try to see if using a non-heuristic method will do better.
                     */
                    // Method #4 - Non-heuristic
                    if (pos % 50 == 0 && records_to_remove_heuristic != lcp_queue.size() - 1) {
                        size_t curr_pos = 0;
                        bool non_heuristic_done = false;
                        non_heuristic_used = true;
                        while (curr_pos < lcp_queue.size() && !non_heuristic_done) {
                            uint8_t curr_ch = lcp_queue[curr_pos].bwt_ch;
                            size_t curr_doc = lcp_queue[curr_pos].doc_num;
                            assert(ch_doc_counters[curr_ch][curr_doc] >= 1);

                            size_t current_run_num = lcp_queue[curr_pos].run_num;
                            size_t count = 0;

                            for (size_t i = 0; i < num_docs; i++) {
                                if (i != curr_doc && ch_doc_counters[curr_ch][i] >= 1)
                                    count++;
                            }
                           /*
                            * IMPORTANT: I added the decrement statements below since we need
                            * an updated ch_doc_counters in order to make correct decisions. The
                            * problem was that sometimes it would remove too many entries
                            */
                            if (count == (num_docs-1)) {
                                records_to_remove_non_heuristic++;
                                ch_doc_counters[curr_ch][curr_doc] -= 1;
                            // TODO: generalize this to take into account characters that only occur once
                            } else if (curr_ch == EndOfDict || (curr_ch != 'A' && curr_ch != 'C'
                                    && curr_ch != 'G' && curr_ch != 'T' && curr_ch != 'U')) {
                                records_to_remove_non_heuristic++;
                                ch_doc_counters[curr_ch][curr_doc] -= 1;
                            } else {
                                non_heuristic_done = true;
                            }
                            curr_pos++;
                        }
                    }

                    // Take the maximum value from the two methods above, and 
                    // depending on which method is used. Make sure to update
                    // the character/doc table.
                    if (records_to_remove_heuristic > records_to_remove_non_heuristic) 
                        update_lcp_queue(std::min(records_to_remove_heuristic, lcp_queue.size()-1), true);
                    else if (records_to_remove_non_heuristic > 0)
                        update_lcp_queue(std::min(records_to_remove_non_heuristic, lcp_queue.size()-1), false);

                    // if (curr_run_num == 10) {
                    //     std::cout << "\n\nFINISHED\npos = " << pos << std::endl;
                    //     std::exit(1);
                    // }
                    
                    // std::cout << tax_sdap_overflow_ptr << " " << tax_edap_overflow_ptr << std::endl;
                    
                    /* End of the DA Profiles code (except for some update statements below) */

                    // Update prevs
                    prev_occ = *curr_occ.first;
                    prev_bwt_ch = curr_bwt_ch;

                    // Update pq
                    curr_occ.first++;
                    if (curr_occ.first != curr_occ.second.first)
                        pq.push(curr_occ);

                    j += 1;
                    pos += 1;
                }

                prev = same_suffix.back();
                curr = next;
            }
            else {
                inc(curr);
            }
        }

        // print last BWT char and SA sample
        print_sa();
        print_bwt();
        print_doc_profiles();
        total_num_runs = curr_run_num;

        // Close output files
        fclose(ssa_file); fclose(esa_file);
        fclose(bwt_file);
        fclose(lcp_file);

        if (rle)
            fclose(bwt_file_len);

        if (taxcomp) {
            fclose(sdap_tax); fclose(edap_tax);
            fclose(sdap_overtax); fclose(edap_overtax);
        } else if (topk) {
            fclose(sdap_topk); fclose(edap_topk);
        } else {
            fclose(sdap_file); fclose(edap_file);
        }
    }

    // Original constructor: The main difference is that it will iterate through 
    // the ALL the queue positions to compute max LCP with predecessor and successor
    // while the new constructor, will compute predecessor using a table and iterate
    // through queue to find the max lcp with successor rows.
    pfp_lcp(pf_parsing &pfp_, std::string filename, RefBuilder* ref_build, bool rle_ = true) : 
                pf(pfp_),
                min_s(1, pf.n),
                pos_s(1,0),
                head(0),
                num_docs(ref_build->num_docs),
                ch_doc_counters(256, std::vector<size_t>(ref_build->num_docs, 0)),
                rle(rle_)
                // heads(1, 0)
    {
        // Opening output files
        std::string outfile = filename + std::string(".lcp");
        if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".ssa");
        if ((ssa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".esa");
        if ((esa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        if (rle) {
            outfile = filename + std::string(".bwt.heads");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

            outfile = filename + std::string(".bwt.len");
            if ((bwt_file_len = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
        } else {
            outfile = filename + std::string(".bwt");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
        }

        outfile = filename + std::string(".sdap");
        if ((sdap_file = fopen(outfile.c_str(), "w")) == nullptr) {
            error("open() file " + outfile + " failed");
        }
        if (fwrite(&ref_build->num_docs, sizeof(size_t), 1, sdap_file) != 1)
            error("SDAP write error: number of documents");

        outfile = filename + std::string(".edap");
        if ((edap_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");
        if (fwrite(&ref_build->num_docs, sizeof(size_t), 1, edap_file) != 1)
            error("SDAP write error: number of documents");

        assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict);

        // variables for bwt/lcp/sa construction
        phrase_suffix_t curr;
        phrase_suffix_t prev;

        // variables for doc profile construction 
        uint8_t prev_bwt_ch = 0;
        size_t curr_run_num = 0;
        size_t pos = 0;
        std::vector<size_t> curr_da_profile (ref_build->num_docs, 0);
        std::vector<bool> docs_to_collect (ref_build->num_docs, false);

        // reserve some room for the lcp queue to avoid reallocation/copying
        //lcp_queue.reserve(100);
        //lcp_queue_profiles.reserve(100 * ref_build->num_docs);

        inc(curr);
        while (curr.i < pf.dict.saD.size())
        {
            // Make sure current suffix is a valid proper phrase suffix (at least w characters but not whole phrase)
            if(is_valid(curr)){

                // Compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);
                phrase_suffix_t next = curr;

                // Go through suffix array of dictionary and store all phrase ids with same suffix
                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {
                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_suffix.push_back(next);
                    }
                }

                // Hard case: phrases with different BWT characters precediing them
                int_t lcp_suffix = compute_lcp_suffix(curr, prev);

                typedef std::pair<int_t *, std::pair<int_t *, uint8_t>> pq_t;

                // using lambda to compare elements.
                auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                    return *lhs.first > *rhs.first;
                };
                
                // Merge a list of occurrences of each phrase in the BWT of the parse
                std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
                for (auto s: same_suffix)
                {
                    size_t begin = pf.pars.select_ilist_s(s.phrase + 1);
                    size_t end = pf.pars.select_ilist_s(s.phrase + 2);
                    pq.push({&pf.pars.ilist[begin], {&pf.pars.ilist[end], s.bwt_char}});
                }

                size_t prev_occ;
                bool first = true;
                while (!pq.empty())
                {
                    auto curr_occ = pq.top();
                    pq.pop();

                    if (!first)
                    {
                        // Compute the minimum s_lcpP of the the current and previous occurrence of the phrase in BWT_P
                        lcp_suffix = curr.suffix_length + min_s_lcp_T(*curr_occ.first, prev_occ);
                    }
                    first = false;

                    // Update min_s
                    print_lcp(lcp_suffix, j);
                    update_ssa(curr, *curr_occ.first);
                    update_bwt(curr_occ.second.second, 1);
                    update_esa(curr, *curr_occ.first);

                    ssa = (pf.pos_T[*curr_occ.first] - curr.suffix_length) % (pf.n - pf.w + 1ULL);
                    esa = (pf.pos_T[*curr_occ.first] - curr.suffix_length) % (pf.n - pf.w + 1ULL);


                    /* Start of the DA Profiles code */

                    uint8_t curr_bwt_ch = curr_occ.second.second;
                    size_t lcp_i = lcp_suffix;
                    size_t sa_i = ssa;
                    size_t doc_i = ref_build->doc_ends_rank(ssa);

                    // Determine whether current suffix is a run boundary
                    bool is_start = (pos == 0 || curr_bwt_ch != prev_bwt_ch) ? 1 : 0;
                    bool is_end = (pos == ref_build->total_length-1); // only special case, common case is below

                    // Handle scenario where the previous suffix was a end of a run (and now we know 
                    // because we see the next character). So we need to reach into queue.
                    if (pos != 0 && prev_bwt_ch != curr_bwt_ch) 
                        lcp_queue.back().is_end = 1; 
                    
                    if (is_start) {curr_run_num++;}
                    size_t pos_of_LF_i = (sa_i > 0) ? (sa_i - 1) : (ref_build->total_length-1);
                    size_t doc_of_LF_i = ref_build->doc_ends_rank(pos_of_LF_i);

                    // Add the current suffix data to LCP queue 
                    queue_entry_t curr_entry = {curr_run_num-1, curr_bwt_ch, doc_of_LF_i, is_start, is_end, lcp_i};
                    lcp_queue.push_back(curr_entry);
                    ch_doc_counters[curr_bwt_ch][doc_of_LF_i] += 1;

                    // Prepare variables to use during the lcp queue traversal
                    size_t min_lcp = lcp_i; // this is lcp with previous suffix
                    bool passed_same_document = false;
                    std::fill(docs_to_collect.begin(), docs_to_collect.end(), false);
                    docs_to_collect[doc_of_LF_i] = true;

                    // Re-initialize doc profiles and max lcp with current document (itself)
                    std::fill(curr_da_profile.begin(), curr_da_profile.end(), 0);
                    curr_da_profile[doc_of_LF_i] = ref_build->total_length - pos_of_LF_i;

                    // lambda to check if we haven't found a certain document
                    auto all = [&](std::vector<bool> docs_collected) {
                            bool found_all = true;
                            for (auto elem: docs_collected)
                                found_all &= elem;
                            return found_all;
                    };

                    // Update the vector containing all the lcp values in queue
                    lcp_vals_in_queue.push_back(lcp_i);

                    // get index before current suffix, use int for signedness
                    int queue_pos = lcp_queue.size() - 2; 
                    while (queue_pos >= 0 && (!all(docs_to_collect) || !passed_same_document)) 
                    {  
                        // Case 1: Cross same character and document pair
                        if (lcp_queue[queue_pos].bwt_ch == curr_bwt_ch && 
                            lcp_queue[queue_pos].doc_num == doc_of_LF_i) 
                        {
                                passed_same_document = true;
                        } 
                        // Case 2: Cross same character, but a different document
                        else if (lcp_queue[queue_pos].bwt_ch == curr_bwt_ch &&
                                 lcp_queue[queue_pos].doc_num != doc_of_LF_i) 
                        {
                            // Case 2a: Update the current DA profile we haven't seen this document yet
                            if (!docs_to_collect[lcp_queue[queue_pos].doc_num]) 
                            {
                                curr_da_profile[lcp_queue[queue_pos].doc_num] = min_lcp + 1;
                                docs_to_collect[lcp_queue[queue_pos].doc_num] = true;
                            }

                            // Case 2b: Update the predecessor DA profile if its run boundary
                            // and we haven't passed the <ch, doc> pair as current suffix
                            if (!passed_same_document && (lcp_queue[queue_pos].is_start or lcp_queue[queue_pos].is_end))
                            {
                                size_t start_pos = ref_build->num_docs * queue_pos;
                                size_t curr_max_lcp = lcp_queue_profiles[start_pos + doc_of_LF_i];
                                lcp_queue_profiles[start_pos + doc_of_LF_i] = std::max(curr_max_lcp, min_lcp+1);
                            }
                        }
                        
                        min_lcp = std::min(min_lcp, lcp_queue[queue_pos].lcp_with_prev_suffix);
                        queue_pos--;
                    }

                    // Add the current profile to the vector (should always be multiple of # of docs)
                    for (auto elem: curr_da_profile)
                        lcp_queue_profiles.push_back(elem);

                    assert(lcp_queue_profiles.size() % ref_build->num_docs == 0);
                    assert(lcp_queue_profiles.size() == (lcp_queue.size() * ref_build->num_docs));

                    // Try to trim the LCP queue and adjust the count matrix ...
                    // Method #1: non-heuristic by removing entries with multiples 
                    size_t curr_pos = 0;
                    size_t records_to_remove = 0;
                    
                    while (curr_pos < lcp_queue.size()) {
                        uint8_t curr_ch = lcp_queue[curr_pos].bwt_ch;
                        size_t curr_doc = lcp_queue[curr_pos].doc_num;
                        assert(ch_doc_counters[curr_ch][curr_doc] >= 1);

                        size_t count = 0;
                        for (size_t i = 0; i < num_docs; i++) {
                            if (i != curr_doc && ch_doc_counters[curr_ch][i] >= 1)
                                count++;
                        }

                        // Means we cannot remove it from LCP queue
                        if ((ch_doc_counters[curr_ch][curr_doc] > 1 && (count == num_docs - 1))  || (curr_ch == EndOfDict))
                            records_to_remove++;
                        else    
                            break;
                        curr_pos++;
                    }

                    // Method #2: heuristic-based, remove records before a 
                    // a small LCP value (7) since they are insignificant
                    /*
                    size_t records_to_remove_for_lcp = 0;
                    curr_pos = 0;
                    while (records_to_remove_for_lcp < lcp_queue.size()) {
                        if (lcp_queue[records_to_remove_for_lcp++].lcp_with_prev_suffix <= 7)
                            break;
                    }
                    records_to_remove_for_lcp = (records_to_remove_for_lcp == lcp_queue.size()) ? (0) : records_to_remove_for_lcp;
                    */

                    // Method #3: non-heuristic, remove all the records starting
                    // from the the front of the queue where the smallest doc array
                    // profile value is larger than the current smallest lcp value in 
                    // queue, meaning that no successors rows will have the max lcp with
                    // these records.
                    /*
                    size_t records_to_remove_method2 = 0;
                    curr_pos = 0;
                    while (curr_pos < lcp_queue.size()) {
                        // Get the smallest LCP in the current profile we are looking at
                        size_t start_pos = ref_build->num_docs * curr_pos;
                        size_t curr_record_min_lcp = ref_build->total_length;
                        for (size_t i = 0; i < num_docs; i++) 
                            curr_record_min_lcp = std::min(lcp_queue_profiles[start_pos+i], curr_record_min_lcp);
                        
                        // Smallest LCP from current record until end of queue
                        size_t curr_min = *std::min_element(lcp_vals_in_queue.begin()+curr_pos+1, lcp_vals_in_queue.end());

                        // If the smallest value less than or equal to curr mininimum,
                        // then we cannot write it out.
                        if (curr_record_min_lcp <= curr_min && lcp_queue[curr_pos].bwt_ch != EndOfDict)
                            break;
                        else 
                            records_to_remove_method2++;
                        curr_pos++; 
                    } */ 


                    // Take the maximum value from the two methods above, BUT
                    // we cannot reduce the queue to empty because we need to 
                    // make sure to update records that are the end of runs
                    //records_to_remove = std::max(records_to_remove, records_to_remove_method2);
                    records_to_remove = std::min(records_to_remove, lcp_queue.size()-1);

                    // Remove those top n elements, and update counts
                    update_lcp_queue(records_to_remove, true);
                    
                    /* End of the DA Profiles code (except for some update statements below) */

                    // Update prevs
                    prev_occ = *curr_occ.first;
                    prev_bwt_ch = curr_bwt_ch;

                    // Update pq
                    curr_occ.first++;
                    if (curr_occ.first != curr_occ.second.first)
                        pq.push(curr_occ);

                    j += 1;
                    pos += 1;
                }

                prev = same_suffix.back();
                curr = next;
            }
            else {
                inc(curr);
            }
        }
        // print last BWT char and SA sample
        print_sa();
        print_bwt();
        print_doc_profiles();

        // Close output files
        fclose(ssa_file);
        fclose(esa_file);
        fclose(bwt_file);
        fclose(lcp_file);

        if (rle)
            fclose(bwt_file_len);
    }



private:
    typedef struct
    {
        size_t i = 0; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t phrase = 0;
        size_t suffix_length = 0;
        int_da sn = 0;
        uint8_t bwt_char = 0;
    } phrase_suffix_t;

    // Contains the positions in the lcp_queue for each
    // <ch, doc> pair. This allows us to only traverse 
    // entries that are relevant to the current row.
    std::vector<std::vector<std::deque<size_t>>> queue_pos_per_tuple;

    // Data-structure to maintain all the lcp values currently in queue
    std::deque<size_t> lcp_vals_in_queue;

    // Data-structure to store the max LCP with respect to all previous <ch, doc> pairs
    std::vector<std::vector<size_t>> predecessor_max_lcp;

    // Each entry in this will represent a row in BWM
    std::deque<queue_entry_t> lcp_queue;

    // This structure contains the raw document array profile values. The 
    // ith continguous group of num_doc integers is a profile for lcp_queue[i]
    std::deque<size_t> lcp_queue_profiles;

    // matrix of counters (alphabet size * number of documents) for the 
    // <ch, doc> pairs in lcp_queue, used for trimming method #1
    std::vector<std::vector<size_t>> ch_doc_counters;

    // matrix of <ch, doc> pairs that keep track of which pairs 
    // we have seen so far. This is used when we initialize
    // the curr_da_profile with the max lcps with previous lcps
    std::vector<std::vector<bool>> ch_doc_encountered;

    size_t j = 0;
    size_t ssa = 0;
    size_t esa = 0;
    size_t num_docs = 0;
    size_t num_records_ejected = 0;
    bool use_taxcomp = false;
    bool use_topk = false;
    size_t tax_sdap_overflow_ptr = 0, tax_edap_overflow_ptr = 0;

    FILE *sdap_file; // start of document array profiles
    FILE *edap_file; // end of document array profiles
    FILE *sdap_tax, *edap_tax;
    FILE *sdap_overtax, *edap_overtax;
    FILE *sdap_topk, *edap_topk;
    FILE *lcp_file; // LCP array
    FILE *bwt_file; // BWT (run characters if using rle)
    FILE *bwt_file_len; // lengths file is using rle
    FILE *ssa_file; // start of suffix array sample
    FILE *esa_file; // end of suffix array sample

    std::ofstream leaveout_fd;

    void print_doc_profiles() {
        /* Go through the leftover lcp queue, and print all the profiles for run boundaries */
        
        std::vector<size_t> curr_profile;
        std::vector<size_t> left_increases, right_increases;
        curr_profile.reserve(num_docs);
        left_increases.reserve(1000);
        right_increases.reserve(1000);

        for (size_t i = 0; i < lcp_queue.size(); i++) {
            uint8_t curr_ch = lcp_queue[i].bwt_ch;
            bool is_start = lcp_queue[i].is_start;
            bool is_end = lcp_queue[i].is_end;

            num_records_ejected++;

            /*
             * Write the data-structure to disk based on the 
             * command-line option chosen.
             */
            if (!use_taxcomp && !use_topk) {
                // remove the DA profile, and print if it's a boundary.
                // The format for each entry is print out the character of this run
                // followed by all the DA entries in DOCWIDTH bytes.
                if (is_start && fwrite(&curr_ch, 1, 1, sdap_file) != 1)
                    FATAL_ERROR("issue occurred while writing to *.sdap file");
                if (is_end && fwrite(&curr_ch, 1, 1, edap_file) != 1)
                    FATAL_ERROR("issue occurred while writing to *.edap file");

                for (size_t j = 0; j < num_docs; j++) {
                    size_t prof_val = std::min(lcp_queue_profiles.front(), (size_t) MAXLCPVALUE);
                    lcp_queue_profiles.pop_front();

                    if (is_start && fwrite(&prof_val, DOCWIDTH, 1, sdap_file) != 1)
                        error("SA write error 1");
                    if (is_end && fwrite(&prof_val, DOCWIDTH, 1, edap_file) != 1)
                        error("SA write error 1");
                }
            } else if (use_taxcomp) {
                size_t curr_val = 0, prev_max = std::min(lcp_queue_profiles.front(), (size_t) MAXLCPVALUE);
                left_increases.push_back(0);
                left_increases.push_back(prev_max);

                // Read the full profile, and keep track of increases from left
                for (size_t j = 0; j < num_docs; j++) {
                    curr_val = std::min(lcp_queue_profiles.front(), (size_t) MAXLCPVALUE);
                    lcp_queue_profiles.pop_front();
                    curr_profile[j] = curr_val;

                    if (curr_val > prev_max) {
                        left_increases.push_back(j);
                        left_increases.push_back(curr_val);
                        prev_max = curr_val;
                    }
                }
                ASSERT((left_increases.size() % 2 == 0), "issue occurred during the compression of profiles.");
                size_t num_left_inc = left_increases.size();
                size_t num_right_inc = 0;
                
                if (is_start || is_end) {
                    // Read the full profile, and keep track of increases from right
                    curr_val = 0; prev_max = curr_profile[num_docs-1];
                    right_increases.push_back(num_docs-1);
                    right_increases.push_back(prev_max);

                    for (int j = num_docs-1; j >= 0; j--) {
                        if (curr_profile[j] > prev_max) {
                            right_increases.push_back(j);
                            right_increases.push_back(curr_profile[j]);
                            prev_max = curr_profile[j];
                        }
                    }
                    ASSERT((right_increases.size() % 2 == 0), "issue occurred during the compression of profiles.");
                    num_right_inc = right_increases.size();

                    if (is_start) {
                        if (fwrite(&curr_ch, 1, 1, sdap_tax) != 1)
                            FATAL_ERROR("issue occurred while writing char to *taxcomp.sdap file");
                        for (size_t i = 0; i < NUMCOLSFORTABLE; i++) 
                            write_to_taxcomp_dap(sdap_tax, left_increases, right_increases, i);
                        append_overflow_pointer(sdap_tax, tax_sdap_overflow_ptr, num_left_inc, num_right_inc);
                        if (num_left_inc > NUMCOLSFORTABLE || num_right_inc > NUMCOLSFORTABLE) {
                            tax_sdap_overflow_ptr = write_remaining_pairs_to_overflow(sdap_overtax, left_increases, tax_sdap_overflow_ptr);
                            tax_sdap_overflow_ptr = write_remaining_pairs_to_overflow(sdap_overtax, right_increases, tax_sdap_overflow_ptr);
                        }
                    }
                    if (is_end) {
                        if (fwrite(&curr_ch, 1, 1, edap_tax) != 1)
                            FATAL_ERROR("issue occurred while writing char to *taxcomp.edap file");
                        for (size_t i = 0; i < NUMCOLSFORTABLE; i++) 
                            write_to_taxcomp_dap(edap_tax, left_increases, right_increases, i);
                        append_overflow_pointer(edap_tax, tax_edap_overflow_ptr, num_left_inc, num_right_inc);
                        if (num_left_inc > NUMCOLSFORTABLE || num_right_inc > NUMCOLSFORTABLE) {
                            tax_edap_overflow_ptr = write_remaining_pairs_to_overflow(edap_overtax, left_increases, tax_edap_overflow_ptr);
                            tax_edap_overflow_ptr = write_remaining_pairs_to_overflow(edap_overtax, right_increases, tax_edap_overflow_ptr);
                        }
                    }
                }                
                left_increases.clear();
                right_increases.clear();
            } else {
                // Read the full profile
                size_t curr_val = 0;
                for (size_t j = 0; j < num_docs; j++) {
                    curr_val = std::min(lcp_queue_profiles.front(), (size_t) MAXLCPVALUE);
                    lcp_queue_profiles.pop_front();
                    curr_profile[j] = curr_val;
                }

                if (is_start || is_end) {
                    // Sort documents by their LCP
                    std::vector<size_t> idx(num_docs);
                    std::iota(idx.begin(), idx.end(), 0);
                    std::stable_sort(idx.begin(), idx.end(),
                                     [&curr_profile](size_t v1, size_t v2) {return curr_profile[v1] > curr_profile[v2];});

                    // Write the topk documents to the index files
                    for (size_t j = 0; j < NUMCOLSFORTABLE; j++) {
                        size_t curr_doc = idx[j];
                        if (is_start && fwrite(&curr_doc, DOCWIDTH, 1, sdap_topk) != 1)
                            FATAL_ERROR("error occurred while writing to *topk.sdap file");
                        if (is_start && fwrite(&curr_profile[curr_doc], DOCWIDTH, 1, sdap_topk) != 1)
                            FATAL_ERROR("error occurred while writing to *topk.sdap file");
                        if (is_end && fwrite(&curr_doc, DOCWIDTH, 1, edap_topk) != 1)
                            FATAL_ERROR("error occurred while writing to *topk.edap file");
                        if (is_end && fwrite(&curr_profile[curr_doc], DOCWIDTH, 1, edap_topk) != 1)
                            FATAL_ERROR("error occurred while writing to *topk.edap file");
                    }
                }
            }
        }
    }

    void update_lcp_queue(size_t num_to_remove, bool update_table) {
        /* remove the first n records from the lcp queue and update count matrix */

        std::vector<size_t> curr_profile;
        std::vector<size_t> left_increases, right_increases;
        curr_profile.reserve(num_docs);
        left_increases.reserve(1000);
        right_increases.reserve(1000);

        for (size_t i = 0; i < num_to_remove; i++) {
            uint8_t curr_ch = lcp_queue[0].bwt_ch;
            size_t curr_doc = lcp_queue[0].doc_num;  
            bool is_start = lcp_queue[0].is_start;
            bool is_end = lcp_queue[0].is_end;

            // Update <ch, doc> count matrix
            if (update_table)
                ch_doc_counters[curr_ch][curr_doc] -= 1;

            // Update the queue position lists
            num_records_ejected++;
            queue_pos_per_tuple[curr_ch][curr_doc].pop_front();
            
            // Update the lcp value in the queue vector, and the queue_entry
            lcp_vals_in_queue.pop_front();
            lcp_queue.pop_front();

            /*
             * Write the data-structure to disk based on the 
             * command-line option chosen.
             */
            if (!use_taxcomp && !use_topk) {
                // remove the DA profile, and print if it's a boundary.
                // The format for each entry is print out the character of this run
                // followed by all the DA entries in DOCWIDTH bytes.
                if (is_start && fwrite(&curr_ch, 1, 1, sdap_file) != 1)
                    FATAL_ERROR("issue occurred while writing to *.sdap file");
                if (is_end && fwrite(&curr_ch, 1, 1, edap_file) != 1)
                    FATAL_ERROR("issue occurred while writing to *.edap file");
                    
                for (size_t j = 0; j < num_docs; j++) {
                    size_t prof_val = std::min(lcp_queue_profiles.front(), (size_t) MAXLCPVALUE);
                    lcp_queue_profiles.pop_front();

                    if (is_start && fwrite(&prof_val, DOCWIDTH, 1, sdap_file) != 1)
                        error("SA write error 1");
                    if (is_end && fwrite(&prof_val, DOCWIDTH, 1, edap_file) != 1)
                        error("SA write error 1");
                }
            } else if (use_taxcomp) {
                size_t curr_val = 0, prev_max = std::min(lcp_queue_profiles.front(), (size_t) MAXLCPVALUE);
                left_increases.push_back(0);
                left_increases.push_back(prev_max);

                // Read the full profile, and keep track of increases from left
                for (size_t j = 0; j < num_docs; j++) {
                    curr_val = std::min(lcp_queue_profiles.front(), (size_t) MAXLCPVALUE);
                    lcp_queue_profiles.pop_front();
                    curr_profile[j] = curr_val;

                    if (curr_val > prev_max) {
                        left_increases.push_back(j);
                        left_increases.push_back(curr_val);
                        prev_max = curr_val;
                    }
                }
                ASSERT((left_increases.size() % 2 == 0), "issue occurred during the compression of profiles.");
                size_t num_left_inc = left_increases.size();
                size_t num_right_inc = 0;
                
                if (is_start || is_end) {
                    // Read the full profile, and keep track of increases from right
                    curr_val = 0; prev_max = curr_profile[num_docs-1];
                    right_increases.push_back(num_docs-1);
                    right_increases.push_back(prev_max);

                    for (int j = num_docs-1; j >= 0; j--) {
                        if (curr_profile[j] > prev_max) {
                            right_increases.push_back(j);
                            right_increases.push_back(curr_profile[j]);
                            prev_max = curr_profile[j];
                        }
                    }
                    ASSERT((right_increases.size() % 2 == 0), "issue occurred during the compression of profiles.");
                    num_right_inc = right_increases.size();

                    if (is_start) {
                        if (fwrite(&curr_ch, 1, 1, sdap_tax) != 1)
                            FATAL_ERROR("issue occurred while writing char to *taxcomp.sdap file");
                        for (size_t i = 0; i < NUMCOLSFORTABLE; i++) 
                            write_to_taxcomp_dap(sdap_tax, left_increases, right_increases, i);
                        append_overflow_pointer(sdap_tax, tax_sdap_overflow_ptr, num_left_inc, num_right_inc);
                        if (num_left_inc > NUMCOLSFORTABLE || num_right_inc > NUMCOLSFORTABLE) {
                            tax_sdap_overflow_ptr = write_remaining_pairs_to_overflow(sdap_overtax, left_increases, tax_sdap_overflow_ptr);
                            tax_sdap_overflow_ptr = write_remaining_pairs_to_overflow(sdap_overtax, right_increases, tax_sdap_overflow_ptr);
                        }
                    }
                    if (is_end) {
                        if (fwrite(&curr_ch, 1, 1, edap_tax) != 1)
                            FATAL_ERROR("issue occurred while writing char to *taxcomp.edap file");
                        for (size_t i = 0; i < NUMCOLSFORTABLE; i++) 
                            write_to_taxcomp_dap(edap_tax, left_increases, right_increases, i);
                        append_overflow_pointer(edap_tax, tax_edap_overflow_ptr, num_left_inc, num_right_inc);
                        if (num_left_inc > NUMCOLSFORTABLE || num_right_inc > NUMCOLSFORTABLE) {
                            tax_edap_overflow_ptr = write_remaining_pairs_to_overflow(edap_overtax, left_increases, tax_edap_overflow_ptr);
                            tax_edap_overflow_ptr = write_remaining_pairs_to_overflow(edap_overtax, right_increases, tax_edap_overflow_ptr);
                        }
                    }
                }            
                left_increases.clear();
                right_increases.clear();
            } else {
                // Read the full profile
                size_t curr_val = 0;
                for (size_t j = 0; j < num_docs; j++) {
                    curr_val = std::min(lcp_queue_profiles.front(), (size_t) MAXLCPVALUE);
                    lcp_queue_profiles.pop_front();
                    curr_profile[j] = curr_val;
                }

                if (is_start || is_end) {
                    // Sort documents by their LCP
                    std::vector<size_t> idx(num_docs);
                    std::iota(idx.begin(), idx.end(), 0);
                    std::stable_sort(idx.begin(), idx.end(),
                                     [&curr_profile](size_t v1, size_t v2) {return curr_profile[v1] > curr_profile[v2];});

                    // Write the topk documents to the index files
                    for (size_t j = 0; j < NUMCOLSFORTABLE; j++) {
                        size_t curr_doc = idx[j];
                        if (is_start && fwrite(&curr_doc, DOCWIDTH, 1, sdap_topk) != 1)
                            FATAL_ERROR("error occurred while writing to *topk.sdap file");
                        if (is_start && fwrite(&curr_profile[curr_doc], DOCWIDTH, 1, sdap_topk) != 1)
                            FATAL_ERROR("error occurred while writing to *topk.sdap file");
                        if (is_end && fwrite(&curr_doc, DOCWIDTH, 1, edap_topk) != 1)
                            FATAL_ERROR("error occurred while writing to *topk.edap file");
                        if (is_end && fwrite(&curr_profile[curr_doc], DOCWIDTH, 1, edap_topk) != 1)
                            FATAL_ERROR("error occurred while writing to *topk.edap file");
                    }
                }
            }
            
        }
    }

    size_t write_remaining_pairs_to_overflow(FILE* outfile, std::vector<size_t> inc_pairs, size_t curr_ptr_pos) {
        // Check if we need to write anything in the first place
        if (inc_pairs.size()/2 > NUMCOLSFORTABLE) {
            size_t num_to_write = inc_pairs.size()/2 - NUMCOLSFORTABLE;
            ASSERT((num_to_write < 256), "we need more than 255 values to compress the profiles.");

            if (fwrite(&num_to_write, 1, 1, outfile) != 1)
                FATAL_ERROR("issue occurred while writing to overflow table.");
            curr_ptr_pos += 1;

            for (size_t j = NUMCOLSFORTABLE; j < inc_pairs.size()/2; j++) {
                bool success = (fwrite(&inc_pairs[j*2], DOCWIDTH, 1, outfile) == 1);
                success &= fwrite(&inc_pairs[j*2+1], DOCWIDTH, 1, outfile) == 1;
                if (!success)
                    FATAL_ERROR("issue occurred while writing the left data to overflow table.");
                curr_ptr_pos += (DOCWIDTH * 2);
            }
        } else {
            // Do not write anything if both left and right are already
            // included in the table...
            // size_t val = 0;
            // if (fwrite(&val, 1, 1, outfile) != 1)
            //     FATAL_ERROR("issue occurred while writing to overflow table.");
            // curr_ptr_pos += 1;
        }
        return curr_ptr_pos;
    }

    void write_to_taxcomp_dap(FILE* outfile, std::vector<size_t> left_inc, std::vector<size_t> right_inc, size_t col_num) {
        size_t index = col_num * 2;
        size_t left_pos = (col_num < left_inc.size()/2) ? left_inc[index]: MAXLCPVALUE;
        size_t left_lcp = (col_num < left_inc.size()/2) ? left_inc[index+1]: 0;
        size_t right_pos = (col_num < right_inc.size()/2) ? right_inc[index]: MAXLCPVALUE;
        size_t right_lcp = (col_num < right_inc.size()/2) ? right_inc[index+1]: 0;

        bool success = fwrite(&left_pos, DOCWIDTH, 1, outfile) == 1;
        success &= fwrite(&left_lcp, DOCWIDTH, 1, outfile) == 1;
        success &= fwrite(&right_pos, DOCWIDTH, 1, outfile) == 1;
        success &= fwrite(&right_lcp, DOCWIDTH, 1, outfile) == 1;
        if (!success)
            FATAL_ERROR("issue occurred during the writing to the *.taxcomp.dap file");
    }

    void append_overflow_pointer(FILE* outfile, size_t curr_pos, size_t num_left_inc, size_t num_right_inc) {
        if (num_left_inc > NUMCOLSFORTABLE || num_right_inc > NUMCOLSFORTABLE) {
            if (fwrite(&curr_pos, sizeof(size_t), 1, outfile) != 1)
                FATAL_ERROR("issue occurred when writing the overflow pointer.");
        } else {
            size_t overflow_pos = 0;
            if (fwrite(&overflow_pos, sizeof(size_t), 1, outfile) != 1)
                FATAL_ERROR("issue occurred when writing the overflow pointer.");
        }
    }

    inline bool inc(phrase_suffix_t& s)
    {
        s.i++;
        if (s.i >= pf.dict.saD.size())
            return false;
        s.sn = pf.dict.saD[s.i];
        s.phrase = pf.dict.rank_b_d(s.sn);
        // s.phrase = pf.dict.daD[s.i] + 1; // + 1 because daD is 0-based
        assert(!is_valid(s) || (s.phrase > 0 && s.phrase < pf.pars.ilist.size()));
        s.suffix_length = pf.dict.select_b_d(pf.dict.rank_b_d(s.sn + 1) + 1) - s.sn - 1;
        if(is_valid(s))
            s.bwt_char = (s.sn == pf.w ? 0 : pf.dict.d[s.sn - 1]);
        return true;
    }

    inline bool is_valid(phrase_suffix_t& s)
    {
        // avoid the extra w # at the beginning of the text
        if (s.sn < pf.w)
            return false;
        // Check if the suffix has length at least w and is not the complete phrase.
        if (pf.dict.b_d[s.sn] != 0 || s.suffix_length < pf.w)
            return false;
        
        return true;
    }
    
    inline int_t min_s_lcp_T(size_t left, size_t right)
    {
        // assume left < right
        if (left > right)
            std::swap(left, right);

        assert(pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] >= pf.w);

        return (pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] - pf.w);
    }

    inline int_t compute_lcp_suffix(phrase_suffix_t& curr, phrase_suffix_t& prev)
    {
        int_t lcp_suffix = 0;

        if (j > 0)
        {
            // Compute phrase boundary lcp
            lcp_suffix = pf.dict.lcpD[curr.i];
            for (size_t k = prev.i + 1; k < curr.i; ++k)
            {
                lcp_suffix = std::min(lcp_suffix, pf.dict.lcpD[k]);
            }

            if (lcp_suffix >= curr.suffix_length && curr.suffix_length == prev.suffix_length)
            {
                // Compute the minimum s_lcpP of the phrases following the two phrases
                // we take the first occurrence of the phrase in BWT_P
                size_t left = pf.pars.ilist[pf.pars.select_ilist_s(curr.phrase + 1)]; //size_t left = first_P_BWT_P[phrase];
                // and the last occurrence of the previous phrase in BWT_P
                size_t right = pf.pars.ilist[pf.pars.select_ilist_s(prev.phrase + 2) - 1]; //last_P_BWT_P[prev_phrase];
                
                lcp_suffix += min_s_lcp_T(left,right);
            }
        }

        return lcp_suffix;
    }

    inline void print_lcp(int_t val, size_t pos)
    {
        size_t tmp_val = val;
        if (fwrite(&tmp_val, THRBYTES, 1, lcp_file) != 1)
            error("LCP write error 1");
    }

    // We can put here the check if we want to store the LCP or stream it out
    inline void new_min_s(int_t val, size_t pos)
    {
        min_s.push_back(val);
        pos_s.push_back(j);
    }

    inline void update_ssa(phrase_suffix_t &curr, size_t pos)
    {
        ssa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
        assert(ssa < (pf.n - pf.w + 1ULL));
    }

    inline void update_esa(phrase_suffix_t &curr, size_t pos)
    {
        esa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
        assert(esa < (pf.n - pf.w + 1ULL));
    }

    inline void print_sa()
    {
        if (j < (pf.n - pf.w + 1ULL))
        {
            size_t pos = j;
            if (fwrite(&pos, SSABYTES, 1, ssa_file) != 1)
                error("SA write error 1");
            if (fwrite(&ssa, SSABYTES, 1, ssa_file) != 1)
                error("SA write error 2");
        }

        if (j > 0)
        {
            size_t pos = j - 1;
            if (fwrite(&pos, SSABYTES, 1, esa_file) != 1)
                error("SA write error 1");
            if (fwrite(&esa, SSABYTES, 1, esa_file) != 1)
                error("SA write error 2");
        }
    }

    inline void print_bwt()
    {   
        if (length > 0)
        {
            if (rle) {
                // write the head character
                if (fputc(head, bwt_file) == EOF)
                    error("BWT write error 1");

                // write the length of that run
                if (fwrite(&length, BWTBYTES, 1, bwt_file_len) != 1)
                    error("BWT write error 2");
            } else {
                for (size_t i = 0; i < length; ++i)
                {
                    if (fputc(head, bwt_file) == EOF)
                        error("BWT write error 1");
                }
            }
        }
    }

    inline void update_bwt(uint8_t next_char, size_t length_)
    {
        if (head != next_char)
        {
            print_sa();
            print_bwt();

            head = next_char;

            length = 0;
        }
        length += length_;

    }
};

#endif /* end of include guard: _LCP_PFP_HH */
