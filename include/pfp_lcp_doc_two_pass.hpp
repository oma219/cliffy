/*
 * File: pfp_lcp_doc_two_pass.hpp
 * Description: Header file containing the two-pass 
 *              construction algorithm for the document
 *              array profiles.
 * Date: September 17th, 2022
 */

#ifndef _LCP_PFP_DOC_TWO_PASS_HH
#define _LCP_PFP_DOC_TWO_PASS_HH

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
#include <cstdio>

#define TEMPDATA_RECORD 6

/* MACROS for reading from file */
#define GET_IS_START(fd, num) (0x2 & fd[num * TEMPDATA_RECORD]) >> 1
#define GET_IS_END(fd, num) (0x1 & fd[num * TEMPDATA_RECORD])
#define GET_BWT_CH(fd, num) (fd[num * TEMPDATA_RECORD + 1])
#define GET_DOC_OF_LF(fd, num) ((0xFF & fd[num * TEMPDATA_RECORD + 2]) | ((0xFF & fd[num * TEMPDATA_RECORD + 3]) << 8))
#define GET_LCP(fd, num) ((0xFF & fd[num * TEMPDATA_RECORD + 4]) | ((0xFF & fd[num * TEMPDATA_RECORD + 5]) << 8))
#define GET_DAP(fd, num) ((0xFF & fd[num * DOCWIDTH]) | ((0xFF & fd[num * DOCWIDTH + 1]) << 8))

/* struct for temp lcp queue data */
typedef struct
{
    bool is_start = false;
    bool is_end = false;
    uint8_t bwt_ch = 0;
    size_t doc_num = 0;
    size_t lcp_i = 0;
} temp_data_entry_t;

/* struct for each phrase suffix */
typedef struct
{
    // this should be safe since the first entry of sa is 
    // always the dollarsign used to compute the sa
    size_t i = 0; 
    size_t phrase = 0;
    size_t suffix_length = 0;
    int_da sn = 0;
    uint8_t bwt_char = 0;
} phrase_suffix_t;

class pfp_lcp_doc_two_pass {
    public:
        size_t total_num_runs = 0;

    pfp_lcp_doc_two_pass(pf_parsing &pfp_, std::string filename, RefBuilder* ref_build, 
                            std::string temp_prefix, size_t tmp_size, bool taxcomp, bool topk, 
                            size_t num_cols, bool rle_ = true) : 
                pf(pfp_),
                min_s(1, pf.n),
                pos_s(1,0),
                head(0),
                NUMCOLSFORTABLE(num_cols),
                num_docs(ref_build->num_docs),
                ch_doc_encountered(256, std::vector<bool>(ref_build->num_docs, false)),
                rle(rle_),
                tmp_file_size(tmp_size),
                use_taxcomp(taxcomp),
                use_topk(topk)
    {   
        /* construction algorithm for document array profiles */
        STATUS_LOG("build_main", "building bwt and doc profiles based on pfp (1st-pass)");
        auto start = std::chrono::system_clock::now();

        initialize_index_files(filename, temp_prefix);  
        initialize_tmp_size_variables();
        assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict);

        // variables for bwt/lcp/sa construction
        phrase_suffix_t curr;
        phrase_suffix_t prev;

        // variables for doc profile construction 
        uint8_t prev_bwt_ch = 0;
        size_t curr_run_num = 0;
        size_t pos = 0;
        std::vector<size_t> curr_da_profile (ref_build->num_docs, 0);

        // create a predecessor max lcp table, and re-initialize with max_lcp
        num_blocks_of_32 = num_docs/32;
        num_blocks_of_32++;
        uint16_t max_lcp_init = 0;
        
        // create a separate block of memory for each character
        predecessor_max_lcp = new uint16_t*[256];
        for (size_t i = 0; i < 256; i++) {
            predecessor_max_lcp[i] =  new uint16_t[num_blocks_of_32 * 32];
        }
        for (size_t i = 0; i < 256; i++)
            for (size_t j = 0; j < num_blocks_of_32 * 32; j++)
                predecessor_max_lcp[i][j] = max_lcp_init;

        // DEBUG: START ------------------------------------------
        // DON'T FORGET TO RELEASE THIS MEMORY
        predecessor_max_lcp2 = new uint16_t*[256];
        for (size_t i = 0; i < 256; i++) {
            predecessor_max_lcp2[i] =  new uint16_t[num_blocks_of_32 * 32];
        }
        for (size_t i = 0; i < 256; i++)
            for (size_t j = 0; j < num_blocks_of_32 * 32; j++)
                predecessor_max_lcp2[i][j] = max_lcp_init;
        
        last_updated_row = 256;
        dirty_lcp_cache = new uint16_t[256];
        for (size_t i = 0; i < 256; i++)
            dirty_lcp_cache[i] = max_lcp_init;
        // DEBUG: END ---------------------------------------------

        // create a struct to store data for temp file
        temp_data_entry_t curr_data_entry;

        // start of construction ... this loop generates the SA, LCP and BWT
        inc(curr);
        while (curr.i < pf.dict.saD.size())
        {
            // make sure current suffix is a valid proper phrase suffix 
            // (at least w characters but not whole phrase)
            if(is_valid(curr))
            {
                // compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);
                phrase_suffix_t next = curr;

                // go through suffix array of dictionary and store all phrase ids with same suffix
                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {
                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_suffix.push_back(next);
                    }
                }

                // hard case: phrases with different BWT characters precediing them
                int_t lcp_suffix = compute_lcp_suffix(curr, prev);
                typedef std::pair<int_t *, std::pair<int_t *, uint8_t>> pq_t;

                // using lambda to compare elements.
                auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                    return *lhs.first > *rhs.first;
                };
                
                // merge a list of occurrences of each phrase in the BWT of the parse
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
                        // compute the minimum s_lcpP of the the current and previous 
                        // occurrence of the phrase in BWT_P
                        lcp_suffix = curr.suffix_length + min_s_lcp_T(*curr_occ.first, prev_occ);
                    }
                    first = false;

                    // update min_s
                    print_lcp(lcp_suffix, j);
                    update_ssa(curr, *curr_occ.first);
                    update_bwt(curr_occ.second.second, 1);
                    update_esa(curr, *curr_occ.first);

                    ssa = (pf.pos_T[*curr_occ.first] - curr.suffix_length) % (pf.n - pf.w + 1ULL);
                    esa = (pf.pos_T[*curr_occ.first] - curr.suffix_length) % (pf.n - pf.w + 1ULL);

                    /* Start of DA Code */
                    curr_bwt_ch = curr_occ.second.second;
                    size_t lcp_i = lcp_suffix;
                    size_t sa_i = ssa;
                    size_t doc_i = ref_build->doc_ends_rank(ssa);

                    // determine whether current suffix is a run boundary   
                    bool is_start = (pos == 0 || curr_bwt_ch != prev_bwt_ch) ? 1 : 0;
                    bool is_end = (pos == ref_build->total_length-1); // only special case, common case is below

                    // handle scenario where the previous suffix was a end of a run
                    if (pos > 0 && prev_bwt_ch != curr_bwt_ch)
                        curr_data_entry.is_end = true;
                    
                    // push previous suffix data into temp data file
                    if (pos > 0) {
                        write_data_to_temp_file(curr_data_entry);
                        if (curr_data_entry.is_start || curr_data_entry.is_end)
                            write_profile_to_temp_file(curr_da_profile);
                    }
                    
                    if (is_start) {curr_run_num++;}
                    size_t pos_of_LF_i = (sa_i > 0) ? (sa_i - 1) : (ref_build->total_length-1);
                    size_t doc_of_LF_i = ref_build->doc_ends_rank(pos_of_LF_i);

                    DEBUG_MSG("pos = " << pos << 
                              ", run_num = " << curr_run_num <<  
                              ", ch = " << curr_bwt_ch << 
                              ", doc = " << doc_of_LF_i << 
                              ", lcp = " << lcp_i << 
                              ", suffix_length = " << (ref_build->total_length - pos_of_LF_i));

                    // gather current suffix data
                    curr_data_entry = {is_start, is_end, curr_bwt_ch, doc_of_LF_i, lcp_i};

                    // update matrices based on current BWT character, and doc
                    ch_doc_encountered[curr_bwt_ch][doc_of_LF_i] = true;

                    // re-initialize doc profiles and max lcp with current document (itself)
                    std::fill(curr_da_profile.begin(), curr_da_profile.end(), 0);
                    curr_da_profile[doc_of_LF_i] = ref_build->total_length - pos_of_LF_i;

                    // ** OLD VERSION **
                    // update_predecessor_max_lcp_table(lcp_i, ref_build->total_length, pos_of_LF_i, doc_of_LF_i, curr_bwt_ch);
                    // initialize_current_row_profile(doc_of_LF_i, curr_da_profile, curr_bwt_ch);
                    
                    // update the predecessor table, and initialize current document array profile
                    update_predecessor_max_lcp_table_lazy_version(lcp_i, ref_build->total_length, pos_of_LF_i, doc_of_LF_i, curr_bwt_ch);
                    initialize_current_row_profile_lazy_version(doc_of_LF_i, curr_da_profile, curr_bwt_ch);

                    /* End of DA Code */

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
        DONE_LOG((std::chrono::system_clock::now() - start));

        // make sure to write the last suffix to temp data
        write_data_to_temp_file(curr_data_entry);
        if (curr_data_entry.is_start || curr_data_entry.is_end)
            write_profile_to_temp_file(curr_da_profile);

        // re-initialize the predecessor table, and other structures
        for (size_t i = 0; i < 256; i++) {
            for (size_t j = 0; j < num_blocks_of_32 * 32; j++)
                predecessor_max_lcp[i][j] = max_lcp_init;
            for (size_t j = 0; j < num_docs; j++)
                ch_doc_encountered[i][j] = false;
        }

        // DEBUG: START ------------------------------------------
        for (size_t i = 0; i < 256; i++)
            for (size_t j = 0; j < num_blocks_of_32 * 32; j++)
                predecessor_max_lcp2[i][j] = max_lcp_init;
        
        last_updated_row = 256;
        for (size_t i = 0; i < 256; i++)
            dirty_lcp_cache[i] = max_lcp_init;
        // DEBUG: END ---------------------------------------------


        // perform 2nd pass which focuses only on document array profiles...
        STATUS_LOG("build_main", "performing 2nd pass to update profiles");
        start = std::chrono::system_clock::now();

        size_t dap_ptr = num_dap_temp_data - num_docs;
        for (size_t i = num_lcp_temp_data; i > 0; i--) {
            // re-initialize at each position
            std::fill(curr_da_profile.begin(), curr_da_profile.end(), 0);

            // grab the data for current suffix
            bool is_start = GET_IS_START(mmap_lcp_inter, (i-1));
            bool is_end = GET_IS_END(mmap_lcp_inter, (i-1));
            uint8_t bwt_ch = GET_BWT_CH(mmap_lcp_inter, (i-1));
            size_t doc_of_LF_i = GET_DOC_OF_LF(mmap_lcp_inter, (i-1));

            flush_row_of_lcp_table_up(bwt_ch);

            DEBUG_MSG("pos = " << pos << 
                      ", run_num = " << curr_run_num <<  
                      ", bwt_ch = " << bwt_ch << 
                      ", doc = " << doc_of_LF_i << 
                      ", lcp = " << lcp_i << 
                      ", suffix_length = " << (ref_build->total_length - pos_of_LF_i));

            // if suffix is a run boundary, update its profile if necessary
            if (is_start || is_end) {
                //initialize_current_row_profile(doc_of_LF_i, curr_da_profile, bwt_ch);
                initialize_current_row_profile_lazy_version(doc_of_LF_i, curr_da_profile, bwt_ch);
                for (size_t j = 0; j < num_docs; j++) {
                    if (GET_DAP(mmap_dap_inter, (dap_ptr+j)) < curr_da_profile[j]){
                        size_t new_lcp_i = std::min((size_t) MAXLCPVALUE, curr_da_profile[j]);
                        mmap_dap_inter[(dap_ptr+j) * DOCWIDTH] = (0xFF & new_lcp_i);
                        mmap_dap_inter[(dap_ptr+j) * DOCWIDTH + 1] = ((0xFF << 8) & new_lcp_i) >> 8;
                    }
                }
                dap_ptr -= num_docs;
            }
            ch_doc_encountered[bwt_ch][doc_of_LF_i] = true;

            // update the predecessor table for next suffix
            size_t lcp_i = GET_LCP(mmap_lcp_inter, (i-1));
            
            //update_predecessor_max_lcp_table_up(lcp_i, doc_of_LF_i, bwt_ch);
            update_predecessor_max_lcp_table_up_lazy_version(lcp_i, doc_of_LF_i, bwt_ch);

        }
        // print out build time
        DONE_LOG((std::chrono::system_clock::now() - start));

        // print the document array profiles
        print_dap(num_lcp_temp_data);

        // print last BWT char and SA sample
        print_sa();
        print_bwt();
        total_num_runs = curr_run_num;

        // close output files
        fclose(ssa_file); fclose(esa_file);
        fclose(bwt_file);
        fclose(bwt_file_len);
        fclose(lcp_file);

        // close the approriate dap files
        if (use_taxcomp) {
            fclose(sdap_tax); fclose(edap_tax);
            fclose(sdap_overtax); fclose(edap_overtax);    
        } else {   
            fclose(sdap_file); 
            fclose(edap_file);
        }

        // delete temporary files
        delete_temp_files(temp_prefix);

        // print out statistics
        size_t total_tmp_used = (num_lcp_temp_data * TEMPDATA_RECORD) + (num_dap_temp_data * DOCWIDTH);
        STATS_LOG("build_main", "stats: n = %ld, r = %ld, total_tmp_used = %ld (%ld + %ld)", 
                  j, total_num_runs, total_tmp_used, (num_lcp_temp_data * TEMPDATA_RECORD), (num_dap_temp_data * DOCWIDTH)); 
    }

    ~pfp_lcp_doc_two_pass() {
        // deallocate memory for predecessor table
        for (size_t i = 0; i < 256; i++)
            delete[] predecessor_max_lcp[i];
        delete[] predecessor_max_lcp;

        // unmap the temp data file
        munmap(mmap_lcp_inter, tmp_file_size);
        munmap(mmap_dap_inter, tmp_file_size);
        close(dap_inter_fd); close(lcp_inter_fd);
    }

    private:
        /*********************************/ 
        /* Private instance variables
        /*********************************/
        pf_parsing& pf;
        std::vector<size_t> min_s; // Value of the minimum lcp_T in each run of BWT_T
        std::vector<size_t> pos_s; // Position of the minimum lcp_T in each run of BWT_T

        uint8_t head; // character of current BWT run
        size_t length = 0; // length of the current BWT run

        bool rle; // run-length encode the BWT
        bool use_taxcomp = false;
        bool use_topk = false; 

        size_t NUMCOLSFORTABLE = 0; 

        size_t tmp_file_size = 0;
        size_t max_dap_records = 0;
        size_t max_lcp_records = 0;

        FILE *lcp_file; // LCP array
        FILE *bwt_file; // BWT (run characters if using rle)
        FILE *bwt_file_len; // lengths file is using rle
        FILE *ssa_file; // start of run: suffix array sample
        FILE *esa_file; // end of run: suffix array sample

        FILE *sdap_file; // start of run: document array profiles
        FILE *edap_file; // end of run: document array profiles
        FILE *sdap_tax, *edap_tax; // start or end when taxonomically compressing DAP
        FILE *sdap_overtax, *edap_overtax; // overflow file when taxonomically compressing DAP

        size_t num_docs = 0;
        size_t j = 0;
        size_t ssa = 0;
        size_t esa = 0;
        size_t num_blocks_of_32 = 0;
        uint8_t curr_bwt_ch = 0;

        size_t tax_sdap_overflow_ptr = 0;
        size_t tax_edap_overflow_ptr = 0;

        size_t num_dap_temp_data = 0;
        size_t num_lcp_temp_data = 0;
        
        /***********************************************************/
        /* matrix of <ch, doc> pairs that keep track of which pairs 
        /* we have seen so far. This is used when we initialize
        /* the curr_da_profile with the max lcps with previous lcps
        /***********************************************************/
        std::vector<std::vector<bool>> ch_doc_encountered;

        /* data-structure to store the max LCP with respect to all previous <ch, doc> pairs */
        uint16_t** predecessor_max_lcp;
        uint16_t** predecessor_max_lcp2;
        uint16_t* dirty_lcp_cache;
        uint16_t last_updated_row;

        /***********************************************************/
        /* This structure contains the raw document array profile 
        /* values. The ith continguous group of num_doc integers is 
        /* a profile for lcp_queue[i]
        /***********************************************************/
        std::deque<size_t> lcp_queue_profiles;

        /* variables are used for the intermediate file */
        int dap_inter_fd, lcp_inter_fd;
        struct stat file_info_dap_inter;
        struct stat file_info_lcp_inter;
        char* mmap_dap_inter;
        char* mmap_lcp_inter;

        /***********************************************************/
        /* Section 1: Document Array Profiles related methods
        /***********************************************************/
        void initialize_index_files(std::string filename, std::string temp_prefix) {
            /* opening output files for data-structures like  LCP, SA, BWT */

            // LCP data-structure
            std::string outfile = filename + std::string(".lcp");
            if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
            
            // Suffix array samples at start of runs
            outfile = filename + std::string(".ssa");
            if ((ssa_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
            
            // Suffix array samples at end of runs
            outfile = filename + std::string(".esa");
            if ((esa_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
            
            // BWT index files
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

            // Document Array Profiles file
            if (use_taxcomp) {
                outfile = filename + std::string(".taxcomp.sdap");
                if ((sdap_tax = fopen(outfile.c_str(), "w")) == nullptr) 
                    error("open() file " + outfile + " failed");
                outfile = filename + std::string(".taxcomp.edap");
                if ((edap_tax = fopen(outfile.c_str(), "w")) == nullptr) 
                    error("open() file " + outfile + " failed");

                outfile = filename + std::string(".taxcomp.of.sdap");
                if ((sdap_overtax = fopen(outfile.c_str(), "w")) == nullptr) 
                    error("open() file " + outfile + " failed");
                if (fwrite(&num_docs, sizeof(size_t), 1, sdap_overtax) != 1)
                    error("SDAP write error: number of documents");
                outfile = filename + std::string(".taxcomp.of.edap");
                if ((edap_overtax = fopen(outfile.c_str(), "w")) == nullptr) 
                    error("open() file " + outfile + " failed");
                if (fwrite(&num_docs, sizeof(size_t), 1, edap_overtax) != 1)
                    error("EDAP write error: number of documents");

                tax_sdap_overflow_ptr += sizeof(size_t);
                tax_edap_overflow_ptr += sizeof(size_t);
            } else {
                outfile = filename + std::string(".sdap");
                if ((sdap_file = fopen(outfile.c_str(), "w")) == nullptr) 
                    error("open() file " + outfile + " failed");
                if (fwrite(&num_docs, sizeof(size_t), 1, sdap_file) != 1)
                    error("SDAP write error: number of documents");

                outfile = filename + std::string(".edap");
                if ((edap_file = fopen(outfile.c_str(), "w")) == nullptr)
                    error("open() file " + outfile + " failed");
                if (fwrite(&num_docs, sizeof(size_t), 1, edap_file) != 1)
                    error("SDAP write error: number of documents"); 
            }

            // Open intermediate files for storing data
            mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;

            // 1. File for LCP Queue data
            outfile = temp_prefix + std::string(".tmp_lcp_data");
            lcp_inter_fd = open(outfile.c_str(), O_RDWR | O_CREAT | O_TRUNC, mode);
            ftruncate(lcp_inter_fd, tmp_file_size);

            if (lcp_inter_fd == -1)
                FATAL_ERROR("Issue occurred when opening file for intermediate file.");
            if (fstat(lcp_inter_fd, &file_info_lcp_inter) == -1)
                FATAL_ERROR("Issue with checking file size on overflow file.");

            mmap_lcp_inter = (char*) mmap(NULL, 
                                         tmp_file_size, 
                                         PROT_READ | PROT_WRITE, 
                                         MAP_SHARED, 
                                         lcp_inter_fd, 0);
            assert(mmap_lcp_inter != MAP_FAILED);
            
            // 2. File for Document Array Profiles
            outfile = temp_prefix + std::string(".tmp_dap_data");
            dap_inter_fd = open(outfile.c_str(), O_RDWR | O_CREAT | O_TRUNC, mode);
            ftruncate(dap_inter_fd, tmp_file_size);

            if (dap_inter_fd == -1)
                FATAL_ERROR("Issue occurred when opening file for intermediate file.");
            if (fstat(dap_inter_fd, &file_info_dap_inter) == -1)
                FATAL_ERROR("Issue with checking file size on overflow file.");

            mmap_dap_inter = (char*) mmap(NULL, 
                                         tmp_file_size, 
                                         PROT_READ | PROT_WRITE, 
                                         MAP_SHARED, 
                                         dap_inter_fd, 0);
            assert(mmap_dap_inter != MAP_FAILED);
        }

        void initialize_tmp_size_variables() {
            /* determine how many records we can store in each temporary file */
            max_lcp_records = tmp_file_size/TEMPDATA_RECORD; 
            max_dap_records = tmp_file_size/DOCWIDTH;
        }

        void update_predecessor_max_lcp_table(size_t lcp_i, size_t total_length, size_t pos_of_LF_i, size_t doc_of_LF_i, uint8_t bwt_ch) {
            /* 
             * Update the predecessor lcp table, this allows us to compute the maximum
             * lcp with respect to all the predecessor occurrences of other documents.
             * For example: if we are at <A, 2>, we will find the maximum lcp with the 
             * the previous occurrences of <A, 0> and <A, 1> for 3 documents.
             */
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

                //std::vector<size_t> dna_chars = {65, 67, 71, 78, 84, 85, 89}; // A, C, G, N, T, U, Y
                //for (size_t ch_num: dna_chars) { // Optimization for DNA

                // This loop only iterates over possible uint8_t chars, and 
                // 0, 1, and 2 are reserved for PFP so the range is [3, 255]
                for (size_t ch_num = 3; ch_num < 256; ch_num++) {
                    // use SIMD for all groups of 32
                    for (size_t i = 0; i < (num_blocks_of_32 * 32); i+=32) {
                        // zero-mask, all the set bit positions are loaded
                        arr1 = _mm512_maskz_loadu_epi16(~0, (const __m512i*) &predecessor_max_lcp[ch_num][i]); 

                        arr3 = _mm512_min_epu16(arr1, arr2);
                        _mm512_mask_storeu_epi16((__m512i*) &predecessor_max_lcp[ch_num][i], k, arr3); 
                    }
                }
                // reset the LCP with respect to the current <ch, doc> pair
                predecessor_max_lcp[bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, total_length - pos_of_LF_i);

                // DEBUG: START --------------------------------------------

                // // avoid overflow issues
                // uint16_t truncated_lcp_i = std::min((size_t) MAXLCPVALUE, lcp_i);

                // // update the dirty lcp table for each character
                // for (size_t i = 0; i < 256; i++) {
                //     dirty_lcp_cache[i] = std::min(dirty_lcp_cache[i], truncated_lcp_i);
                // }

                // // update the position corresponding the previous bwt ch since it was
                // // set to zero 
                // if (last_updated_row < 256)
                //     dirty_lcp_cache[last_updated_row] = truncated_lcp_i;
                
                // // update the bwt_ch row of table with minimum lcp since last update
                // uint16_t min_lcp_to_flush = dirty_lcp_cache[bwt_ch];
                // //std::cout << "min_lcp_to_flush = " << min_lcp_to_flush << std::endl;
                // for (size_t i = 0; i < num_docs; i++) {
                    
                //     predecessor_max_lcp2[bwt_ch][i] = std::min(predecessor_max_lcp2[bwt_ch][i], min_lcp_to_flush);
                //     //std::cout << "lcp2[" << bwt_ch << "][" << i << "] = " << predecessor_max_lcp2[bwt_ch][i] << std::endl;
                // }
                
                // // update table with length of current suffix
                // predecessor_max_lcp2[bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, total_length-pos_of_LF_i);
                
                // // update dirty lcp cache since we have already flush lcp_i to table
                // dirty_lcp_cache[bwt_ch] = 0;
                // last_updated_row = bwt_ch;
                // DEBUG: END ----------------------------------------------

            #else                        
                for (size_t ch_num = 0; ch_num < 256; ch_num++) {
                    for (size_t doc_num = 0; doc_num < num_docs; doc_num++) {
                        predecessor_max_lcp[ch_num][doc_num] = std::min(predecessor_max_lcp[ch_num][doc_num], (uint16_t) std::min(lcp_i, (size_t) MAXLCPVALUE));
                    }
                }
                // reset the LCP with respect to the current <ch, doc> pair
                predecessor_max_lcp[bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, total_length - pos_of_LF_i);


                // DEBUG: START --------------------------------------------

                // // avoid overflow issues
                // uint16_t truncated_lcp_i = std::min((size_t) MAXLCPVALUE, lcp_i);

                // // update the dirty lcp table for each character
                // for (size_t i = 0; i < 256; i++) {
                //     dirty_lcp_cache[i] = std::min(dirty_lcp_cache[i], truncated_lcp_i);
                // }

                // // update the position corresponding the previous bwt ch since it was
                // // set to zero 
                // if (last_updated_row < 256)
                //     dirty_lcp_cache[last_updated_row] = truncated_lcp_i;
                
                // // update the bwt_ch row of table with minimum lcp since last update
                // uint16_t min_lcp_to_flush = dirty_lcp_cache[bwt_ch];
                // for (size_t i = 0; i < (num_blocks_of_32 * 32); i+= 32) {
                //     predecessor_max_lcp2[bwt_ch][i] = std::min(predecessor_max_lcp2[bwt_ch][i], min_lcp_to_flush);
                // }
                
                // // update table with length of current suffix
                // predecessor_max_lcp2[bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, total_length-pos_of_LF_i);
                
                // // update dirty lcp cache since we have already flush lcp_i to table
                // dirty_lcp_cache[bwt_ch] = 0;
                // last_updated_row = bwt_ch;
                // DEBUG: END ----------------------------------------------



            #endif
        }

        void update_predecessor_max_lcp_table_lazy_version(size_t lcp_i, size_t total_length, size_t pos_of_LF_i, size_t doc_of_LF_i, uint8_t bwt_ch) {
            /* 
             * Update the predecessor lcp table, this allows us to compute the maximum
             * lcp with respect to all the predecessor occurrences of other documents.
             * For example: if we are at <A, 2>, we will find the maximum lcp with the 
             * the previous occurrences of <A, 0> and <A, 1> for 3 documents.
             */

            // avoid overflow issues
            uint16_t truncated_lcp_i = std::min((size_t) MAXLCPVALUE, lcp_i);

            // update the dirty lcp table for each character
            for (size_t i = 0; i < 256; i++) {
                dirty_lcp_cache[i] = std::min(dirty_lcp_cache[i], truncated_lcp_i);
            }

            // update the position corresponding the previous bwt_ch since
            // it was set to zero after flushing the previous lcp_i value
            if (last_updated_row < 256)
                dirty_lcp_cache[last_updated_row] = truncated_lcp_i;
            
            // use this lcp value update current bwt char
            uint16_t min_lcp_to_flush = dirty_lcp_cache[bwt_ch];

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
                    lcp_i_vector[i] = min_lcp_to_flush;

                // load the array of constants (lcp_i)
                __m512i arr1, arr2, arr3; 
                __mmask32 k = ~0; // all 32 bits on, means all 32 values will be written
                arr2 = _mm512_maskz_loadu_epi16(~0, (const __m512i*) &lcp_i_vector[0]);
                
                // // update the bwt_ch row of table with minimum lcp since last update
                // uint16_t min_lcp_to_flush = dirty_lcp_cache[bwt_ch];
                // for (size_t i = 0; i < num_docs; i++) {
                    
                //     predecessor_max_lcp2[bwt_ch][i] = std::min(predecessor_max_lcp2[bwt_ch][i], min_lcp_to_flush);
                //     //std::cout << "lcp2[" << bwt_ch << "][" << i << "] = " << predecessor_max_lcp2[bwt_ch][i] << std::endl;
                // }

                // use SIMD for all groups of 32
                for (size_t i = 0; i < (num_blocks_of_32 * 32); i+=32) {
                    // zero-mask, all the set bit positions are loaded
                    arr1 = _mm512_maskz_loadu_epi16(~0, (const __m512i*) &predecessor_max_lcp2[bwt_ch][i]); 
                    arr3 = _mm512_min_epu16(arr1, arr2);
                    _mm512_mask_storeu_epi16((__m512i*) &predecessor_max_lcp2[bwt_ch][i], k, arr3); 
                }
                
            #else                                   
                // update the bwt_ch row of table with minimum lcp since last update
                for (size_t i = 0; i < num_docs; i++) {
                    predecessor_max_lcp2[bwt_ch][i] = std::min(predecessor_max_lcp2[bwt_ch][i], min_lcp_to_flush);
                }
                
            #endif   

            // update table with length of current suffix
            predecessor_max_lcp2[bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, total_length-pos_of_LF_i);
            
            // update dirty lcp cache since we have already flush lcp_i to table
            dirty_lcp_cache[bwt_ch] = 0;
            last_updated_row = bwt_ch;
        }

        void update_predecessor_max_lcp_table_up(size_t lcp_i, size_t doc_of_LF_i, uint8_t bwt_ch) {
            /* 
             * Update the predecessor lcp table, this allows us to compute the maximum
             * lcp with respect to all the predecessor occurrences of other documents.
             * For example: if we are at <A, 2>, we will find the maximum lcp with the 
             * the previous occurrences of <A, 0> and <A, 1> for 3 documents.
             */
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

                // Reset the LCP with respect to the current <ch, doc> pair
                predecessor_max_lcp[bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, lcp_i);

                //std::vector<size_t> dna_chars = {65, 67, 71, 78, 84, 85, 89}; // A, C, G, N, T, U, Y
                //for (size_t ch_num: dna_chars) { // Optimization for DNA
                
                // This loop only iterates over possible uint8_t chars, and 
                // 0, 1, and 2 are reserved for PFP so the range is [3, 255]
                for (size_t ch_num = 3; ch_num < 256; ch_num++) {
                    // use SIMD for all groups of 32
                    for (size_t i = 0; i < (num_blocks_of_32 * 32); i+=32) {
                        // zero-mask, all the set bit positions are loaded
                        arr1 = _mm512_maskz_loadu_epi16(~0, (const __m512i*) &predecessor_max_lcp[ch_num][i]); 

                        arr3 = _mm512_min_epu16(arr1, arr2);
                        _mm512_mask_storeu_epi16((__m512i*) &predecessor_max_lcp[ch_num][i], k, arr3); 
                    }
                }    

                // DEBUG: START --------------------------------------------
                // uint16_t truncated_lcp_i = std::min((size_t) MAXLCPVALUE, lcp_i);
                
                // // update the dirty lcp table for each character
                // for (size_t i = 0; i < 256; i++) {
                //     dirty_lcp_cache[i] = std::min(dirty_lcp_cache[i], truncated_lcp_i);
                // }
                // dirty_lcp_cache[bwt_ch] = truncated_lcp_i;


                // // update the position corresponding the previous bwt ch since it was
                // // set to zero 
                // // if (last_updated_row < 256)
                // //     dirty_lcp_cache[last_updated_row] = truncated_lcp_i;

                // // update table with length of current suffix, because we are 
                // // going up so we want to update the current suffix with lcp_i
                // predecessor_max_lcp2[bwt_ch][doc_of_LF_i] = truncated_lcp_i;

                // //std::cout << "initialize: lcp2[" << bwt_ch << "][" << doc_of_LF_i << "] = " << predecessor_max_lcp2[bwt_ch][doc_of_LF_i] << std::endl;

                // // update the bwt_ch row of table with minimum lcp since last update
                // // uint16_t min_lcp_to_flush = dirty_lcp_cache[bwt_ch];
                // // for (size_t i = 0; i < num_docs; i++) {
                // //     predecessor_max_lcp2[bwt_ch][i] = std::min(predecessor_max_lcp2[bwt_ch][i], min_lcp_to_flush);
                // // }
                // //std::cout << "CHECK: lcp2[" << "T" << "][" << "4" << "] = " << predecessor_max_lcp2['T'][4] << std::endl;
                // //std::cout << "CHECK: dirty[T] = " << dirty_lcp_cache['T'] << std::endl;
                       
                // // update dirty lcp cache since we have already flush lcp_i to table
                // //dirty_lcp_cache[bwt_ch] = truncated_lcp_i;
                // //last_updated_row = bwt_ch;

                // DEBUG: END --------------------------------------------        
            #else        
                // Reset the LCP with respect to the current <ch, doc> pair
                predecessor_max_lcp[bwt_ch][doc_of_LF_i] = std::min((size_t) MAXLCPVALUE, lcp_i);

                for (size_t ch_num = 0; ch_num < 256; ch_num++) {
                    for (size_t doc_num = 0; doc_num < num_docs; doc_num++) {
                        predecessor_max_lcp[ch_num][doc_num] = std::min(predecessor_max_lcp[ch_num][doc_num], (uint16_t) std::min(lcp_i, (size_t) MAXLCPVALUE));
                    }
                }

                // DEBUG: START --------------------------------------------
                // uint16_t truncated_lcp_i = std::min((size_t) MAXLCPVALUE, lcp_i);
                
                // // update the dirty lcp table for each character
                // for (size_t i = 0; i < 256; i++) {
                //     dirty_lcp_cache[i] = std::min(dirty_lcp_cache[i], truncated_lcp_i);
                // }
                
                // // update the position corresponding the previous bwt ch since it was
                // // set to zero 
                // if (last_updated_row < 256)
                //     dirty_lcp_cache[last_updated_row] = truncated_lcp_i;

                // // update table with length of current suffix, because we are 
                // // going up so we want to update the current suffix with lcp_i
                // predecessor_max_lcp2[bwt_ch][doc_of_LF_i] = truncated_lcp_i;

                // // update the bwt_ch row of table with minimum lcp since last update
                // uint16_t min_lcp_to_flush = dirty_lcp_cache[bwt_ch];
                // for (size_t i = 0; i < num_docs; i++) {
                //     predecessor_max_lcp2[bwt_ch][i] = std::min(predecessor_max_lcp2[bwt_ch][i], min_lcp_to_flush);
                // }
                       
                // // update dirty lcp cache since we have already flush lcp_i to table
                // //dirty_lcp_cache[bwt_ch] = 0;
                // last_updated_row = bwt_ch;

                // DEBUG: END --------------------------------------------   

            #endif    
        }

        void update_predecessor_max_lcp_table_up_lazy_version(size_t lcp_i, size_t doc_of_LF_i, uint8_t bwt_ch) {
            /* 
             * Update the predecessor lcp table, this allows us to compute the maximum
             * lcp with respect to all the predecessor occurrences of other documents.
             * For example: if we are at <A, 2>, we will find the maximum lcp with the 
             * the previous occurrences of <A, 0> and <A, 1> for 3 documents.
             */
            uint16_t truncated_lcp_i = std::min((size_t) MAXLCPVALUE, lcp_i);
            
            // update the dirty lcp table for each character
            for (size_t i = 0; i < 256; i++) {
                dirty_lcp_cache[i] = std::min(dirty_lcp_cache[i], truncated_lcp_i);
            }
            dirty_lcp_cache[bwt_ch] = truncated_lcp_i;
            predecessor_max_lcp2[bwt_ch][doc_of_LF_i] = truncated_lcp_i;
        }

        void flush_row_of_lcp_table_up(uint8_t bwt_ch) {
            uint16_t min_lcp_to_flush = dirty_lcp_cache[bwt_ch];
            #if AVX512BW_PRESENT
                // initialize an lcp_i vector
                uint16_t lcp_i_vector[32];
                for (size_t i = 0; i < 32; i++)
                    lcp_i_vector[i] = min_lcp_to_flush;

                // load the array of constants (lcp_i)
                __m512i arr1, arr2, arr3; 
                __mmask32 k = ~0; // all 32 bits on, means all 32 values will be written
                arr2 = _mm512_maskz_loadu_epi16(~0, (const __m512i*) &lcp_i_vector[0]);          

                // use SIMD for all groups of 32
                for (size_t i = 0; i < (num_blocks_of_32 * 32); i+=32) {
                    // zero-mask, all the set bit positions are loaded
                    arr1 = _mm512_maskz_loadu_epi16(~0, (const __m512i*) &predecessor_max_lcp2[bwt_ch][i]); 
                    arr3 = _mm512_min_epu16(arr1, arr2);
                    _mm512_mask_storeu_epi16((__m512i*) &predecessor_max_lcp2[bwt_ch][i], k, arr3); 
                }

            #else
                for (size_t i = 0; i < num_docs; i++) {
                    predecessor_max_lcp2[bwt_ch][i] = std::min(predecessor_max_lcp2[bwt_ch][i], min_lcp_to_flush);
                }
            #endif
        }

        void initialize_current_row_profile(size_t doc_of_LF_i, std::vector<size_t>& curr_da_profile, uint8_t bwt_ch){
            /* 
             * Initialize the curr_da_profile with max lcp for predecessor 
             * occurrences of the same BWT character from another document, and
             * we check and make sure they occurred to avoid initializing it
             * with 1 (0 + 1 = 1)
             */
             //std::cout << "initialize: bwt_ch = " << bwt_ch << ", ";
            for (size_t i = 0; i < num_docs; i++) {
                if (i != doc_of_LF_i && ch_doc_encountered[bwt_ch][i]) {
                    curr_da_profile[i] = predecessor_max_lcp[bwt_ch][i] + 1;
                    
                    // if (predecessor_max_lcp[bwt_ch][i] != predecessor_max_lcp2[bwt_ch][i]) {
                    //     std::cout << "lcp[" << bwt_ch << "][" << i << "] = " << predecessor_max_lcp[bwt_ch][i] << ", ";
                    //     std::cout << "lcp2[" << bwt_ch << "][" << i << "] = " << predecessor_max_lcp2[bwt_ch][i] << ", \n";
                    //     std::cout << "j = " << j << std::endl;
                    // }
                    
                    // ASSERT((predecessor_max_lcp[bwt_ch][i] == predecessor_max_lcp2[bwt_ch][i]),
                    //         "difference found.");
                }
            }
        }

        void initialize_current_row_profile_lazy_version(size_t doc_of_LF_i, std::vector<size_t>& curr_da_profile, uint8_t bwt_ch){
            /* 
             * Initialize the curr_da_profile with max lcp for predecessor 
             * occurrences of the same BWT character from another document, and
             * we check and make sure they occurred to avoid initializing it
             * with 1 (0 + 1 = 1)
             */
            for (size_t i = 0; i < num_docs; i++) {
                if (i != doc_of_LF_i && ch_doc_encountered[bwt_ch][i]) {
                    curr_da_profile[i] = predecessor_max_lcp2[bwt_ch][i] + 1;
                }
            }
        }

        void write_data_to_temp_file(temp_data_entry_t data_entry) {
            /* write to temp lcp queue data to the file */
            size_t start_pos = num_lcp_temp_data * TEMPDATA_RECORD;
            mmap_lcp_inter[start_pos] = 0x00 | (data_entry.is_start << 1)  | (data_entry.is_end);
            mmap_lcp_inter[start_pos + 1] = data_entry.bwt_ch;
            
            // little-endian orientation
            ASSERT((data_entry.doc_num < MAXDOCS), "invalid document number encountered when writing temp file.");
            mmap_lcp_inter[start_pos + 2] = (data_entry.doc_num & 0xFF);
            mmap_lcp_inter[start_pos + 3] = (data_entry.doc_num & (0xFF << 8)) >> 8;

            size_t new_lcp_i = std::min((size_t) MAXLCPVALUE, data_entry.lcp_i);
            mmap_lcp_inter[start_pos + 4] = (new_lcp_i & 0xFF);
            mmap_lcp_inter[start_pos + 5] = (new_lcp_i & (0xFF << 8)) >> 8;
            
            num_lcp_temp_data += 1;
            ASSERT((max_lcp_records > num_lcp_temp_data), "the space in the temporary file has been used up.");
        }

        void write_profile_to_temp_file(std::vector<size_t>& curr_prof) {
            /* write the current profile to the temporary storage */
            size_t start_pos = num_dap_temp_data * DOCWIDTH;
            for (size_t i = 0; i < curr_prof.size(); i++) {
                uint16_t curr_lcp = std::min((size_t) MAXLCPVALUE, curr_prof[i]);
                mmap_dap_inter[start_pos + (i*DOCWIDTH)] = (0xFF & curr_lcp);
                mmap_dap_inter[start_pos + (i*DOCWIDTH) + 1] = ((0xFF << 8) & curr_lcp) >> 8;
            } 
            num_dap_temp_data += curr_prof.size();
            
            ASSERT((max_dap_records > (num_dap_temp_data+num_docs)), 
                    "the space in the document array temp file is not large enough.");
        }

        void delete_temp_files(std::string filename) {
            /* delete the temporary files */
            if (std::remove((filename + ".tmp_dap_data").data()))
                FATAL_ERROR("issue occurred while deleting temporary file."); 
            if (std::remove((filename + ".tmp_lcp_data").data()))
                FATAL_ERROR("issue occurred while deleting temporary file."); 
        }

        void print_dap(size_t num_entries) {
            /* writes the document-array profiles to file */
            std::vector<size_t> curr_profile;
            std::vector<size_t> left_increases, right_increases;
            curr_profile.reserve(num_docs);

            size_t dap_ptr = 0;
            for (size_t i = 0; i < num_entries; i++) {
                // Grab the data for current suffix
                bool is_start = GET_IS_START(mmap_lcp_inter, i);
                bool is_end = GET_IS_END(mmap_lcp_inter, i);
                uint8_t bwt_ch = GET_BWT_CH(mmap_lcp_inter, i);
                size_t doc_of_LF_i = GET_DOC_OF_LF(mmap_lcp_inter, i);

                // Load profiles if its a run-boundary
                if (is_start || is_end) {
                    
                    // option 1: write full document array profile
                    if (!use_taxcomp && !use_topk) {

                        // write out the BWT character
                        if (is_start && fwrite(&bwt_ch, 1, 1, sdap_file) != 1)
                            FATAL_ERROR("issue occurred while writing to *.sdap file");
                        if (is_end && fwrite(&bwt_ch, 1, 1, edap_file) != 1)
                            FATAL_ERROR("issue occurred while writing to *.sdap file");
                        
                        // write out the lcps
                        for (size_t j = 0; j < num_docs; j++) {
                            size_t prof_val = GET_DAP(mmap_dap_inter, (dap_ptr+j));
                            assert(prof_val <= MAXLCPVALUE);

                            if (is_start && fwrite(&prof_val, DOCWIDTH, 1, sdap_file) != 1)
                                FATAL_ERROR("issue occurred while writing to *.sdap file");
                            if (is_end && fwrite(&prof_val, DOCWIDTH, 1, edap_file) != 1)
                                FATAL_ERROR("issue occurred while writing to *.edap file");
                        }
                    // option 2: use taxonomic compression
                    } else if (use_taxcomp) {
                        size_t curr_val = 0; 
                        size_t prev_max = std::min(GET_DAP(mmap_dap_inter, (dap_ptr)), MAXLCPVALUE);
                        left_increases.push_back(0);
                        left_increases.push_back(prev_max);

                        // keep track of increases from left
                        for (size_t j = 0; j < num_docs; j++) {
                            curr_val = GET_DAP(mmap_dap_inter, (dap_ptr+j));
                            ASSERT((curr_val <= MAXLCPVALUE), "issue occurred in writing of DAP. (taxcomp - 1)");
                            curr_profile[j] = curr_val;

                            if (curr_val > prev_max) {
                                left_increases.push_back(j);
                                left_increases.push_back(curr_val);
                                prev_max = curr_val;
                            }
                        }
                        ASSERT((left_increases.size() % 2 == 0), "issue occurred in writing of DAP. (taxcomp - 2)");
                        size_t num_left_inc = left_increases.size()/2;
                        size_t num_right_inc = 0;

                        // gather the increases from right
                        prev_max = curr_profile[num_docs-1];
                        right_increases.push_back(num_docs-1);
                        right_increases.push_back(prev_max);

                        for (int j = num_docs-1; j >= 0; j--) {
                            if (curr_profile[j] > prev_max) {
                                right_increases.push_back(j);
                                right_increases.push_back(curr_profile[j]);
                                prev_max = curr_profile[j];
                            }
                        }
                        ASSERT((right_increases.size() % 2 == 0), "issue occurred in writing of DAP. (taxcomp - 3)");
                        num_right_inc = right_increases.size() / 2;

                        size_t old_sdap_overflow_pos = tax_sdap_overflow_ptr;
                        size_t old_edap_overflow_pos = tax_edap_overflow_ptr;

                        // number of overflow pairs
                        size_t num_of_pairs = (num_left_inc > NUMCOLSFORTABLE) ? (num_left_inc-NUMCOLSFORTABLE) : 0;
                        num_of_pairs += (num_right_inc > NUMCOLSFORTABLE) ? (num_right_inc-NUMCOLSFORTABLE) : 0;

                        // write to start-runs file
                        if (is_start) {
                            // Step 1: Write the BWT char
                            if (fwrite(&bwt_ch, 1, 1, sdap_tax) != 1)
                                FATAL_ERROR("issue occurred while writing char to *taxcomp.sdap file");
                            
                            // Step 2: Writes the ordered pairs to the file
                            for (size_t i = 0; i < NUMCOLSFORTABLE; i++) 
                                write_to_taxcomp_dap(sdap_tax, left_increases, right_increases, i);
                            
                            // Step 3: Writes the overflow pointer to a file in order to be able to go to overflow file
                            append_overflow_pointer(sdap_tax, tax_sdap_overflow_ptr, num_left_inc, num_right_inc);
                            
                            // Step 4: Write any leftover pairs to overflow table
                            if (num_left_inc > NUMCOLSFORTABLE || num_right_inc > NUMCOLSFORTABLE) {
                                tax_sdap_overflow_ptr = write_remaining_pairs_to_overflow(sdap_overtax, left_increases, tax_sdap_overflow_ptr);
                                tax_sdap_overflow_ptr = write_remaining_pairs_to_overflow(sdap_overtax, right_increases, tax_sdap_overflow_ptr);
                            }

                            // Check that overflow ptr is being incremented as expected
                            if (num_of_pairs > 0)
                                ASSERT((tax_sdap_overflow_ptr == (old_sdap_overflow_pos + 2 + (DOCWIDTH * 2 * num_of_pairs))),  "issue occurred in writing of DAP. (taxcomp - 4)");
                            else
                                ASSERT((tax_sdap_overflow_ptr == old_sdap_overflow_pos), "issue occurred in writing of DAP. (taxcomp - 5)");
                        }
                        // write to end-runs file
                        if (is_end) {
                            // Step 1: Write the BWT char
                            if (fwrite(&bwt_ch, 1, 1, edap_tax) != 1)
                                FATAL_ERROR("issue occurred while writing char to *taxcomp.edap file");

                            // Step 2: Writes the ordered pairs to the file
                            for (size_t i = 0; i < NUMCOLSFORTABLE; i++) 
                                write_to_taxcomp_dap(edap_tax, left_increases, right_increases, i);
                            
                            // Step 3: Writes the overflow pointer to a file in order to be able to go to overflow file
                            append_overflow_pointer(edap_tax, tax_edap_overflow_ptr, num_left_inc, num_right_inc);
                            
                            // Step 4: Write any leftover pairs to overflow table
                            if (num_left_inc > NUMCOLSFORTABLE || num_right_inc > NUMCOLSFORTABLE) {
                                tax_edap_overflow_ptr = write_remaining_pairs_to_overflow(edap_overtax, left_increases, tax_edap_overflow_ptr);
                                tax_edap_overflow_ptr = write_remaining_pairs_to_overflow(edap_overtax, right_increases, tax_edap_overflow_ptr);
                            }

                            // Check that overflow ptr is being incremented as expected
                            if (num_of_pairs > 0)
                                ASSERT ((tax_edap_overflow_ptr == (old_edap_overflow_pos + 2 + (DOCWIDTH * 2 * num_of_pairs))), "issue occurred in writing of DAP. (taxcomp - 6)");
                            else
                                ASSERT ((tax_edap_overflow_ptr == old_edap_overflow_pos), "issue occurred in writing of DAP. (taxcomp - 7)");
                        }

                        // empty the vectors
                        left_increases.clear();
                        right_increases.clear();
                    // option 3: top-k compression
                    } else {
                        FATAL_ERROR("Not implemented yet ...");
                    }

                    // increment after writing the profile
                    dap_ptr += num_docs;
                }
            }
        }

        void write_to_taxcomp_dap(FILE* outfile, std::vector<size_t>& left_inc, std::vector<size_t>& right_inc, size_t col_num){
            size_t index = col_num * 2;
            size_t left_pos = (col_num < left_inc.size()/2) ? left_inc[index]: MAXLCPVALUE;
            size_t left_lcp = (col_num < left_inc.size()/2) ? left_inc[index+1]: MAXLCPVALUE;
            size_t right_pos = (col_num < right_inc.size()/2) ? right_inc[index]: MAXLCPVALUE;
            size_t right_lcp = (col_num < right_inc.size()/2) ? right_inc[index+1]: MAXLCPVALUE;

            bool success = fwrite(&left_pos, DOCWIDTH, 1, outfile) == 1;
            success &= fwrite(&left_lcp, DOCWIDTH, 1, outfile) == 1;
            success &= fwrite(&right_pos, DOCWIDTH, 1, outfile) == 1;
            success &= fwrite(&right_lcp, DOCWIDTH, 1, outfile) == 1;
            if (!success)
                FATAL_ERROR("issue occurred during the writing to the *.taxcomp.dap file");
        }

        void append_overflow_pointer(FILE* outfile, size_t curr_pos, size_t num_left_inc, size_t num_right_inc){
            // When there is overflow, we write the pointer to the starting position
            // in overflow file, otherwise it will just be 0. Take note that 0 would
            // never be a valid position since we initialize the overflow document
            // with number of documents so first possible position is 8.
            if (num_left_inc > NUMCOLSFORTABLE || num_right_inc > NUMCOLSFORTABLE) {
                if (fwrite(&curr_pos, sizeof(size_t), 1, outfile) != 1)
                    FATAL_ERROR("issue occurred when writing the overflow pointer.");
            } else {
                size_t overflow_pos = 0;
                if (fwrite(&overflow_pos, sizeof(size_t), 1, outfile) != 1)
                    FATAL_ERROR("issue occurred when writing the overflow pointer.");
            }
        }

        size_t write_remaining_pairs_to_overflow(FILE* outfile, std::vector<size_t>& inc_pairs, size_t curr_ptr_pos){
           /* 
            * Note: 
            * This method is only called when the document array profile has an 
            * overflow either in the left to right or right to left direction so
            * we always write something to the overflow even if it just 1 byte saying
            * there are no leftover pairs in one of the direction in order to guarantee
            * we can find all the data.
            */
            int diff = (inc_pairs.size()/2 - NUMCOLSFORTABLE);
            size_t num_to_write = (diff > 0) ? diff : 0; 
            ASSERT((num_to_write < 256), "we need less than 255 values to compress the profiles.");

            if (fwrite(&num_to_write, 1, 1, outfile) != 1)
                FATAL_ERROR("issue occurred while writing to overflow table.");
            curr_ptr_pos += 1;

            if (num_to_write > 0) {
                for (size_t j = NUMCOLSFORTABLE; j < inc_pairs.size()/2; j++) {
                    bool success = (fwrite(&inc_pairs[j*2], DOCWIDTH, 1, outfile) == 1);
                    success &= fwrite(&inc_pairs[j*2+1], DOCWIDTH, 1, outfile) == 1;
                    if (!success)
                        FATAL_ERROR("issue occurred while writing the left data to overflow table.");
                    curr_ptr_pos += (DOCWIDTH * 2);
                }
            }
            return curr_ptr_pos;
        }

        /***********************************************************/
        /* Section 2: BWT, LCP, and SA construction related methods
        /***********************************************************/
        inline bool inc(phrase_suffix_t& s){
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

        inline bool is_valid(phrase_suffix_t& s){
            // avoid the extra w # at the beginning of the text
            if (s.sn < pf.w)
                return false;
            // Check if the suffix has length at least w and is not the complete phrase.
            if (pf.dict.b_d[s.sn] != 0 || s.suffix_length < pf.w)
                return false;
            
            return true;
        }
        
        inline int_t min_s_lcp_T(size_t left, size_t right){
            // assume left < right
            if (left > right)
                std::swap(left, right);

            assert(pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] >= pf.w);

            return (pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] - pf.w);
        }

        inline int_t compute_lcp_suffix(phrase_suffix_t& curr, phrase_suffix_t& prev){
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

        inline void print_lcp(int_t val, size_t pos){
            size_t tmp_val = val;
            if (fwrite(&tmp_val, THRBYTES, 1, lcp_file) != 1)
                error("LCP write error 1");
        }

        inline void new_min_s(int_t val, size_t pos){
            // We can put here the check if we want to store the LCP or stream it out
            min_s.push_back(val);
            pos_s.push_back(j);
        }

        inline void update_ssa(phrase_suffix_t &curr, size_t pos){
            ssa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
            assert(ssa < (pf.n - pf.w + 1ULL));
        }

        inline void update_esa(phrase_suffix_t &curr, size_t pos){
            esa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
            assert(esa < (pf.n - pf.w + 1ULL));
        }

        inline void print_sa(){
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

        inline void print_bwt(){   
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

        inline void update_bwt(uint8_t next_char, size_t length_){
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