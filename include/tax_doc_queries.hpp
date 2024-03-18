/*
 * File: tax_doc_queries.hpp
 * Description: ....
 * Date: Aug. 30th, 2023
 */

#include <minimizer_digest.hpp>

#ifndef _TAX_DOC_QUERIES_H
#define _TAX_DOC_QUERIES_H

#define READ_NUM_DOCS(fd) (0xFF & fd[0]) | ((0xFF & fd[1]) << 8) \
                        | ((0xFF & fd[2]) << 16) | ((0xFF & fd[3]) << 24) \
                        | ((0xFF & fd[4]) << 32) | ((0xFF & fd[5]) << 40) \
                        | ((0xFF & fd[6]) << 48) | ((0xFF & fd[7]) << 56) 
#define READ_NUM_PAIRS(fd, pos) (fd[pos] & 0xFF) 
#define READ_DOC_ID_OR_LCP_VAL(fd, pos) (0xFF & fd[pos]) | ((0xFF & fd[pos+1]) << 8)

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd>
class tax_doc_queries : ri::r_index<sparse_bv_type, rle_string_t>
{
    public:
        bool print_profiles = false;
        std::string output_csv_path = "";
        size_t profiles_to_print = 0;

        std::ofstream csv_sdap_output;
        std::ofstream csv_edap_output;

        int sdap_of_fd;
        int edap_of_fd;
        struct stat file_info_sdap;
        struct stat file_info_edap;

        char* mmap_sdap_of;
        char* mmap_edap_of;

        // This vectors has the following dimensions: [256][num of ith char][num_docs]
        // This structure stores the DA profiles for each character separately.
        std::vector<std::vector<std::vector<uint64_t>>> start_doc_profiles;
        std::vector<std::vector<std::vector<uint64_t>>> end_doc_profiles;

        tax_doc_queries(std::string filename,
                        size_t num_cols,
                        size_t num_profiles = 0,
                        std::string output_path = "",
                        bool rle = true):
                        profiles_to_print(num_profiles),
                        output_ref(filename),
                        rle (rle),
                        num_cols (num_cols),
                        ri::r_index<sparse_bv_type, rle_string_t>(),
                        start_doc_profiles(256, std::vector<std::vector<uint64_t>>(0, std::vector<uint64_t>(0))),
                        end_doc_profiles(256, std::vector<std::vector<uint64_t>>(0, std::vector<uint64_t>(0)))
        {
            // load the BWT 
            STATUS_LOG("query_main", "loading the bwt of the input text");
            auto start = std::chrono::system_clock::now();

            std::string bwt_fname = filename + ".bwt";
            load_bwt_structure(bwt_fname);
            DONE_LOG((std::chrono::system_clock::now() - start));

            // gather some statistics on the BWT
            this->r = this->bwt.number_of_runs();
            ri::ulint n = this->bwt.size();
            size_t log_r = bitsize(uint64_t(this->r));
            size_t log_n = bitsize(uint64_t(this->bwt.size()));
            FORCE_LOG("query_main", "bwt statistics: n = %ld, r = %ld\n" , this->bwt.size(), this->r);

            // check file size and make sure it is valid
            check_doc_array_files(filename + ".taxcomp.sdap");
            check_doc_array_files(filename + ".taxcomp.sdap");

            // check if we are trying to print out profiles
            output_csv_path.assign(output_path);
            if (output_csv_path.length()) {
                print_profiles = true;
                this->profiles_to_print = (num_profiles == 0 || num_profiles > this->r) ? this->r : num_profiles;
                csv_sdap_output.open(output_csv_path + ".sdap.csv");
                csv_edap_output.open(output_csv_path + ".edap.csv");
            }

            // read the document array main table
            STATUS_LOG("build_profiles", "loading the document array profiles");
            start = std::chrono::system_clock::now();

            read_doc_profiles_main_table(start_doc_profiles, filename + ".taxcomp.sdap", &csv_sdap_output);
            read_doc_profiles_main_table(end_doc_profiles, filename + ".taxcomp.edap", &csv_edap_output);

            if (print_profiles) {
                csv_sdap_output.close();
                csv_edap_output.close();
            }
            DONE_LOG((std::chrono::system_clock::now() - start));

            // load the overflow files
            std::string sdap_of_path = filename + ".taxcomp.of.sdap";
            std::string edap_of_path = filename + ".taxcomp.of.edap";

            sdap_of_fd = open(sdap_of_path.data(), O_RDONLY);
            edap_of_fd = open(edap_of_path.data(), O_RDONLY);
            if (sdap_of_fd == -1 || edap_of_fd == -1)
                FATAL_ERROR("Issue with opening up the overflow files.");

            if (fstat(sdap_of_fd, &file_info_sdap) == -1)
                FATAL_ERROR("Issue with checking file size on overflow file.");
            if (fstat (edap_of_fd, &file_info_edap) == -1)
                FATAL_ERROR("Issue with checking file size on overflow file.");
        
            mmap_sdap_of = (char*) mmap(NULL, file_info_sdap.st_size, PROT_READ, MAP_PRIVATE, sdap_of_fd, 0);
            mmap_edap_of = (char*) mmap(NULL, file_info_edap.st_size, PROT_READ, MAP_PRIVATE, edap_of_fd, 0);
            assert(mmap_sdap_of != MAP_FAILED);
            assert(mmap_edap_of != MAP_FAILED);
        }

        ~tax_doc_queries()
        {
            /* Destructor - unmap files and close files */
            munmap(mmap_sdap_of, file_info_sdap.st_size);
            munmap(mmap_edap_of, file_info_edap.st_size);
            close(sdap_of_fd); close(edap_of_fd);
            csv_sdap_output.close(); csv_edap_output.close();
        }

        void query_profiles(std::string pattern_file, 
                            std::string output_prefix,
                            ref_type database_type,
                            size_t small_window_l,
                            size_t large_window_l){
            /* Go through and query the reads and print the leftmost, rightmost occurrence of string */

            STATUS_LOG("query_main", "processing the patterns");
            auto start_time = std::chrono::system_clock::now();

            // open output/input files
            std::ofstream listings_fd (output_prefix + ".listings");
            gzFile fp; kseq_t* seq;

            fp = gzopen(pattern_file.data(), "r"); 
            if(fp == 0) {std::exit(1);}
            seq = kseq_init(fp);

            // variable used keep track of position in lamdbas
            size_t idx = 0;

            // alternate lambda to print out the document listing, that
            // will be tested to see if it improves classification
            auto process_profile_with_subtree = [&](std::vector<uint64_t> profile, uint16_t length, bool use_end) {
                    bool left_done = false, right_done = false;
                    std::vector<uint64_t> left_docs;
                    std::vector<uint64_t> right_docs;
                    size_t tuple_num = 0;

                    // std::cout << "length: " << length << std::endl;
                    // std::cout << "profile: ";
                    // for (auto x: profile)
                    //     std::cout << x << " ";
                    // std::cout << "\n";

                    // Iterates through the main table ...
                    while ((!left_done || !right_done) && (tuple_num < (this->num_cols * 4))) {   

                        // Checks if the LCP value has reached the max, and
                        // the monotonic increases have ended
                        if (profile[tuple_num] == MAXLCPVALUE)
                            left_done = true;
                        if (profile[tuple_num+2] == MAXLCPVALUE)
                            right_done = true;

                        // Check if increases in either direction are still going, and 
                        // lcp is greater than query length
                        if (!left_done && profile[tuple_num+1] >= length) {
                            left_docs.push_back(profile[tuple_num]);
                        }
                        if (!right_done && profile[tuple_num+3] >= length) {
                            right_docs.push_back(profile[tuple_num+2]);
                        }
                        tuple_num += 4;
                    }

                    // Double-check that if we haven't found the leftmost or rightmost
                    // document that there are still some documents to look at in overflow file
                    if ((left_docs.size() == 0 && left_done) || (right_docs.size() == 0 && right_done))
                        FATAL_ERROR("Issue occurred when querying the taxonomic document array.");

                    // Grab the overflow pointer based on which type of profiles was used
                    char* of_ptr = nullptr;
                    if (use_end)
                        of_ptr = mmap_edap_of;
                    else
                        of_ptr = mmap_sdap_of;
                    size_t of_pos = profile[this->num_cols * 4];

                    // Make sure that both directions are done based on overflow ptr
                    if (of_pos == 0) {
                        left_done = true;
                        right_done = true;
                        assert((left_done && right_done));
                    }
                    // std::cout << leftmost_found << " " << rightmost_found << std::endl;
                    // std::cout << of_pos << std::endl;
                    
                    // Use overflow file for left to right direction ...
                    if (!left_done) {
                        // std::cout << "left is not done!" << std::endl;
                        assert(of_pos != 0);

                        size_t num_left_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_left_pairs < 256);
                        of_pos += 1;
                        // std::cout << "num_pairs_left = " << num_left_pairs << "\n";

                        // Go through all the L2R pairs
                        for (size_t i = 0; i < num_left_pairs; i++){
                            size_t doc_id = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            size_t lcp_val = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            // std::cout << "doc = " << doc_id << ", lcp = " << lcp_val << std::endl;

                            if (lcp_val >= length) {
                                left_docs.push_back(doc_id);
                                //leftmost_doc = doc_id;
                                //leftmost_found = true;
                                //break;
                            }
                        }
                    }
                    // Use overflow file for right to left direction ...
                    if (!right_done) {
                        // std::cout << "right is not done!" << std::endl;
                        of_pos = profile[this->num_cols * 4];
                        assert(of_pos != 0);

                        // Move past the left pairs to right pairs
                        size_t num_left_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_left_pairs < 256);
                        of_pos += 1 + (DOCWIDTH * 2 * num_left_pairs);

                        size_t num_right_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_right_pairs < 256);
                        of_pos += 1;
                        // std::cout << "num_pairs_right = " << num_right_pairs << "\n";

                        // Go through all the R2L pairs
                        for (size_t i = 0; i < num_right_pairs; i++){
                            size_t doc_id = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            size_t lcp_val = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            // std::cout << "doc = " << doc_id << ", lcp = " << lcp_val << std::endl;
                            
                            if (lcp_val >= length) {
                                right_docs.push_back(doc_id);
                                //rightmost_doc = doc_id;
                                //rightmost_found = true;
                                //break;
                            }
                        }
                    }
                    // Make sure we have found both documents in 
                    // both directions, and they end with same doc
                    assert((left_docs.size() && right_docs.size()));
                    assert ((left_docs.back() == right_docs.back()));

                    // Remove the last document from one of the vectors
                    left_docs.pop_back();

                    // std::cout << "left docs: ";
                    // for (auto x: left_docs)
                    //     std::cout << x << " ";
                    // std::cout << "\n";
                    // std::cout << "right docs: ";
                    // for (auto x: right_docs)
                    //     std::cout << x << " ";
                    // std::cout << "\n";

                    // Output the leftmost/rightmost nodes along with
                    // nodes underneath with hit
                    std::string output_str = "{";
                    for (auto x: left_docs)
                        output_str += std::to_string(x) + ",";
                    for (int i = right_docs.size()-1; i >= 0; i--)
                        output_str += std::to_string(right_docs[i]) + ",";

                    // Remove comma and close bracket
                    output_str.pop_back();
                    output_str += "} ";
                    listings_fd << output_str;
            };

            // lambda to print out the document listing
            auto process_profile = [&](std::vector<uint64_t> profile, uint16_t length, bool use_end) {
                    bool leftmost_found = false, rightmost_found = false;
                    bool left_done = false, right_done = false;

                    size_t leftmost_doc = 0, rightmost_doc = 0;
                    size_t tuple_num = 0;

                    // std::cout << "\nlength = " << length << std::endl;
                    // std::cout << "\nprofile: ";
                    // for (auto x: profile)
                    //     std::cout << x << " ";
                    // std::cout << "\n";

                    // Iterates through the main table ...
                    while ((!leftmost_found || !rightmost_found) && (tuple_num < (this->num_cols * 4))) {
                        // Check if the increases in each direction is done or not
                        if (profile[tuple_num] == MAXLCPVALUE)
                            left_done = true;
                        if (profile[tuple_num+2] == MAXLCPVALUE)
                            right_done = true;
                        
                        // Check if increases in each direction is not done, and we haven't
                        // found an matching document, and then compare lengths
                        if (!left_done && !leftmost_found && profile[tuple_num+1] >= length) {
                            leftmost_doc = profile[tuple_num];
                            leftmost_found = true;
                        }
                        if (!right_done && !rightmost_found && profile[tuple_num+3] >= length) {
                            rightmost_doc = profile[tuple_num+2];
                            rightmost_found = true;
                        }
                        tuple_num += 4;
                    }

                    // Double-check that if we haven't found the leftmost or rightmost
                    // document that there are still some documents to look at in overflow file
                    if ((!leftmost_found && left_done) || (!rightmost_found && right_done))
                        FATAL_ERROR("Issue occurred when querying the taxonomic document array.");
                    
                    // Grab the overflow pointer based on which type of profiles was used
                    char* of_ptr = nullptr;
                    if (use_end)
                        of_ptr = mmap_edap_of;
                    else
                        of_ptr = mmap_sdap_of;
                    size_t of_pos = profile[this->num_cols * 4];

                    // std::cout << leftmost_found << " " << rightmost_found << std::endl;
                    // std::cout << of_pos << std::endl;
                    
                    // Use overflow file for left to right direction ...
                    if (!leftmost_found && !left_done) {
                        assert(of_pos != 0);

                        size_t num_left_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_left_pairs < 256);
                        of_pos += 1;
                        // std::cout << "num_pairs_left = " << num_left_pairs << "\n";

                        // Go through all the L2R pairs
                        for (size_t i = 0; i < num_left_pairs; i++){
                            size_t doc_id = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            size_t lcp_val = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            // std::cout << "doc = " << doc_id << ", lcp = " << lcp_val << std::endl;

                            if (lcp_val >= length) {
                                leftmost_doc = doc_id;
                                leftmost_found = true;
                                break;
                            }
                        }
                    }
                    // Use overflow file for right to left direction ...
                    if (!rightmost_found && !right_done) {
                        of_pos = profile[this->num_cols * 4];
                        assert(of_pos != 0);

                        // Move past the left pairs to right pairs
                        size_t num_left_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_left_pairs < 256);
                        of_pos += 1 + (DOCWIDTH * 2 * num_left_pairs);

                        size_t num_right_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_right_pairs < 256);
                        of_pos += 1;
                        // std::cout << "num_pairs_right = " << num_right_pairs << "\n";

                        // Go through all the R2L pairs
                        for (size_t i = 0; i < num_right_pairs; i++){
                            size_t doc_id = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            size_t lcp_val = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            // std::cout << "doc = " << doc_id << ", lcp = " << lcp_val << std::endl;
                            
                            if (lcp_val >= length) {
                                rightmost_doc = doc_id;
                                rightmost_found = true;
                                break;
                            }
                        }
                    }
                    // Make sure we have found both directions
                    assert((rightmost_found && leftmost_found));

                    // Output the leftmost/rightmost nodes
                    std::string output_str = "";
                    output_str += "{" + std::to_string(leftmost_doc) + 
                                  "," + std::to_string(rightmost_doc) + "} ";
                    listings_fd << output_str;

                    // std::cout << output_str << "\n";
                    // std::exit(1);
            };

            // build minimizer digest object (only used if needed)
            MinimizerDigest digester; 
            digester.set_lexorder(false);
            digester.set_windows(small_window_l, large_window_l);
            if (database_type == MINIMIZER) digester.set_minimizer_alp(true);

            // process each read, and print out the document lists
            while (kseq_read(seq)>=0) {

                // uppercase every letter in read
                for (int j = 0; j < seq->seq.l; j++) {seq->seq.s[j] = up_tab[(int) seq->seq.s[j]];}

                // digest read if needed
                std::string input_read = seq->seq.s;
                size_t read_length = 0;

                if (database_type == DNA_MINIMIZER) {input_read = digester.compute_digest(input_read);}
                else if (database_type == MINIMIZER) {input_read = digester.compute_digest(input_read);}
                read_length = input_read.size();

                listings_fd << ">" << seq->name.s << "\n";

                // initialize variables to use for backward search
                size_t start = 0, end = this->bwt.size(), end_pos_of_match = read_length-1;
                size_t start_run = 0, end_run = this->r - 1;
                size_t num_ch_before_start = 0, num_ch_before_end = 0;
                std::vector<uint64_t> curr_profile (this->num_cols * 4, 0);
                uint16_t length = 0;

                // pointer variables to current profile being "used"
                uint8_t curr_prof_ch = 0;
                size_t curr_prof_pos = 0, num_LF_steps = 0;

                // tell us what type of profile to grab based on pointer variables
                bool use_start = false, use_end = false;
                bool pointer_set = false;

                // perform backward search and report document listings when
                // range goes empty or we reach the end
                for (int i = (read_length-1); i >= 0; i--) {
                    
                    // extract the next char
                    uint8_t next_ch = input_read[i];

                    // special case: character is not present in text
                    if (this->bwt.number_of_letter(next_ch) == 0) {
                        // step 1: print listing for [i+1,end]
                        if (pointer_set) {
                            if (use_end)
                                curr_profile = end_doc_profiles[curr_prof_ch][curr_prof_pos];
                            else
                                curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];

                            idx = 0;
                            std::for_each(curr_profile.begin(), 
                                            curr_profile.end(), 
                                            [&](uint64_t &x){
                                            if (idx % 2 == 1) 
                                            {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                                            idx++;});

                            listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";

                            length = std::min((size_t) MAXLCPVALUE, (end_pos_of_match-i));
                            process_profile_with_subtree(curr_profile, length, use_end);
                            end_pos_of_match = i;
                        }
                        // step 2: print an empty listing for current character
                        listings_fd << "[" << i << "," << i << "] {} ";
                        end_pos_of_match = i - 1;

                        // step 3: reset variables to full bwt range
                        start = 0; end = this->bwt.size();
                        start_run = 0; end_run = this->r - 1;
                        num_ch_before_start = 0;
                        num_ch_before_end = 0;

                        // step 4: reset pointer and go back to top of loop
                        pointer_set = false;
                        continue;
                    }

                    // identify the run # that start and end are in
                    num_ch_before_start = this->bwt.rank(start, next_ch); 
                    num_ch_before_end = this->bwt.rank(end, next_ch);
                    start_run = this->bwt.run_of_position(start);
                    end_run = this->bwt.run_of_position(end-1);
                    
                    // range spans runs, so there are different BWT characters
                    if (start_run != end_run) 
                    {       
                        // bwt range is empty, so will reset start and end
                        if (num_ch_before_end == num_ch_before_start) {

                            // grab the correct current profile, and update with steps
                            if (use_end)
                                curr_profile = end_doc_profiles[curr_prof_ch][curr_prof_pos];
                            else
                                curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];

                            idx = 0;
                            std::for_each(curr_profile.begin(), 
                                          curr_profile.end(), 
                                          [&](uint64_t &x){
                                          if (idx % 2 == 1) 
                                          {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                                          idx++;});
                            listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";

                            length = std::min((size_t) MAXLCPVALUE, (end_pos_of_match-i));
                            process_profile_with_subtree(curr_profile, length, use_end);
                            end_pos_of_match = i;

                            start = 0; end = this->bwt.size();
                            start_run = 0; end_run = this->r - 1;
                            num_ch_before_start = 0;
                            num_ch_before_end = this->bwt.number_of_letter(next_ch);
                        }

                        // grab any profile at run boundary of next_ch (choose first one)
                        curr_prof_ch = next_ch;
                        curr_prof_pos = this->bwt.run_head_rank(start_run, next_ch);
                        num_LF_steps = 0;
                        use_start = false; use_end = false;
                        pointer_set = true;

                        // if the start position run is the same as query
                        // ch, then we can guarantee that end of run is in the range
                        // otherwise, we can guarantee the start of run is in range.
                        if (this->bwt[start] == next_ch)
                            use_end = true;
                        else
                            use_start = true;
                    } 
                    // range is within BWT run, but wrong character 
                    else if (this->bwt[start] != next_ch) 
                    {
                        // grab the correct current profile, and update with steps
                        if (use_end)
                            curr_profile = end_doc_profiles[curr_prof_ch][curr_prof_pos];
                        else
                            curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];

                        idx = 0;
                        std::for_each(curr_profile.begin(), 
                                        curr_profile.end(), 
                                        [&](uint64_t &x){
                                        if (idx % 2 == 1) 
                                        {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                                        idx++;});
                        listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";

                        length = std::min((size_t) MAXLCPVALUE, (end_pos_of_match-i));
                        process_profile_with_subtree(curr_profile, length, use_end);
                        end_pos_of_match = i;

                        start = 0; end = this->bwt.size();
                        start_run = 0; end_run = this->r - 1;
                        num_ch_before_start = 0;
                        num_ch_before_end = this->bwt.number_of_letter(next_ch);
                        pointer_set = true;

                        // grab any profile at run boundary of next_ch (choose first one)
                        curr_prof_ch = next_ch;
                        curr_prof_pos = this->bwt.run_head_rank(start_run, next_ch);
                        num_LF_steps = 0;
                        use_start = false; use_end = false;

                        // if the start position run is the same as query
                        // ch, then we can guarantee that end of run is in the range
                        // otherwise, we can guarantee the start of run is in range.
                        if (this->bwt[start] == next_ch)
                            use_end = true;
                        else
                            use_start = true;
                    }
                    // range is within BWT run, and is the correct character
                    else 
                    {
                        num_LF_steps++;
                    }

                    // Perform an LF step
                    start = num_ch_before_start + this->F[next_ch]; 
                    end = num_ch_before_end + this->F[next_ch];
                }
                // grab the current profile, and update with steps
                if (use_end)
                    curr_profile = end_doc_profiles[curr_prof_ch][curr_prof_pos];
                else
                    curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];

                // Update profile based on LF steps
                idx = 0;
                std::for_each(curr_profile.begin(), 
                              curr_profile.end(), 
                              [&](uint64_t &x){
                              if (idx % 2 == 1) 
                              {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                              idx++;});

                listings_fd << "[" << 0 << "," << end_pos_of_match << "] ";
                length = std::min((size_t) MAXLCPVALUE, end_pos_of_match+1);
                process_profile_with_subtree(curr_profile, length, use_end);
                listings_fd << "\n";
            }
            DONE_LOG((std::chrono::system_clock::now() - start_time));
        }

        void query_profiles_optimized(std::string pattern_file, 
                                      std::string output_prefix,
                                      ref_type database_type,
                                      size_t small_window_l,
                                      size_t large_window_l) {
            /* query reads using the optimized algorithm to accelerate queries */

            STATUS_LOG("query_main", "processing the patterns");
            auto start_time = std::chrono::system_clock::now();

            // open output/input files
            std::ofstream listings_fd (output_prefix + ".listings");
            gzFile fp; kseq_t* seq;
            fp = gzopen(pattern_file.data(), "r"); 

            if(fp == 0) {std::exit(1);}
            seq = kseq_init(fp);

            // initialize variables for backward search
            size_t start = 0, end = this->bwt.size(), end_pos_of_match = 0;
            size_t start_run = 0, end_run = 0;
            size_t num_ch_before_start = 0, num_ch_before_end = 0;

            std::vector<uint64_t> curr_profile (this->num_cols * 4, 0);
            uint16_t length = 0; uint8_t next_ch = 0;
            int i = 0; size_t read_length = 0;
            bool pointer_set = false;

            // variable used keep track of position in lamdbas
            size_t idx = 0;

            // initialize variables for identifying correct document array profile
            uint8_t curr_prof_ch = 0;
            size_t curr_prof_pos = 0, num_LF_steps = 0;

            // initialize variables for which profile to use
            bool use_start = false, use_end = false;
        
            // alternate lambda to print out the document listing, that
            // will be tested to see if it improves classification
            auto process_profile_with_subtree = [&](std::vector<uint64_t> profile, uint16_t length, bool use_end) {
                    bool left_done = false, right_done = false;
                    std::vector<uint64_t> left_docs;
                    std::vector<uint64_t> right_docs;
                    size_t tuple_num = 0;

                    // std::cout << "length: " << length << std::endl;
                    // std::cout << "profile: ";
                    // for (auto x: profile)
                    //     std::cout << x << " ";
                    // std::cout << "\n";

                    // Iterates through the main table ...
                    while ((!left_done || !right_done) && (tuple_num < (this->num_cols * 4))) {   

                        // Checks if the LCP value has reached the max, and
                        // the monotonic increases have ended
                        if (profile[tuple_num] == MAXLCPVALUE)
                            left_done = true;
                        if (profile[tuple_num+2] == MAXLCPVALUE)
                            right_done = true;

                        // Check if increases in either direction are still going, and 
                        // lcp is greater than query length
                        if (!left_done && profile[tuple_num+1] >= length) {
                            left_docs.push_back(profile[tuple_num]);
                        }
                        if (!right_done && profile[tuple_num+3] >= length) {
                            right_docs.push_back(profile[tuple_num+2]);
                        }
                        tuple_num += 4;
                    }

                    // Double-check that if we haven't found the leftmost or rightmost
                    // document that there are still some documents to look at in overflow file
                    if ((left_docs.size() == 0 && left_done) || (right_docs.size() == 0 && right_done))
                        FATAL_ERROR("Issue occurred when querying the taxonomic document array.");

                    // Grab the overflow pointer based on which type of profiles was used
                    char* of_ptr = nullptr;
                    if (use_end)
                        of_ptr = mmap_edap_of;
                    else
                        of_ptr = mmap_sdap_of;
                    size_t of_pos = profile[this->num_cols * 4];

                    // Make sure that both directions are done based on overflow ptr
                    if (of_pos == 0) {
                        left_done = true;
                        right_done = true;
                        assert((left_done && right_done));
                    }
                    // std::cout << leftmost_found << " " << rightmost_found << std::endl;
                    // std::cout << of_pos << std::endl;
                    
                    // Use overflow file for left to right direction ...
                    if (!left_done) {
                        // std::cout << "left is not done!" << std::endl;
                        assert(of_pos != 0);

                        size_t num_left_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_left_pairs < 256);
                        of_pos += 1;
                        // std::cout << "num_pairs_left = " << num_left_pairs << "\n";

                        // Go through all the L2R pairs
                        for (size_t i = 0; i < num_left_pairs; i++){
                            size_t doc_id = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            size_t lcp_val = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            // std::cout << "doc = " << doc_id << ", lcp = " << lcp_val << std::endl;

                            if (lcp_val >= length) {
                                left_docs.push_back(doc_id);
                                //leftmost_doc = doc_id;
                                //leftmost_found = true;
                                //break;
                            }
                        }
                    }
                    // Use overflow file for right to left direction ...
                    if (!right_done) {
                        // std::cout << "right is not done!" << std::endl;
                        of_pos = profile[this->num_cols * 4];
                        assert(of_pos != 0);

                        // Move past the left pairs to right pairs
                        size_t num_left_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_left_pairs < 256);
                        of_pos += 1 + (DOCWIDTH * 2 * num_left_pairs);

                        size_t num_right_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_right_pairs < 256);
                        of_pos += 1;
                        // std::cout << "num_pairs_right = " << num_right_pairs << "\n";

                        // Go through all the R2L pairs
                        for (size_t i = 0; i < num_right_pairs; i++){
                            size_t doc_id = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            size_t lcp_val = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            // std::cout << "doc = " << doc_id << ", lcp = " << lcp_val << std::endl;
                            
                            if (lcp_val >= length) {
                                right_docs.push_back(doc_id);
                                //rightmost_doc = doc_id;
                                //rightmost_found = true;
                                //break;
                            }
                        }
                    }
                    // Make sure we have found both documents in 
                    // both directions, and they end with same doc
                    assert((left_docs.size() && right_docs.size()));
                    assert ((left_docs.back() == right_docs.back()));

                    // Remove the last document from one of the vectors
                    left_docs.pop_back();

                    // std::cout << "left docs: ";
                    // for (auto x: left_docs)
                    //     std::cout << x << " ";
                    // std::cout << "\n";
                    // std::cout << "right docs: ";
                    // for (auto x: right_docs)
                    //     std::cout << x << " ";
                    // std::cout << "\n";

                    // Output the leftmost/rightmost nodes along with
                    // nodes underneath with hit
                    std::string output_str = "{";
                    for (auto x: left_docs)
                        output_str += std::to_string(x) + ",";
                    for (int i = right_docs.size()-1; i >= 0; i--)
                        output_str += std::to_string(right_docs[i]) + ",";

                    // Remove comma and close bracket
                    output_str.pop_back();
                    output_str += "} ";
                    listings_fd << output_str;
            };
        
            // lambda to output the document listing for current exact match 
            auto generate_listing = [this, &process_profile_with_subtree, &listings_fd,
                                         &curr_prof_ch, &curr_prof_pos, &num_LF_steps,
                                         &i, &end_pos_of_match, &use_end, &length, &pointer_set, &idx] 
                                        () {
                // only grab pointer if it has been set already
                if (pointer_set) {

                    // grab the correct current profile
                    std::vector<uint64_t> curr_profile (this->num_cols * 4, 0);
                    if (use_end)
                        curr_profile = end_doc_profiles[curr_prof_ch][curr_prof_pos];
                    else
                        curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];

                    // update profile based on LF steps
                    idx = 0;
                    std::for_each(curr_profile.begin(), 
                                curr_profile.end(), 
                                [&](uint64_t &x){
                                if (idx % 2 == 1) 
                                {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                                idx++;});

                    listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";

                    length = std::min((size_t) MAXLCPVALUE, (end_pos_of_match-i));
                    process_profile_with_subtree(curr_profile, length, use_end);
                    end_pos_of_match = i;
                }
            };

            // lambda to output empty listing for current position
            auto generate_empty_listing = [this, &i, &end_pos_of_match, &listings_fd] {
                listings_fd << "[" << i << "," << end_pos_of_match << "] {} ";
                end_pos_of_match = i - 1;
            };

            // lambda to reset bwt variables to full range
            auto initialize_to_full_range = [this, &start, &end, &start_run, &end_run,
                                             &num_ch_before_start, &num_ch_before_end, &next_ch] () {
                start = 0; end = this->bwt.size();
                start_run = 0; end_run = this->r - 1;
                num_ch_before_start = 0;
                num_ch_before_end = this->bwt.number_of_letter(next_ch);
            };

            // lambda to reset bwt variables to full range after empty listing
            auto initialize_to_full_range_after_empty = [this, &start, &end, &start_run, &end_run,
                                             &num_ch_before_start, &num_ch_before_end, &next_ch] () {
                start = 0; end = this->bwt.size();
                start_run = 0; end_run = this->r - 1;
                num_ch_before_start = 0;
                num_ch_before_end = 0;
            };

            // lambda to choose a new document array profile to use
            auto update_profile_pointer = [this, &curr_prof_ch, &next_ch, &curr_prof_pos,
                                           &start_run, &num_LF_steps, &use_start, &use_end, 
                                           &start, &pointer_set] () {
                
                // turn on variable to show we have found set it
                pointer_set = true;

                // grab any profile at run boundary of next_ch (choose first one)
                curr_prof_ch = next_ch;
                curr_prof_pos = this->bwt.run_head_rank(start_run, next_ch);
                num_LF_steps = 0;
                use_start = false; use_end = false;

                // if the start position run is the same as query
                // ch, then we can guarantee that end of run is in the range
                // otherwise, we can guarantee the start of run is in range.
                if (this->bwt[start] == next_ch)
                    use_end = true;
                else
                    use_start = true;
            };

            // build minimizer digest object (only used if needed)
            MinimizerDigest digester; 
            digester.set_lexorder(false);
            digester.set_windows(small_window_l, large_window_l);
            if (database_type == MINIMIZER) digester.set_minimizer_alp(true);

            // process each read, and print out the document lists
            while (kseq_read(seq)>=0) {
                // uppercase every letter in read
                for (int j = 0; j < seq->seq.l; j++) {seq->seq.s[j] = up_tab[(int) seq->seq.s[j]];}

                // digest read if needed
                std::string input_read = seq->seq.s;
                read_length = 0;

                if (database_type == DNA_MINIMIZER) {input_read = digester.compute_digest(input_read);}
                else if (database_type == MINIMIZER) {input_read = digester.compute_digest(input_read);}
                read_length = input_read.size();

                // include read name in output file
                listings_fd << ">" << seq->name.s << "\n";

                // re-initialize variables to use for backward search
                start = 0; end = this->bwt.size(); 
                end_pos_of_match = read_length-1; length = 0;

                curr_prof_ch = 0; curr_prof_pos = 0, num_LF_steps = 0;
                use_start = false, use_end = false; pointer_set = false;

                // initialize position variable 
                i = (read_length-1);
                // perform backward search and report document listings
                while (i >= 0) {
                    // already upper-cased
                    next_ch = input_read[i];

                    // special case: next_ch does not occur in text
                    if (this->bwt.number_of_letter(next_ch) == 0) {
                        generate_listing(); // listing for [i+1,end] ...
                        generate_empty_listing(); // empty listing for [i,i] ...

                        i--; 
                        initialize_to_full_range_after_empty();
                        pointer_set = false;
                        continue;
                    }

                    // identify the run # that start and end are in
                    num_ch_before_start = this->bwt.rank(start, next_ch); 
                    num_ch_before_end = this->bwt.rank(end, next_ch);
                    start_run = 0; end_run = 0;

                    // calculate key variables
                    size_t range_length = end - start;
                    size_t num_next_ch_in_range = num_ch_before_end - num_ch_before_start;

                    // range is within BWT run, and it is the next character
                    if (range_length == num_next_ch_in_range) {
                        num_LF_steps++;
                    } 
                    // range does not contain the next character 
                    else if (num_next_ch_in_range == 0) {
                        generate_listing();
                        initialize_to_full_range(); 
                        update_profile_pointer();
                    }
                    // range spans a run boundary, and contains the next character
                    else {
                        // update start_run variable here since it is used by following lambda
                        start_run = this->bwt.run_of_position(start);
                        update_profile_pointer();
                    }

                    // Perform an LF step
                    start = num_ch_before_start + this->F[next_ch]; 
                    end = num_ch_before_end + this->F[next_ch];

                    // move to next character
                    i--;
                }
                generate_listing(); listings_fd << "\n";
            }
            listings_fd.close();
            DONE_LOG((std::chrono::system_clock::now() - start_time));
        }

        void query_profiles_with_ftab(std::string pattern_file, 
                                      std::string output_prefix, 
                                      ref_type database_type,
                                      size_t small_window_l,
                                      size_t large_window_l) {
            /* query reads using the ftab to accelerate queries */

            STATUS_LOG("query_main", "querying the patterns");
            auto start_time = std::chrono::system_clock::now();

            // open output/input files
            std::ofstream listings_fd (output_prefix + ".listings");
            gzFile fp; kseq_t* seq;
            fp = gzopen(pattern_file.data(), "r"); 

            if(fp == 0) {std::exit(1);}
            seq = kseq_init(fp);

            // initialize variables for tracking length
            size_t length_sum = 0, num_matches = 0, num_reads = 0;

            // initialize ftab-related variables
            size_t curr_ftab_length = FTAB_ENTRY_LENGTH;
            if (database_type == MINIMIZER) curr_ftab_length = FTAB_ENTRY_LENGTH_MIN;

            char curr_ftab_char[curr_ftab_length+1] = "";
            std::string curr_ftab_str = "";
            size_t curr_ftab_num = 0;
            int ftab_pos = 0;

            // initialize variables for backward search
            size_t start = 0, end = this->bwt.size(), end_pos_of_match = 0;
            size_t start_run = 0, end_run = 0;
            size_t num_ch_before_start = 0, num_ch_before_end = 0;

            std::vector<uint64_t> curr_profile (this->num_cols * 4, 0);
            uint16_t length = 0; uint8_t next_ch = 0;
            int i = 0; size_t read_length = 0;
            bool pointer_set = false;

            // variable used keep track of position in lamdbas
            size_t idx = 0;

            // initialize variables for identifying correct document array profile
            uint8_t curr_prof_ch = 0;
            size_t curr_prof_pos = 0, num_LF_steps = 0;

            // initialize variables for which profile to use
            bool use_start = false, use_end = false;
        
            // alternate lambda to print out the document listing, that
            // will be tested to see if it improves classification
            auto process_profile_with_subtree = [&](std::vector<uint64_t> profile, uint16_t length, bool use_end) {
                    bool left_done = false, right_done = false;
                    std::vector<uint64_t> left_docs;
                    std::vector<uint64_t> right_docs;
                    size_t tuple_num = 0;

                    // std::cout << "length: " << length << std::endl;
                    // std::cout << "profile: ";
                    // for (auto x: profile)
                    //     std::cout << x << " ";
                    // std::cout << "\n";

                    // Iterates through the main table ...
                    while ((!left_done || !right_done) && (tuple_num < (this->num_cols * 4))) {   

                        // Checks if the LCP value has reached the max, and
                        // the monotonic increases have ended
                        if (profile[tuple_num] == MAXLCPVALUE)
                            left_done = true;
                        if (profile[tuple_num+2] == MAXLCPVALUE)
                            right_done = true;

                        // Check if increases in either direction are still going, and 
                        // lcp is greater than query length
                        if (!left_done && profile[tuple_num+1] >= length) {
                            left_docs.push_back(profile[tuple_num]);
                        }
                        if (!right_done && profile[tuple_num+3] >= length) {
                            right_docs.push_back(profile[tuple_num+2]);
                        }
                        tuple_num += 4;
                    }

                    // Double-check that if we haven't found the leftmost or rightmost
                    // document that there are still some documents to look at in overflow file
                    if ((left_docs.size() == 0 && left_done) || (right_docs.size() == 0 && right_done))
                        FATAL_ERROR("Issue occurred when querying the taxonomic document array.");

                    // Grab the overflow pointer based on which type of profiles was used
                    char* of_ptr = nullptr;
                    if (use_end)
                        of_ptr = mmap_edap_of;
                    else
                        of_ptr = mmap_sdap_of;
                    size_t of_pos = profile[this->num_cols * 4];

                    // Make sure that both directions are done based on overflow ptr
                    if (of_pos == 0) {
                        left_done = true;
                        right_done = true;
                        assert((left_done && right_done));
                    }
                    // std::cout << leftmost_found << " " << rightmost_found << std::endl;
                    // std::cout << of_pos << std::endl;
                    
                    // Use overflow file for left to right direction ...
                    if (!left_done) {
                        // std::cout << "left is not done!" << std::endl;
                        assert(of_pos != 0);

                        size_t num_left_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_left_pairs < 256);
                        of_pos += 1;
                        // std::cout << "num_pairs_left = " << num_left_pairs << "\n";

                        // Go through all the L2R pairs
                        for (size_t i = 0; i < num_left_pairs; i++){
                            size_t doc_id = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            size_t lcp_val = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            // std::cout << "doc = " << doc_id << ", lcp = " << lcp_val << std::endl;

                            if (lcp_val >= length) {
                                left_docs.push_back(doc_id);
                                //leftmost_doc = doc_id;
                                //leftmost_found = true;
                                //break;
                            }
                        }
                    }
                    // Use overflow file for right to left direction ...
                    if (!right_done) {
                        // std::cout << "right is not done!" << std::endl;
                        of_pos = profile[this->num_cols * 4];
                        assert(of_pos != 0);

                        // Move past the left pairs to right pairs
                        size_t num_left_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_left_pairs < 256);
                        of_pos += 1 + (DOCWIDTH * 2 * num_left_pairs);

                        size_t num_right_pairs = READ_NUM_PAIRS(of_ptr, of_pos);
                        assert(num_right_pairs < 256);
                        of_pos += 1;
                        // std::cout << "num_pairs_right = " << num_right_pairs << "\n";

                        // Go through all the R2L pairs
                        for (size_t i = 0; i < num_right_pairs; i++){
                            size_t doc_id = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            size_t lcp_val = READ_DOC_ID_OR_LCP_VAL(of_ptr, of_pos); of_pos += 2;
                            // std::cout << "doc = " << doc_id << ", lcp = " << lcp_val << std::endl;
                            
                            if (lcp_val >= length) {
                                right_docs.push_back(doc_id);
                                //rightmost_doc = doc_id;
                                //rightmost_found = true;
                                //break;
                            }
                        }
                    }
                    // Make sure we have found both documents in 
                    // both directions, and they end with same doc
                    assert((left_docs.size() && right_docs.size()));
                    assert ((left_docs.back() == right_docs.back()));

                    // Remove the last document from one of the vectors
                    left_docs.pop_back();

                    // std::cout << "left docs: ";
                    // for (auto x: left_docs)
                    //     std::cout << x << " ";
                    // std::cout << "\n";
                    // std::cout << "right docs: ";
                    // for (auto x: right_docs)
                    //     std::cout << x << " ";
                    // std::cout << "\n";

                    // Output the leftmost/rightmost nodes along with
                    // nodes underneath with hit
                    std::string output_str = "{";
                    for (auto x: left_docs)
                        output_str += std::to_string(x) + ",";
                    for (int i = right_docs.size()-1; i >= 0; i--)
                        output_str += std::to_string(right_docs[i]) + ",";

                    // Remove comma and close bracket
                    output_str.pop_back();
                    output_str += "} ";
                    listings_fd << output_str;
            };
        
            // lambda to output the document listing for current exact match 
            auto generate_listing = [this, &process_profile_with_subtree, &listings_fd,
                                         &curr_prof_ch, &curr_prof_pos, &num_LF_steps,
                                         &i, &end_pos_of_match, &use_end, &length, &pointer_set, &idx,
                                         &length_sum, &num_matches] 
                                        () {
                // only grab pointer if it has been set already
                if (pointer_set) {

                    // grab the correct current profile
                    std::vector<uint64_t> curr_profile (this->num_cols * 4, 0);
                    if (use_end)
                        curr_profile = end_doc_profiles[curr_prof_ch][curr_prof_pos];
                    else
                        curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];

                    // update profile based on LF steps
                    idx = 0;
                    std::for_each(curr_profile.begin(), 
                                curr_profile.end(), 
                                [&](uint64_t &x){
                                if (idx % 2 == 1) 
                                {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                                idx++;});

                    listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";

                    length = std::min((size_t) MAXLCPVALUE, (end_pos_of_match-i));
                    process_profile_with_subtree(curr_profile, length, use_end);
                    end_pos_of_match = i;

                    length_sum += length;
                    num_matches++;
                }
            };

            // lambda to output empty listing for current position
            auto generate_empty_listing = [this, &i, &end_pos_of_match, &listings_fd] {
                listings_fd << "[" << i << "," << end_pos_of_match << "] {} ";
                end_pos_of_match = i - 1;
            };

            // lambda to reset bwt variables to full range
            auto initialize_to_full_range = [this, &start, &end, &start_run, &end_run,
                                             &num_ch_before_start, &num_ch_before_end, &next_ch] () {
                start = 0; end = this->bwt.size();
                start_run = 0; end_run = this->r - 1;
                num_ch_before_start = 0;
                num_ch_before_end = this->bwt.number_of_letter(next_ch);
            };

            // lambda to reset bwt variables to full range after empty listing
            auto initialize_to_full_range_after_empty = [this, &start, &end, &start_run, &end_run,
                                             &num_ch_before_start, &num_ch_before_end, &next_ch] () {
                start = 0; end = this->bwt.size();
                start_run = 0; end_run = this->r - 1;
                num_ch_before_start = 0;
                num_ch_before_end = 0;
            };

            // lambda to choose a new document array profile to use
            auto update_profile_pointer = [this, &curr_prof_ch, &next_ch, &curr_prof_pos,
                                           &start_run, &num_LF_steps, &use_start, &use_end, 
                                           &start, &pointer_set] () {
                
                // turn on variable to show we have found set it
                pointer_set = true;

                // grab any profile at run boundary of next_ch (choose first one)
                curr_prof_ch = next_ch;
                curr_prof_pos = this->bwt.run_head_rank(start_run, next_ch);
                num_LF_steps = 0;
                use_start = false; use_end = false;

                // if the start position run is the same as query
                // ch, then we can guarantee that end of run is in the range
                // otherwise, we can guarantee the start of run is in range.
                if (this->bwt[start] == next_ch)
                    use_end = true;
                else
                    use_start = true;
            };

            // lambda to extract ftab entry and load variables
            auto load_ftab_entry = [this, &ftab_pos, &start, &end, &curr_prof_pos, &num_LF_steps,
                                    &curr_prof_ch, &use_start, &use_end, &i, &curr_ftab_length, &pointer_set] () {
                start = ftab[ftab_pos][0]; 
                end = ftab[ftab_pos][1];
                curr_prof_pos = ftab[ftab_pos][2]; 
                num_LF_steps = ftab[ftab_pos][3]; 
                curr_prof_ch = ftab[ftab_pos][4]; 
                use_start = ftab[ftab_pos][5];
                use_end = ftab[ftab_pos][6];
                i -= curr_ftab_length;

                // turn on variable to show we have found set it
                pointer_set = true;
            };

            // build minimizer digest object (only used if needed)
            MinimizerDigest digester; 
            digester.set_lexorder(false);
            digester.set_windows(small_window_l, large_window_l);
            if (database_type == MINIMIZER) digester.set_minimizer_alp(true);

            // process each read, and print out the document lists
            while (kseq_read(seq)>=0) {
                // uppercase every letter in read
                for (int j = 0; j < seq->seq.l; j++) {seq->seq.s[j] = up_tab[(int) seq->seq.s[j]];}
                num_reads++;

                // digest read if needed
                std::string input_read = seq->seq.s;
                read_length = 0;

                if (database_type == DNA_MINIMIZER) {input_read = digester.compute_digest(input_read);}
                else if (database_type == MINIMIZER) {input_read = digester.compute_digest(input_read);}
                read_length = input_read.size();

                // include read name in output file
                listings_fd << ">" << seq->name.s << "\n";

                // re-initialize variables to use for backward search
                start = 0; end = this->bwt.size(); 
                end_pos_of_match = read_length-1; length = 0;

                curr_prof_ch = 0; curr_prof_pos = 0, num_LF_steps = 0;
                use_start = false, use_end = false; pointer_set = false;

                // re-initialize ftab-related variables
                curr_ftab_str = ""; curr_ftab_num = 0;

                // initialize position variable 
                i = (read_length-1);

                // determine if we can use ftab to start out
                ftab_pos = check_if_ftab_can_be_used(input_read.data(), i, (database_type == MINIMIZER)); 
                if (ftab_pos >= 0) {load_ftab_entry();}

                // perform backward search and report document listings
                while (i >= 0) {
                    bool ftab_used = false; 
                    // already upper-cased
                    next_ch = input_read[i];

                    // special case: next_ch does not occur in text
                    if (this->bwt.number_of_letter(next_ch) == 0) {
                        generate_listing(); // listing for [i+1,end] ...
                        generate_empty_listing(); // empty listing for [i,i] ...

                        i--; 
                        initialize_to_full_range_after_empty();
                        pointer_set = false;
                        continue;
                    }

                    // identify the run # that start and end are in
                    num_ch_before_start = this->bwt.rank(start, next_ch); 
                    num_ch_before_end = this->bwt.rank(end, next_ch);
                    start_run = 0; end_run = 0;

                    // calculate key variables
                    size_t range_length = end - start;
                    size_t num_next_ch_in_range = num_ch_before_end - num_ch_before_start;

                    // range is within BWT run, and it is the next character
                    if (range_length == num_next_ch_in_range) {
                        num_LF_steps++;
                    } 
                    // range does not contain the next character 
                    else if (num_next_ch_in_range == 0) {
                        generate_listing();
                        initialize_to_full_range(); 

                        ftab_pos = check_if_ftab_can_be_used(input_read.data(), i, (database_type == MINIMIZER));
                        if (ftab_pos >= 0) {
                            load_ftab_entry();
                            ftab_used = true;
                        } else {
                            update_profile_pointer();
                        }
                    }
                    // range spans a run boundary, and contains the next character
                    else {
                        // update start_run variable here since it is used by following lambda
                        start_run = this->bwt.run_of_position(start);
                        update_profile_pointer();
                    }

                    if (!ftab_used) {
                        // perform an LF step
                        start = num_ch_before_start + this->F[next_ch]; 
                        end = num_ch_before_end + this->F[next_ch];

                        // move to next character
                        i--;
                    }
                }
                generate_listing(); listings_fd << "\n";
            }
            listings_fd.close();
            DONE_LOG((std::chrono::system_clock::now() - start_time));
        
                  // report overall metrics on querying
            FORCE_LOG("query_main", 
                      "processed %ld reads (%ld total chars / %ld exact matches = %.2f avg length)",
                      num_reads, length_sum, num_matches, (length_sum/(num_matches+0.0)));  
        }

        std::pair<size_t, size_t> build_ftab(bool minimizer_alphabet) {
            /* build the f_tab and store it in the *.fna.ftab file */

            // table to convert 2-bits into a character (00, 01, 10, 11)
            char nuc_tab[] = {'A', 'C', 'T', 'G'}; 

            // initialize some variables
            std::vector<size_t> index_vec(FTAB_ENTRY_LENGTH) ; 
            std::iota(index_vec.begin(), index_vec.end(), 0); 

            size_t num_entries = std::pow(FTAB_ALPHABET_SIZE, FTAB_ENTRY_LENGTH);
            size_t ftab_file_size = num_entries * FTAB_ENTRY_SIZE;
            size_t length_of_ftab_entry = FTAB_ENTRY_LENGTH;

            // change variables if using minimizer alphabet
            if (minimizer_alphabet) {
                index_vec.resize(FTAB_ENTRY_LENGTH_MIN);
                std::iota(index_vec.begin(), index_vec.end(), 0);

                num_entries = std::pow(FTAB_ALPHABET_SIZE_MIN, FTAB_ENTRY_LENGTH_MIN);
                ftab_file_size = num_entries * FTAB_ENTRY_SIZE;
                length_of_ftab_entry = FTAB_ENTRY_LENGTH_MIN;
            }

            // create a memory-mapped file for writing the ftab
            mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
            struct stat ftab_file_info;
            char* mmap_ftab;

            std::string outfile = output_ref + std::string(".ftab");
            int ftab_fd = open(outfile.c_str(), O_RDWR | O_CREAT | O_TRUNC, mode);
            ftruncate(ftab_fd, ftab_file_size);

            if (ftab_fd == -1)
                FATAL_ERROR("issue occurred when opening up ftab file.");
            if (fstat(ftab_fd, &ftab_file_info) == -1)
                FATAL_ERROR("issue occurred when checking file size on ftab file.");

            mmap_ftab = (char*) mmap(NULL,
                                     ftab_file_size,
                                     PROT_READ | PROT_WRITE, 
                                     MAP_SHARED, 
                                     ftab_fd, 0);
            ASSERT((mmap_ftab != MAP_FAILED),
                    "issue occurred when mapping file for ftab.");

            // DEBUG: open file and store results for checking
            std::ofstream debug_fd(output_ref + ".ftab.debug.txt");
                      
            // iterate through all possible k-mers in dictionary
            size_t num_not_found = 0, num_found = 0;
            for (size_t loop_index = 0; loop_index < num_entries; loop_index++) {
                std::string curr_seq(length_of_ftab_entry, '*');

                if (!minimizer_alphabet) {
                    // case DNA-seq: generates current 10-mer by converting 
                    // every 2-bits into a character
                    std::transform(index_vec.begin(), index_vec.end(), curr_seq.begin(),
                                [&] (size_t pos) {
                                        uint8_t code = FTAB_GRAB_CODE(loop_index, pos);
                                        ASSERT((code <= 3), "invalid DNA code encounted in ftab generation.");
                                        return nuc_tab[code];});
                } else {
                    // case minimizer-alphabet: generate 3-mer by converting the decimal
                    // number into a 253-base
                    size_t curr_quotient = loop_index;
                    size_t pos = 0;
                    while (curr_quotient > 0) {
                        size_t remainder = curr_quotient % FTAB_ALPHABET_SIZE_MIN;
                        curr_quotient = curr_quotient / FTAB_ALPHABET_SIZE_MIN;
                        ASSERT((remainder < 253), "invalid minimizer character encountered in ftab generation.");
                        curr_seq[pos++] = (uint8_t) (remainder+3); // avoid 0, 1, 2 due to PFP
                    }
                    for (;pos < length_of_ftab_entry; pos++) {curr_seq[pos] = (uint8_t) 3;} // smallest char is 3
                }

                // we reverse here since we grab bits from LSB to MSB and 
                // build string from left to right, so we reverse it to make
                // the two-bits correspond to letters                  
                std::reverse(curr_seq.begin(), curr_seq.end());

                // initialize variables for backward search
                size_t start = 0, end = this->bwt.size();
                bool not_found = false;

                // initialize variables that are need to grab correct profile 
                uint8_t curr_prof_ch = 0;
                size_t curr_prof_pos = 0, num_LF_steps = 0;
                bool use_start = false, use_end = false;

                // perform backward search for current 10-mer
                int j = (length_of_ftab_entry-1);
                for (; j>=0; j--) {
                    // get character, keep in mind it is already capitalized
                    uint8_t next_ch = curr_seq[j];

                    // identify the run # for start and end
                    size_t num_ch_before_start = this->bwt.rank(start, next_ch);
                    size_t num_ch_before_end = this->bwt.rank(end, next_ch);
                    size_t start_run = this->bwt.run_of_position(start);
                    size_t end_run = this->bwt.run_of_position(end-1);

                    // range spans run boundary, so there are different BWT chars
                    if (start_run != end_run) {
                        // bwt range is going to be empty, so we didn't find it ...
                        if (num_ch_before_start == num_ch_before_end) {
                            not_found = true;
                            break;
                        }

                        // grab any profile with next_ch in BWT (policy: choose first one)
                        curr_prof_ch = next_ch;
                        curr_prof_pos = this->bwt.run_head_rank(start_run, next_ch);
                        num_LF_steps = 0;
                        use_start = false; use_end = false;

                        // if current start position is a run of next_ch characters,
                        // then we want to use the profile at the of run, otherwise,
                        // use the profile at the beginning of a run
                        if (this->bwt[start] == next_ch)
                            use_end = true;
                        else
                            use_start = true;
                    }
                    // range is within a BWT run, but wrong character
                    else if(this->bwt[start] != next_ch) {
                        not_found = true;
                        break;
                    }
                    // range is within BWT run, and is correct character
                    else {
                        num_LF_steps++;
                    }

                    // perform backward step
                    start = num_ch_before_start + this->F[next_ch];
                    end = num_ch_before_end + this->F[next_ch];
                }

                // keep track of how many are not found
                if (not_found) {
                    num_not_found++;  
                    start = 1; end = 0; // signal that not found
                    curr_prof_ch = 0; curr_prof_pos = 0;
                    num_LF_steps = 0;
                } else {
                    num_found++;
                }

                // write out the first three elements of ftab entry
                size_t start_pos = loop_index * FTAB_ENTRY_SIZE;
                ASSERT((sizeof(size_t) == 8), "issue occurred when checking size_t size");

                for (int k = 0; k < sizeof(size_t); k++) {
                    /* 1. BWT range start position (little-endian) */
                    mmap_ftab[start_pos + k] = ((0xFF << (8*k)) & start) >> (8*k);
                    /* 2. BWT range end position (little-endian) */
                    mmap_ftab[start_pos + sizeof(size_t) + k] = ((0xFF << (8*k)) & end) >> (8*k);
                    /* 3. Current profile position (little-endian) */
                    mmap_ftab[start_pos + (2*sizeof(size_t)) + k] = ((0xFF << (8*k)) & curr_prof_pos) >> (8*k);
                }

                // write the remaining two elements of ftab entry
                ASSERT((num_LF_steps < length_of_ftab_entry), "Issue with ftab build. (1)");

                /* 4. Number of LF Steps - USE_START (top bit) - USE_END (2nd top bit) */
                ASSERT((use_start || use_end || not_found), "Issue with ftab build. (2)");
                mmap_ftab[start_pos + (3*sizeof(size_t))] = (0xFF & num_LF_steps) | (use_start << 7) | (use_end << 6);

                /* 5. BWT character */
                mmap_ftab[start_pos + (3*sizeof(size_t)) + 1] = (0xFF & curr_prof_ch);

                // create readable string for debug file
                if (minimizer_alphabet) {
                    std::string temp_str = "";
                    uint8_t curr_ch = 0;
                    for (size_t i = 0; i < curr_seq.size()-1; i++) {
                        curr_ch = curr_seq[i];
                        temp_str += std::to_string(curr_ch) + "-";
                    }
                    curr_ch = curr_seq[curr_seq.size()-1];
                    temp_str += std::to_string(curr_ch);
                    curr_seq = temp_str;
                }

                debug_fd << loop_index
                         << "," << curr_seq 
                         << "," << start 
                         << "," << end 
                         << "," << curr_prof_pos 
                         << "," << num_LF_steps 
                         << "," << unsigned(curr_prof_ch)
                         << "," << use_start 
                         << "," << use_end << std::endl;
            }

            // unmap the ftab and debug files
            munmap(mmap_ftab, ftab_file_size);
            close(ftab_fd); 
            debug_fd.close();

            return std::make_pair(num_found, num_not_found);
        }

        void load_ftab_from_file(bool minimizer_alp) {
            /* loads ftab from a file into a vector structure that can used during query */
            check_ftab_file(minimizer_alp);

            STATUS_LOG("query_main", "loading the ftab structure");
            auto start = std::chrono::system_clock::now();
            read_ftab_file(minimizer_alp);

            DONE_LOG((std::chrono::system_clock::now() - start));
            FORCE_LOG("query_main", "number of entries in ftab: %ld" , ftab.size());
            std::cerr << "\n";
        }

    private:
        /*********************************/ 
        /* Private instance variables
        /*********************************/
        bool rle = true;
        size_t num_cols = 0;
        std::string output_ref = "";

        // This represents the ftab, it is represented as a vector
        // of vectors, indexed by the k-mer and the vectors stores
        // start, end, profile pos, LF steps, and BWT char.
        std::vector<std::vector<size_t>> ftab;

        /*********************************/ 
        /* Private instance methods
        /*********************************/

        void load_bwt_structure(std::string bwt_fname) {
            /* loads the BWT data-structure from the build sub-command */
            if(rle)
            {
                std::string bwt_heads_fname = bwt_fname + ".heads";
                std::ifstream ifs_heads(bwt_heads_fname);

                std::string bwt_len_fname = bwt_fname + ".len";
                std::ifstream ifs_len(bwt_len_fname);

                this->bwt = rle_string_t(ifs_heads, ifs_len);

                ifs_heads.seekg(0);
                ifs_len.seekg(0);
                this->build_F_(ifs_heads, ifs_len);

                ifs_heads.close();
                ifs_len.close();
            }
            else
            {
                std::ifstream ifs(bwt_fname);
                this->bwt = rle_string_t(ifs);

                ifs.seekg(0);
                this->build_F(ifs);

                ifs.close();
            }
        }

        void check_doc_array_files(std::string fname) {
            /* Examines the file size and make sure it is the correct size */
            struct stat filestat;
            FILE *fd;

            if ((fd = fopen(fname.c_str(), "r")) == nullptr)
                error("open() file " + fname + " failed");
            if (fstat(fileno(fd), &filestat) < 0)
                error("stat() file " + fname + " failed");

            // calculate the size of each row of table and check file is correct size
            size_t record_size = 1 + (8 * this->num_cols) + 8;
            if (filestat.st_size != (record_size * this->r))
                FATAL_ERROR("the file size for the document array is not valid.");
            fclose(fd);
        }

        void read_doc_profiles_main_table(std::vector<std::vector<std::vector<uint64_t>>>& prof_matrix, 
                                          std::string input_file,
                                          std::ofstream* out_fd = nullptr) {
            /* reads the main document array file that is taxonomic compressed */
            FILE *fd;
            if ((fd = fopen(input_file.c_str(), "r")) == nullptr)
                FATAL_ERROR("open() file failed");
            
            // Go through each record, it starts with the BWT ch and then
            // the left and right monotonic increase pairs and lastly
            // followed by index to the overflow file.
            size_t curr_val = 0;
            uint8_t curr_bwt_ch = 0;
            for (size_t i = 0; i < this->r; i++){
                // Step 1: reading bwt character
                if (fread(&curr_bwt_ch, 1, 1, fd) != 1)
                    FATAL_ERROR("issue occurred while reading in bwt character from doc profiles file.");

                // Step 2: reading all the left and right increase pairs
                std::vector<uint64_t> curr_record ((num_cols * 4) + 1, 0);
                for (size_t j = 0; j < (num_cols*4); j++) {
                    if ((fread(&curr_val, DOCWIDTH, 1, fd)) != 1)
                        FATAL_ERROR("fread() failed"); 
                    curr_record[j] = curr_val;
                }

                // Step 3: read the overflow pointer value
                size_t of_ptr = 0;
                if (fread(&of_ptr, 8, 1, fd) != 1)
                    FATAL_ERROR("fread() failed");
                curr_record[(num_cols * 4)] = of_ptr;

                // If we are print out the profiles ...
                if (print_profiles && i < profiles_to_print) {
                    *out_fd << curr_bwt_ch << ",";
                    for (auto x: curr_record)
                        *out_fd << x << ",";
                    *out_fd << of_ptr << "\n";
                }
                prof_matrix[curr_bwt_ch].push_back(curr_record);
            }
        }

        void check_ftab_file(bool minimizer_alp) {
            /* makes sure the ftab file is the correct length */
            struct stat filestat;
            FILE *fd;

            std::string file_path = output_ref + std::string(".ftab");
            if ((fd = fopen(file_path.c_str(), "r")) == nullptr)
                FATAL_ERROR(("open() file " + file_path + " failed").data());
            if (fstat(fileno(fd), &filestat) < 0)
                FATAL_ERROR(("stat() file " + file_path + " failed").data());
            fclose(fd);

            // compute the expected size of the ftab file
            size_t expected_size = FTAB_ENTRY_SIZE * std::pow(FTAB_ALPHABET_SIZE, FTAB_ENTRY_LENGTH);
            if (minimizer_alp)
                expected_size = FTAB_ENTRY_SIZE * std::pow(FTAB_ALPHABET_SIZE_MIN, FTAB_ENTRY_LENGTH_MIN);

            if (filestat.st_size != (expected_size))
                FATAL_ERROR("invalid file size for *.ftab files");  
        }

        void read_ftab_file(bool minimizer_alp) {
            /* load ftab into vector for usage */

            // open the ftab file
            std::string ftab_file = output_ref + std::string(".ftab");

            FILE *fd;
            if ((fd = fopen(ftab_file.c_str(), "r")) == nullptr)
                error("open() file " + ftab_file + " failed");

            size_t start = 0, end = 0, pos = 0;
            uint8_t num_LF_steps = 0, bwt_ch = 0;
            bool use_start = false, use_end = false;

            // go through reach entry in the ftab, and add to vector
            size_t num_entries = std::pow(FTAB_ALPHABET_SIZE, FTAB_ENTRY_LENGTH);
            if (minimizer_alp) {num_entries = std::pow(FTAB_ALPHABET_SIZE_MIN, FTAB_ENTRY_LENGTH_MIN);}

            for (int i = 0; i < num_entries; i++) {
                // read the current record
                if (fread(&start, sizeof(size_t), 1, fd) != 1)
                    FATAL_ERROR("issue occurred when reading from ftab file. (1)");
                if (fread(&end, sizeof(size_t), 1, fd) != 1)
                    FATAL_ERROR("issue occurred when reading from ftab file. (2)");
                if (fread(&pos, sizeof(size_t), 1, fd) != 1)
                    FATAL_ERROR("issue occurred when reading from ftab file. (3)");
                if (fread(&num_LF_steps, 1, 1, fd) != 1)
                    FATAL_ERROR("issue occurred when reading from ftab file. (4)");
                if (fread(&bwt_ch, 1, 1, fd) != 1)
                    FATAL_ERROR("issue occurred when reading from ftab file. (5)");
                
                use_start = num_LF_steps & (1 << 7);
                use_end = num_LF_steps & (1 << 6);
                num_LF_steps = num_LF_steps & 0x3F; // ignore top two bits

                ASSERT(((unsigned(num_LF_steps) < FTAB_ENTRY_LENGTH)), 
                        "issue occurred when loading ftab. (6)");
                ASSERT((use_start || use_end || (start > end)),
                        "issue occurred when loading ftab. (7)");

                // add current entry to table
                std::vector<size_t> curr_entry {start, end, pos, num_LF_steps, bwt_ch, use_start, use_end};
                ftab.push_back(curr_entry);
            }
            fclose(fd);
        }

        int check_if_ftab_can_be_used(const char* seq, int i, bool minimizer_alp) {
            /* extracts the current query from read, and checks table */

            // initialize variable for ftab query string
            size_t ftab_entry_length = FTAB_ENTRY_LENGTH;
            if (minimizer_alp) {ftab_entry_length = FTAB_ENTRY_LENGTH_MIN;}

            char curr_ftab_char[ftab_entry_length+1] = "";
            std::string curr_ftab_str = "";
            size_t curr_ftab_num = 0;

            // check if there is room left to take a k-mer
            if (i >= ftab_entry_length-1) {
                std::memcpy(curr_ftab_char, 
                            (char *) &seq[i-ftab_entry_length+1], 
                            ftab_entry_length);
                curr_ftab_char[ftab_entry_length] = '\0'; 
                curr_ftab_str = curr_ftab_char;

                // generate the pointer for the current k-mer
                curr_ftab_num = 0;
                if (!minimizer_alp) {
                    for (int j = 0; j < ftab_entry_length; j++){
                        if (curr_ftab_str[j] == 'A')
                            curr_ftab_num = curr_ftab_num | (0x0 << (2*(ftab_entry_length-j-1)));
                        else if (curr_ftab_str[j] == 'C')
                            curr_ftab_num = curr_ftab_num | (0x1 << (2*(ftab_entry_length-j-1)));
                        else if (curr_ftab_str[j] == 'T')
                            curr_ftab_num = curr_ftab_num | (0x2 << (2*(ftab_entry_length-j-1)));
                        else if (curr_ftab_str[j] == 'G')
                            curr_ftab_num = curr_ftab_num | (0x3 << (2*(ftab_entry_length-j-1)));
                        else
                            return -1;
                    } 
                } else {
                    // subtract 3 from each char and the convert back to 
                    // decimal from 253-base
                    for (size_t j = 0; j < curr_ftab_str.size(); j++) {
                        curr_ftab_str[j] -= 3;
                        uint8_t curr_ch = unsigned(curr_ftab_str[j]);
                        ASSERT((curr_ch >= 0 && curr_ch <= 252), 
                                "unexpected string found in check_if_ftab().");

                        curr_ftab_num += std::pow(253, ftab_entry_length-1-j) * unsigned(curr_ch);
                    }
                }

                // use ftab entry to see if it is present in text
                size_t start = ftab[curr_ftab_num][0];
                size_t end = ftab[curr_ftab_num][1];

                if (start > end) 
                    return -1;
                else 
                    return curr_ftab_num;
            } else {
                return -1;
            }
        }

        vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths){
            heads.clear();
            heads.seekg(0);
            lengths.clear();
            lengths.seekg(0);

            this->F = vector<ulint>(256, 0);
            int c;
            ulint i = 0;
            while ((c = heads.get()) != EOF)
            {
                size_t length;
                lengths.read((char *)&length, 5);
                if (c > TERMINATOR)
                    this->F[c]+=length;
                else
                {
                    this->F[TERMINATOR]+=length;
                    this->terminator_position = i;
                }
                i++;
            }
            for (ulint i = 255; i > 0; --i)
                this->F[i] = this->F[i - 1];
            this->F[0] = 0;
            for (ulint i = 1; i < 256; ++i)
                this->F[i] += this->F[i - 1];
            return this->F;
        }
};

#endif /* end of include guard: _TAX_DOC_QUERIES_H */