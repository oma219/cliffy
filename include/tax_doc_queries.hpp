/*
 * File: tax_doc_queries.hpp
 * Description: ....
 * Date: Aug. 30th, 2023
 */

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
        size_t num_cols = 0;
        bool rle = true;
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
                        rle (rle),
                        num_cols (num_cols),
                        ri::r_index<sparse_bv_type, rle_string_t>(),
                        start_doc_profiles(256, std::vector<std::vector<uint64_t>>(0, std::vector<uint64_t>(0))),
                        end_doc_profiles(256, std::vector<std::vector<uint64_t>>(0, std::vector<uint64_t>(0)))
        {
            // load the BWT 
            STATUS_LOG("build_profiles", "loading the bwt of the input text");
            auto start = std::chrono::system_clock::now();

            std::string bwt_fname = filename + ".bwt";
            load_bwt_structure(bwt_fname);
            DONE_LOG((std::chrono::system_clock::now() - start));

            // gather some statistics on the BWT
            this->r = this->bwt.number_of_runs();
            ri::ulint n = this->bwt.size();
            size_t log_r = bitsize(uint64_t(this->r));
            size_t log_n = bitsize(uint64_t(this->bwt.size()));
            FORCE_LOG("build_profiles", "bwt statistics: n = %ld, r = %ld\n" , this->bwt.size(), this->r);

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

        void query_profiles(std::string pattern_file){
            /* Go through and query the reads and print the leftmost, rightmost occurrence of string */

            // Open output/input files
            std::ofstream listings_fd (pattern_file + ".listings");
            gzFile fp; kseq_t* seq;

            fp = gzopen(pattern_file.data(), "r"); 
            if(fp == 0) {std::exit(1);}
            seq = kseq_init(fp);

            // Variable used keep track of position in lamdbas
            size_t idx = 0;

            // Alternate lambda to print out the document listing, that
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

            // Process each read, and print out the document lists
            while (kseq_read(seq)>=0) {

                // Uppercase every character in read
                for (size_t i = 0; i < seq->seq.l; ++i) 
                {
                    seq->seq.s[i] = static_cast<char>(std::toupper(seq->seq.s[i]));
                }

                size_t start = 0, end = this->bwt.size(), end_pos_of_match = seq->seq.l-1;
                std::vector<uint64_t> curr_profile (this->num_cols * 4, 0);
                uint16_t length = 0;
                listings_fd << ">" << seq->name.s << "\n";

                // Pointer variables to current profile being "used"
                uint8_t curr_prof_ch = 0;
                size_t curr_prof_pos = 0, num_LF_steps = 0;

                // Tell us what type of profile to grab based on pointer variables
                bool use_start = false, use_end = false;

                // Perform backward search and report document listings when
                // range goes empty or we reach the end
                for (int i = (seq->seq.l-1); i >= 0; i--) {
                    uint8_t next_ch = seq->seq.s[i];

                    size_t num_ch_before_start = this->bwt.rank(start, next_ch); 
                    size_t num_ch_before_end = this->bwt.rank(end, next_ch);
                    size_t start_run = this->bwt.run_of_position(start);
                    size_t end_run = this->bwt.run_of_position(end-1);
                    
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

                        // If the start position run is the same as query
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

                        // grab any profile at run boundary of next_ch (choose first one)
                        curr_prof_ch = next_ch;
                        curr_prof_pos = this->bwt.run_head_rank(start_run, next_ch);
                        num_LF_steps = 0;
                        use_start = false; use_end = false;

                        // If the start position run is the same as query
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

                // std::exit(1);
            }

        }

    private:
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
                                          std::ofstream* out_fd = nullptr) 
        {
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