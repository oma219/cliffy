/*
 * File: topk_doc_queries.hpp
 * Description: ....
 * Date: Sept. 14th, 2023
 */

#ifndef _TOPK_DOC_QUERIES_H
#define _TOPK_DOC_QUERIES_H

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd>
class topk_doc_queries : ri::r_index<sparse_bv_type, rle_string_t>
{
    public:
        size_t num_cols = 0;
        bool rle = true;
        bool print_profiles = false;
        std::string output_csv_path = "";
        size_t profiles_to_print = 0;

        // Used for printing out data-structure
        std::ofstream csv_sdap_output;
        std::ofstream csv_edap_output;

        // This vectors has the following dimensions: [256][num of ith char][num_docs]
        // This structure stores the DA profiles for each character separately.
        std::vector<std::vector<std::vector<uint16_t>>> start_doc_profiles;
        std::vector<std::vector<std::vector<uint16_t>>> end_doc_profiles;

        topk_doc_queries(std::string filename,
                         size_t num_cols,
                         size_t num_profiles = 0,
                         std::string output_path = "",
                         bool rle = true):
                         profiles_to_print(num_profiles),
                         rle (rle),
                         num_cols (num_cols),
                         ri::r_index<sparse_bv_type, rle_string_t>(),
                         start_doc_profiles(256, std::vector<std::vector<uint16_t>>(0, std::vector<uint16_t>(0))),
                         end_doc_profiles(256, std::vector<std::vector<uint16_t>>(0, std::vector<uint16_t>(0)))
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
            check_doc_array_files(filename + ".topk.sdap");
            check_doc_array_files(filename + ".topk.sdap");

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

            read_doc_profiles_topk_table(start_doc_profiles, filename + ".topk.sdap", &csv_sdap_output);
            read_doc_profiles_topk_table(end_doc_profiles, filename + ".topk.edap", &csv_edap_output);

            if (print_profiles) {
                csv_sdap_output.close();
                csv_edap_output.close();
            }
            DONE_LOG((std::chrono::system_clock::now() - start));
        }

        ~topk_doc_queries()
        {
            /* Destructor - close files */
            csv_sdap_output.close(); csv_edap_output.close();
        }

        void query_profiles(std::string pattern_file){
            /* Go through and query the reads and print all occurrences among top-k */

            // Open output/input files
            std::ofstream listings_fd (pattern_file + ".listings");
            gzFile fp; kseq_t* seq;

            fp = gzopen(pattern_file.data(), "r"); 
            if(fp == 0) {std::exit(1);}
            seq = kseq_init(fp);

            // Variable used keep track of position in lamdbas
            size_t idx = 0;

            // lambda to print out the document listing
            auto process_profile = [&](std::vector<uint16_t> profile, uint16_t length) {
                std::vector<size_t> found_docs = {};

                // Go from largest lcp to smallest lcp and check for membership
                assert(profile.size() == this->num_cols);
                for (int i = profile.size()-1; i>=0; i-=2) {
                    if (profile[i] >= length) {
                        found_docs.push_back(profile[i-1]);
                    }
                }

                // Print the listings
                assert(found_docs.size() > 0);
                listings_fd << "{";
                for (size_t i = 0; i < found_docs.size(); i++) {
                    if (i != found_docs.size()-1) 
                        listings_fd << std::to_string(found_docs[i]) << ",";
                    else 
                        listings_fd << std::to_string(found_docs[i]) << "} ";
                }
            };

            // Process each read, and print out the document lists
            while (kseq_read(seq)>=0) {

                // Uppercase every character in read
                for (size_t i = 0; i < seq->seq.l; ++i) 
                {
                    seq->seq.s[i] = static_cast<char>(std::toupper(seq->seq.s[i]));
                }

                size_t start = 0, end = this->bwt.size(), end_pos_of_match = seq->seq.l-1;
                std::vector<uint16_t> curr_profile (this->num_cols * 2, 0);
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
                                          [&](uint16_t &x){
                                          if (idx % 2 == 1) 
                                          {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                                          idx++;});
                            listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";

                            length = std::min((size_t) MAXLCPVALUE, (end_pos_of_match-i));
                            process_profile(curr_profile, length);
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
                                        [&](uint16_t &x){
                                        if (idx % 2 == 1) 
                                        {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                                        idx++;});
                        listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";

                        length = std::min((size_t) MAXLCPVALUE, (end_pos_of_match-i));
                        process_profile(curr_profile, length);
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
                              [&](uint16_t &x){
                              if (idx % 2 == 1) 
                              {x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);} 
                              idx++;});

                listings_fd << "[" << 0 << "," << end_pos_of_match << "] ";
                length = std::min((size_t) MAXLCPVALUE, end_pos_of_match+1);
                process_profile(curr_profile, length);
                listings_fd << "\n";
            }
            listings_fd.close();
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

        void check_doc_array_files(std::string fname) {
            /* Examines the file size and make sure it is the correct size */
            struct stat filestat;
            FILE *fd;

            if ((fd = fopen(fname.c_str(), "r")) == nullptr)
                error("open() file " + fname + " failed");
            if (fstat(fileno(fd), &filestat) < 0)
                error("stat() file " + fname + " failed");

            // calculate the size of each row of table and check file is correct size
            size_t record_size = (4 * this->num_cols) + 1;
            if (filestat.st_size != (record_size * this->r))
                FATAL_ERROR("the file size for the document array is not valid.");
            fclose(fd);
        }

        void read_doc_profiles_topk_table(std::vector<std::vector<std::vector<uint16_t>>>& prof_matrix, 
                                          std::string input_file,
                                          std::ofstream* out_fd = nullptr) 
        {
            /* reads the main document array file that is topk */
            FILE *fd;
            if ((fd = fopen(input_file.c_str(), "r")) == nullptr)
                FATAL_ERROR("open() file failed");
            
            // Go through each record, it starts with the BWT ch and then
            // the a # of pairs of document number and lcp value. Representing
            // the top-k lcps for that suffix.
            size_t curr_val = 0;
            uint8_t curr_bwt_ch = 0;
            for (size_t i = 0; i < this->r; i++){
                // Step 1: reading bwt character
                if (fread(&curr_bwt_ch, 1, 1, fd) != 1)
                    FATAL_ERROR("issue occurred while reading in bwt character from doc profiles file.");

                // Step 2: reading all the left and right increase pairs
                std::vector<uint16_t> curr_record ((this->num_cols * 2), 0);
                for (size_t j = 0; j < (this->num_cols*2); j++) {
                    if ((fread(&curr_val, DOCWIDTH, 1, fd)) != 1)
                        FATAL_ERROR("fread() failed"); 
                    curr_record[j] = curr_val;
                }

                // If we are printing out the profiles ...
                if (this->print_profiles && i < this->profiles_to_print) {
                    *out_fd << curr_bwt_ch << ",";
                    for (auto x: curr_record)
                        *out_fd << x << ",";
                    *out_fd << "\n";
                }
                prof_matrix[curr_bwt_ch].push_back(curr_record);
            }
        }

};

#endif /* end of include guard: _TOPK_DOC_QUERIES_H */