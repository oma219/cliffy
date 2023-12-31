/*
 * File: doc_queries.hpp
 * Description: Definition of doc_queries objects that performs queries
 *              on the document array profiles data-structure. It is based
 *              ms_pointers.hpp file that was written by Massimiliano Rossi.
 * Date: Sept. 10th, 2022
 */

#ifndef _DOC_QUERIES_H
#define _DOC_QUERIES_H

#include <kseq.h>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <r_index.hpp>
#include <ms_rle_string.hpp>
#include <pfp_doc.hpp>
#include <string>

char up_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
	'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',	91,	 92,  93,  94,	95,
	 96, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
	'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 123, 124, 125, 126, 127
};

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd>
class doc_queries : ri::r_index<sparse_bv_type, rle_string_t>
{
    public:
        doc_queries(std::string filename,
                    std::string output_path="", 
                    size_t num_profiles=0, 
                    bool rle = true): 
                    ri::r_index<sparse_bv_type, rle_string_t>(),
                    rle(rle),
                    start_doc_profiles(256, std::vector<std::vector<uint16_t>>(0, std::vector<uint16_t>(0))),
                    end_doc_profiles(256, std::vector<std::vector<uint16_t>>(0, std::vector<uint16_t>(0)))
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

            // check file sizes and make sure it is valid
            check_doc_array_files(filename + ".sdap");
            check_doc_array_files(filename + ".edap");
            FORCE_LOG("query_main", "number of documents: d = %ld" , num_docs);
           
            // load the profiles for starts and ends
            STATUS_LOG("query_main", "loading the document array profiles");
            start = std::chrono::system_clock::now();
            
            read_doc_profiles(start_doc_profiles, filename + ".sdap", 
                              this->num_docs, this->r, 
                              output_path + ".sdap.csv", num_profiles);
            read_doc_profiles(end_doc_profiles, filename + ".edap", 
                              this->num_docs, this->r, 
                              output_path + ".edap.csv", num_profiles);
            DONE_LOG((std::chrono::system_clock::now() - start));
        }

        void query_profiles(std::string pattern_file, std::string output_prefix){
            /* Takes in a file of reads, and lists all the documents containing the read */

            // open output/input files
            std::ofstream listings_fd (output_prefix + ".listings");
            gzFile fp; kseq_t* seq;
            fp = gzopen(pattern_file.data(), "r"); 

            if(fp == 0) {std::exit(1);}
            seq = kseq_init(fp);

            // lambda to print out the document listing
            auto process_profile = [&](std::vector<uint16_t> profile, uint16_t length) {
                    std::string output_str = "{";
                    bool found_one = false;

                    for (size_t i = 0; i < profile.size(); i++) {
                        //std::cout << "length = " << length << "  profile_val = " << profile[i] << std::endl;
                        if (profile[i] >= length) {
                            output_str += "," + std::to_string(i);
                            found_one = true;
                        }
                    }

                    if (found_one)
                        output_str = "{" + output_str.substr(2) + "} ";
                    else {   
                        output_str = "{} ";
                    }                
                    listings_fd << output_str;
            };

            // process each read, and print out the document lists
            while (kseq_read(seq)>=0) {
                
                // initialize variables to use for backward search
                size_t start = 0, end = this->bwt.size(), end_pos_of_match = seq->seq.l-1;
                std::vector<uint16_t> curr_profile (this->num_docs, 0);
                uint16_t length = 0;

                listings_fd << ">" << seq->name.s << "\n";

                // pointer variables to current profile being "used"
                uint8_t curr_prof_ch = 0;
                size_t curr_prof_pos = 0, num_LF_steps = 0;

                // tell us what type of profile to grab based on pointer variables
                bool use_start = false, use_end = false;

                // perform backward search and report document listings when
                // range goes empty or we reach the end
                for (int i = (seq->seq.l-1); i >= 0; i--) {

                    // Uppercase the letter prior to using it
                    uint8_t next_ch = up_tab[(int) seq->seq.s[i]];

                    // Identify the run # that start and end are in
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
                            std::for_each(curr_profile.begin(), curr_profile.end(), [&](uint16_t &x){x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);});

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
                        std::for_each(curr_profile.begin(), curr_profile.end(), [&](uint16_t &x){x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);});

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
                std::for_each(curr_profile.begin(), curr_profile.end(), [&](uint16_t &x){x = std::min((size_t) MAXLCPVALUE, x+num_LF_steps);});

                listings_fd << "[" << 0 << "," << end_pos_of_match << "] ";
                length = std::min((size_t) MAXLCPVALUE, end_pos_of_match+1);

                process_profile(curr_profile, length);
                listings_fd << "\n";
            }
        }

        size_t serialize(std::ostream &out, std::ostream &out_bwt, int read_length, 
                        sdsl::structure_tree_node *v = nullptr, std::string name = "") {
            /* serializes data-structures to disk for seeing index size */

            // Build int-vectors to write to disk ...
            size_t max_width = (read_length == 0) ? (DOCWIDTH * 8) : std::ceil(std::log2(read_length));
            max_width = std::min(static_cast<size_t>(DOCWIDTH * 8), max_width);

            std::cout << "\n";
            FORCE_LOG("serialize", "writing doc profiles to disk where each entry is %d bits", max_width);

            sdsl::int_vector<> start_da_profiles (this->r * num_docs, 0, max_width);
            sdsl::int_vector<> end_da_profiles (this->r * num_docs, 0, max_width);

            // Fill in the integer vectors in order BWT runs
            bool is_start = false, is_end = false;
            uint8_t prev_ch = 0, curr_ch = 0;

            size_t n = this->bwt.size();
            size_t curr_run_num  = 0;
            std::vector<size_t> ch_pos (256, 0);
            
            std::vector<uint16_t> curr_start_profile(num_docs, 0);
            std::vector<uint16_t> curr_end_profile(num_docs, 0);

            uint16_t max_poss_lcp = (read_length == 0) ? (std::pow(2, DOCWIDTH*8)-1) : read_length;

            for (size_t i = 0; i < n; i++){
                // Determine if we have started a new run
                curr_ch = this->bwt[i];
                is_start = (curr_ch != prev_ch);
                prev_ch = curr_ch;

                // Grab the profiles at the start and end of this run (if it is)
                if (is_start) {
                    curr_ch = (curr_ch == 1) ? 0 : curr_ch; // TODO: figure out why $ is 1, but placed in 0
                    size_t curr_pos = ch_pos[curr_ch];

                    curr_start_profile = start_doc_profiles[curr_ch][curr_pos];
                    curr_end_profile = end_doc_profiles[curr_ch][curr_pos];

                    // Write the profiles to the int vectors
                    size_t start_pos = curr_run_num * num_docs;
                    for (size_t j = 0; j < num_docs; j++) {
                        start_da_profiles[start_pos + j] = std::min(curr_start_profile[j], max_poss_lcp);
                        end_da_profiles[start_pos + j] = std::min(curr_end_profile[j], max_poss_lcp); 
                    }
                    curr_run_num++; ch_pos[curr_ch]++;
                }
                is_start = false;
            }

            // Bit-compress (will mostly affect the size when not specifying the max read length)
            sdsl::util::bit_compress(start_da_profiles);
            sdsl::util::bit_compress(end_da_profiles);

            // Start to write data to files ...
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_t written_bytes = 0;

            out.write((char *)&this->num_docs, sizeof(this->num_docs));
            written_bytes += sizeof(this->num_docs);
            written_bytes += start_da_profiles.serialize(out);
            written_bytes += end_da_profiles.serialize(out);
            written_bytes += this->bwt.serialize(out_bwt);

            FORCE_LOG("serialize", "finished writing index to *.docprofiles and *.docprofiles.bwt files");
            
            return written_bytes;
        }

        std::pair<size_t, size_t> build_ftab(std::string output_ref) {
            /* build the f_tab and store it in the *.fna.ftab file */

            // table to convert 2-bits into a character (00, 01, 10, 11)
            char nuc_tab[] = {'A', 'C', 'T', 'G'}; 

            // initialize some variables
            size_t num_entries = std::pow(FTAB_ALPHABET_SIZE, FTAB_ENTRY_LENGTH);
            std::vector<size_t> index_vec(FTAB_ENTRY_LENGTH) ; 
            std::iota(index_vec.begin(), index_vec.end(), 0); 
            size_t ftab_file_size = num_entries * FTAB_ENTRY_SIZE;

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
            assert(mmap_ftab != MAP_FAILED);

            // DEBUG: open file and store results for checking
            std::ofstream debug_fd(output_ref + ".ftab.debug.txt");
                      
            // iterate through all possible 10-mer in dictionary
            size_t num_not_found = 0, num_found = 0;
            for (size_t loop_index = 0; loop_index < num_entries; loop_index++) {
                std::string curr_seq(FTAB_ENTRY_LENGTH, '*');

                // generates current 10-mer by converting every 2-bits into a character
                std::transform(index_vec.begin(), index_vec.end(), curr_seq.begin(),
                               [&] (size_t pos) {
                                    uint8_t code = FTAB_GRAB_CODE(loop_index, pos);
                                    assert(code <= 3);
                                    return nuc_tab[code];});
                
                // initialize variables for backward search
                size_t start = 0, end = this->bwt.size();
                bool not_found = false;

                // initialize variables that are need to grab correct profile 
                uint8_t curr_prof_ch = 0;
                size_t curr_prof_pos = 0, num_LF_steps = 0;
                bool use_start = false, use_end = false;

                // perform backward search for current 10-mer
                int j = (FTAB_ENTRY_LENGTH-1);
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
                assert(sizeof(size_t) == 8);

                for (int k = 0; k < sizeof(size_t); k++) {
                    /* 1. BWT range start position (little-endian) */
                    mmap_ftab[start_pos + k] = ((0xFF << (8*k)) & start) >> (8*k);
                    /* 2. BWT range end position (little-endian) */
                    mmap_ftab[start_pos + sizeof(size_t) + k] = ((0xFF << (8*k)) & end) >> (8*k);
                    /* 3. Current profile position (little-endian) */
                    mmap_ftab[start_pos + (2*sizeof(size_t)) + k] = ((0xFF << (8*k)) & curr_prof_pos) >> (8*k);
                }

                // write the remaining two elements of ftab entry
                assert(num_LF_steps < 10);

                /* 4. Number of LF Steps */
                mmap_ftab[start_pos + (3*sizeof(size_t))] = (0xFF & num_LF_steps);

                /* 5. BWT character */
                mmap_ftab[start_pos + (3*sizeof(size_t)) + 1] = (0xFF & curr_prof_ch);

                debug_fd << curr_seq << "," << start 
                                     << "," << end 
                                     << "," << curr_prof_pos 
                                     << "," << num_LF_steps 
                                     << "," << curr_prof_ch << std::endl;
            }

            // unmap the ftab and debug files
            munmap(mmap_ftab, ftab_file_size);
            close(ftab_fd); 
            debug_fd.close();

            return std::make_pair(num_found, num_not_found);
        }

    private:
        /*********************************/ 
        /* Private instance variables
        /*********************************/
        size_t num_docs = 0;
        bool rle = true;

        // This vectors has the following dimensions: [256][num of ith char][num_docs]
        // This structure stores the DA profiles for each
        // character separately.
        std::vector<std::vector<std::vector<uint16_t>>> start_doc_profiles;
        std::vector<std::vector<std::vector<uint16_t>>> end_doc_profiles;

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
                FATAL_ERROR(("open() file " + fname + " failed").data());
            if (fstat(fileno(fd), &filestat) < 0)
                FATAL_ERROR(("stat() file " + fname + " failed").data());

            num_docs = 0;
            if ((fread(&num_docs, sizeof(size_t), 1, fd)) != 1)
                FATAL_ERROR(("fread() file " + fname + " failed").data()); 
            fclose(fd);

            if (filestat.st_size != ((num_docs * this->r * DOCWIDTH) + sizeof(size_t) + this->r))
                FATAL_ERROR("invalid file size for *dap files");      
        }

        static void read_doc_profiles(std::vector<std::vector<std::vector<uint16_t>>>& prof_matrix, std::string input_file, size_t num_docs, 
                                    size_t num_runs, std::string output_path, size_t num_profiles) {
            /* loads a set of document array profiles into their respective matrix */

            // First, lets open the file and verify the size/# of docs are valid
            struct stat filestat; FILE *fd;

            if ((fd = fopen(input_file.c_str(), "r")) == nullptr)
                error("open() file " + input_file + " failed");

            if (fstat(fileno(fd), &filestat) < 0)
                error("stat() file " + input_file + " failed");

            size_t num_docs_found = 0;
            if ((fread(&num_docs_found, sizeof(size_t), 1, fd)) != 1)
                error("fread() file " + input_file + " failed"); 

            ASSERT((num_docs_found == num_docs), "mismatch in the number of documents.");
            ASSERT((filestat.st_size == ((num_docs * num_runs * DOCWIDTH) + sizeof(size_t) + num_runs)), "invalid file size.");

            // Second, lets open files if we want to write out some number of the profiles
            std::ofstream dap_csv_file;
            bool print_to_file = (output_path.size() > 9); // Not just .sdap.csv
            size_t profiles_to_print = 0;

            if (print_to_file){ 
                dap_csv_file.open(output_path);
                profiles_to_print = (num_profiles == 0 || num_profiles > num_runs) ? num_runs : num_profiles; 
            }

            // Thirdly, go through the rest of file and fill in the profiles. Each 
            // profile will start with the BWT character which we will use figure out which
            // list to put it in.
            size_t curr_val = 0;
            for (size_t i = 0; i < num_runs; i++) {
                uint8_t curr_bwt_ch = 0;

                if (fread(&curr_bwt_ch, 1, 1, fd) != 1)
                    FATAL_ERROR("issue occurred while reading in bwt character from doc profiles file.");

                std::vector<uint16_t> curr_profile (num_docs, 0);
                for (size_t j = 0; j < num_docs; j++) {
                    if ((fread(&curr_val, DOCWIDTH, 1, fd)) != 1)
                        error("fread() file " + input_file + " failed"); 
                    curr_profile[j] = curr_val;
                    
                    // Check if we want to print
                    if (print_to_file) {
                        if (j < (num_docs-1))
                            dap_csv_file << curr_val << ",";
                        else
                            dap_csv_file << curr_val << "\n";
                    }
                }
                prof_matrix[curr_bwt_ch].push_back(curr_profile);
            }
            fclose(fd);

            if (print_to_file)
                dap_csv_file.close();
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

        ulint LF(ri::ulint i, ri::uchar c){
            // number of c before the interval
            ri::ulint c_before = this->bwt.rank(i, c);
            // number of c inside the interval rn
            ri::ulint l = this->F[c] + c_before;
            return l;
        }
};

#endif /* end of include guard:_DOC_QUERIES_H */
