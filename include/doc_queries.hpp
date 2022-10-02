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

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd>
class doc_queries : ri::r_index<sparse_bv_type, rle_string_t>
{
    public:

    size_t num_docs = 0;

    // This vectors has the following dimensions: [256][num of ith char][num_docs]
    // This structure stores the DA profiles for each
    // character separately.
    std::vector<std::vector<std::vector<size_t>>> start_doc_profiles;
    std::vector<std::vector<std::vector<size_t>>> end_doc_profiles;

    typedef size_t size_type;

    doc_queries(std::string filename, bool rle = true): 
            ri::r_index<sparse_bv_type, rle_string_t>(),
            start_doc_profiles(256, std::vector<std::vector<size_t>>(0, std::vector<size_t>(0))),
            end_doc_profiles(256, std::vector<std::vector<size_t>>(0, std::vector<size_t>(0)))
    {
        std::string bwt_fname = filename + ".bwt";

        STATUS_LOG("build_profiles", "loading the bwt of the input text");
        auto start = std::chrono::system_clock::now();

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
        DONE_LOG((std::chrono::system_clock::now() - start));

        // gather some statistics on the BWT
        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        size_t log_r = bitsize(uint64_t(this->r));
        size_t log_n = bitsize(uint64_t(this->bwt.size()));

        FORCE_LOG("build_profiles", "bwt statistics: n = %ld, r = %ld\n" , this->bwt.size(), this->r);
        
        // determine the number of documents and verify the that file
        // sizes are correct.
        std::string tmp_filename = filename + ".sdap";

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        if (fstat(fileno(fd), &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        num_docs = 0;
        if ((fread(&num_docs, sizeof(size_t), 1, fd)) != 1)
            error("fread() file " + tmp_filename + " failed"); 
        
        ASSERT((filestat.st_size == ((num_docs * this->r * DOCWIDTH) + sizeof(size_t) + this->r)), "invalid file size.");
        fclose(fd);
        
        // load the profiles for starts and ends
        STATUS_LOG("build_profiles", "loading the document array profiles");
        start = std::chrono::system_clock::now();

        read_doc_profiles(start_doc_profiles, filename + ".sdap", this->num_docs, this->r);
        read_doc_profiles(end_doc_profiles, filename + ".edap", this->num_docs, this->r);

        DONE_LOG((std::chrono::system_clock::now() - start));
    }

    static void read_doc_profiles(std::vector<std::vector<std::vector<size_t>>>& prof_matrix, std::string input_file, size_t num_docs, size_t num_runs) {
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

        // Secondly, go through the rest of file and fill in the profiles. Each 
        // profile will start with the BWT character which we will use figure out which
        // list to put it in.
        size_t curr_val = 0;
        for (size_t i = 0; i < num_runs; i++) {
            uint8_t curr_bwt_ch = 0;

            if (fread(&curr_bwt_ch, 1, 1, fd) != 1)
                FATAL_ERROR("issue occurred while reading in bwt character from doc profiles file.");
            
            std::vector<size_t> curr_profile (num_docs, 0);
            for (size_t j = 0; j < num_docs; j++) {
                if ((fread(&curr_val, DOCWIDTH, 1, fd)) != 1)
                    error("fread() file " + input_file + " failed"); 
                curr_profile[j] = curr_val;
            }
            prof_matrix[curr_bwt_ch].push_back(curr_profile);
        }
        fclose(fd);
    }

    void query_profiles(std::string pattern_file) {
        /* Takes in a file of reads, and lists all the documents containing the read */
    
        // Open output/input files
        std::ofstream listings_fd (pattern_file + ".listings");

        gzFile fp; kseq_t* seq;
        fp = gzopen(pattern_file.data(), "r"); 

        if(fp == 0) {std::exit(1);}
        seq = kseq_init(fp);

        // lambda to print out the document listing
        auto process_profile = [&](std::vector<size_t> profile, size_t length) {
                std::vector<size_t> docs_found;
                std::string output_str = "{";
                bool found_one = false;

                for (size_t i = 0; i < profile.size(); i++) {
                    if (profile[i] >= length) {
                        output_str += "," + std::to_string(i);
                        found_one = true;
                    }
                }
                if (found_one)
                    output_str = "{" + output_str.substr(2) + "} ";
                else    
                    output_str = "{} ";
                listings_fd << output_str;
        };

    
        // Process each read, and print out the document lists
        while (kseq_read(seq)>=0) {
            
            // Uppercase every character in read
			for (size_t i = 0; i < seq->seq.l; ++i) {
				seq->seq.s[i] = static_cast<char>(std::toupper(seq->seq.s[i]));
            }

            size_t start = 0, end = this->bwt.size(), end_pos_of_match = seq->seq.l-1;
            std::vector<size_t> curr_profile (this->num_docs, 0);

            listings_fd << ">" << seq->name.s << "\n";

            // Pointer variables to current profile being "used"
            uint8_t curr_prof_ch = 0;
            size_t curr_prof_pos = 0, num_LF_steps = 0;

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
                        // grab the current profile, and update with steps
                        curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];
                        std::for_each(curr_profile.begin(), curr_profile.end(), [&](size_t &x){x+=num_LF_steps;});

                        listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";
                        process_profile(curr_profile, (end_pos_of_match-i));
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
                } 
                // range is within BWT run, but wrong character 
                else if (this->bwt[start] != next_ch) 
                {
                    // grab the current profile, and update with steps
                    curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];
                    std::for_each(curr_profile.begin(), curr_profile.end(), [&](size_t &x){x+=num_LF_steps;});

                    listings_fd << "[" << (i+1) << "," << end_pos_of_match << "] ";
                    process_profile(curr_profile, (end_pos_of_match-i));
                    end_pos_of_match = i;

                    start = 0; end = this->bwt.size();
                    start_run = 0; end_run = this->r - 1;
                    num_ch_before_start = 0;
                    num_ch_before_end = this->bwt.number_of_letter(next_ch);

                    // grab any profile at run boundary of next_ch (choose first one)
                    curr_prof_ch = next_ch;
                    curr_prof_pos = this->bwt.run_head_rank(start_run, next_ch);
                    num_LF_steps = 0;
                }
                // range is within BWT run, and is the correct character
                else 
                {
                    num_LF_steps++;
                    //std::transform(curr_profile.begin(), curr_profile.end(), curr_profile.begin(), 
                    //                [](size_t x) { return (++x); });   
                }

                // Perform an LF step
                start = num_ch_before_start + this->F[next_ch]; 
                end = num_ch_before_end + this->F[next_ch];
            }
            // grab the current profile, and update with steps
            curr_profile = start_doc_profiles[curr_prof_ch][curr_prof_pos];
            std::for_each(curr_profile.begin(), curr_profile.end(), [&](size_t &x){x+=num_LF_steps;});

            listings_fd << "[" << 0 << "," << end_pos_of_match << "] ";
            process_profile(curr_profile, end_pos_of_match+1);
            listings_fd << "\n";
        }
    }

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths)
    {
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

    ulint LF(ri::ulint i, ri::uchar c)
    {
        // //if character does not appear in the text, return empty pair
        // if ((c == 255 and this->F[c] == this->bwt_size()) || this->F[c] >= this->F[c + 1])
        //     return {1, 0};
        //number of c before the interval
        ri::ulint c_before = this->bwt.rank(i, c);
        // number of c inside the interval rn
        ri::ulint l = this->F[c] + c_before;
        return l;
    }

    size_t serialize(std::ostream &out, std::ostream &out_bwt, int read_length, 
                     sdsl::structure_tree_node *v = nullptr, std::string name = "") {
        /* serializes data-structures to disk for seeing index size */

        /*
        // Build int-vectors to write to disk ...
        size_t max_width = (read_length == 0) ? (DOCWIDTH * 8) : std::ceil(std::log2(read_length));
        max_width = std::min(static_cast<size_t>(DOCWIDTH * 8), max_width);

        std::cout << "\n";
        FORCE_LOG("serialize", "writing doc profiles to disk where each entry is %d bits", max_width);

        sdsl::int_vector<> start_da_profiles (this->r * num_docs, 0, max_width);
        sdsl::int_vector<> end_da_profiles (this->r * num_docs, 0, max_width);

        size_t iter_pos = 0;
        for (size_t i = 0; i < start_doc_profiles.size(); i++) {
            for (size_t j = 0; j < num_docs; j++) {
                start_da_profiles[iter_pos] = start_doc_profiles[i][j];
                end_da_profiles[iter_pos] = end_doc_profiles[i][j];
                iter_pos++;
            }
        }

        // Start to write data to files ...
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        out.write((char *)&this->num_docs, sizeof(this->num_docs));
        written_bytes += sizeof(this->num_docs);
        written_bytes += start_da_profiles.serialize(out);
        written_bytes += end_da_profiles.serialize(out);

        written_bytes += this->bwt.serialize(out_bwt);

        FORCE_LOG("serialize", "finished writing index to *.docprofiles and *.docprofiles.bwt files");
        */
        return 0;
        //return written_bytes;

    }

};

#endif /* end of include guard:_DOC_QUERIES_H */
