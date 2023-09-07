/*
 * File: tax_doc_queries.hpp
 * Description: ....
 * Date: Aug. 30th, 2023
 */

#ifndef _TAX_DOC_QUERIES_H
#define _TAX_DOC_QUERIES_H

//#define READ_NUM_DOCS(fd) unsigned(fd[0])
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
        std::vector<std::vector<std::vector<uint16_t>>> start_doc_profiles;
        std::vector<std::vector<std::vector<uint16_t>>> end_doc_profiles;

        tax_doc_queries(std::string filename,
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


            // size_t num = READ_NUM_DOCS(mmap_sdap_of);
            // std::cout << "num_docs = " << num << std::endl;

            // num = READ_NUM_PAIRS(mmap_sdap_of, 8);
            // std::cout << num << std::endl;

            // num = READ_NUM_PAIRS(mmap_sdap_of, 9);
            // std::cout << num << std::endl;

            // num = READ_DOC_ID_OR_LCP_VAL(mmap_sdap_of, 10);
            // std::cout << num << std::endl;

            // num = READ_DOC_ID_OR_LCP_VAL(mmap_sdap_of, 12);
            // std::cout << num << std::endl;

            // num = READ_NUM_PAIRS(mmap_sdap_of, 14);
            // std::cout << num << std::endl;

            // num = READ_NUM_PAIRS(mmap_sdap_of, 15);
            // std::cout << num << std::endl;

            // num = READ_DOC_ID_OR_LCP_VAL(mmap_sdap_of, 16);
            // std::cout << num << std::endl;

            // num = READ_DOC_ID_OR_LCP_VAL(mmap_sdap_of, 18);
            // std::cout << num << std::endl;

            // num = READ_NUM_PAIRS(mmap_sdap_of, 20);
            // std::cout << num << std::endl;

            // num = READ_NUM_PAIRS(mmap_sdap_of, 21);
            // std::cout << num << std::endl;

            // num = READ_DOC_ID_OR_LCP_VAL(mmap_sdap_of, 22);
            // std::cout << num << std::endl;

            // num = READ_DOC_ID_OR_LCP_VAL(mmap_sdap_of, 24);
            // std::cout << num << std::endl;
        }

        ~tax_doc_queries()
        {
            /* Destructor - unmap files and close files */
            munmap(mmap_sdap_of, file_info_sdap.st_size);
            munmap(mmap_edap_of, file_info_edap.st_size);
            close(sdap_of_fd); close(edap_of_fd);
            csv_sdap_output.close(); csv_edap_output.close();
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

        void read_doc_profiles_main_table(std::vector<std::vector<std::vector<uint16_t>>>& prof_matrix, 
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
                std::vector<uint16_t> curr_record (num_cols * 4, 0);
                for (size_t j = 0; j < (num_cols*4); j++) {
                    if ((fread(&curr_val, DOCWIDTH, 1, fd)) != 1)
                        FATAL_ERROR("fread() failed"); 
                    curr_record[j] = curr_val;
                }

                // Step 3: read the overflow pointer value
                size_t of_ptr = 0;
                if (fread(&of_ptr, 8, 1, fd) != 1)
                    FATAL_ERROR("fread() failed");

                // If we are print out the profiles ...
                if (print_profiles && i < profiles_to_print) {
                    *out_fd << curr_bwt_ch << ",";
                    for (auto x: curr_record)
                        *out_fd << x << ",";
                    *out_fd << of_ptr << "\n";
                }
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