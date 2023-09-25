/*
 * File: pfp_doc.hpp
 * Description: Header file for pfp_doc data-structure, contains the
 *              construction of the document profile data-structure
 *              and other needed structs for pfp_doc.cpp
 * Date: August 30th, 2022
 */

#ifndef PFP_DOC_H
#define PFP_DOC_H

#include <string>
#include <iostream>
#include <filesystem>
#include <vector>

/* Useful MACROs */
#define FATAL_ERROR(...) do {std::fprintf(stderr, "\n\033[31mError: \033[m"); std::fprintf(stderr, __VA_ARGS__);\
                              std::fprintf(stderr, "\n\n"); std::exit(1);} while(0)
#define ASSERT(condition, msg) do {if (!condition){std::fprintf(stderr, "\n\n\033[31mAssertion Failed:\033[m %s\n\n", msg); \
                                                   std::exit(1);}} while(0)
#define STATUS_LOG(x, ...) do {std::fprintf(stderr, "\033[32m[%s] \033[0m", x); std::fprintf(stderr, __VA_ARGS__ ); \
                               std::fprintf(stderr, " ... ");} while(0)
#define DONE_LOG(x) do {auto sec = std::chrono::duration<double>(x); \
                        std::fprintf(stderr, "done.  (%.3f sec)\n", sec.count());} while(0)
#define FORCE_LOG(func, ...)  do {std::fprintf(stderr, "\033[32m[%s] \033[m", func); \
                                  std::fprintf(stderr, __VA_ARGS__); \
                                  std::fprintf(stderr, "\n");} while (0)

// Defintions
#define PFPDOC_VERSION "1.0.8"

#define DOCWIDTH 2
#define MAXQUEUELENGTH 1000000
#define MAXLCPVALUE 65535 // 2^16 - 1
#define MAXDOCS 65535

#define AVX2_PRESENT __AVX2__ 
#define AVX512BW_PRESENT __AVX512BW__ 


/* Function declations */
int pfpdoc_usage();
int build_main(int argc, char** argv);
int run_main(int argc, char** argv);
int info_main(int argc, char** argv);
int pfpdoc_build_usage();
int pfpdoc_run_usage();
int pfpdoc_info_usage();
int is_file(std::string path);
int is_dir(std::string path);
std::vector<std::string> split(std::string input, char delim);
bool is_integer(const std::string& str);
bool endsWith(const std::string& str, const std::string& suffix);
std::string execute_cmd(const char* cmd);

struct PFPDocInfoOptions {
    public:
        std::string ref_file = "";
        std::string output_path = "";
        size_t num_profiles = 0;
        bool use_taxcomp = false;
        bool use_topk = false;
        size_t num_cols =  0;

        void validate() {
            /* checks the arguments and makes sure they are valid */
            if (!ref_file.size())
                FATAL_ERROR("An index prefix must be provided using the -r option.");
            
            // check the input files
            ref_file += ".fna";
            if (!is_file(ref_file))
                FATAL_ERROR("The index prefix provided is not valid. %s does not exist.", ref_file.data());
            std::filesystem::path p (output_path);
            if (output_path.size() && !is_dir(p.parent_path().string()))
                FATAL_ERROR("The output path for profiles is not valid.");  

            // check the core index files (bwt)
            if (!is_file(ref_file + ".bwt.heads") || !is_file(ref_file + ".bwt.len"))
                FATAL_ERROR("At least one of the index files is not present.");
            
            // check for the document array profiles files
            if ((!use_taxcomp && !use_topk) && (!is_file(ref_file + ".sdap") || !is_file(ref_file + ".edap")))
                FATAL_ERROR("At least one of the index files is not present.");
            if (use_taxcomp && (!is_file(ref_file + ".taxcomp.sdap") || !is_file(ref_file + ".taxcomp.edap") ||
                                !is_file(ref_file + ".taxcomp.of.sdap") || !is_file(ref_file + ".taxcomp.of.edap")))
                FATAL_ERROR("At least one of the index files is not present.");
            if (use_topk && (!is_file(ref_file + ".topk.sdap") || !is_file(ref_file + ".topk.edap")))
                FATAL_ERROR("At least one of the index files is not present.");

            // check the structure type, can only be 1
            if (use_taxcomp && use_topk)
                FATAL_ERROR("Cannot use both taxonomic and top-k document array structures.");
            
            // checks the number of column argument
            if ((use_taxcomp || use_topk) && (num_cols < 2 || num_cols > 20))
                FATAL_ERROR("Invalid number of columns in compressed table, make sure to set it with -c, --num-col");
        }
};

struct PFPDocBuildOptions {
    public:
        std::string input_list = "";
        std::string output_prefix = "";
        std::string output_ref = "";
        std::string temp_prefix = "";
        std::string tmp_size_str = "";
        bool use_rcomp = false;
        size_t pfp_w = 10;
        size_t hash_mod = 100;
        size_t threads = 0;
        bool is_fasta = true;
        bool use_taxcomp = false;
        bool use_topk = false;
        size_t numcolsintable = 7;
        size_t doc_to_extract = 0;
        size_t use_heuristics = true;
        bool use_two_pass = false;
        size_t tmp_size = 0;

        void validate() {
            /* checks the arguments and make sure they are valid */
            if (input_list.length() && !is_file(input_list)) // provided a file-list
                FATAL_ERROR("The provided file-list is not valid.");
            else if (input_list.length() == 0)
                FATAL_ERROR("need to provide a file-list for processing.");

            std::filesystem::path p (output_prefix);
            if (!is_dir(p.parent_path().string()))
                FATAL_ERROR("output path prefix is not in a valid directory."); 
            
            if (use_two_pass) {
                std::filesystem::path p (temp_prefix);
                if (!is_dir(p.parent_path().string()))
                    FATAL_ERROR("output path prefix for temporary file is not valid.");
                if (tmp_size_str.length() == 0 || tmp_size_str.find("GB") == std::string::npos)
                    FATAL_ERROR("temporary file size argument needs to be in this form (e.g. 4GB)");

                tmp_size_str.erase(tmp_size_str.find("GB"), tmp_size_str.length());
                tmp_size = std::atoi(tmp_size_str.data());
                tmp_size *= 1073741824;
            }

            if (use_taxcomp && use_topk)
                FATAL_ERROR("taxonomic and top-k compression cannot be used together.");   

            if (doc_to_extract > 0 && (use_taxcomp || use_topk)) 
                FATAL_ERROR("cannot extract document array when using compression.");       
        }
};

struct PFPDocRunOptions {
    public:
        std::string ref_file = "";
        std::string pattern_file = "";
        int read_length = 0;
        int num_cols = -1;
        bool write_to_file = false;
        bool use_taxcomp = false;
        bool use_topk = false;
    
    void validate() {
        /* checks arguments and makes sure they are valid files */
        ref_file += ".fna";

        // check the input files
        if (!is_file(ref_file) || !is_file(pattern_file))
            FATAL_ERROR("At least one of the input files is not valid.");
        
        // check that the upper bound on the read length is positive
        if (read_length < 0)
            FATAL_ERROR("Length of read must be positive.");

        // check the index files
        if (!is_file(ref_file + ".bwt.heads") || !is_file(ref_file + ".bwt.len"))
            FATAL_ERROR("At least one of the index files is not present.");

        // make sure that we didn't try to turn on both types of doc array
        if (use_taxcomp && use_topk)
            FATAL_ERROR("Cannot use both a 'taxonomic' and 'top-k' compressed doc array.");
        
        // depending on the compression scheme used, lets
        // check that the files are present
        if (!use_taxcomp && !use_topk){
            if (!is_file(ref_file + ".sdap") || !is_file(ref_file + ".edap"))
                FATAL_ERROR("One or more of the document array files (*.sdap, *.edap) is not present.");

        } else if (use_taxcomp) {
            if (!is_file(ref_file + ".taxcomp.sdap") || !is_file(ref_file + ".taxcomp.edap")
               || !is_file(ref_file + ".taxcomp.of.sdap") || !is_file(ref_file + ".taxcomp.of.edap"))
               FATAL_ERROR("One or more of the 'taxonomic' document array files is not present.");

        } else if (use_topk) {
            if (!is_file(ref_file + ".topk.sdap") || !is_file(ref_file + ".topk.edap"))
                FATAL_ERROR("One or more of the 'top-k' document array files is not present.");
        }

        // check if the user provided a number of columns for table
        if ((use_taxcomp || use_topk) && (num_cols < 1 || num_cols > 20))
            FATAL_ERROR("the number of columns in the document array is not valid or not provided.");
    }
};

struct HelperPrograms {
  /* Contains paths to run helper programs */
  std::string base_path = "";
  std::string parseNT_bin = "newscanNT.x";
  std::string parse_fasta_bin = "newscan.x";
  std::string parse_bin = "pscan.x";
  
public:
  void build_paths(std::string base) {
      /* Takes the base path, and combines it with names to generate executable paths */
      base_path.assign(base);
      parseNT_bin.assign(base_path + parseNT_bin);
      parse_fasta_bin.assign(base_path + parse_fasta_bin);
      parse_bin.assign(base_path + parse_bin);
  }

  void validate() const {
      /* Makes sure that each path for an executable is valid */
      bool invalid_path = !is_file(parseNT_bin) | !is_file(parse_fasta_bin) | !is_file(parse_bin);
      if (invalid_path) {FATAL_ERROR("One or more of helper program paths are invalid.");}
  }
};

/* Function Declartions involving structs */
void parse_build_options(int argc, char** argv, PFPDocBuildOptions* opts);
void parse_run_options(int argc, char** argv, PFPDocRunOptions* opts);
void parse_info_options(int argc, char** argv, PFPDocInfoOptions* opts);
void print_build_status_info(PFPDocBuildOptions* opts);
void run_build_parse_cmd(PFPDocBuildOptions* build_opts, HelperPrograms* helper_bins);

#endif /* End of PFP_DOC_H */