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
#define FATAL_ERROR(...) do {std::fprintf(stderr, "\nError: "); std::fprintf(stderr, __VA_ARGS__);\
                              std::fprintf(stderr, "\n\n"); std::exit(1);} while(0)
#define ASSERT(condition, msg) do {if (!condition){std::fprintf(stderr, "Assertion Failed: %s\n", msg); \
                                                   std::exit(1);}} while(0)
#define STATUS_LOG(x, ...) do {std::fprintf(stderr, "[%s] ", x); std::fprintf(stderr, __VA_ARGS__ ); \
                               std::fprintf(stderr, " ... ");} while(0)
#define DONE_LOG(x) do {auto sec = std::chrono::duration<double>(x); \
                        std::fprintf(stderr, "done.  (%.3f sec)\n", sec.count());} while(0)
#define FORCE_LOG(func, ...)  do {std::fprintf(stderr, "[%s] ", func); \
                                  std::fprintf(stderr, __VA_ARGS__); \
                                  std::fprintf(stderr, "\n");} while (0)

// Defintions
#define DOCWIDTH 1 // 5
#define MAXQUEUELENGTH 1000000
#define PFPDOC_VERSION "1.0.1"

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

            // check the index files
            if (!is_file(ref_file + ".bwt.heads") || !is_file(ref_file + ".bwt.len")
                || !is_file(ref_file + ".sdap") || !is_file(ref_file + ".edap"))
                FATAL_ERROR("At least one of the index files is not present.");
        }
};

struct PFPDocBuildOptions {
    public:
        std::string input_list = "";
        std::string output_prefix = "";
        std::string output_ref = "";
        bool use_rcomp = false;
        size_t pfp_w = 10;
        size_t hash_mod = 100;
        size_t threads = 0;
        bool is_fasta = true;

        void validate() {
            /* checks the arguments and make sure they are valid */
            if (input_list.length() && !is_file(input_list)) // provided a file-list
                FATAL_ERROR("The provided file-list is not valid.");
            else if (input_list.length() == 0)
                FATAL_ERROR("Need to provide a file-list for processing.");

            std::filesystem::path p (output_prefix);
            if (!is_dir(p.parent_path().string()))
                FATAL_ERROR("Output path prefix is not in a valid directory.");                 
        }
};

struct PFPDocRunOptions {
    public:
        std::string ref_file = "";
        std::string pattern_file = "";
        int read_length = 0;
        bool write_to_file = false;
    
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
        if (!is_file(ref_file + ".bwt.heads") || !is_file(ref_file + ".bwt.len")
            || !is_file(ref_file + ".sdap") || !is_file(ref_file + ".edap"))
            FATAL_ERROR("At least one of the index files is not present.");
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