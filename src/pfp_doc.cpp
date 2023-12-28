/*
 * File: pfp_doc.cpp
 * Description: Main file for building the document array profiles
 *              based on the prefix-free parse. This workflow is 
 *              based on the pfp_lcp.hpp developed by Massimiliano
 *              Rossi.
 *
 * Date: August 30th, 2022
 */

#include <pfp_doc.hpp> 
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <filesystem>
#include <ref_builder.hpp>
#include <vector>
#include <pfp.hpp>
#include <pfp_lcp_doc.hpp>
#include <pfp_lcp_doc_two_pass.hpp>
#include <doc_queries.hpp>
#include <tax_doc_queries.hpp>
#include <topk_doc_queries.hpp>
#include <immintrin.h>
#include <getopt.h>
#include <queue>

int build_main(int argc, char** argv) {
    /* main method for build the document profiles */
    if (argc == 1) return pfpdoc_build_usage();
    auto build_start = std::chrono::system_clock::now();

    // grab the command-line options, and validate them
    PFPDocBuildOptions build_opts;
    parse_build_options(argc, argv, &build_opts);
    build_opts.validate();

    // determine output path for reference, and print all info
    build_opts.output_ref.assign(build_opts.output_prefix + ".fna");
    print_build_status_info(&build_opts);

    // build the input reference file, and bitvector labeling the end for each doc
    STATUS_LOG("build_main", "building the reference file based on file-list");
    auto start = std::chrono::system_clock::now();
    RefBuilder ref_build(build_opts.input_list, build_opts.output_prefix, build_opts.use_rcomp);
    DONE_LOG((std::chrono::system_clock::now() - start));

    // make sure that document numbers can be store in 2 bytes, and 
    // it makes sense with respect to k
    if (ref_build.num_docs >= MAXDOCS)
        FATAL_ERROR("An index cannot be build over %ld documents, "
                    "please reduce to a max of 65,535 docs.", ref_build.num_docs);
    if ((build_opts.use_taxcomp || build_opts.use_topk) 
        && ref_build.num_docs < build_opts.numcolsintable)
        FATAL_ERROR("the k provided is larger than the number of documents.");

    // check and make sure the document to extract is a valid id
    if (build_opts.doc_to_extract > ref_build.num_docs) {
        FATAL_ERROR("Document #%d was requested to be extracted,"
                    " but there are only %d documents", build_opts.doc_to_extract, ref_build.num_docs);
    }

    // determine the paths to the BigBWT executables
    HelperPrograms helper_bins;
    if (!std::getenv("PFPDOC_BUILD_DIR")) {FATAL_ERROR("Need to set PFPDOC_BUILD_DIR environment variable.");}
    helper_bins.build_paths((std::string(std::getenv("PFPDOC_BUILD_DIR")) + "/bin/").data());
    helper_bins.validate();

    // parse the input text with BigBWT
    STATUS_LOG("build_main", "generating the prefix-free parse for given reference");
    start = std::chrono::system_clock::now();
    run_build_parse_cmd(&build_opts, &helper_bins);
    DONE_LOG((std::chrono::system_clock::now() - start));

    // load the parse and dictionary into pf object
    STATUS_LOG("build_main", "building the parse and dictionary objects");
    start = std::chrono::system_clock::now();
    pf_parsing pf(build_opts.output_ref, build_opts.pfp_w);
    DONE_LOG((std::chrono::system_clock::now() - start));

    // print info regarding the compression scheme being used
    std::cerr << "\n";
    if (build_opts.use_taxcomp)
        FORCE_LOG("build_main", "taxonomic compression of the doc profiles will be used");
    else if (build_opts.use_topk)
        FORCE_LOG("build_main", "top-k compression of the doc profile will be used");
    else   
        FORCE_LOG("build_main", "no compression scheme will be used for the doc profiles");

    // builds the BWT, SA, LCP, and document array profiles and writes to a file
    size_t num_runs = 0;
    if (!build_opts.use_two_pass) {
        pfp_lcp lcp(build_opts.use_heuristics, build_opts.doc_to_extract, build_opts.use_taxcomp, build_opts.use_topk, 
                build_opts.numcolsintable, pf, build_opts.output_ref, &ref_build);
        num_runs = lcp.total_num_runs;
    } else {
        pfp_lcp_doc_two_pass lcp(pf, build_opts.output_ref, &ref_build, 
                                 build_opts.temp_prefix, build_opts.tmp_size,
                                 build_opts.use_taxcomp, build_opts.use_topk,
                                 build_opts.numcolsintable);
        num_runs = lcp.total_num_runs;
    }

    // print stats before closing out
    auto build_time = std::chrono::duration<double>((std::chrono::system_clock::now() - build_start));
    FORCE_LOG("build_main", "finished: build time (s) = %.2f", build_time.count());
    std::cerr << "\n";
    
    return 0;
}

int run_main(int argc, char** argv) {
    /* main method for querying the data-structure */
    if (argc == 1) return pfpdoc_run_usage();
    std::cerr << "\n";

    // grab the command-line options, and validate them
    PFPDocRunOptions run_opts;
    parse_run_options(argc, argv, &run_opts);
    run_opts.validate();

    if (!run_opts.use_taxcomp && !run_opts.use_topk) {
        // build the doc_queries object (load data-structures)
        doc_queries doc_queries_obj (run_opts.ref_file);

        // query the doc_profiles with the given reads
        std::cerr << "\n";
        STATUS_LOG("query_main", "processing the patterns");

        auto start = std::chrono::system_clock::now();
        doc_queries_obj.query_profiles(run_opts.pattern_file, run_opts.output_prefix);
        DONE_LOG((std::chrono::system_clock::now() - start));
        
        // write index to disk
        if (run_opts.write_to_file) {
            std::string outfile = run_opts.ref_file + ".docprofiles";
            std::string outfile_bwt = run_opts.ref_file + ".docprofiles.bwt";
            std::ofstream out(outfile), out_bwt(outfile_bwt);

            doc_queries_obj.serialize(out, out_bwt, run_opts.read_length);
            out.close(); out_bwt.close();
        }
    } else if (run_opts.use_taxcomp) {
        // build the tax_doc_queries object (load data-structures)
        tax_doc_queries tax_doc_queries_obj(run_opts.ref_file,
                                            run_opts.num_cols);
        // query the doc_profiles with the given reads
        STATUS_LOG("run_main", "processing the patterns");

        auto start = std::chrono::system_clock::now();
        tax_doc_queries_obj.query_profiles(run_opts.pattern_file);
        DONE_LOG((std::chrono::system_clock::now() - start));
    } else if (run_opts.use_topk) {
        // build the topk_doc_queries object (load data-structures)
        topk_doc_queries topk_doc_queries_obj(run_opts.ref_file,
                                              run_opts.num_cols);
        // query the doc_profiles with the given reads
        STATUS_LOG("run_main", "processing the patterns");

        auto start = std::chrono::system_clock::now();
        topk_doc_queries_obj.query_profiles(run_opts.pattern_file);
        DONE_LOG((std::chrono::system_clock::now() - start));
    }
    std::cerr << "\n";
    return 0;
}

int info_main(int argc, char** argv) {
    /* main method for info sub-command */
    if (argc == 1) return pfpdoc_info_usage();
    std::cerr << "\n";

    // grab the command-line options, and validate them
    PFPDocInfoOptions info_opts;
    parse_info_options(argc, argv, &info_opts);
    info_opts.validate();

    if (!info_opts.use_taxcomp && !info_opts.use_topk) {
        doc_queries doc_queries_obj (info_opts.ref_file, 
                                     info_opts.output_path, 
                                     info_opts.num_profiles);
    } else if (info_opts.use_taxcomp) {
        tax_doc_queries tax_doc_queries_obj(info_opts.ref_file,
                                            info_opts.num_cols,
                                            info_opts.num_profiles,
                                            info_opts.output_path);
    } else if (info_opts.use_topk) {
        topk_doc_queries topk_doc_queries_obj(info_opts.ref_file,
                                              info_opts.num_cols,
                                              info_opts.num_profiles,
                                              info_opts.output_path);
    }
    std::cerr << "\n";
    return 0;
}

void run_build_parse_cmd(PFPDocBuildOptions* build_opts, HelperPrograms* helper_bins) {
    // Generates and runs the command-line for executing the PFP of the reference 
    std::ostringstream command_stream;
    if (build_opts->threads > 0) {
        std::string curr_exe = "";
        if (build_opts->is_fasta) {curr_exe.assign(helper_bins->parse_fasta_bin);}
        else {curr_exe.assign(helper_bins->parse_bin);}

        command_stream << curr_exe << " "; //<< " -i ";
        command_stream << build_opts->output_ref << " ";
        command_stream << "-w " << build_opts->pfp_w;
        command_stream << " -p " << build_opts->hash_mod;
        command_stream << " -t " << build_opts->threads;
    }
    else {
        std::string curr_exe = "";
        command_stream << helper_bins->parseNT_bin << " "; // << " -i ";
        command_stream << build_opts->output_ref << " ";
        command_stream << "-w " << build_opts->pfp_w;
        command_stream << " -p " << build_opts->hash_mod;
    }
    if (build_opts->is_fasta) {command_stream << " -f";}

    // std::cout << command_stream.str() << std::endl;
    //LOG(build_opts->verbose, "build_parse", ("Executing this command: " + command_stream.str()).data());
    auto parse_log = execute_cmd(command_stream.str().c_str());
    //OTHER_LOG(parse_log.data());
}

std::string execute_cmd(const char* cmd) {
    std::array<char, 256> buffer{};
    std::string output = "";

    std::string cmd_plus_stderr = std::string(cmd) + " 2>&1";
    FILE* pipe = popen(cmd_plus_stderr.data(), "r"); // Extract stderr as well
    if (!pipe) {FATAL_ERROR("popen() failed!");}

    try {
        std::size_t bytes;
        while ((bytes = fread(buffer.data(), sizeof(char), sizeof(buffer), pipe))) {
            output += std::string(buffer.data(), bytes);
        }
    } catch (...) {
        pclose(pipe);
        FATAL_ERROR("Error occurred while reading popen() stream.");
    }

    // Check if the command failed
    size_t exit_code = WEXITSTATUS(pclose(pipe));
    if (exit_code != 0) {
        std::cout << "\n";
        FATAL_ERROR("external command failed ... here is the error message:\n%s", output.data());
    }

    return output;
}

bool endsWith(const std::string& str, const std::string& suffix) {
    // Checks if the string ends the suffix
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

bool is_integer(const std::string& str) {
    /* Checks if string passed is an integer */
    std::string::const_iterator iter = str.begin();
    while (iter != str.end() && std::isdigit(*iter)) {++iter;}
    return !str.empty() && iter == str.end();
}

std::vector<std::string> split(std::string input, char delim) {
    /* Takes in a string, and splits it based on delimiters */
    std::vector<std::string> word_list;
    std::string curr_word = "";

    for (char ch: input) {
        if (ch == delim && curr_word.length()) {word_list.push_back(curr_word); curr_word = "";}
        else if (ch != delim) {curr_word += ch;}
    }
    if (curr_word.length()) {word_list.push_back(curr_word);}
    return word_list;
}

int is_file(std::string path) {
    /* Checks if the path is a valid file-path */
    std::ifstream test_file(path.data());
    if (test_file.fail()) {return 0;}
    
    test_file.close();
    return 1;
}

int is_dir(std::string path) {
    /* Checks if the directory is a valid path */
    return std::filesystem::exists(path);
}

void print_build_status_info(PFPDocBuildOptions* opts) {
    /* prints out the information being used in the current run */
    std::fprintf(stderr, "\nOverview of Parameters:\n");
    std::fprintf(stderr, "\tInput file-list: %s\n", opts->input_list.data());
    std::fprintf(stderr, "\tOutput ref path: %s\n", opts->output_ref.data());
    std::fprintf(stderr, "\tPFP window size: %d\n", opts->pfp_w);
    std::fprintf(stderr, "\tInclude rev-comp?: %d\n", opts->use_rcomp);
    if (opts->use_two_pass)
        std::fprintf(stderr, "\tBuild Algorithm: Two-Pass\n\n");
    else {
        std::fprintf(stderr, "\tBuild Algorithm: One-Pass\n");
        std::fprintf(stderr, "\tUse heuristics?: %d\n\n", opts->use_heuristics);
    }
}

void parse_build_options(int argc, char** argv, PFPDocBuildOptions* opts) {
    /* parses the arguments for the build sub-command, and returns a struct with arguments */

    static struct option long_options[] = {
        {"help",      no_argument, NULL,  'h'},
        {"filelist",   required_argument, NULL,  'f'},
        {"output",       required_argument, NULL,  'o'},
        {"revcomp",   no_argument, NULL,  'r'},
        {"taxcomp",   no_argument, NULL,  't'},
        {"num-col",   required_argument, NULL,  'k'},
        {"top-k",   no_argument, NULL,  'p'},
        {"print-doc", required_argument, NULL, 'e'},
        {"no-heuristic", no_argument, NULL, 'n'},
        {"modulus", required_argument, NULL, 'm'},
        {"two-pass", required_argument, NULL, 'a'},
        {"tmp-size", required_argument, NULL, 's'},
        {0, 0, 0,  0}
    };

    int c = 0;
    int long_index = 0;
    while ((c = getopt_long(argc, argv, "hf:o:w:rtk:pe:nm:a:s:", long_options, &long_index)) >= 0) {
        switch(c) {
            case 'h': pfpdoc_build_usage(); std::exit(1);
            case 'f': opts->input_list.assign(optarg); break;
            case 'o': opts->output_prefix.assign(optarg); break;
            case 'w': opts->pfp_w = std::atoi(optarg); break;
            case 'r': opts->use_rcomp = true; break;
            case 't': opts->use_taxcomp = true; break;
            case 'p': opts->use_topk = true; break;
            case 'k': opts->numcolsintable = std::max(std::atoi(optarg), 2); break;
            case 'e': opts->doc_to_extract = std::atoi(optarg); break;
            case 'n': opts->use_heuristics = false; break;
            case 'm': opts->hash_mod = std::atoi(optarg); break;
            case 'a': opts->use_two_pass = true; opts->temp_prefix.assign(optarg); break;
            case 's': opts->tmp_size_str.assign(optarg); break;
            default: pfpdoc_build_usage(); std::exit(1);
        }
    }
}

void parse_run_options(int argc, char** argv, PFPDocRunOptions* opts) {
    /* parses the arguments for the run sub-command, and returns struct */
    static struct option long_options[] = {
        {"help",      no_argument, NULL,  'h'},
        {"ref",   required_argument, NULL,  'r'},
        {"pattern",       required_argument, NULL,  'p'},
        {"output",       required_argument, NULL,  'o'},
        {"write",   no_argument, NULL,  's'},
        {"length",   required_argument, NULL,  'l'},
        {"taxcomp",   no_argument, NULL,  't'},
        {"top-k",   no_argument, NULL,  'k'},
        {"num-col",   required_argument, NULL,  'c'},
        {0, 0, 0,  0}
    };

    int c = 0;
    int long_index = 0;
    while ((c = getopt_long(argc, argv, "hr:p:o:sl:tkc:", long_options, &long_index)) >= 0) {
        switch(c) {
            case 'h': pfpdoc_run_usage(); std::exit(1);
            case 'r': opts->ref_file.assign(optarg); break;
            case 'p': opts->pattern_file.assign(optarg); break;
            case 'o': opts->output_prefix.assign(optarg); break;
            case 'l': opts->read_length = std::atoi(optarg); break;
            case 's': opts->write_to_file = true; break;
            case 't': opts->use_taxcomp = true; break;
            case 'k': opts->use_topk = true; break;
            case 'c': opts->num_cols = std::atoi(optarg); break;
            default: pfpdoc_run_usage(); std::exit(1);
        }
    }
}

void parse_info_options(int argc, char** argv, PFPDocInfoOptions* opts) {
    /* parses the arguments for the info sub-command, and returns struct */
    static struct option long_options[] = {
        {"help",      no_argument, NULL,  'h'},
        {"ref",   required_argument, NULL,  'r'},
        {"output",       required_argument, NULL,  'o'},
        {"num",       required_argument, NULL,  'n'},
        {"taxcomp",       no_argument, NULL,  't'},
        {"top-k",       no_argument, NULL,  'k'},
        {"num-col",       required_argument, NULL,  'c'},
        {0, 0, 0,  0}
    };

    int c = 0;
    int long_index = 0;
    while ((c = getopt_long(argc, argv, "hr:o:n:tkc:", long_options, &long_index)) >= 0) {
        switch(c) {
            case 'h': pfpdoc_run_usage(); std::exit(1);
            case 'r': opts->ref_file.assign(optarg); break;
            case 'o': opts->output_path.assign(optarg); break;
            case 'n': opts->num_profiles = std::max(0, std::atoi(optarg)); break;
            case 't': opts->use_taxcomp = true; break;
            case 'k': opts->use_topk = true; break;
            case 'c': opts->num_cols = std::atoi(optarg); break;
            default: pfpdoc_info_usage(); std::exit(1);
        }
    }
}

int pfpdoc_build_usage() {
    /* prints out the usage information for the build method */
    std::fprintf(stderr, "\npfp_doc build - builds the document array profiles using PFP.\n");
    std::fprintf(stderr, "Usage: pfp_doc build [options]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-28sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-18s%-10spath to a file-list of genomes to use\n", "-f, --filelist", "[FILE]");
    std::fprintf(stderr, "\t%-18s%-10soutput prefix path if using -f option\n", "-o, --output", "[PREFIX]");
    std::fprintf(stderr, "\t%-28sinclude the reverse-complement of sequence (default: false)\n\n", "-r, --revcomp");

    std::fprintf(stderr, "\t%-18s%-10suse the 2-pass construction algorithm (default: false)\n", "-a, --two-pass", "[PREFIX]");
    std::fprintf(stderr, "\t%-18s%-10ssize of temporary storage files in GB (e.g. 4GB)\n\n", "-s, --tmp-size", "[ARG]");
    
    std::fprintf(stderr, "\t%-28suse taxonomic compression of the document array (default: false)\n", "-t, --taxcomp");
    std::fprintf(stderr, "\t%-28suse top-k compression of the document array (default: false)\n", "-p, --top-k");
    std::fprintf(stderr, "\t%-18s%-10snumber of columns to include in the main table (default: 7)\n\n", "-k, --num-col", "[INT]");
    
    std::fprintf(stderr, "\t%-18s%-10sdocument number whose profiles to extract\n", "-e, --print-doc", "[INT]");
    std::fprintf(stderr, "\t%-28sturn off any heuristics used to build profiles\n\n", "-n, --no-heuristic");

    std::fprintf(stderr, "\t%-18s%-10swindow size used for pfp (default: 10)\n", "-w, --window", "[INT]");
    std::fprintf(stderr, "\t%-18s%-10shash-modulus used for pfp (default: 100)\n\n", "-m, --modulus", "[INT]");

    return 0;
}

int pfpdoc_run_usage() {
    /* prints out the usage information for the run method */
    std::fprintf(stderr, "\npfp_doc run - processes a set of reads using the document array profiles.\n");
    std::fprintf(stderr, "Usage: pfp_doc run [options]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-28sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-18s%-10spath to the input reference file\n", "-r, --ref", "[arg]");
    std::fprintf(stderr, "\t%-18s%-10spath to the pattern file\n", "-p, --pattern", "[arg]");
    std::fprintf(stderr, "\t%-18s%-10soutput path prefix\n\n", "-o, --output", "[arg]");

    std::fprintf(stderr, "\t%-28swrite data-structures to disk\n", "-s, --write");
    std::fprintf(stderr, "\t%-18s%-10supper-bound on read length (used to shrink size of index using -s)\n\n", "-l, --length", "[arg]");

    std::fprintf(stderr, "\t%-28sload document array that is 'taxonomic' compressed (default: false)\n", "-t, --taxcomp");
    std::fprintf(stderr, "\t%-28sload document array that is 'top-k' compressed (default: false)\n", "-k, --top-k");
    std::fprintf(stderr, "\t%-18s%-10snumber of columns used in main table\n\n", "-c, --num-col", "[arg]");
    
    return 0;
}

int pfpdoc_info_usage() {
    /* prints out the usage information for the info sub-command */
    std::fprintf(stderr, "\npfp_doc info - read in document array profiles and print information.\n");
    std::fprintf(stderr, "Usage: pfp_doc info [options]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-28sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-18s%-10soutput prefix used for index\n", "-r, --ref", "[arg]");
    std::fprintf(stderr, "\t%-18s%-10soutput file path for printing document array profiles (csv format)\n", "-o, --output", "[arg]");
    std::fprintf(stderr, "\t%-18s%-10snumber of profiles to print if using -o\n\n", "-n, --num", "[arg]");

    std::fprintf(stderr, "\t%-18s%-10sdocument array is compressed using taxonomy\n", "-t, --taxcomp", "[arg]");
    std::fprintf(stderr, "\t%-18s%-10sdocument array is compressed using top-k\n", "-k, --top-k", "[arg]");
    std::fprintf(stderr, "\t%-18s%-10snumber of columns in compressed table\n\n", "-c, --num-col", "[arg]");
    
    return 0;
}

int pfpdoc_usage() {
    /* Prints the usage information for pfp_doc */
    std::fprintf(stderr, "\npfp_doc has different sub-commands to run:\n");
    std::fprintf(stderr, "Usage: pfp_doc <sub-command> [options]\n\n");

    std::fprintf(stderr, "Commands:\n");
    std::fprintf(stderr, "\tbuild\tbuilds the document profile data-structure\n");
    std::fprintf(stderr, "\trun\truns queries with respect to the document array structure\n");
    std::fprintf(stderr, "\tinfo\tprint out information regarding this index and document array\n\n");
    return 0;
}

int main(int argc, char** argv) {
    /* main method for pfp_doc */
    std::fprintf(stderr, "\033[1m\033[31m\npfp-doc version: %s\033[m\033[0m\n", PFPDOC_VERSION);

    if (argc > 1) {
        if (std::strcmp(argv[1], "build") == 0) 
            return build_main(argc-1, argv+1);
        else if (std::strcmp(argv[1], "run") == 0)
            return run_main(argc-1, argv+1);
        else if (std::strcmp(argv[1], "info") == 0)
            return info_main(argc-1, argv+1);
    }
    return pfpdoc_usage();
}