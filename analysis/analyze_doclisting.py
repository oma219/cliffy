#!/usr/bin/env python3

# Name: analyze_doclisting.py
# Description: script that will take the output from Cliffy when
#              when performing document listing (not taxonomic mode)
#              and report the top hits for each read.

import argparse
import os
import sys
import subprocess
import copy

########################################################
# extract sub-command #1: extract the data for each
# example and prepare the data in order process it.
########################################################
 
def extract_main(args):
    # make sure both files are present
    check_if_file_exists(args.input_fasta)
    check_if_file_exists(args.input_data)
    
    # check if output directory is present
    check_if_dir_exists(args.output_dir)
    if args.output_dir[-1] != '/':
        args.output_dir += '/'
    
    # extract the *.tar.gz file to the specified output directory
    try:
        log_message(f"starting to extract: {args.input_data}")
        subprocess.run(f"tar -xvf {args.input_data} -C {args.output_dir}", shell=True, check=True)
        log_message(f"successfully extracted {args.input_data}\n")
        
        log_message(f"starting to extract: {args.input_fasta}")
        subprocess.run(f"tar -xvf {args.input_fasta} -C {args.output_dir}", shell=True, check=True)
        log_message(f"successfully extracted {args.input_fasta}\n")
    except subprocess.CalledProcessError:
        error_message(f"issue occurred when extracting {args.input_data} or {args.input_fasta}")

    # iterate through the filelist and make sure each document is present
    orig_filelist = args.output_dir + os.path.basename(args.input_data).split(".tar.gz")[0] + "/orig_filelist.txt"
    check_if_file_exists(orig_filelist)
    
    fasta_dir = args.output_dir + os.path.basename(args.input_fasta).split(".tar.gz")[0]
    avail_fasta_files = [path for path in os.listdir(fasta_dir) if os.path.isfile(os.path.join(fasta_dir, path))]
    
    output_filelist = args.output_dir + os.path.basename(args.input_fasta).split(".tar.gz")[0] + "/input_filelist.txt"
    num_docs = 0
    
    log_message(f"writing new filelist to: {output_filelist}")
    with open(orig_filelist, "r") as in_fd, open(output_filelist, "w") as out_fd:
        for line in in_fd:
            line_split = line.split()
            assert len(line_split) == 2 
            assert line_split[0] in avail_fasta_files
            
            if int(line_split[1]) > num_docs:
                num_docs = int(line_split[1])
            elif int(line_split) < num_docs:
                error_message("unexpected order of document numbers.")
            
            out_fd.write(f"{os.path.join(fasta_dir, line_split[0])} {line_split[1]}\n")
    
    log_message(f"finished writing filelist: found {num_docs} documents.")        

########################################################
# analyze sub-command #1: take the output from cliffy
# and output a file that summarizes the results and
# makes it more human readable.
########################################################

def analyze_main(args):
    """ process the listings file and report what document has the highest weight """
    # check if files and output directory is valid
    paired_end = False
    check_if_file_exists(args.mate1_listings)
    
    if args.mate2_listings != "":
        check_if_file_exists(args.mate2_listings)
        paired_end = True
    
    check_if_file_exists(args.doc_labels_path)
    check_if_dir_exists(os.path.dirname(args.output_path))
    log_message("all files have been verified")
    
    # load the document labels
    doc_to_label = {}
    with open(args.doc_labels_path, "r") as in_fd:
        line_i = 0
        for line in in_fd:
            doc_to_label[line_i] = line.strip()
            line_i += 1  
    num_docs = len(doc_to_label)
    log_message("loaded a dictionary mapping document id to label")
    
    # load the reads and classify them
    if paired_end == True:
        read_list = load_paired_end_reads(args.mate1_listings, args.mate2_listings)
    else:
        read_list = load_single_end_reads(args.mate1_listings)
    log_message(f"loaded {len(read_list)} reads for classification\n")
    
    # process each read and find classification
    with open(args.output_path, "w") as out_fd:
        for curr_read in read_list:
            curr_name = curr_read[0]
            weights = [0.0 for x in range(num_docs)]
            
            # add weight for each mate
            for listings in curr_read[1:]:
                weights = add_mate_listings_to_weight_vector(listings, weights)
            
            # identify the top n taxa
            sorted_weights = sorted(list(enumerate(weights)), key=lambda x: x[1], reverse=True)
            top_n = 3
            
            out_fd.write(f"{curr_name}")
            for i in range(min(top_n, len(sorted_weights))):
                curr_doc_id = sorted_weights[i][0]
                curr_weight = sorted_weights[i][1]
                
                out_fd.write(f"\t{doc_to_label[curr_doc_id]}\t{curr_weight}")
            out_fd.write("\n")
    log_message("finished writing classifications to output file.")
        
def add_mate_listings_to_weight_vector(mate_listings, weight_dict):
    mate_listings_list = mate_listings.split()
    assert len(mate_listings_list) % 2 == 0

    for i in range(0, len(mate_listings_list), 2):
        # get the length of exact match
        start_and_stop = mate_listings_list[i][1:-1].split(",")
        stop = int(start_and_stop[1])
        start = int(start_and_stop[0])
        length = stop - start + 1
        assert length > 0

        # get the traversal of documents that are hit
        doc_listing = [int(x) for x in mate_listings_list[i+1][1:-1].split(",")]
        for doc_id in doc_listing:
            assert doc_id < len(weight_dict), f"doc_id = {doc_id}"
            weight_dict[doc_id] += length
    return weight_dict
    
def load_paired_end_reads(mate1_path, mate2_path):
    """ load the results into a list """
    with open(mate1_path, "r") as mate1_fd, open(mate2_path, "r") as mate2_fd:
        mate1_lines = [x.strip() for x in  mate1_fd.readlines()]
        mate2_lines = [x.strip() for x in  mate2_fd.readlines()]
    
    header_m1 = ""; listing_m1 = ""; 
    header_m2 = ""; listing_m2 = ""; pos = 0
    final_read_list = []
    
    for mate1_line, mate2_line in zip(mate1_lines, mate2_lines):
        if mate1_line.startswith(">") and mate2_line.startswith(">"):
            header_m1 = mate1_line[1:]; header_m2 = mate2_line[1:]
            pos += 1
        elif pos == 1:
            listing_m1 = mate1_line; listing_m2 = mate2_line 
            pos = 0

            # Save reads to list ...               
            final_read_list.append((header_m1, listing_m1, listing_m2))
    return final_read_list    

def load_single_end_reads(mate1_path):
    """ load the results into Read objects """
    with open(mate1_path, "r") as mate1_fd:
        mate1_lines = [x.strip() for x in  mate1_fd.readlines()]
    
    header_m1 = ""; listing_m1 = ""; pos = 0
    final_read_list = []
    
    for mate1_line in mate1_lines:
        if mate1_line.startswith(">"):
            header_m1 = mate1_line[1:]
            pos += 1
        elif pos == 1:
            listing_m1 = mate1_line
            pos = 0

            # Save read object ...
            read_group_name = int(header_m1.split("_")[1])                
            final_read_list.append((header_m1, listing_m1))
    return final_read_list
    
########################################################
# helper method: argument parsing, file checking, etc.
########################################################

def parse_arguments():
    main_parser = argparse.ArgumentParser(description="analyze the output files from Cliffy classifier when doing document listing.")
    sub_parser = main_parser.add_subparsers(dest="command", 
                                            help="available sub-commands",
                                            required=True)
    
    # sub-command 1: extract data from compressed files
    extract_parser = sub_parser.add_parser("extract", help="prepare the data in order to build index.")
    extract_parser.add_argument("--input-text", 
                                dest="input_fasta", 
                                type=str, 
                                help="path to compressed document files. (e.g. FASTA files)",
                                required=True)
    extract_parser.add_argument("--input-data", 
                                dest="input_data", 
                                type=str, 
                                help="path to compressed metadata files",
                                required=True)
    extract_parser.add_argument("--output-dir",
                                dest="output_dir",
                                help="path to output directory for FASTA files",
                                type=str,
                                required=True)
    
    # sub-command 2: analyze cliffy output results
    analyze_parser = sub_parser.add_parser("analyze", help="analyze the results from Cliffy and summarize it.")
    analyze_parser.add_argument("--mate-1", 
                                dest="mate1_listings", 
                                type=str, 
                                help="path to mate1 listings file (or listing file if not paired-end)",
                                required=True)
    analyze_parser.add_argument("--mate-2", 
                                dest="mate2_listings", 
                                type=str, default="",
                                help="path to mate2 listings file (leave empty if not paired-end)",
                                required=False)
    analyze_parser.add_argument("--doc-labels",
                                dest="doc_labels_path",
                                type=str, default="",
                                help="path to file describing what group each document corresponds to.",
                                required=True)
    analyze_parser.add_argument("--output-file",
                                dest="output_path",
                                type=str, default="",
                                help="path to file to write out output.",
                                required=True)

    args = main_parser.parse_args()
    return args

def check_if_file_exists(path):
    if not os.path.exists(path):    
        error_message(f"file {path} doesn't exist.")

def check_if_dir_exists(dir_path):
    if not os.path.exists(dir_path):
        error_message(f"directory {dir_path} does not exist, please make it.")

def error_message(msg):
    red_start = "\033[91m"; red_end = "\033[0m"
    print(f"\n{red_start}[Error]{red_end} {msg}\n")
    exit(1)

def log_message(msg):
    print(f"\033[92m[cliffy::log]\033[0m {msg}")

if __name__ == "__main__":
    args = parse_arguments()

    if args == None or args.command not in ["extract", "analyze"]:
        error_message("need to specify a sub-command.")
    elif args.command == "extract":
        print(); extract_main(args)
    elif args.command == "analyze": 
        print(); analyze_main(args)
    print()