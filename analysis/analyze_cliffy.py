#!/usr/bin/env python3

# Name: analyze_cliffy.py
# Description: script that will take the output from Cliffy when
#              when performing taxonomic classification with respect
#              to the SILVA database, it is capable of reporting
#              classifications at read-level or genera.

import argparse
import sys
import os
import re
import subprocess
import time
import multiprocess as mp
import copy
from collections import Counter

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
        log_message(f"successfully extracted {args.input_data}")
        
        log_message(f"starting to extract: {args.input_fasta}")
        subprocess.run(f"tar -xvf {args.input_fasta} -C {args.output_dir}", shell=True, check=True)
        log_message(f"successfully extracted {args.input_fasta}")
    except subprocess.CalledProcessError:
        error_message(f"issue occurred when extracting {args.input_data} or {args.input_fasta}")

    # verify the example data is as expected and find number of documents
    fasta_folder_name = args.output_dir + os.path.basename(args.input_fasta).split(".tar.gz")[0]
    pattern = r"doc_[0-9]+_seq\.fa$"
    max_number = 0; num_docs = 0;
    
    for file_name in os.listdir(fasta_folder_name):
        if re.match(pattern, file_name):
            num_docs += 1
            max_number = max(max_number, int(file_name.split("_")[1]))
        else:
            error_message(f"found file that is unexpected: {file_name}")
            
    if num_docs != max_number:
        error_message(f"issue occurred document numbers are not expected.")
    
    print(); log_message(f"found {max_number} documents/different groups in the example.")
    
    # write a filelist file for cliffy
    metadata_folder = args.output_dir + os.path.basename(args.input_data).split(".tar.gz")[0] + "/"
    filelist_path = metadata_folder + "input_filelist.txt"
    
    with open(filelist_path, "w") as out_fd:
        for i in range(1, max_number+1):
            out_fd.write(f"{fasta_folder_name}/doc_{i}_seq.fa {i}\n")
    log_message(f"finished writing filelist: {filelist_path}")

########################################################
# analyze sub-command #1: take the output from cliffy
# and output a file that summarizes the results and
# makes it more human readable.
########################################################

def analyze_main(args):
    # check the files and output directories
    paired_end = False
    if args.mate1_listings == "":
        error_message("need to specify a output listing files.")
    check_if_file_exists(args.mate1_listings)
    
    if args.mate2_listings != "":
        check_if_file_exists(args.mate2_listings)
        paired_end = True
        log_message("found paired-end results to process.")
    else:
        log_message("found single-end results to process.")
    
    read_level = True
    if not args.read_level and not args.abundance_profile:
        error_message("need to specify either --read-level or --abundance-profile.")
    elif args.read_level and args.abundance_profile:
        error_message("cannot specify both --read-level and --abundance-profile")
    else:
        if args.abundance_profile:
            read_level = False
            
    check_if_dir_exists(os.path.dirname(args.output_path))
    check_if_file_exists(args.doc_id_to_taxa_path)
    check_if_file_exists(args.silva_taxa_ranks_path)
    log_message("all files have been verified ... loading reads now\n")
    
    # create a readset object
    if paired_end == True:
        cliffy_rs = CliffyReadSet(args.mate1_listings,
                                  args.mate2_listings,
                                  args.doc_id_to_taxa_path,
                                  args.silva_taxa_ranks_path,
                                  True)
    else:
        cliffy_rs = CliffyReadSet(args.mate1_listings,
                                  "",
                                  args.doc_id_to_taxa_path,
                                  args.silva_taxa_ranks_path,
                                  False)
    
    # classify the reads in the readset object
    if read_level:
        log_message("starting to classify reads at the read-level.")
        cliffy_rs.classify_at_read_level(args.output_path)
    else:
        log_message("starting to compute the abundance profile.")
        cliffy_rs.generate_abundance(args.output_path)

class SILVA_taxonomy:
    def __init__(self,
                 silva_taxonomy_file):
        """
        Input Fields:
            silva_taxonomy_file: SILVA file that maps every taxonomic
                                 traversal to node id
        """
        self.raw_data = []
        self.initialize_raw_data(silva_taxonomy_file)

        self.trav_to_rank = {}
        for _, trav, rank in self.raw_data:
            self.trav_to_rank[trav] = rank
        
        self.id_to_trav = {}
        for node_id, trav, _ in self.raw_data:
            self.id_to_trav[int(node_id)] = trav

    def initialize_raw_data(self, 
                            input_path):
        """
        Initializes the raw_data attribute with
        tuples containing the node_id, traversal, and
        taxonomic rank
        """
        with open(input_path, "r") as in_fd:
            for line in in_fd:
                last_pos_of_traversal, traversal = SILVA_taxonomy.get_traversal_from_full_line(line)
                
                line_split = line.split()
                node_id = int(line_split[last_pos_of_traversal+1])
                rank = line_split[last_pos_of_traversal+2]

                self.raw_data.append([node_id, traversal, rank])

    def get_traversal_from_full_line(line):
        """
        Get the traversal part of line, the presence of spaces in
        in traversal string makes it slightly less than trivial

        line: string with entire lines contents from silva taxonomy file
        """
        line_split = line.split()
        last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
        traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
        
        assert (len(line_split) - last_pos_of_traversal) <= 4   
        return last_pos_of_traversal, traversal

    def get_rank_from_traversal(self,
                                traversal):
        """ 
        Returns taxonomic rank of a traversal path if
        it exists, otherwise None
        """
        if traversal not in self.trav_to_rank:
            return None 
        else:
            return self.trav_to_rank[traversal]

    def check_traversal_presence(self, 
                                 traversal):
        """
        Returns True if traversal is present in SILVA
        file, and False otherwise
        """
        if traversal in self.trav_to_rank:
            return True 
        else:
            return False

    def get_silva_dict_for_level_to_list(self):
        """
        Creates a dictionary thats takes a taxonomic rank like
        'domain', 'phylum' and returns a list of the names in the
        SILVA taxonomy. 
        """
        tree_levels = ["domain", "phylum", "class", "order", "family", "genus"]
        silva_ranks_dict = {}

        for level in tree_levels:
            silva_level_names = []
            for tup in self.raw_data:
                if tup[2] == level:
                    name = tup[1].split(";")[-2]
                    silva_level_names.append(TaxonomyTraversal(tup[1], self))
            silva_ranks_dict[level] = silva_level_names
        return silva_ranks_dict

    def get_traversal_from_id(self,
                              node_id):
        """ 
        Returns traversal path based on node id if
        it exists, otherwise None
        """
        if node_id not in self.id_to_trav:
            return None 
        else:
            return self.id_to_trav[int(node_id)]

class TaxonomyTraversal:
    def __init__(self,
                traversal,
                silva_taxonomy):
        """
        Input Fields:
            traversal: a string from the doc_to_traversal file
            silva_taxonomy_file: SILVA file that maps every traversal to rank
        """
        self.trav_str = traversal
        self._domain = None
        self._phylum = None
        self._class = None
        self._order = None
        self._family = None
        self._genus = None

        self.initialize_ranks(traversal, silva_taxonomy)
    
    def initialize_ranks(self,
                         traversal,
                         silva_taxonomy):
        """ parse the traversal and fill in taxonomic ranks """
        assert silva_taxonomy.check_traversal_presence(traversal)
        curr_str = ""
        for ch in traversal:
            curr_str += ch
            if ch == ";":
                rank = silva_taxonomy.get_rank_from_traversal(curr_str)
                if rank != None:
                    if rank == "domain":
                        self._domain = curr_str
                    elif rank == "phylum":
                        self._phylum = curr_str
                    elif rank == "class":
                        self._class = curr_str
                    elif rank == "order":
                        self._order = curr_str
                    elif rank == "family":
                        self._family = curr_str
                    elif rank == "genus":
                        self._genus = curr_str

    def get_certain_level(self,
                         level_name):
        """
        Return the specific level requested if it 
        is present, if not present return None
        """
        assert level_name in ["domain", "phylum", "class", "order", "family", "genus"]
        curr_trav = ""

        if level_name == "domain":
            curr_trav = self._domain
        elif level_name == "phylum":
            curr_trav = self._phylum
        elif level_name == "class":
            curr_trav = self._class
        elif level_name == "order":
            curr_trav = self._order
        elif level_name == "family":
            curr_trav = self._family
        elif level_name == "genus":
            curr_trav = self._genus
        
        if curr_trav is None:
            return None
        else:
            return curr_trav.split(";")[-2]
       
class CliffyRead:
    def __init__(self,
                 name,
                 mate1_listings,
                 mate2_listings,
                 paired_end=True):
        """
        Input Fields:
            name: name of the read
            mate{1,2}_listings: output from PFP_DOC64
            correct_taxa: ReadsetTraversal Object that represents the correct classification for this read
        """
        self.name = name 
        self.mate1_listings = mate1_listings 
        self.mate2_listings = mate2_listings
        self.paired_end = True

    def classify_with_lca(self,
                          level,
                          doc_id_to_clade_input,
                          top_n=3):
        """ 
        Idea: Take leftmost and rightmost distribute the weight
              to all the nodes equally between those two.

        level: string, that tells you what level of tree we are classifying at
        list_of_classes: list of TaxonomyTraversal objects that are all possible classes
        doc_id_to_traversal: dictionary from node ids to TaxonomyTraversal objects

        Returns True if Correct, and False Otherwise
        """
        # create a dictionary of names to weights
        doc_id_to_weight = copy.deepcopy(doc_id_to_clade_input)

        # helper method to add to dictionary
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
                left_doc = doc_listing[0]; right_doc = doc_listing[-1]
                doc_listing_new = list(set([left_doc, right_doc]))

                # iterate through all documents between left and rightmost and add weight
                weight_to_distribute = (length)/(right_doc+1-left_doc)
                for doc_id in range(left_doc, right_doc+1):
                    assert doc_id in weight_dict
                    weight_dict[doc_id][0] += weight_to_distribute #length

            return weight_dict

        # add each mate results to list
        name_to_weight = add_mate_listings_to_weight_vector(self.mate1_listings, doc_id_to_weight)
        if self.paired_end == True:
            name_to_weight = add_mate_listings_to_weight_vector(self.mate2_listings, doc_id_to_weight)

        # get the predicted taxa, and return results
        top_taxa = []
        sorted_taxa = sorted(doc_id_to_weight.items(), key=lambda x: x[1][0], reverse=True)
        for doc_id, (weight, taxa) in sorted_taxa:
            top_taxa.append((taxa, round(weight, 3)))
            if len(top_taxa) == top_n:
                break
        return top_taxa
    
    def classify_with_approx_doc_list(self,
                                      level,
                                      doc_id_to_clade_input,
                                      top_n=3):
        """ 
        Idea: Takes all the occurrences in the document listing
              and adds the full weight to each of the nodes.

        level: string, that tells you what level of tree we are classifying at
        list_of_classes: list of TaxonomyTraversal objects that are all possible classes
        doc_id_to_traversal: dictionary from node ids to TaxonomyTraversal objects

        Returns True if Correct, and False Otherwise
        """
        # create a dictionary of names to weights
        doc_id_to_weight = copy.deepcopy(doc_id_to_clade_input)

        # helper method to add to dictionary
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
                    assert doc_id in weight_dict
                    weight_dict[doc_id][0] += length
            return weight_dict
        
        # add each mate1 results to list
        name_to_weight = add_mate_listings_to_weight_vector(self.mate1_listings, doc_id_to_weight)
        if self.paired_end == True:
            name_to_weight = add_mate_listings_to_weight_vector(self.mate2_listings, doc_id_to_weight)


        # get the predicted taxa, and return results
        max_count = max(name_to_weight.values(), key=lambda x: x[0])[0]
        max_name_taxa = [k for (k, v) in name_to_weight.items() if v[0] == max_count]
        
        # get the predicted taxa, and return results
        top_taxa = []
        sorted_taxa = sorted(doc_id_to_weight.items(), key=lambda x: x[1][0], reverse=True)
        for doc_id, (weight, taxa) in sorted_taxa:
            top_taxa.append((taxa, round(weight, 5)))
            if len(top_taxa) == top_n:
                break
        return top_taxa

class CliffyReadSet:
    def __init__(self,
                mate1_listings,
                mate2_listings,
                doc_to_traveral_file,
                silva_traversal_file,
                paired_end=True):
        """
        Input Fields:
            mate{1,2}_listings: output files from cliffy
            doc_to_traversal_file: index file from PFP_DOC64
            silva_traversal_file: taxonomy file from SILVA database
        """
        self.paths = {}
        self.initialize_paths_variable(mate1_listings, 
                                       mate2_listings, 
                                       doc_to_traveral_file, 
                                       silva_traversal_file)

        self.silva_taxonomy = SILVA_taxonomy(self.paths["silva_traversal"])
        log_message("loaded the SILVA taxonomy.")

        self.doc_id_to_full_traversal = {}
        self.initialize_pfpdoc_doc_to_traversal_dict(self.paths["doc_id_to_traversal"])
        log_message("loaded the mapping from document id to taxa.")

        self.silva_nodes_at_each_level = self.silva_taxonomy.get_silva_dict_for_level_to_list()

        self.read_list = []
        if paired_end == True:
            self.load_paired_end_reads(self.paths["mate1_listings"],
                                       self.paths["mate2_listings"])
        else:
            self.load_single_end_reads(self.paths["mate1_listings"])
        log_message(f"finished loading {len(self.read_list)} reads (or read-pairs)\n")
        
    def initialize_paths_variable(self,
                                  mate1_listings,
                                  mate2_listings,
                                  doc_to_traveral_file,
                                  silva_traversal_file):
        """ initialize a dictionary with all the file paths """
        self.paths["mate1_listings"] = mate1_listings
        self.paths["mate2_listings"] = mate2_listings
        self.paths["doc_id_to_traversal"] = doc_to_traveral_file
        self.paths["silva_traversal"] = silva_traversal_file

    def initialize_pfpdoc_doc_to_traversal_dict(self,
                                                file_path):
        """ parse the two-column file that maps doc # to genus """
        with open(file_path, "r") as in_fd:
            for line in in_fd:
                line_split = [x.strip() for x in line.split()]
                curr_traversal = " ".join(line_split[1:])
                self.doc_id_to_full_traversal[int(line_split[0])-1] = TaxonomyTraversal(curr_traversal,
                                                                                        self.silva_taxonomy)

    def initialize_trav_to_length_dict(self,
                                       file_path):
        """ parse two-column file that maps reference sequences to seq length """
        with open(file_path, "r") as in_fd:
            for line in in_fd:
                line_split = [x.strip() for x in line.split()]
                curr_traversal = " ".join(line_split[:len(line_split)-1])
                self.trav_to_length[TaxonomyTraversal(curr_traversal, self.silva_taxonomy)] = int(line_split[-1])

    def initialize_readset_truthset_dict(self,
                                         input_path):
        """ parse the two-column file into dictionary mapping read name to traversal """
        with open(input_path, "r") as in_fd:
            for line in in_fd:
                line_split = [x.strip() for x in line.split()]

                # parse out the traversal and genus num
                last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])  
                traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
                curr_genus_num = int(line_split[last_pos_of_traversal+1])

                self.read_name_to_traversal[curr_genus_num] = TaxonomyTraversal(traversal, self.silva_taxonomy)

    def load_paired_end_reads(self,
                  mate1_path,
                  mate2_path):
        """ load the results into Read objects """
        with open(mate1_path, "r") as mate1_fd, open(mate2_path, "r") as mate2_fd:
            mate1_lines = [x.strip() for x in  mate1_fd.readlines()]
            mate2_lines = [x.strip() for x in  mate2_fd.readlines()]
        
        header_m1 = ""; listing_m1 = ""; 
        header_m2 = ""; listing_m2 = ""; pos = 0

        for mate1_line, mate2_line in zip(mate1_lines, mate2_lines):
            if mate1_line.startswith(">") and mate2_line.startswith(">"):
                header_m1 = mate1_line[1:]; header_m2 = mate2_line[1:]
                pos += 1
            elif pos == 1:
                listing_m1 = mate1_line; listing_m2 = mate2_line 
                pos = 0

                # Save read object ...
                assert header_m1[:-1] == header_m2[:-1]                
                self.read_list.append(CliffyRead(header_m1, 
                                                 listing_m1, 
                                                 listing_m2,
                                                 True))

    def load_single_end_reads(self,
                              mate1_path):
        """ load the results into Read objects """
        with open(mate1_path, "r") as mate1_fd:
            mate1_lines = [x.strip() for x in  mate1_fd.readlines()]
        
        header_m1 = ""; listing_m1 = ""; pos = 0
        for mate1_line in mate1_lines:
            if mate1_line.startswith(">"):
                header_m1 = mate1_line[1:]
                pos += 1
            elif pos == 1:
                listing_m1 = mate1_line
                pos = 0

                # Save read object ...
                read_group_name = int(header_m1.split("_")[1])                
                self.read_list.append(CliffyRead(header_m1, 
                                                 listing_m1, 
                                                 "",
                                                 False))

    def classify_at_read_level(self,
                               output_path):
        """ generate an output file with stats for each read """

        def batch(iterable, batch_size=5000):
            l = len(iterable)
            for ndx in range(0, l, batch_size):
                yield iterable[ndx:min(ndx + batch_size, l)]

        full_results = []
        with open(output_path, "w") as out_fd:
            # define level of interest
            clade = "genus"
            
            # create a dict from doc_id to weight and taxa name
            doc_id_to_clade = {i: [0, x.trav_str] for i, x in self.doc_id_to_full_traversal.items()}
                        
            # create a list of tuples to pass to worker threads
            curr_read_list = [(i, 
                               curr_read, 
                               clade, 
                               doc_id_to_clade ,
                               True) for i, curr_read in enumerate(self.read_list)]

            # convert to a list of batches
            batches = list(batch(curr_read_list, 1000))

            # classify all the reads at this level
            log_message(f"starting to classify {len(curr_read_list)} reads at the genus level.")
            
            start = time.time()
            with mp.Pool() as pool:
                full_results = pool.map(classify_cliffy_read_worker, batches)
                
            log_message(f"finished classifying: {round(time.time()-start,4)} seconds")
            flattened_results = [item for batch_results in full_results for item in batch_results]
            print()
            
            # write out results to file
            log_message(f"writing results to: {output_path}")
            for read_name, results_tup in flattened_results:
                out_fd.write(read_name)
                for taxa_path, weight in results_tup:
                    out_fd.write(f"\t{taxa_path}\t{weight}")
                out_fd.write("\n")
            log_message("finished.")
              
    def generate_abundance(self,
                           output_path):
        """ generate an abundance plot at genus level """

        # helper method: split reads into batches
        def batch(iterable, batch_size=5000):
            l = len(iterable)
            for ndx in range(0, l, batch_size):
                yield iterable[ndx:min(ndx + batch_size, l)]

        # create a dict from doc_id to weight and taxa name
        doc_id_to_clade = {i: [0, x.trav_str] for i, x in self.doc_id_to_full_traversal.items()}
                    
        # create a list of tuples to pass to worker threads
        curr_read_list = [(i, 
                            curr_read, 
                            "genus", 
                            doc_id_to_clade ,
                            True) for i, curr_read in enumerate(self.read_list)]

        # convert to a list of batches
        batches = list(batch(curr_read_list, 1000))

        # classify all the reads at genus level
        with mp.Pool() as pool:
            full_results = pool.map(classify_cliffy_read_worker, batches)
        flattened_results = [item for batch_results in full_results for item in batch_results]

        # extract the genus call for each read
        classified_genera = [x[1][0][0] for x in flattened_results]
        genera_counts = {k: v for k, v in sorted(dict(Counter(classified_genera)).items(), key=lambda item: item[1], reverse=True)}

        # write to file
        log_message(f"writing results to: {output_path}")
        with open(output_path, "w") as out_fd:
            for taxa, count in genera_counts.items():
                out_fd.write(f"{taxa}\t{count}\n")
        log_message("finished writing.")
        
def classify_cliffy_read_worker(input_list):
    """ function that is called in parallel to classify all the reads """
    results = []
    for input_tuple in input_list:
        i, curr_read_obj, clade, doc_id_to_clade, use_approx_doc = input_tuple
        if use_approx_doc == True:
            result_tup = curr_read_obj.classify_with_approx_doc_list(clade, doc_id_to_clade)
        else:
            result_tup = curr_read_obj.classify_with_lca(clade, doc_id_to_clade)
        results.append((curr_read_obj.name, result_tup))
    return results

########################################################
# helper method: argument parsing, file checking, etc.
########################################################

def parse_arguments():
    main_parser = argparse.ArgumentParser(description="analyze the output files from Cliffy classifier.")
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
    analyze_parser.add_argument("--doc-id-to-taxa",
                                dest="doc_id_to_taxa_path",
                                type=str, default="",
                                help="path to file describing what group each document corresponds to.",
                                required=True)
    analyze_parser.add_argument("--silva-taxa-names",
                                dest="silva_taxa_ranks_path",
                                type=str, default="",
                                help="path to file mapping SILVA taxas to ranks",
                                required=True)
    analyze_parser.add_argument("--output-file",
                                dest="output_path",
                                type=str, default="",
                                help="path to file to write out output.",
                                required=True)
    analyze_parser.add_argument("--read-level",
                                dest="read_level",
                                action="store_true",
                                default=False,
                                help="report classifications at the read-level (default: false)")
    analyze_parser.add_argument("--abundance-profile",
                                dest="abundance_profile",
                                action="store_true",
                                default=False,
                                help="report classifications as abundances at genus level (default: false)")

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