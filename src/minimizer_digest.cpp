/*
 * File: minimizer_digest.cpp
 * Description: Implementation of MinimizerDigest object
 *              that constructs the minimizer digest.
 *
 *              This code is based on MinimizerScanner 
 *              object defined in the Kraken2 classification
 *              system: https://github.com/DerrickWood/kraken2/blob/master/src/mmscanner.cc
 *
 * Date: January 20th, 2024
 */

#include <minimizer_digest.hpp>
#include <pfp_doc.hpp>
#include <hash_func.hpp>
#include <string.h>

MinimizerDigest::MinimizerDigest(uint64_t k, uint64_t w, bool lex_order, bool minimizer_alp): 
        k(k), 
        w(w), 
        loaded_kmers(0), 
        lex_order(lex_order), 
        minimizer_alp(minimizer_alp) 
{
    // make sure the small window size is not too large
    uint64_t limit = ((sizeof(uint64_t) * 8 - 1)/BITS_PER_CHAR);
    if (k > limit) {FATAL_ERROR("the small-window size provided is too large.");}

    // assert that small window is 4 if using minimizer alphabet
    if (minimizer_alp && k != 4) {FATAL_ERROR("in order to use minimizer alphabet, k must be 4.");}

    // make sure the large window size is larger
    ASSERT((w >= k), "the large-window size cannot be smaller than small-window size.");

    // initialize the lookup table 
    for (size_t i = 0; i < UINT8_MAX; i++){
        lookup_table[i] = UINT8_MAX;
    }
    update_lookup_table('A', 0x00);
    update_lookup_table('C', 0x01);
    update_lookup_table('G', 0x02);
    update_lookup_table('T', 0x03);
}

void MinimizerDigest::update_lookup_table(char ch, uint8_t val) {
    lookup_table[ch] = val;
    lookup_table[std::tolower(ch)] = val;
}

std::string MinimizerDigest::compute_digest(std::string input_seq) {
    // return full sequnece if it less than large window
    if (input_seq.size() < w)
        return input_seq;
    
    // define key variables
    uint64_t curr_kmer_val = 0x00, str_pos = 0, loaded_chs = 0;
    uint64_t last_n_bits = ((uint64_t) 1 << (2*k))-1;
    uint64_t last_minimizer = UINT64_MAX;
    std::string digest = "";

    // lambda to load up the first kmer of string, or first kmer after non-ACGT char
    auto load_first_kmer_to_queue = [&] () {
        // perform operations to reset variable
        queue.clear();
        loaded_chs = 0; 
        curr_kmer_val = 0x00; 
        loaded_kmers = 0;

        // keep going until we get our first k-mer
        while (loaded_chs < k && str_pos < input_seq.size()) {
            if (lookup_table[input_seq[str_pos]] == UINT8_MAX) {
                queue.clear();
                loaded_chs = 0; loaded_kmers = 0;
                curr_kmer_val = 0x00;
                str_pos++;
                continue;
            }
        
            curr_kmer_val <<= BITS_PER_CHAR;
            curr_kmer_val |= lookup_table[input_seq[str_pos++]];
            loaded_chs++;
        }
        
        // add first k-mer to the queue
        uint64_t hash_val = (lex_order) ? (curr_kmer_val) : MurmurHash3(curr_kmer_val); 
        MinimizerData curr_kmer_data {(str_pos-k), curr_kmer_val, hash_val};
        queue.push_back(curr_kmer_data);
        loaded_kmers++;
    };

    // lambda to add current k-mer to queue
    auto append_kmer_to_queue = [&] (uint64_t new_pos, uint64_t new_val) {
        // compute hash value based on whether using lex order or not
        uint64_t hash_val = (lex_order) ? new_val : MurmurHash3(new_val);

        // go from the back until we see a value less then the new one
        int it = queue.size()-1;
        for (; it >= 0; it--) {
            if (queue[it].hash_val < hash_val) {break;}
        }
        if (queue.size() > 0){
            queue.erase(queue.begin()+it+1, queue.end());
        }

        // append new kmer data
        MinimizerData new_kmer {new_pos, new_val, hash_val};
        queue.push_back(new_kmer);
        loaded_kmers++;
    };

    // load the first kmer, and then continue char by char
    load_first_kmer_to_queue();

    // first, we check if have enough k-mers to report a minimizer,
    // this would happen if k = w
    if (loaded_kmers >= (w-k+1) && queue[0].val != last_minimizer) { 
        DEBUG_MSG("MINIMIZER FOUND: pos = " << queue[0].pos << 
                                 ", val = " << queue[0].val << 
                                 ", hash_val = " << queue[0].hash_val << 
                                 ", string = " << input_seq.substr(queue[0].pos,k));   

        if (minimizer_alp) {
            ASSERT((queue[0].val <= 255), "unexpected value seen for k-mer");
            uint8_t min_char = KMER_TO_UINT8_CHAR(queue[0].val);
            digest += min_char;
        } else {
            digest += input_seq.substr(queue[0].pos, k);
        }
        //digest += (!minimizer_alp) ? input_seq.substr(queue[0].pos, k) : KMER_TO_UINT8_CHAR(queue[0].val);
        last_minimizer = queue[0].val;
    }
    // remove first element if out of range
    if (str_pos >= w && queue[0].pos < (str_pos-w+1)){
        queue.erase(queue.begin());
        loaded_kmers--;
    }

    while (str_pos <= (input_seq.size()-1)) {
        // for non-ACGT chars, reset large window
        if (lookup_table[input_seq[str_pos]] == UINT8_MAX) {
            if ((input_seq.size() - (++str_pos)) >= w)
                load_first_kmer_to_queue();
            else
                return digest;
        // for ACGT chars, add kmer to window
        } else {
            curr_kmer_val <<= BITS_PER_CHAR;
            curr_kmer_val = (curr_kmer_val | lookup_table[input_seq[str_pos++]]) & last_n_bits;
            append_kmer_to_queue((str_pos-k), curr_kmer_val);
        }

        // report minimizer if we have a full range
        if (loaded_kmers >= (w-k+1) && queue[0].val != last_minimizer) { 
            DEBUG_MSG("MINIMIZER FOUND: pos = " << queue[0].pos << 
                                     ", val = " << queue[0].val << 
                                     ", hash_val = " << queue[0].hash_val << 
                                     ", string = " << input_seq.substr(queue[0].pos,k));
            
            if (minimizer_alp) {
                ASSERT((queue[0].val <= 255), "unexpected value seen for k-mer");
                uint8_t min_char = KMER_TO_UINT8_CHAR(queue[0].val);
                digest += min_char;
            } else {
                digest += input_seq.substr(queue[0].pos, k);
            }
            //digest += (!minimizer_alp) ? input_seq.substr(queue[0].pos, k) : KMER_TO_UINT8_CHAR(queue[0].val);
            last_minimizer = queue[0].val;
        }

        // remove first element if out of range
        if (str_pos >= w && queue[0].pos < (str_pos-w+1)) {
            queue.erase(queue.begin());
            loaded_kmers--;
        }
    }
    return digest;
}