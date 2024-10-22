/*
 * File: minimizer_digest.hpp
 * Description: Definition of MinimizerDigest object that can
 *              construct a concatenation of minimizers
 *              which are called a digest.
 *
 *              This code is based on MinimizerScanner 
 *              object defined in the Kraken2 classification
 *              system: https://github.com/DerrickWood/kraken2/blob/master/src/mmscanner.h
 *
 * Date: January 20th, 2024
 */

#ifndef _MINIMIZER_DIGEST_H
#define _MINIMIZER_DIGEST_H

#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include <pfp_doc.hpp>

#define BITS_PER_CHAR 2
#define KMER_TO_UINT8_CHAR(x) (char) ((x <= 2) ? (x+3) : x)

struct MinimizerData {
    uint64_t pos;
    uint64_t val;
    uint64_t hash_val;
};

class MinimizerDigest {
    public:
        MinimizerDigest();
        MinimizerDigest(uint64_t k, uint64_t w, bool lex_order=true, bool minimizer_alp=false);

        std::string compute_digest(std::string input_seq);
        uint64_t get_k() {return k;}
        uint64_t get_w() {return w;}
        void set_windows(uint64_t k, uint64_t w);
        void set_lexorder(bool status) {this->lex_order = status;}
        void set_minimizer_alp(bool status) {this->minimizer_alp = status;}
    private:
        uint64_t k = 0;
        uint64_t w = 0;
        uint8_t lookup_table[UINT8_MAX+1];
        uint64_t loaded_kmers;
        bool lex_order = true;
        bool minimizer_alp = false;
        std::vector<MinimizerData> queue;

        void update_lookup_table(char ch, uint8_t val);
        void inititalize_lookup_table();
};

#endif /* end of include guard: _MINIMIZER_DIGEST_H */
