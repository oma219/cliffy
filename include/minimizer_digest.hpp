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

#define BITS_PER_CHAR 2

struct MinimizerData {
    uint64_t pos;
    uint64_t val;
};

class MinimizerDigest {
    public:
        MinimizerDigest(uint64_t k, uint64_t w);

        std::string compute_digest(std::string input_seq);
        uint64_t get_k() {return k;}
        uint64_t get_w() {return w;}
    private:
        uint64_t k;
        uint64_t w;
        uint8_t lookup_table[UINT8_MAX+1];
        uint64_t loaded_kmers;
        std::vector<MinimizerData> queue;

        void update_lookup_table(char ch, uint8_t val);
};

#endif /* end of include guard: _MINIMIZER_DIGEST_H */
