/* ms_rle_string - Extension of the r-index rle_string to compute matching statistics
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file ms_rle_string.hpp
   \brief ms_rle_string.hpp Extension of the r-index rle_string to compute matching statistics.
   \author Massimiliano Rossi
   \date 10/07/2020
*/

#ifndef _MS_RLE_STRING_HH
#define _MS_RLE_STRING_HH

#include <common.hpp>
#include <rle_string.hpp>
#include <bitset>

template <
    class sparse_bitvector_t = ri::sparse_sd_vector, //predecessor structure storing run length
    class string_t = ri::huff_string                 //run heads
    >
class ms_rle_string : public ri::rle_string<sparse_bitvector_t, string_t>
{
    public:
        static const u_char TERMINATOR = 1;

    ms_rle_string():ri::rle_string<sparse_bitvector_t, string_t>() {}

    /*
     * constructor: build structure on the input string
     * \param input the input string without 0x0 bytes in it.
     * \param B block size. The main sparse bitvector has R/B bits set (R being number of runs)
     *
     */
    ms_rle_string(string &input, ulint B = 2): ri::rle_string<sparse_bitvector_t, string_t>(input, B) {}

    ms_rle_string(std::ifstream &ifs, ulint B = 2): ri::rle_string<sparse_bitvector_t, string_t>(ifs, B) {}

    // Construction from run-length encoded BWT
    ms_rle_string(std::ifstream& heads, std::ifstream& lengths, ulint B = 2) {
        heads.clear(); heads.seekg(0);
        lengths.clear(); lengths.seekg(0);

        this->B = B;
        auto runs_per_letter_bv = vector<vector<bool> >(256);

        //runs in main bitvector
        vector<bool> runs_bv;

        // Reads the run heads file
        string run_heads_s;
        heads.seekg(0, heads.end);
        run_heads_s.resize(heads.tellg());
        heads.seekg(0, heads.beg);
        heads.read(&run_heads_s[0], run_heads_s.size());

        size_t pos = 0;
        this->n = 0;
        this->R = run_heads_s.size();

        // Compute runs_bv and runs_per_letter_bv
        for (size_t i = 0; i < run_heads_s.size(); ++i)
        {
            size_t length = 0;
            lengths.read((char*)&length, 5);
            
            // fixed old bug since previous version did not
            // have the unsigned(), otherwise it will be 
            // treated as signed char
            uint8_t curr_ch = unsigned(run_heads_s[i]);
            if(curr_ch <= TERMINATOR) { // change 0 to 1
                run_heads_s[i]=TERMINATOR;
            }

            std::fill_n(std::back_inserter(runs_bv), length-1, false);
            runs_bv.push_back(i%B==B-1);

            std::fill_n(std::back_inserter(runs_per_letter_bv[curr_ch]), length-1, false);
            runs_per_letter_bv[curr_ch].push_back(true);

            this->n += length;
        }
        runs_bv.push_back(false);

        if (runs_bv.size() != (this->n + 1)) {
            std::cout << "Assertion Error: runs_bv.size() should equal n." << std::endl;
            std::exit(1);
        }

        //now compact structures
        ulint t = 0;
        for(ulint i=0;i<256;++i)
            t += runs_per_letter_bv[i].size();
        
        if (t != this->n) {
            std::cout << "Assertion Error: t != n" << std::endl;
            std::exit(1);
        }

        // construct the runs data-structure
        this->runs = sparse_bitvector_t(runs_bv);

        // fast direct array: char -> bitvector.
        this->runs_per_letter = vector<sparse_bitvector_t>(256);
        for(ulint i=0;i<256;++i)
            this->runs_per_letter[i] = sparse_bitvector_t(runs_per_letter_bv[i]);

        // construct the huffman wavelet tree
        this->run_heads = string_t(run_heads_s);

        if (this->run_heads.size() != this->R) {
            std::cout << "Assertion Error: run_heads.size() != R" << std::endl;
            std::exit(1);
        }
        

        // debugging:
        // for (size_t i = 0; i < 256; i++) {
        //     std::cout << "ms_rle_string: i = " << i << ", length = " << this->runs_per_letter[i].size() << std::endl;
        // }   
        // std::cout << "wt-tree size: " << this->run_heads.size() << std::endl;
        // for (size_t i = 0; i < 256; i++) {
        //     std::cout << "ms_rle_string: i = " << i << ", huff.rank() = " << this->run_heads.rank(this->run_heads.size(), i) << std::endl;
        // }
    }

    size_t number_of_runs_of_letter(uint8_t c)
    {
        return this->runs_per_letter[c].number_of_1();
    }

    size_t number_of_letter(uint8_t c)
    {
        return this->runs_per_letter[c].size();
    }

    /*
     * The next three methods were added in order to support
     * the new version of the thresholds data-structure that
     * is bitvector-based.
     */

    // number of runs of character c in in position i
    size_t run_head_rank(const size_t i, const uint8_t c) {
        assert(i < this->R);
        size_t j = this->run_heads.rank(i, c);
        return j;
    }

    inline std::pair<size_t,size_t> run_and_head_rank(const size_t i, const uint8_t c) {
        assert(i < this->R);
        const size_t j = this->run_heads.rank(i, c);
        if( j < 1)
            return make_pair(j,j);
        const size_t k = this->runs_per_letter[c].select(j - 1) + 1; // j-1 because the select is 0 based
        return make_pair(j, k);
    }

    // Select the i-th run of c
    size_t run_head_select(const size_t i, const uint8_t c) {
        assert(i < this->R and i > 0);
        return this->run_heads.select(i - 1, c);
    }


    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    ulint serialize(std::ostream &out) 
    {
        return ri::rle_string<sparse_bitvector_t, string_t>::serialize(out);
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in)
    {
        ri::rle_string<sparse_bitvector_t, string_t>::load(in);
    }

    private :
};

typedef ms_rle_string<ri::sparse_sd_vector> ms_rle_string_sd;
typedef ms_rle_string<ri::sparse_hyb_vector> ms_rle_string_hyb;

#endif /* end of include guard: _MS_RLE_STRING_HH */
