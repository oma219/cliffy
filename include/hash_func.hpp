/*
 * File: hash_func.hpp
 * Description: Contains definitions of hash functions
 *              used in the docprofiles repository.
 *
 * Date: January 28th, 2024
 */

#ifndef _HASH_FUNC_H
#define _HASH_FUNC_H

/* 
 * MurmurHash3 developed by Austin Appleby 
 * (github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp)
 */

uint64_t inline MurmurHash3(uint64_t key) {
    uint64_t hash = key;
    hash ^= hash >> 33;
    hash *= 0xff51afd7ed558ccd;
    hash ^= hash >> 33;
    hash *= 0xc4ceb9fe1a85ec53;
    hash ^= hash >> 33;
    return hash;
}

#endif