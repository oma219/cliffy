#include <iostream>
#include <minimizer_digest.hpp>
#include <pfp_doc.hpp>
#include <random>
#include <algorithm>
#include <chrono>

std::string generate_random_dna_string(size_t length) {
    // declare needed variables
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d({100, 100, 100, 100, 0});
    const char* nucs[5] = { "A", "C", "G", "T", "N"};

    // generate the sequence of requested length
    std::string dna_seq = "";
    for (size_t i = 0; i < length; i++) {dna_seq += nucs[d(gen)];}
    return dna_seq;
}

int main() {

    // declare a digest object, and input/output
    MinimizerDigest obj(5, 10);
    MinimizerDigest obj2 (31, 35);
    std::string dna = "";
    std::string digest = "";

    // test 1: simple example with all ACGT characters
    dna = "ACATAGCCATACATATAGATACGATTACCA";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAAGATAACGAT"), "digestion #1 is not correct.");
    std::cout << "test 1: passed" << std::endl;

    // test 2: small example where sequence is smaller than large window size
    dna = "ACAT";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "ACAT"), "digestion #2 is not correct.");
    std::cout << "test 2: passed" << std::endl;

    // test 3: similar to test 1, but inserted an N near th middle of sequence
    dna = "ACATAGCCATACATANTAGATACGATTACCA";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAACGAT"), "digestion #3 is not correct.");
    std::cout << "test 3: passed" << std::endl;

    // test 4: similar dna sequence as test 1, but within first w of sequence
    dna = "ANCATAGCCATACATATAGATACGATTACCA";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "AGCCAACATAAGATAACGAT"), "digestion #4 is not correct.");
    std::cout << "test 4: passed" << std::endl;

    // test 5: similar dna sequence as test 1, but with N within last w of sequence
    dna = "ACATAGCCATACATATAGATACGATTANCCA";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAAGATAACGAT"), "digestion #5 is not correct.");
    std::cout << "test 5: passed" << std::endl;

    // test 6: using the same parameters as Kraken2
    dna = "ACATAGCCATACATATAGATACGATTACCAAAAAA";
    digest = obj2.compute_digest(dna);
    ASSERT((digest == "ACATAGCCATACATATAGATACGATTACCAA"), "digestion #6 is not correct.");

    // test 7: try example for k = w
    MinimizerDigest obj3 (5, 5);
    dna = "AAAAACCCCCGGGGGTTTTT";
    digest = obj3.compute_digest(dna);
    ASSERT((digest == "AAAAAAAAACAAACCAACCCACCCCCCCCCCCCCGCCCGGCCGGGCGGGGGGGGGGGGGTGGGTTGGTTTGTTTT"), "digestion #7 is not correct.");
    std::cout << "\n";

    // experiment 1: try using k=31 with various large windows
    size_t input_length = 1000000;
    std::string full_dna_seq = generate_random_dna_string(input_length);
    size_t input_size = full_dna_seq.size();

    std::vector<size_t> list_of_w (49, 0); // vector from 31 to 199, in steps of 4
    std::generate(list_of_w.begin(), list_of_w.begin()+49, [curr_w=3] () mutable { return curr_w+=4;});

    size_t k = 7;
    std::cout << "k,w,inputlength,digestlength,time" << std::endl;

    for (auto curr_w: list_of_w) {
        // create object with specified parameters
        MinimizerDigest new_obj (k, curr_w);

        // compute digest and measure time
        auto start_time = std::chrono::high_resolution_clock::now();
        std::string curr_digest = new_obj.compute_digest(full_dna_seq);

        auto end_time = std::chrono::high_resolution_clock::now();
        double elapsed_time_ms = std::chrono::duration<double, std::milli>(end_time-start_time).count();

        // output results in csv format
        std::cout << k << "," << curr_w << ","  << input_length << "," << curr_digest.size() << "," << elapsed_time_ms << std::endl;
    }

    return 0;
}