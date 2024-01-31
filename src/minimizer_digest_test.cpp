#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
#include <minimizer_digest.hpp>
#include <pfp_doc.hpp>

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

void test1() {
    // test 1: simple example with all ACGT characters
    std::string dna = "ACATAGCCATACATATAGATACGATTACCA";
    MinimizerDigest obj(5, 10);

    std::string digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAAGATAACGAT"), "digestion #1 is not correct.");
    std::cout << "test 1: passed" << std::endl; 
}

void test2() {
    // test 2: small example where sequence is smaller than large window size
    MinimizerDigest obj(5, 10);
    std::string dna = "ACAT";

    std::string digest = obj.compute_digest(dna);
    ASSERT((digest == "ACAT"), "digestion #2 is not correct.");
    std::cout << "test 2: passed" << std::endl;
}

void test3() {
    // test 3: similar to test 1, but inserted an N near th middle of sequence
    MinimizerDigest obj(5, 10);
    std::string dna = "ACATAGCCATACATANTAGATACGATTACCA";

    std::string digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAACGAT"), "digestion #3 is not correct.");
    std::cout << "test 3: passed" << std::endl;
}

void test4() {
    // test 4: similar dna sequence as test 1, but within first w of sequence
    MinimizerDigest obj(5, 10);
    std::string dna = "ANCATAGCCATACATATAGATACGATTACCA";

    std::string digest = obj.compute_digest(dna);
    ASSERT((digest == "AGCCAACATAAGATAACGAT"), "digestion #4 is not correct.");
    std::cout << "test 4: passed" << std::endl;
}

void test5() {
    // test 5: similar dna sequence as test 1, but with N within last w of sequence
    MinimizerDigest obj(5, 10);
    std::string dna = "ACATAGCCATACATATAGATACGATTANCCA";

    std::string digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAAGATAACGAT"), "digestion #5 is not correct.");
    std::cout << "test 5: passed" << std::endl;
}

void test6() {
    // test 6: using the same parameters as Kraken2
    MinimizerDigest obj2 (31, 35);
    std::string dna = "ACATAGCCATACATATAGATACGATTACCAAAAAA";

    std::string digest = obj2.compute_digest(dna);
    ASSERT((digest == "ACATAGCCATACATATAGATACGATTACCAA"), "digestion #6 is not correct.");
    std::cout << "test 6: passed" << std::endl;
}

void test7() {
    // test 7: try example for k = w
    MinimizerDigest obj3 (5, 5);
    std::string dna = "AAAAACCCCCGGGGGTTTTT";

    std::string digest = obj3.compute_digest(dna);
    ASSERT((digest == "AAAAAAAAACAAACCAACCCACCCCCCCCCCCCCGCCCGGCCGGGCGGGGGGGGGGGGGTGGGTTGGTTTGTTTTTTTTT"), "digestion #7 is not correct.");
    std::cout << "test 7: passed" << std::endl;
}

void test8() {
    // test 8: try example where k=4 to test for minimizer alphabet
    MinimizerDigest obj(4, 4);
    std::string dna = "AAAAAGCCATACATATAGATACGATTTTTTA";

    std::string digest = obj.compute_digest(dna);
    ASSERT((digest=="AAAAAAAGAAGCAGCCGCCACCATCATAATACTACAACATCATAATATTATAATAGTAGAAGATGATAATACTACGACGACGATGATTATTTTTTTTTTA"),
            "digestion #8 is not correct.");
    std::cout << "test 8: passed" << std::endl;
}

void test9() {
    // test 9: try example where k=4 and using hash order for first time
    MinimizerDigest obj(4, 10, false);
    std::string dna = "AAAAAGCCATACATATAGATACGATTTTTTA";

    std::string digest = obj.compute_digest(dna);
    ASSERT((digest=="AAAAAAAGAGCCTACAAGATGATAATTT"),
            "digestion #9 is not correct.");
    std::cout << "test 9: passed" << std::endl;
}

void test10() {
    // test 10: try the minimizer alp digestion and make sure it is the right length
    MinimizerDigest obj(4, 10, false, true);
    std::string dna = "CCCCCCCCCCCCCCCCCC";

    std::string digest = obj.compute_digest(dna);
    ASSERT((digest.size() == 1), "digestion #10 is not correct.");
    std::cout << "test 10: passed" << std::endl;
}

void run_exp1(size_t k, bool use_minimizer_alp) {
    // experiment 1: try using k with various large windows
    size_t input_length = 1000000;
    std::string full_dna_seq = generate_random_dna_string(input_length);
    size_t input_size = full_dna_seq.size();

    std::vector<size_t> list_of_w (49, 0); // vector from 31 to 199, in steps of 4
    std::generate(list_of_w.begin(), list_of_w.begin()+49, [curr_w=(k-4)] () mutable { return curr_w+=4;});

    // iterate through different window sizes
    std::cout << "k,w,inputlength,digestlength,time" << std::endl;
    for (auto curr_w: list_of_w) {
        // create object with specified parameters
        MinimizerDigest new_obj (k, curr_w, false, use_minimizer_alp);

        // compute digest and measure time
        auto start_time = std::chrono::high_resolution_clock::now();
        std::string curr_digest = new_obj.compute_digest(full_dna_seq);

        auto end_time = std::chrono::high_resolution_clock::now();
        double elapsed_time_ms = std::chrono::duration<double, std::milli>(end_time-start_time).count();

        // output results in csv format
        std::cout << k << "," << curr_w << ","  << input_length << "," << curr_digest.size() << "," << elapsed_time_ms << std::endl;
    }
}

int main() {


    // unit tests
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();
    test9();
    test10();
    std::cout << "\n";

    // experiments
    run_exp1(7, false); std::cout << "\n";
    run_exp1(31, false); std::cout << "\n";
    run_exp1(4, true); std::cout << "\n";

    return 0;
}