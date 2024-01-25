#include <iostream>
#include <minimizer_digest.hpp>
#include <pfp_doc.hpp>

int main() {

    // declare a digest object, and input/output
    MinimizerDigest obj(5, 10);
    std::string dna = "";
    std::string digest = "";

    // test 1: simple example with all ACGT characters
    dna = "ACATAGCCATACATATAGATACGATTACCA";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAAGATAACGAT"), "digestion #1 is not correct.");

    // test 2: small example where sequence is smaller than large window size
    dna = "ACAT";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "ACAT"), "digestion #2 is not correct.");

    // test 3: similar to test 1, but inserted an N near th middle of sequence
    dna = "ACATAGCCATACATANTAGATACGATTACCA";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAACGAT"), "digestion #3 is not correct.");

    // test 4: similar dna sequence as test 1, but within first w of sequence
    dna = "ANCATAGCCATACATATAGATACGATTACCA";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "AGCCAACATAAGATAACGAT"), "digestion #4 is not correct.");

    // test 5: similar dna sequence as test 1, but with N within last w of sequence
    dna = "ACATAGCCATACATATAGATACGATTANCCA";
    digest = obj.compute_digest(dna);
    ASSERT((digest == "ACATAAGCCAACATAAGATAACGAT"), "digestion #5 is not correct.");

    return 0;
}