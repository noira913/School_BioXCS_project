#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cctype>


// DNA to RNA conversion
string dnaToRna(const std::string& dna) {
    std::string rna = dna;
    for (char& nucleotide : rna) {
        if (nucleotide == 'T') nucleotide = 'U';
    }
    return rna;
}

// Map codons to amino acids
unordered_map<std::string, std::string> getCodonTable() {
    return {
        {"AUG", "Methionine"}, {"UUU", "Phenylalanine"}, {"UUC", "Phenylalanine"},
        {"UUA", "Leucine"},    {"UUG", "Leucine"},       {"UCU", "Serine"},
        {"UCC", "Serine"},     {"UCA", "Serine"},        {"UCG", "Serine"},
        {"UAU", "Tyrosine"},   {"UAC", "Tyrosine"},      {"UGU", "Cysteine"},
        {"UGC", "Cysteine"},   {"UGG", "Tryptophan"},    {"UAA", "STOP"},
        {"UAG", "STOP"},       {"UGA", "STOP"}
        // Add more codons as needed
    };
}

// Translate RNA sequence to codons
std::vector<std::string> rnaToCodons(const std::string& rna) {
    std::vector<std::string> codons;
    for (size_t i = 0; i < rna.length(); i += 3) {
        if (i + 2 < rna.length()) {
            codons.push_back(rna.substr(i, 3));
        }
    }
    return codons;
}

// Translate codons to amino acids
std::vector<std::string> codonsToAminoAcids(const std::vector<std::string>& codons) {
    std::vector<std::string> aminoAcids;
    auto codonTable = getCodonTable();
    for (const std::string& codon : codons) {
        if (codonTable.find(codon) != codonTable.end()) {
            std::string aminoAcid = codonTable[codon];
            if (aminoAcid == "STOP") break; // Stop translation at stop codon
            aminoAcids.push_back(aminoAcid);
        } else {
            aminoAcids.push_back("Unknown");
        }
    }
    return aminoAcids;
}

// Simulate protein identification
std::string identifyProtein(const std::vector<std::string>& aminoAcids) {
    // Simulated hexokinase sequence
    std::vector<std::string> catalase = {"Methionine", "Phenylalanine", "Leucine", "Serine"};
    if (aminoAcids.size() >= catalase.size() &&
        equal(catalase.begin(), catalse.end(), catalase.begin())) {
        return "catalase (Homo Sapiens)";
    }
    return "Unknown Protein";
}

string cleanSequence(const string& sequence) {
    std::string cleaned;
    for (char c : sequence) {
        if (isalpha(c)) {
            cleaned += toupper(c);
        }
    }
    return cleaned;
}

int main() {
    // Open the input file
    std::ifstream inputFile("human_catalase_dna_sequence.txt");
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file." << std::endl;
        return 1;
    }

    // Read the file and construct the DNA sequence
    std::string line, dnaSequence;
    while (getline(inputFile, line)) {
        std::stringstream ss(line);
        std::string segment;
        while (ss >> segment) {
            dnaSequence += cleanSequence(segment);
        }
    }
    inputFile.close();

    // Convert DNA to RNA
    std::string rnaSequence = dnaToRna(dnaSequence);
    cout << "RNA Sequence: " << rnaSequence << endl;

    // Convert RNA to codons
    std::vector<std::string> codons = rnaToCodons(rnaSequence);
    cout << "Codons: ";
    for (const string& codon : codons) {
        cout << codon << " ";
    }
    cout << endl;

    // Translate codons to amino acids
    vector<string> aminoAcids = codonsToAminoAcids(codons);
    cout << "Amino Acids: ";
    for (const string& aminoAcid : aminoAcids) {
        cout << aminoAcid << " ";
    }
    cout << endl;

    // Identify protein
    string protein = identifyProtein(aminoAcids);
    cout << "Protein: " << protein << endl;

    return 0;
}