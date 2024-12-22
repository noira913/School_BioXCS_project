#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cctype>




// DNA to RNA conversion
std::string dnaToRna(const std::string& dna) {
    std::string rna = dna;
    for (char& nucleotide : rna) {
        if (nucleotide == 'T') nucleotide = 'U';
    }
    return rna;
}

// Map codons to amino acids
std::unordered_map<std::string, char> getCodonTable() {
    return {
        {"AUG", 'M'}, {"UUG", 'L'}, {"GUG", 'V'},

        {"UUU", 'F'}, {"UUC", 'F'}, {"UUA", 'L'},  
        {"UCU", 'S'}, {"UCC", 'S'}, {"UCA", 'S'}, {"UCG", 'S'},
        {"UAU", 'Y'}, {"UAC", 'Y'}, 
        {"UGU", 'C'}, {"UGC", 'C'}, {"UGG", 'W'}, 
        {"CUU", 'L'}, {"CUC", 'L'}, {"CUA", 'L'}, {"CUG", 'L'},
        {"CCU", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAU", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGU", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"AUU", 'I'}, {"AUC", 'I'}, {"AUA", 'I'},
        {"ACU", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAU", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGU", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GUU", 'V'}, {"GUC", 'V'}, {"GUA", 'V'},
        {"GCU", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAU", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGU", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'},

        {"UAA", 'STOP'}, {"UAG", 'STOP'}, {"UGA", 'STOP'},
        
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
std::vector<char> codonsToAminoAcids(const std::vector<std::string>& codons) {
    std::vector<char> aminoAcids;
    auto codonTable = getCodonTable();
    for (const std::string& codon : codons) {
        if (codonTable.find(codon) != codonTable.end()) {
            char aminoAcid = codonTable[codon];
            if (aminoAcid == 'STOP') break; // Stop translation at stop codon
            aminoAcids.push_back(aminoAcid);
        } else {
            aminoAcids.push_back('NULL');
        }
    }
    return aminoAcids;
}

// Simulate protein identification
std::string identifyProtein(const std::vector<char>& aminoAcids) {

    std::ifstream controlFile("human_catalase_amino_acid_sequence-control.txt");

    std::string line, amino_acid_sequence;
    while (getline(controlFile, line)) {
        std::stringstream ss(line);
        std::string segment;
        while (ss >> segment) {
            amino_acid_sequence += segment;
        }
    }
    controlFile.close();


    std::vector<char> catalase;

    for (char ch : amino_acid_sequence) {
        if (ch != ' ' && ch != '\0') {
            catalase.push_back(ch);
        }
    }

    if (aminoAcids.size() == catalase.size() && equal(catalase.begin(), catalase.end(), aminoAcids.begin())) {
        return "catalase (Homo Sapiens)";
    }
    else {
        return "Unknown Protein";
    }
}

std::string cleanSequence(const std::string& sequence) {
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
    std::cout << "RNA Sequence: " << rnaSequence << '\n' << std::endl;

    // Convert RNA to codons
    std::vector<std::string> codons = rnaToCodons(rnaSequence);
    std::cout << "Codons: ";
    for (const std::string& codon : codons) {
        std::cout << codon << " ";
    }
    std::cout << '\n' << std::endl;

    // Translate codons to amino acids
    std::vector<char> aminoAcids = codonsToAminoAcids(codons);
    std::cout << "Amino Acids: ";
    for (const char& aminoAcid : aminoAcids) {
        std::cout << aminoAcid;
    }
    std::cout <<  '\n' <<std::endl;

    // Identify protein
    std::string protein = identifyProtein(aminoAcids);
    std::cout << "Protein: " << protein << std::endl;

    return 0;
}