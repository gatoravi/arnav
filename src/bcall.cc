/*  bcall.cc -- main

    Copyright (c) 2016, Avinash Ramu

    Author: Avinash Ramu <avinash3003@yahoo.co.in>

    This file is part of bcall.

    bcall is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bcall is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bcall.  If not, see <http://www.gnu.org/licenses/>. */

#include <bitset>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include "cereal/types/unordered_map.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/archives/binary.hpp"
#include "gzstream/gzstream.h"
#include "bh-fdr.h"
#include "common.h"
#include "Rmath.h"

using namespace std;

//Map from chromosome to integer
std::map<string, int> chr_to_int = {
    {"1", 0}, {"2", 1}, {"3", 2}, {"4", 3},
    {"5", 4}, {"6", 5}, {"7", 6}, {"8", 7},
    {"9", 8}, {"10", 9}, {"11", 10}, {"12", 11},
    {"13", 12}, {"14", 13}, {"15", 14}, {"16", 15},
    {"17", 16}, {"18", 17}, {"19", 18}, {"20", 19},
    {"21", 20}, {"22", 21}, {"X", 22}, {"Y", 23},
    {"MT", 24}
};

//Map from chromosome to integer
std::map<int, string> int_to_chr = {
    {0, "1"}, {1, "2"}, {2, "3"}, {3, "4"},
    {4, "5"}, {5, "6"}, {6, "7"}, {7, "8"},
    {8, "9"}, {9, "10"}, {10, "11"}, {11, "12"},
    {12, "13"}, {13, "14"}, {14, "15"}, {15, "16"},
    {16, "17"}, {17, "18"}, {18, "19"}, {19, "20"},
    {20, "21"}, {21, "22"}, {22, "X"}, {23, "Y"},
    {24, "MT"}
};

//key is sampleID, value is path to readcount.gz file
std::unordered_map<string, string> sample_to_readcountfile;

//struct to store accumulated read counts
struct readcounts {
    uint64_t total_ref_count;
    uint64_t total_alt_count;
    template <class Archive>
    void serialize(Archive& ar) {
        ar(total_ref_count, total_alt_count);
    }
};

//key is pos << 5 | chr_int, value is struct readcounts
std::unordered_map<uint64_t, readcounts> site_readcounts;
//vector of p-values for a sample, apply FDR adjustment to these
std::vector<double> pvalues;
//Map of line => binomial_test_pval
std::unordered_map<string, double> lines_pvalues;
//Number of calls not in the site of interest
uint64_t not_in_map = 0;

//usage message
int usage() {
    cerr << endl << "./bcall ";
    cerr << endl << "\t prior-and-call file_with_mpileupcounts op_variants_file_name";
    cerr << endl << "\t prior-dump file_with_mpileupcounts op_priors_dump_file_name";
    cerr << endl << "\t prior-dump-fixed file_with_mpileupcounts op_priors_dump_file_name fixed-sites.bed.gz";
    cerr << endl << "\t prior-merge priors_dump_file_list merged_dump_file";
    cerr << endl << "\t prior-print priors_dump_file_list";
    cerr << endl << "\t call-using-merged file_with_mpileupcounts merged_dump_file";
    cerr << endl;
    cerr << endl << "The file_with_mpileupcounts has two columns, sample_name and path to "
                    "\nfile with readcounts that have been compressed with "
                    "\nbgzip/gzip, for e.g `SRR1 SRR1_readcounts.gz`";
    cerr << endl << "The prior-dump command creates a file that is a C++ map serialized using cereal "
                    "\nand written to disk. This can then be read by a different "
                    "\nprocess. The prior-dump-fixed only looks at sites specified by the bed.gz file.";
    cerr << endl << "The prior-merge command requires a file that has two columns, "
                    "\ndump_name and path to dump file.";
    cerr << endl << "The prior-print command requires a dump file.";
    cerr << endl << endl;
    return 0;
}

//Split a key into constituent chr and pos
pair<string, uint32_t> decode_key(uint64_t key) {
    int chr_index = static_cast<uint32_t>(key);
    if(chr_index < 0 || chr_index > 24) {
        throw runtime_error("Unable to decode key " + to_string(key));
    }
    //left most 32 bits
    string chr = int_to_chr[chr_index];
    //Right most 32 bits
    uint32_t pos = static_cast<uint32_t>(key >> 32);
    return std::make_pair(chr, pos);
}

//Create a key that is of type double
//The key is unique for a chr:pos combination
//Left shift the position, AND the chr bits
uint64_t create_key(string chr, uint64_t pos) {
    uint32_t chr_int = chr_to_int[chr];
    /* //see binary encoding
    bitset<5> x(chr_int);
    cerr << "chr_x is " << x << endl;
    bitset<32> pos_x(pos);
    cerr << "pos_x is " << pos_x << endl;
    cerr << "unique key is " << unique_key << endl;
    bitset<64> key_x(unique_key);
    cerr << "key_x is " << key_x << endl;
    */
    uint64_t unique_key = (uint64_t) (pos << 32) | chr_int;
    return unique_key;
}

//Print output header
void print_header(ostream& out = cout) {
    out << "chr" << "\t"
        << "pos" << "\t"
        << "depth" << "\t"
        << "ref_base" << "\t"
        << "refcount" << "\t"
        << "altcount" << "\t"
        << "acount" << "\t"
        << "ccount" << "\t"
        << "gcount" << "\t"
        << "tcount" << "\t"
        << "ncount" << "\t"
        << "indelcount" << "\t"
        << "site_total_ref_count" << "\t"
        << "site_total_alt_count" << "\t"
        << "region" << "\t"
        << "p_value" << "\t"
        << "sample"
        << endl;
}

//Print output line
inline void print_significant_lines(double pval_cutoff, ostream& out = cout) {
    for (auto& kv : lines_pvalues) {
        //Print the line if significant
        if (kv.second <= pval_cutoff) {
            out << kv.first << endl;
        }
    }
}

//Takes a line of the readcount file as input and applies
//the binomial test, the `p` is calculated using all the
//readcounts at this site across samples
void apply_model_readcount_line(string sample, string line, bool fixed_sites = false) {
    stringstream ss(line);
    string chr, ref;
    uint32_t pos, depth, ref_count, all_alt_count, alt_count,
             acount, ccount, gcount, tcount;
    ss >> chr >> pos >> depth >> ref >> ref_count;
    //Get counts for specific nucleotides
    ss >> all_alt_count >> acount >> ccount >> gcount >> tcount;
    switch(ref[0]) {
        case 'A':
            alt_count = std::max({ccount, gcount, tcount});
            break;
        case 'C':
            alt_count = std::max({acount, gcount, tcount});
            break;
        case 'G':
            alt_count = std::max({acount, ccount, tcount});
            break;
        case 'T':
            alt_count = std::max({acount, ccount, gcount});
            break;
        default://Don't call at this position for N etc
            return;
    }
    if(chr_to_int.find(chr) != chr_to_int.end()) {
        uint64_t key = create_key(chr, pos);
        if(site_readcounts.find(key) == site_readcounts.end()) {
            //throw runtime_error("Unable to find chr/pos " + chr + " " + to_string(pos));
            //Not in the merged-map
            //cerr << "not in map" << endl;
            not_in_map++;
            return;
        }
        //Subtract this sample's counts from the prior
        uint64_t total_alt_count = site_readcounts[key].total_alt_count - all_alt_count;
        uint64_t total_ref_count = site_readcounts[key].total_ref_count - ref_count;
        uint64_t total_rc =
            site_readcounts[key].total_ref_count + site_readcounts[key].total_alt_count - ref_count - all_alt_count;
        double prior_p =
            (double)total_alt_count / (double) total_rc;
        if (prior_p == 0) { //Set lower bound on the sequencing error rate
            prior_p = 0.001;
        }
        //This should help with pulling out from tabix
        string region = chr + ":" + common::num_to_str(pos);
        if (alt_count >= 2) { //only look at sites with 2 or more variant supporting reads
            //(1 - pbinom(8, 10, 0.5))   == binom.test(9, 10, 0.5, alternative="greater")
            double p_value = 1 - pbinom(alt_count - 1, ref_count + alt_count, prior_p, true, false);
            pvalues.push_back(p_value);
            if (p_value <= 0.05) { //Only store lines <= 0.05, these will be further filtered out while printing
                line = line + "\t" +
                       common::num_to_str(total_ref_count) + "\t" + common::num_to_str(total_alt_count) +
                       "\t" + region + "\t" + common::num_to_str(p_value) + "\t" + sample;
                //Store the line with its pvalue in the map
                lines_pvalues[line] = p_value;
            }
        } else {
            //cerr << "not greater than 2 variant reads" << endl;
        }
    }
}

//parse a line from the readcount file
void calculate_prior_line(string sample, string line, bool fixed_sites = false) {
    stringstream ss(line);
    string chr, ref;
    uint32_t pos, depth, ref_count, alt_count;
    ss >> chr >> pos >> depth >> ref;
    ss >> ref_count >> alt_count;
    if(chr_to_int.find(chr) != chr_to_int.end()) {
        uint64_t key = create_key(chr, pos);
        if(site_readcounts.find(key) == site_readcounts.end()) {
            //Sites are fixed by the BED file, don't add new sites
            if (fixed_sites) {
                return;
            }
            //if (alt_count == 0) { //Experimental - only look at sites with Non-zero alt
            //    return;
            //}
            site_readcounts[key].total_ref_count = 0;
            site_readcounts[key].total_alt_count = 0;
        }
        site_readcounts[key].total_ref_count += ref_count;
        site_readcounts[key].total_alt_count += alt_count;
    }
}

//iterate through readcount file - And apply model to each line
//first arg is name of the gz file
//second argument is the function to apply to each line of the file
//The third argument specifies if we should look only at fixed sites of interest(passed as a param)
void parse_readcount_file(string sample, string gzfile, function<void(string, string, bool)> func,
                          bool fixed_sites = false) {
    igzstream in(gzfile.c_str());
    cerr << "Opening " << gzfile << endl;
    std::string line;
    int line_count = 0;
    while(std::getline(in, line)){
        if(line.substr(0, 3) == "chr") {
            continue;
        }
        func(sample, line, fixed_sites);
        line_count += 1;
    }
    if(!line_count) {
        throw runtime_error("Readcount file empty - " + gzfile);
    }
    cerr << "Read " << line_count << " lines from " << gzfile << endl;
}

//Iterate through each sample's readcounts
void calculate_priors(bool fixed_sites = false) {
    for (auto& kv : sample_to_readcountfile) {
        cerr << "Processing " << kv.first << endl;
        parse_readcount_file(kv.first, kv.second, calculate_prior_line, fixed_sites);
        cerr << "Size of readcount map is " << site_readcounts.size() << endl;
    }
}

//Iterate through each sample's readcounts and call
void apply_model() {
    for (auto& kv : sample_to_readcountfile) {
        cerr << endl << "Applying model to " << kv.first << endl;
        not_in_map = 0; //Global variable - BAD
        pvalues.clear(); //Remove all p-values from previous sample
        lines_pvalues.clear(); //Remove all lines and p-values from previous sample
        parse_readcount_file(kv.first, kv.second, apply_model_readcount_line);
        cerr << "The number of tests performed for this sample is: " << pvalues.size() << endl;
        cerr << "Applying the BH-procedure to control FDR." << endl;
        double pvalue_cutoff = 0.05/(float)9597; //Bonferroni
        //double pvalue_cutoff = bh_fdr(pvalues, 0.05); //FDR
        cerr << "The p-value cutoff is:" << pvalue_cutoff << endl;
        cerr << "The number of sites not in the region of interest is " << not_in_map;
        print_header();
        print_significant_lines(pvalue_cutoff);
    }
}

//Print the header for priors for each site
void print_priors_header(ostream& fout) {
    fout << "#chr" << "\t" << "pos" << "\t";
    fout << "key" << "\t";
    fout << "total_ref_count" << "\t";
    fout << "total_alt_count";
    fout << endl;
}

//Print the priors for each site
void print_priors(ostream& fout, bool print_zeros = true) {
    print_priors_header(fout);
    for (auto& kv : site_readcounts) {
        if(print_zeros == false &&
           kv.second.total_ref_count == 0 &&
           kv.second.total_alt_count == 0) {
            continue;
        }
        auto decoded = decode_key(kv.first);
        fout << decoded.first << "\t" << decoded.second << "\t";
        fout << kv.first << "\t";
        fout << kv.second.total_ref_count << "\t";
        fout << kv.second.total_alt_count;
        fout << endl;
    }
}

//Write the priors map to a file
void write_priors(const string& output_file) {
    ofstream fout(output_file);
    if (!fout.is_open()) {
        throw runtime_error("unable to open " + output_file +
                            " for writing.");
    }
    cereal::BinaryOutputArchive archive(fout);
    archive(site_readcounts);
    fout.close();
}

//Read the prior from the merged prior-dump file
void read_priors_merged(string merged_prior_dump) {
    cerr << "Reading merged dump " << merged_prior_dump << endl;
    ifstream fin(merged_prior_dump);
    if (!fin.is_open()) {
        throw runtime_error("unable to open " + merged_prior_dump +
                            " for reading priors.");
    }
    cereal::BinaryInputArchive archive(fin);
    archive(site_readcounts);
    fin.close();
}

//Read the priors map from a file
void read_priors() {
    for (auto& kv : sample_to_readcountfile) {
        cerr << "Reading dump " << kv.first << endl;
        string prior_file = kv.second;
        ifstream fin(prior_file);
        if (!fin.is_open()) {
            throw runtime_error("unable to open " + prior_file +
                                " for reading priors.");
        }
        cereal::BinaryInputArchive archive(fin);
        std::unordered_map<uint64_t, readcounts> temp_site_readcounts;
        archive(temp_site_readcounts);
        //Aggregate temp with main
        for (auto& kv : temp_site_readcounts) {
            if (site_readcounts.find(kv.first) == site_readcounts.end()) {
                site_readcounts[kv.first] = kv.second;
            } else {
                site_readcounts[kv.first].total_ref_count += kv.second.total_ref_count;
                site_readcounts[kv.first].total_alt_count += kv.second.total_alt_count;
            }
        }
        fin.close();
    }
}

//Iterate through the samples file and calc site-specific priors
void read_samples(char* samples_file) {
    ifstream sample_fh(samples_file, ios::in);
    string line;
    string sample, readcountfile;
    int line_count = 0;
    while (getline(sample_fh, line)) {
        stringstream iss(line);
        iss >> sample >> readcountfile;
        sample_to_readcountfile[sample] = readcountfile;
        line_count++;
    }
    if(!line_count) {
        throw runtime_error("Sample file empty - " + string(samples_file));
    }
}

void add_bedline_to_map(string line) {
    string chr;
    uint32_t start, end;
    stringstream ss(line);
    ss >> chr >> start >> end;
    //Skip header
    if(chr == "track") {
        return;
    }
    for (uint32_t pos = start + 1; pos <= end; pos++) {
        uint64_t key = create_key(chr, pos);
        //Initialize total_ref and total_alt to zero
        site_readcounts[key] = {0, 0};
    }
}

//Read a BED file and only store those sites in the map
void initialize_fixed_map(string bedFile) {
    igzstream in(bedFile.c_str());
    cerr << "Initializing map with sites in  " << bedFile << endl;
    std::string line;
    int line_count = 0;
    while(std::getline(in, line)){
        add_bedline_to_map(line);
        line_count += 1;
    }
    if(!line_count) {
        throw runtime_error("Bedfile empty - " + bedFile);
    }
    cerr << "Read " << line_count << " lines from " << bedFile << endl;
    cerr << "Size of readcount map is " << site_readcounts.size() << endl;
}

int main(int argc, char* argv[]) {
    try {
        if(argc >= 4) {
            read_samples(argv[2]);
            if (string(argv[1]) == "prior-and-call") {
                    calculate_priors();
                    //print_priors(cout);
                    print_header();
                    apply_model();
                    return 0;
            }
            else if (string(argv[1]) == "prior-dump") {
                    calculate_priors();
                    //print_priors(cout, false);
                    write_priors(string(argv[3]));
                    return 0;
            }
            else if (argc > 4 && string(argv[1]) == "prior-dump-fixed") { //prior dump on fixed sites
                    initialize_fixed_map(string(argv[4]));
                    calculate_priors(true);
                    print_priors(cout, false);
                    write_priors(string(argv[3]));
                    return 0;
            }
            else if (string(argv[1]) == "prior-merge") {
                read_samples(argv[2]);
                read_priors();
                print_priors(cout);
                write_priors(string(argv[3]));
                return 0;
            }
            else if (string(argv[1]) == "call-using-merged") {
                read_samples(argv[2]);
                read_priors_merged(argv[3]);
                apply_model();
                return 0;
            }
        }
        else if (argc == 3 && string(argv[1]) == "prior-print") {
                read_priors_merged(argv[2]);
                print_priors(cout);
                return 0;
        }
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return 1;
    }
    return usage();
}
