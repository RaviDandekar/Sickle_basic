//
//  sickle_processing.h
//  Sickle_basic
//
//  Created by Ravi Dandekar on 9/6/13.
//  Copyright (c) 2013 Korf Lab. All rights reserved.
//

#ifndef __Sickle_basic__sickle_processing__
#define __Sickle_basic__sickle_processing__

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <algorithm>
#include <bitset>
#include <bitset>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <map>
#include <utility>
/*  HMMER 3.0 Headers  */
extern "C" {
#include "p7_config.h"
#include "easel.h"
#include "hmmer.h"
}
/*  StochHMM Headers  */
#include "StochHMMlib.h"
#include "text.h"
#include "track.h"
#include "stochMath.h"




class SEQanalysis{
public:
    enum AMINO_ACIDS {
        alanine = 0,            // A
        arginine = 1,           // R
        asparagine = 2,         // N
        aspartic_acid = 3,      // D
        cysteine = 4,           // C
        glutamic_acid = 5,      // E
        glutamine = 6,          // Q
        glycine = 7,            // G
        histidine = 8,          // H
        isoleucine = 9,         // I
        leucine = 10,           // L
        lysine = 11,            // K
        methionine = 12,        // M
        phenylalanine = 13,     // F
        proline = 14,           // P
        serine = 15,            // S
        threonine = 16,         // T
        tryptophan = 17,        // W
        tyrosine = 18,          // Y
        valine = 19,            // V
        start_stop = 20,        // *
        not_defined = 21,       // X
    };
    
    SEQanalysis(){};
    SEQanalysis(std::string& trans_tbl, std::map<std::string, SEQanalysis::AMINO_ACIDS>& aa_dictionary, const P7_HMM *hmm, const ESL_ALPHABET *abc, std::string& protein_seq, std::string& codon, const std::string): \
    trans_tb_fh(trans_tbl), translation_table(aa_dictionary), hmm_file(hmm), digital_alphabet(abc), protein_sequence(protein_seq){};
    
    
    /*----------------------------------------------*/
    /*FUNCTIONS                                     */
    /*  1. parse translation table to define codons */
    /*  2. translate fasta seq                      */
    /*  3. claculatate HMMER score                  */
    /*----------------------------------------------*/
    static std::map<std::string, AMINO_ACIDS> parse_translation_table(std::string trans_tbl);
    std::string translate (const std::string *my_seq, std::map<std::string, SEQanalysis::AMINO_ACIDS> &aa_dictionary);
    static float hmmer_score (const ESL_ALPHABET *abc, const P7_HMM *hmm, std::string protein_seq);
    
    
private:
    std::map<std::string, SEQanalysis::AMINO_ACIDS> translation_table;
    std::string trans_tb_fh;
    std::string protein_sequence;               //Translated Protein Sequence
    const P7_HMM *hmm_file;                     //Pointer to HMM file
    const ESL_ALPHABET *digital_alphabet;       //Pointer to digital alphabet
};


#endif /* defined(__Sickle_basic__sickle_processing__) */
