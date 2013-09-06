//
//  sickle_processing.cpp
//  Sickle_basic
//
//  Created by Ravi Dandekar on 9/6/13.
//  Copyright (c) 2013 Korf Lab. All rights reserved.
//

#include "sickle_processing.h"


/*---------------------------*/
/*  PARSE TRANSLATION TABLE  */
/*---------------------------*/

std::map<std::string, SEQanalysis::AMINO_ACIDS> SEQanalysis::parse_translation_table (std::string trans_tbl) {
    enum parse_state {
        parse_aa = 0,
        parse_starts = 1,
        parse_base1 = 2,
        parse_base2 = 3,
        parse_base3 = 4,
    };
    parse_state pstate = parse_aa;
    
    std::ifstream file;
	file.open(trans_tbl.c_str());
	if (!file.is_open()){
		std::cerr << "Couldn't open Translation Table file\n";
    }
    
    char peek_char;
    std::string amino_acids("");
    std::string base1("");
    std::string base2("");
    std::string base3("");
    std::map<std::string, SEQanalysis::AMINO_ACIDS> codon_key;
    
    while (!file.eof()) {
        peek_char = file.get();
        
        // if there is a space or tab empty string for each AA and base lines
        if (peek_char == ' ' || peek_char == '\t') {
            switch (pstate){
                case parse_aa:
                    amino_acids = "";
                    break;
                case parse_starts:
                    break;
                case parse_base1:
                    base1 = "";
                    break;
                case parse_base2:
                    base2 = "";
                    break;
                case parse_base3:
                    base3 = "";
                    break;
            }
        }
        //If it reaches the end of the line add 1 to the state in order to move to the next line (refer to enum)
        else if (peek_char == '\n') {
            pstate++;
        }
        //Else will contain the nucleotide and amino acid keys based on translation table
        else {
            switch (pstate){
                case parse_aa:
                    amino_acids.append(1, peek_char);
                    break;
                case parse_starts:
                    break;
                case parse_base1:
                    base1.append(1, peek_char);
                    break;
                case parse_base2:
                    base2.append(1, peek_char);
                    break;
                case parse_base3:
                    base3.append(1, peek_char);
                    break;
            }
        }
    }
    // Creates codon dictionary (similar to hash created in <translation_table_formatter.pl>)
    for(int i=0; i < amino_acids.length(); i++) {
        std::string temp_str("");
        SEQanalysis::AMINO_ACIDS temp_enum;
        switch (amino_acids[i]) {
            case 'A': case 'a':
                temp_enum = SEQanalysis::alanine;
                break;
            case 'R': case 'r':
                temp_enum = SEQanalysis::arginine;
                break;
            case 'N': case 'n':
                temp_enum = SEQanalysis::asparagine;
                break;
            case 'D': case 'd':
                temp_enum = SEQanalysis::aspartic_acid;
                break;
            case 'C': case 'c':
                temp_enum = SEQanalysis::cysteine;
                break;
            case 'E': case 'e':
                temp_enum = SEQanalysis::glutamic_acid;
                break;
            case 'Q': case 'q':
                temp_enum = SEQanalysis::glutamine;
                break;
            case 'G': case 'g':
                temp_enum = SEQanalysis::glycine;
                break;
            case 'H': case 'h':
                temp_enum = SEQanalysis::histidine;
                break;
            case 'I': case 'i':
                temp_enum = SEQanalysis::isoleucine;
                break;
            case 'L': case 'l':
                temp_enum = SEQanalysis::leucine;
                break;
            case 'K': case 'k':
                temp_enum = SEQanalysis::lysine;
                break;
            case 'M': case 'm':
                temp_enum = SEQanalysis::methionine;
                break;
            case 'F': case 'f':
                temp_enum = SEQanalysis::phenylalanine;
                break;
            case 'P': case 'p':
                temp_enum = SEQanalysis::proline;
                break;
            case 'S': case 's':
                temp_enum = SEQanalysis::serine;
                break;
            case 'T': case 't':
                temp_enum = SEQanalysis::threonine;
                break;
            case 'W': case 'w':
                temp_enum = SEQanalysis::tryptophan;
                break;
            case 'Y': case 'y':
                temp_enum = SEQanalysis::tyrosine;
                break;
            case 'V': case 'v':
                temp_enum = SEQanalysis::valine;
                break;
            case '*':
                temp_enum = SEQanalysis::start_stop;
                break;
            case 'X': case 'x':
                temp_enum = SEQanalysis::not_defined;
                break;
        }
        temp_str.append(1, base1[i]);
        temp_str.append(1, base2[i]);
        temp_str.append(1, base3[i]);
        
        if (temp_str.length() == 3) {
            codon_key.insert(make_pair(temp_str, temp_enum));
        }
    }
    file.close();
    return codon_key;
}




/*----------------------*/
/*  TRANSLATE SEQUENCE  */
/*----------------------*/

std::string SEQanalysis::translate(const std::string *my_seq, std::map<std::string, SEQanalysis::AMINO_ACIDS> &aa_dictionary){
    
    // Translate FASTA sequence
    std::string protein_seq("");
    static const char str[] = "ARNDCEQGHILKMFPSTWYV*X";
    
    for (int i=0; i < my_seq->length(); i += 3) {
        std::string codon("");
        codon = my_seq->substr(i,3);
        std::map<std::string, SEQanalysis::AMINO_ACIDS>::iterator codon_to_aa = aa_dictionary.find(codon);
        
        if (codon_to_aa != aa_dictionary.end()) {
            protein_seq.append(1,str[(*codon_to_aa).second]);
        }
    }
    return protein_seq;
}




/*-------------------------*/
/*  CALCULATE HMMER SCORE  */
/*-------------------------*/

float SEQanalysis::hmmer_score(const ESL_ALPHABET *abc, const P7_HMM *hmm, std::string protein_seq) {
    
    //Convert string into char array
    char *seq = (char*)protein_seq.c_str();
    std::cout << (*seq) << std::endl;
    
    P7_BG        *bg     = NULL; /* background */
	P7_PROFILE   *gm     = NULL; /* generic model */
	P7_GMX       *mx     = NULL; /* viterbi matrix */
	ESL_DSQ      *dsq    = NULL; /* digital sequence */
	int           L      = 0;    /* length of sequence */
	float         score;
	
	/* digitize sequence */
	L = strlen(seq);
	esl_abc_CreateDsq(abc, seq, &dsq);
	
	/* background */
	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, L);
	
	/* profile */
	gm = p7_profile_Create(hmm->M, abc);
	p7_ProfileConfig(hmm, bg, gm, L, p7_GLOCAL);
	static int i(0);
    std::cout << i++ << std::endl;
    
	/* viterbi */
	mx = p7_gmx_Create(gm->M, L);
    try{
        p7_GViterbi(dsq, L, gm, mx, &score);
    }
    catch(...){
        std::cerr << "Error\n";
    }
	
	/* local clean up */
	p7_gmx_Destroy(mx);
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);
	free(dsq);
	
	return score;
}
