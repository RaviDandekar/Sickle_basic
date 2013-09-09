//
//  main.cpp
//  Sickle_basic
//
//  Created by Ravi Dandekar on 9/5/13.
//  Copyright (c) 2013 Korf Lab. All rights reserved.
//

#include "sickle_processing.h"

/*  GLOBAL variables    */
std::map<std::string, SEQanalysis::AMINO_ACIDS> aa_dictionary;
         /* stocchmm alphabet */
static ESL_ALPHABET *ALPHABET   = NULL;    /* hmmer digitial alphabet */
static P7_HMM       *PROFILE    = NULL;    /* hmm */


/*-----------------------------*/
/*  STOCHHMM LINKING FUNCTION  */
/*-----------------------------*/
double traceback_model_eval (const std::string *genome, size_t pos, const std::string *mRNA, size_t traceback) {
    SEQanalysis compute;
    
    /*  HMMER score of HMM parsed seq  */
    std::string protein_seq = compute.translate(mRNA, aa_dictionary);
    double score            = compute.hmmer_score(ALPHABET, PROFILE, protein_seq);
    
    std::cout << "HMMER score:  " << score << std::endl;
    
    return 0.01;
}


/*=============*/
/*  MAIN CODE  */
/*=============*/

std::string usage = "usage: Sickle_basic <translation table> <profile> <HMM> <FASTA seq>";
int main(int argc, char* const argv[]){
    
    /*  INPUT   */
    if (argc != 5){
        std::cout << usage << std::endl;
        exit(2);
    }
    std::string trans_tbl     = argv[1];        /* Translation Table File */
    char *profile_file        = argv[2];        /* Hmmer (profile) HMM file */
    std::string stochhmm_file = argv[3];        /* StochHMM HMM file*/
    std::string fasta_file    = argv[4];        /* Sequence input (FASTA file) */
    P7_HMMFILE  *hfp;
    
    
    /*  SETUP HMMER PROFILE    */
    if (p7_hmmfile_Open(profile_file, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", profile_file);
	if (p7_hmmfile_Read(hfp, &ALPHABET, &PROFILE) != eslOK) p7_Fail("Failed to read HMM");
	p7_hmmfile_Close(hfp);
    
    
    /*  CREATE CODON DICTIONARY */
    std::vector<std::string> alphabet;
    alphabet.push_back("A");
    alphabet.push_back("C");
    alphabet.push_back("G");
    alphabet.push_back("T");
    StochHMM::track tr = alphabet;
    aa_dictionary = SEQanalysis::parse_translation_table(trans_tbl);
    
    /*  VITERBI TRACEBACK */
    StochHMM::StateFuncs myTransitions;
    StochHMM::model myModel;
    StochHMM::seqTracks jobs;
    std::vector<StochHMM::gff_feature> seq_coord;
    
    myTransitions.assignTransitionFunction("HMMER", *traceback_model_eval);
    myModel.import(stochhmm_file, &myTransitions);
    jobs.loadSeqs(myModel, fasta_file);
    StochHMM::seqJob *job = jobs.getJob();
    
    while (job != NULL) {
        // Perform Viterbi
        StochHMM::trellis trell(&myModel, job->getSeqs());
        trell.viterbi();
        
        // Traceback
        StochHMM::traceback_path tb(&myModel);
        trell.traceback(tb);
        
        // Print GFF Output
        tb.gff(seq_coord, job->getSeqs()->getHeader());
        for(size_t gff=0;gff<seq_coord.size();gff++){
            std::cout << seq_coord[gff].strand << "\t" << seq_coord[gff].feature << "\t" << seq_coord[gff].start << "\t" << seq_coord[gff].end << std::endl;
        }
        
        job = jobs.getJob();
        seq_coord.clear();
    }
    
    
    /* GLOBAL CLEANUP */
	p7_hmm_Destroy(PROFILE);
	esl_alphabet_Destroy(ALPHABET);

    return 0;
}

