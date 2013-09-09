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
std::vector<std::string> alphabet;
ESL_ALPHABET *abc    = NULL;    /* digitial alphabet */
P7_HMMFILE   *hfp    = NULL;    /* hmm file */
P7_HMM       *hmm    = NULL;    /* hmm */
char         *profile_fh;       /* HMMER Model File */
double       unparsed_score;
double       score;             /* HMMER score of edited(HMM) sequence */
static std::string unparsed_protein_seq;
static std::string protein_seq;



/*-----------------------------*/
/*  STOCHHMM LINKING FUNCTION  */
/*-----------------------------*/
double traceback_model_eval (const std::string *nucleotide_seq, size_t current_position, const std::string *edited_seq, size_t traceback){
    SEQanalysis compute;
    //std::cout << *edited_seq << std::endl;
    
    // HMMER score of HMM parsed seq 'edited_seq'
    protein_seq = compute.translate(edited_seq, aa_dictionary);
    unparsed_protein_seq = compute.translate(nucleotide_seq, aa_dictionary);
    
    score = compute.hmmer_score(abc, hmm, protein_seq);
    unparsed_score = compute.hmmer_score(abc,hmm, unparsed_protein_seq);
    std::cout << "Unparsed HMMER score:  " << unparsed_score << "\tParsed HMMER score: " << score << std::endl;
    
    // Adjust HMMER score by a scaling factor then return Sickle score
    
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
    std::string trans_tbl      (argv[1]);    /* Translation Table File */
    profile_fh =    argv[2];     /* Profile HMM file */
    std::string hmm_fh =        argv[3];     /* HMM files for StochHMM */
    std::string fasta_fh =      argv[4];     /* Sequence input (FASTA file) */
    
    
    /*  SETUP HMMER PROFILE    */
    if (p7_hmmfile_Open(profile_fh, NULL, &hfp) != eslOK)
		p7_Fail("Failed to open HMM file %s", profile_fh);
	
	if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK)
		p7_Fail("Failed to read HMM");
	p7_hmmfile_Close(hfp);
    
    
    /*  CREATE CODON DICTIONARY */
    alphabet.push_back("A");
    alphabet.push_back("C");
    alphabet.push_back("G");
    alphabet.push_back("T");
    StochHMM::track tr(alphabet);
    aa_dictionary = SEQanalysis::parse_translation_table(trans_tbl);
    
    
    /*  VITERBI TRACEBACK */
    StochHMM::StateFuncs myTransitions;
    StochHMM::model myModel;
    StochHMM::seqTracks jobs;
    std::vector<StochHMM::gff_feature> seq_coord;
    
    myTransitions.assignTransitionFunction("HMMER", *traceback_model_eval);
    myModel.import(hmm_fh, &myTransitions);
    jobs.loadSeqs(myModel, fasta_fh);
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
	p7_hmm_Destroy(hmm);
	esl_alphabet_Destroy(abc);

    return 0;
}

