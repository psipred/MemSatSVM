// *********************************************************
// *    MEMSAT-SVM - MEMbrane protein Structure            *
// *    And Topology using Support Vector Machines         *
// * Integral Membrane Protein Topology Prediction Program *
// *    Copyright (C) 2008 David T. Jones and Tim Nugent   *
// *********************************************************
// 
// This program is copyright and may not be distributed without
// permission of the author unless specifically permitted under
// the terms of the license agreement.
// 
// THIS SOFTWARE MAY ONLY BE USED FOR NON-COMMERCIAL PURPOSES. PLEASE CONTACT
// THE AUTHOR IF YOU REQUIRE A LICENSE FOR COMMERCIAL USE.
// 
// Based on code by David T. Jones: August 6th 1992
// Last edit by Tim Nugent: September 2nd 2008
// 
// Description: This program uses dynamic programming to
// predict the secondary structure and topology of integral
// membrane proteins based on multiple sequence profiles.
// 

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

// Verbose
bool VERBOSE = false;

// Topology string format
unsigned int format = 1;

// Maximum number of predicted helices 
const unsigned int MAXNHEL = 50;

// Minimum length of loop 
int MINLLEN = 3;

// Minimum length of helix 
int MINHLEN = 16;

// Maximum length of helix 
int MAXHLEN = 31;

// Minimum length of re-entrant helix 
int MINREHLEN = 8;

// Maximum length of re-entrant helix 
int MAXREHLEN = 34;

// The minimum sequence length required for one predicted helix 
int HELDIV = MINHLEN + MINLLEN;

// Maximum length of "topogenic" loop 
int LIMITLOOP = 60;

// Minimum score for helix
int MINHSC = 220;

// Maximum length of sequence 
const int MAXSEQLEN = 4000;

// Minimum length of sequence 
const int MINSEQLEN = 20;

// Modify re-entrant helix SVM score
float RE_MOD = 0.7;

// Re-helix threshold
int MINRHSC = 178;

// Signal peptide threshold
double SIGNAL_THRESHOLD = 8.5;

// Minimum signal peptide cleavage position
int MIN_CLEAVAGE = 11; 

// Modifier for constrained residues
unsigned int CONSTRAINED_VALUE = 1000;

// Set to 0 with -s 0 if you know the sequence doesn't have a signal
int SIGNAL = 1;

// Big value
const long int BIG = 100000000;

// Containers etc.
float mat[MAXNHEL][MAXSEQLEN];
short length[50][MAXSEQLEN], path[50][MAXSEQLEN];
vector<char> seq;
vector<float> svm_hl,svm_i,svm_o,svm_re,svm_sp;

// Classes

// Topology class
class Topology{

        public:
                Topology(vector<int> top_vec, const string in_out, vector<float> re_vec, const float scr, int iflg, int seqlen){
			
			topology_vector = top_vec;
			re_helix_vector = re_vec;
			n_term = in_out;
			score = scr;
			inflg = iflg;
			length = seqlen;
			FillCompareToConstraintsVector();
			
                }
		
		void AlterScore(float new_value){
			score = new_value;
		}

		void AlterNterm(string new_n_term){
			n_term = new_n_term;
		}
			
		const float Score(){
			return(score);
		}	

		unsigned int FirstHelixBoundary(){
			return(topology_vector[0]);
		}
			
		vector<string> CompareToConstraintsVector(){
			return compare_to_constraints;
		}
		
		unsigned int HelixCount(){
			return(topology_vector.size()/2);
		}

		string Nterm(){
			return(n_term);
		}
			
		void FillCompareToConstraintsVector();	
				
		void Print();

                ~Topology(){
                }

        private:
	
                string topology_string,re_helix_string,n_term;
		vector<int> topology_vector;
		vector<float> re_helix_vector;
		vector<string> compare_to_constraints;
                float score;
		int inflg;
		unsigned int length;
		

};

// Fill a vector to compare to the contraints vector
void Topology::FillCompareToConstraintsVector(){

	for (unsigned int c = 0; c < length; c++){
		compare_to_constraints.push_back("x");
	}

	for (unsigned int i = 0; i < topology_vector.size(); i+= 2){
	
		for (int b = topology_vector[i]; b <= topology_vector[i+1]; b++){
			compare_to_constraints[b-1] = 'm';
		}	
	}	

	char pos;
	int in_helix = 0;
	if (n_term == "in"){
		pos = 'i';
	}else{
		pos = 'o';
	}
	
	for (unsigned int i = 0; i < length; i++){
		if (compare_to_constraints[i] == "x"){
			compare_to_constraints[i] = pos;
			if (in_helix == 1){
				in_helix = 0;
			}
			
		}else if (compare_to_constraints[i] == "m"){

			if (in_helix == 0){
				if (pos == 'i'){
					pos = 'o';
				}else{
					pos = 'i';
				}
				in_helix = 1;
			}
		}	
	}
}

// Print results in different formats
void Topology::Print(){

	char letters[26] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};

	// Create topology string
	for (unsigned int h = 0; h < topology_vector.size(); h += 2){
				
		// Convert int to string
		string s1;
		stringstream b1;
		b1 << topology_vector[h];
		s1 = b1.str();
		string s2;
		stringstream b2;
		b2 << topology_vector[h+1];
		s2 = b2.str();					
					
		// i6-21o40-57i142-172o214-236i276-302o
		if (format == 3){

			if (h != 0){
					
				if (n_term == "in"){
					if ((h/2) % 2){
						topology_string.append("o");
					}else{
						topology_string.append("i");
					}						
				}else{
					if ((h/2) % 2){
						topology_string.append("i");
					}else{
						topology_string.append("o");
					}	
				}
			}else{
					
				if (n_term == "in"){
					topology_string.append("i");
				}else{
					topology_string.append("o");
				}					
						
			}

			topology_string.append(s1);
			topology_string.append("-");
			topology_string.append(s2);

			if (h == topology_vector.size() - 2){
					
				if (n_term == "in"){
					if ((h/2) % 2){
						topology_string.append("i");
					}else{
						topology_string.append("o");
					}						
				}else{
					if ((h/2) % 2){
						topology_string.append("o");
					}else{
						topology_string.append("i");
					}	
				}
			}
				
		// A.6,21;B.40,57;C.142,172;D.214,236;E.276,302	
		}else if (format == 2){
				
			topology_string += letters[h/2];
			topology_string.append(".");
			topology_string.append(s1);
			topology_string.append(",");
			topology_string.append(s2);
			
			if (h != topology_vector.size() - 2) topology_string.append(";");
				
		// 6-21,40-57,142-172,214-236,276-302
		}else{
				
			topology_string.append(s1);
			topology_string.append("-");
			topology_string.append(s2);
			if (h != topology_vector.size() - 2) topology_string.append(",");				
	
		}
	}	

	// Create re-entrant helix string
	unsigned int h = 1;
	for (unsigned int r = 0;r < re_helix_vector.size();r += 4){

		// Convert int to string
		string s1;
		stringstream b1;
		b1 << re_helix_vector[r+1];
		s1 = b1.str();
		string s2;
		stringstream b2;
		b2 << re_helix_vector[r+2];
		s2 = b2.str();	

		
		if (inflg){
		
			if ((int)re_helix_vector[r] % 2){

				// i6-21o40-57i142-172o214-236i276-302o
				if (format == 3){

					re_helix_string.append("o");
					re_helix_string.append(s1);
					re_helix_string.append("-");
					re_helix_string.append(s2);
					re_helix_string.append("o");
				
				// A.6,21;B.40,57;C.142,172;D.214,236;E.276,302	
				}else if (format == 2){

					re_helix_string += letters[h-1];
					re_helix_string.append(".");				
					re_helix_string.append(s1);
					re_helix_string.append(",");
					re_helix_string.append(s2);
					if (r != re_helix_vector.size() - 4) re_helix_string.append(";");

				// 6-21,40-57,142-172,214-236,276-302
				}else{

					re_helix_string.append(s1);
					re_helix_string.append("-");
					re_helix_string.append(s2);
					if (r != re_helix_vector.size() - 4) re_helix_string.append(",");
				
				}
				
			}else{
				if (format == 3){

					re_helix_string.append("i");
					re_helix_string.append(s1);
					re_helix_string.append("-");
					re_helix_string.append(s2);
					re_helix_string.append("i");
				
				}else if (format == 2){

					re_helix_string += letters[h-1];
					re_helix_string.append(".");				
					re_helix_string.append(s1);
					re_helix_string.append(",");
					re_helix_string.append(s2);
					if (r != re_helix_vector.size() - 4) re_helix_string.append(";");

				}else{

					re_helix_string.append(s1);
					re_helix_string.append("-");
					re_helix_string.append(s2);
					if (r != re_helix_vector.size() - 4) re_helix_string.append(",");
				
				}
			}
		}else{
			if ((int)re_helix_vector[r] % 2){
			
				if (format == 3){

					re_helix_string.append("i");
					re_helix_string.append(s1);
					re_helix_string.append("-");
					re_helix_string.append(s2);
					re_helix_string.append("i");
				
				}else if (format == 2){

					re_helix_string += letters[h-1];
					re_helix_string.append(".");				
					re_helix_string.append(s1);
					re_helix_string.append(",");
					re_helix_string.append(s2);
					if (r != re_helix_vector.size() - 4) re_helix_string.append(";");

				}else{

					re_helix_string.append(s1);
					re_helix_string.append("-");
					re_helix_string.append(s2);
					if (r != re_helix_vector.size() - 4) re_helix_string.append(",");
				
				}
				
			}else{
			
				if (format == 3){

					re_helix_string.append("o");
					re_helix_string.append(s1);
					re_helix_string.append("-");
					re_helix_string.append(s2);
					re_helix_string.append("o");
				
				}else if (format == 2){

					re_helix_string += letters[h-1];
					re_helix_string.append(".");				
					re_helix_string.append(s1);
					re_helix_string.append(",");
					re_helix_string.append(s2);
					if (r != re_helix_vector.size() - 4) re_helix_string.append(";");

				}else{

					re_helix_string.append(s1);
					re_helix_string.append("-");
					re_helix_string.append(s2);
					if (r != re_helix_vector.size() - 4) re_helix_string.append(",");
				
				}
			}
		}
		
		h++;
	}
	
	if (re_helix_string == ("")) re_helix_string = ("Not detected.");
		
	cout << "Topology:\t\t" << topology_string << endl;
	cout << "Re-entrant helices:\t" << re_helix_string << endl;
	cout << "Helix count:\t\t" << topology_vector.size()/2 << endl;
	cout << "N-terminal:\t\t" << n_term << endl;
	cout << "Score:\t\t\t" << score << endl;

}

// Subroutines

// Print usage and exit
static void usage(){

	cout << endl << "MEMSAT-SVM: Alpha-helical transmembrane protein topology prediction" << endl;
	cout << "using Support Vector Machines." << endl;
	cout << endl << "Version 1.0" << endl;
	cout << endl << "Usage: memsat-svm [options] <SVM File>" << endl << endl;
	cout << "Options:" << endl << endl;
	cout << "  -c <filename>   Constrained prediction. The constraints file should have the" << endl;
	cout << "                  following format, where s,o,m,i are signal peptide, outside" << endl;
	cout << "                  loop, membrane and inside loop:" << endl;	
	cout << "                  s:   1-15" << endl;
	cout << "                  o:   1-30" << endl;	
	cout << "                  m:   37-59,82-100" << endl;	
	cout << "                  i:   65,80,220-230" << endl;	
	cout << "  -s <0|1>        Run with signal peptide function.        Default: 1" << endl;
	cout << "  -f <1|2|3>      Output format for topology string.       Default: 1" << endl;
	cout << "                  1 = 6-21,40-57,142-172,214-236,276-302" << endl;	
	cout << "                  2 = A.6,21;B.40,57;C.142,172;D.214,236;E.276,302" << endl;	
	cout << "                  3 = i6-21o40-57i142-172o214-236i276-302o" << endl;		
	cout << "  -m <int>        Minimum score for a transmembrane helix. Default: " << MINHSC << endl;
	cout << "  -r <int>        Minimum score for a re-entrant helix.    Default: " << MINRHSC << endl;
	cout << "  -h <int>        Show help." << endl;
	cout << endl;
	exit(1);

}

// Parse file containing raw SVM scores and put into vectors
static void parse_input_file(string& file, int & seqlen){

        // Convert string to file, open stream
        ifstream is(file.c_str());

        // Check to make sure the stream is ok
        if(!is.good()){
                cout << "Cannot open file : "<< file << endl;
                exit(1);
        }

        //cout << "Parsing SVM scores file " << file << " ... ";

        string line;
        char ch[50],s;
        float h,i,r,sig;

        // Get line from stream
        while(getline(is,line)){

                // copy the string to a character array
                strcpy(ch, line.c_str());


                if (sscanf(ch, "%c%f%f%f%f", &s, &h, &i, &r, &sig) != 5){
		
                        cout << "Bad SVM input format!" << endl;
			
                }else{

			r += RE_MOD;

			seq.push_back(s);
			svm_re.push_back(r);
			svm_hl.push_back(h);
			svm_i.push_back(i);
			svm_o.push_back(-i);
			svm_sp.push_back(sig);
			seqlen++;
                }  
        }

	// Close stream
	is.close();

        if (seqlen){
	
		//cout << "Done." << endl << endl;
                if (VERBOSE) cerr << "Length of sequence is " << seqlen << " residues." << endl << endl;
		
		if (seqlen > MAXSEQLEN){

			cout << "Maximum sequence length is " << MAXSEQLEN << "." << endl << endl;
			exit(1);
		
		}
		
		if (seqlen < MINSEQLEN){

			cout << "Minimum sequence length is " << MINSEQLEN << "." << endl << endl;
			exit(1);
		
		}

        }else{
	
		cout << "Couldn't find any SVM data!" << endl << endl;
		exit(1);
	}
	

}

// Functions to split strings and place contents in a vector
static void tokenize(const string& str, vector<string>& tokens, const string& delimiters = " "){
    
	// Skip delimiters at beginning
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	
	// Find first non-delimiter
	string::size_type pos = str.find_first_of(delimiters, lastPos);

 	while (string::npos != pos || string::npos != lastPos){
	
        	// Found a token, add it to the vector
        	tokens.push_back(str.substr(lastPos, pos - lastPos));
		
        	// Skip delimiters.  Note the not_of
        	lastPos = str.find_first_not_of(delimiters, pos);
		
        	// Find next non-delimiter
        	pos = str.find_first_of(delimiters, lastPos);
    	}
	
}

// Parse constraints file and modify SVM scores
static void parse_constraints_file(string& file,vector<string>& constraints_vector,int seqlen){

	// Open stream
        ifstream is(file.c_str());

        if(!is.good()){
                cout << "Cannot open file : "<< file << endl;
                exit(1);
        }
	cerr << file << endl;
        cout << "Parsing constraints file " << file << " ... ";

        string line;

	// Constraints file format is as follows. s = signal, o = outside loop, m/h = membrane, i = inside loop
	// s:   1-10
	// o:   1-15,30-80,120-140
	// m:   22,140-165
	// i:   100,166-192

        // Get line from stream
        while(getline(is,line)){

		vector<string> region;
		tokenize(line, region, ":");
		
		vector<string> constrined_residues;
		tokenize(region[1], constrined_residues, (","));
		
		vector<unsigned int> all_residues;
		
		for (unsigned int i = 0; i < constrined_residues.size(); i++){
		
			if (constrined_residues[i].find("-") < constrined_residues[i].length()){

				vector<string> range;
				tokenize(constrined_residues[i], range, ("-"));
				for (int r = atoi(range[0].c_str()); r <= atoi(range[1].c_str()); r++){
					if (r){
						if (r < seqlen) all_residues.push_back(r);
						//cerr << r << endl;
					}
				}
				
				if (VERBOSE) cerr << region[0] << ": Adding range: " << range[0] << " - " << range[1] << endl;
			
			}else if (atoi(constrined_residues[i].c_str())){

				if (VERBOSE) cerr << region[0] << ": Adding single residue: " << atoi(constrined_residues[i].c_str()) << endl;
				all_residues.push_back(atoi(constrined_residues[i].c_str()));

			}	
		
		}

		// Vector all_residues now contains all residues to be constrained to a particular region
		// Modify four vectors containing SVM scores

		if (region[0] == "s" || region[0] == "S"){
			if (VERBOSE) cerr << "Signal:" << endl;
			for (unsigned int a = 0; a < all_residues.size(); a++){
			
				svm_sp[all_residues[a]-1] += CONSTRAINED_VALUE;					
				if (VERBOSE) cerr << all_residues[a] << endl;	
				constraints_vector[all_residues[a]-1] = "o";			
			}
		}else if (region[0] == "o" || region[0] == "O"){
			if (VERBOSE) cerr << "Outside:" << endl;
			for (unsigned int a = 0; a < all_residues.size(); a++){	
			
				svm_o[all_residues[a]-1] += CONSTRAINED_VALUE;				
				if (VERBOSE) cerr << all_residues[a] << endl;	
				constraints_vector[all_residues[a]-1] = "o";	
						
			}		
		}else if (region[0] == "m" || region[0] == "M" || region[0] == "h" || region[0] == "H"){
			if (VERBOSE) cerr << "Membrane:" << endl;
			for (unsigned int a = 0; a < all_residues.size(); a++){	
					
				svm_hl[all_residues[a]-1] += CONSTRAINED_VALUE;		
				if (VERBOSE) cerr << all_residues[a] << endl;
				constraints_vector[all_residues[a]-1] = "m";	
							
			}		
		}else if (region[0] == "i" || region[0] == "I"){
			if (VERBOSE) cerr << "Inside:" << endl;
			for (unsigned int a = 0; a < all_residues.size(); a++){	
			
				svm_i[all_residues[a]-1] += CONSTRAINED_VALUE;			
				if (VERBOSE) cerr << all_residues[a] << endl;
				constraints_vector[all_residues[a]-1] = "i";	
							
			}		
		}
        }

	// Close stream
	is.close();

	cout << "Done." << endl << endl;

}

// Deal with signal peptides
static void process_signal_peptide(float & sig_total, int & cleavage){

	// Only process the first 30 residues
	for (unsigned int s = 0; s < 40; s++){

		if (svm_sp[s] > 0){
				
			if (VERBOSE) cerr << "Adding SP score " << s << " - " << svm_sp[s] << endl;
			
			// Scale up signal peptide score by 10
			svm_hl[s] -= 10 * svm_sp[s];
			svm_i[s]  -= svm_sp[s];
			svm_o[s]  += svm_sp[s];
			sig_total += svm_sp[s];

			// The residue with the highest signal peptide SVM score > 0 will be the cleavage point
			if (s > 0){
				if (svm_sp[s] > svm_sp[s-1]){
					cleavage = s+1;
					if (VERBOSE) cerr << "Signal peptide cleavage point: " << cleavage << endl;
				}				
			}
		}
	}	
}

// Generate TM helix score
float helscore(int start, int length, int inflg){

	float score = 0,helix_mod = 100,neg_helix_mod = 50;

        for (int i = start; i <= start + length; i++){

		if (svm_hl[i-1] > 0){
			score += (helix_mod * svm_hl[i-1]);
 		}else{
			score += (neg_helix_mod * svm_hl[i-1]);
		}
        }

	for (int i = 0; i < 8; i++){
	
		if (inflg){
	  		score += svm_i[start + i];
			score += svm_o[start + length - i - 1];
		}else{
	   		score += svm_o[start + i];
			score += svm_i[start + length - i - 1];
		}
    	}

    	return(score);
}

// Generate re-entrant helix score
float rehelscore(int start, int length){

	float score = 0,re_helix_mod = 100;

        for (int i = start; i <= start + length; i++){

		if (svm_re[i-1] > 0){
			score += (re_helix_mod * svm_re[i-1]);
		}
        }

    	return(score);
}

// Search for re-entrant helices in the loops
vector<float> find_re_helices(vector<int>& topology_array,int inflg){

	float re_hsc,last_re_hsc;
	int start, stop, h = 0;
	vector<float> re_topology_vector;

	if (topology_array.size() > 1){
	
		for (unsigned int s = 1;s < topology_array.size() - 2;s += 2){

			start = topology_array[s] + 1;
			stop =  topology_array[s+1] - 1;

			for (int i = start + MINLLEN; i < stop - MINLLEN; i++){

				for (int l = MINREHLEN; l <= MAXREHLEN; l++){

					re_hsc = rehelscore(i, l);
					if (i + l > stop - MINLLEN){
						re_hsc = -BIG;
					}
					
					if (re_hsc >= MINRHSC){
		
						if (re_topology_vector.size()){
						
							// If this re-entrant helix overlaps with the previous one
							if (i <= re_topology_vector[re_topology_vector.size()-2]){
							
								//And if the new score is higher than the last
								if (re_hsc > last_re_hsc){
								
									// Replace the last one
									re_topology_vector[re_topology_vector.size()-4] = (s+1)/2;
									re_topology_vector[re_topology_vector.size()-3] = i;
									re_topology_vector[re_topology_vector.size()-2] = i+l;
									re_topology_vector[re_topology_vector.size()-1] = re_hsc;
									last_re_hsc = re_hsc;
								}
															
							}else{
							
								// Otherside add a new one
								re_topology_vector.push_back((s+1)/2);
								re_topology_vector.push_back(i);
								re_topology_vector.push_back(i+l);
								re_topology_vector.push_back(re_hsc);
								last_re_hsc = re_hsc;
							}
						}else{		
						
							// First helix	
							re_topology_vector.push_back((s+1)/2);			
							re_topology_vector.push_back(i);
							re_topology_vector.push_back(i+l);
							re_topology_vector.push_back(re_hsc);
							last_re_hsc = re_hsc;
							
						}

					}else{
						//cout << "Below threshold: " << i << " - " << i+l << " score " << re_hsc << endl;
					}
				}	
			}
		}
	}

	if (VERBOSE) cerr << re_topology_vector.size()/4 << " re-entrant helices found." << endl;
	
	for (unsigned int r = 0;r < re_topology_vector.size();r += 4){

		h++;

		if (inflg){
			if ((int)re_topology_vector[r] % 2){
				cout << "Re-entrant helix " << h << " from " << re_topology_vector[r+1] << " (out) to "
				<< re_topology_vector[r+2] << " (out) :     \t" << re_topology_vector[r+3] << endl;
			}else{
				cout << "Re-entrant helix " << h << " from " << re_topology_vector[r+1] << " (in) to "
				<< re_topology_vector[r+2] << " (in) :     \t" << re_topology_vector[r+3] << endl;
			}
		}else{
			if ((int)re_topology_vector[r] % 2){
				cout << "Re-entrant helix " << h << " from " << re_topology_vector[r+1] << " (in) to "
				<< re_topology_vector[r+2] << " (in) :     \t" << re_topology_vector[r+3] << endl;
			}else{
				cout << "Re-entrant helix " << h << " from " << re_topology_vector[r+1] << " (out) to "
				<< re_topology_vector[r+2] << " (out) :     \t" << re_topology_vector[r+3] << endl;
			}
		}

	}
	
	return(re_topology_vector);
}	
	

// Calculate score/path matrix by dynamic programming 
void calcmat(int inflg, int nhelix, unsigned int seqlen){

	int i, j, h, k, l, maxj, maxl;
	float maxsc, lsc, intog, hsc;

	// Toggle inflag to 1/0
	intog = inflg ^ (nhelix & 1);
	
	// For each helix between 20 & 1
	for (h = nhelix - 1; h >= 0; h--, intog = !intog){
		
		// Fill the matrix for each residue in the seq, for this helix number
		for (i = 0; i < (int)seqlen; i++){
	   		 mat[h][i] = -BIG;
		}
	 
	 	// Work back from end
		for (i = seqlen - MINHLEN - MINLLEN; i >= MINLLEN; i--){
		
			maxsc = -BIG;
			maxj = maxl = 0;
			
			// Generate all possible helices (length l) for this residue i
			for (l = MINHLEN; l <= MAXHLEN; l++){
		
				// As long as we're within the length of the sequence
				if (i + l + MINLLEN <= (int)seqlen){
				
					// Generate a score for this helix
		    			hsc = helscore(i, l, !intog);

					// If this is the C-terminal helix					
		    			if (h == nhelix - 1){
					
						// Generate the C-terminal loop score, from the end of the helix to seqlen
						lsc = 0;
						for (k = i + l; k < (int)seqlen; k++){
			    				if (intog){
								// Add the inside score
								lsc += svm_i[k];
			    				}else{
								// Add the outside score
								lsc += svm_o[k];
							}
						}

						// Loop score = 0 if length > 60
						if (k - i - l + 1 > LIMITLOOP){
							//lsc -= intog ? svm_i[k] : svm_o[k];
			    				lsc = 0;
						}

						// hsc + lsc = maxsc if this is the highest scoring helix found so far
						if (hsc + lsc > maxsc){
			   		 		maxsc = hsc + lsc;
			    				maxl = l;
						}
					
					// For all other helices						
					}else{
	
						// Calculate initial loop score 
						lsc = 0;
						for (k = i + l; k < i + l + MINLLEN - 1; k++){
							lsc += intog ? svm_i[k] : svm_o[k];
						}
						
						// Between helix end + min loop length and seqleng - 1 loop and 1 helix
						// E.g. finds optimal helix and following loop score
						// And records loop end as maxj
						for (j = i + l + MINLLEN; j < (int)seqlen - MINHLEN - MINLLEN; j++){
		
							lsc += intog ? svm_i[j - 1] : svm_o[j - 1];
		
							// extend loop max 60 residues
							if (j - i - l > LIMITLOOP){
								//lsc -= intog ? svm_i[j - 1] : svm_o[j - 1];
								lsc = 0;
							}

							if (hsc + lsc + mat[h + 1][j] > maxsc){
								maxsc = hsc + lsc + mat[h + 1][j];
								maxl = l;
								maxj = j;
			    				}
						}
					}				
				}
			}

			mat[h][i] = maxsc;
    			length[h][i] = maxl;
    			path[h][i] = maxj;
		}
	}
}

// Generate the N-terminal loop score
void firstlp(int starthel, int inflg, int nhelix, unsigned int seqlen){

	float lsc = 0;

	for (int i = 0; i < (int)seqlen - (nhelix - starthel) * (MINHLEN + MINLLEN); i++){	

		lsc += inflg ? svm_i[i] : svm_o[i];

    		if (i+1 > LIMITLOOP){
    	    		lsc = 0;
	    	}

    		mat[starthel][i+1] += lsc;
    	}
}

// Trace back highest scoring path through matrix 
float trace_back(int starthel, int inflg, int nhelix, vector<Topology>& topologies, unsigned int seqlen){

	int i, h, res = 1, intog = inflg;
	bool weak_h = false;
	float maxsc = -BIG, hsc, helix_total, score;

	if (nhelix-starthel  == 1){
		cout << "Processing 1 helix:" << endl;
	}else{
		cout << "Processing " << nhelix-starthel << " helices:" << endl;
	}

	// Find the best helix start position = res + 1
	// The range to search changes as more helices are added
    	for (i = 0; i < (int)seqlen - (nhelix - starthel) * (MINHLEN + MINLLEN); i++){

		//These values generated in firstlp
		if (mat[starthel][i] > maxsc){
	   		maxsc = mat[starthel][i];
	    		res = i;
		}
		
	}
	
	if (VERBOSE) cerr << "maxsc = " << maxsc << endl << "res = " << res << endl;

	string n_term;
	vector<int> topology_vector;
	vector<float> re_helix_vector;

    	for (h = 0; h < nhelix-starthel; h++, intog = !intog){

	    	if (intog){
		
			printf("Transmembrane helix %d from %3d (in) to %3d (out) :\t", h + 1, res+1, res + length[h+starthel][res]);
			hsc = helscore(res, length[h+starthel][res], intog);
			cout << hsc << endl;
			helix_total += hsc;
	    	
			topology_vector.push_back(res+1);
			topology_vector.push_back(res + length[h+starthel][res]);
			
			if (h == 0){
				n_term = "in";
			}	
		
		}else{
		
			printf("Transmembrane helix %d from %3d (out) to %3d (in) :\t", h + 1, res+1, res + length[h+starthel][res]);
			hsc = helscore(res, length[h+starthel][res], intog);
			cout << hsc << endl;
			helix_total += hsc;
			
			topology_vector.push_back(res+1);
			topology_vector.push_back(res + length[h+starthel][res]);

			if (h == 0){
				n_term = "out";
			}
		}

		if (hsc < MINHSC){
	    		weak_h = true;
		}
		
		// go to the next row/res in the matrix
		res = path[h+starthel][res];

    	}

	re_helix_vector = find_re_helices(topology_vector,inflg);
	score = (nhelix == 1) ? maxsc/1000.0 : (weak_h ? -BIG/1000.0 : maxsc/1000.0);
	topologies.push_back(Topology(topology_vector,n_term,re_helix_vector,score,inflg,seqlen));
    	return (score);
	
}

int main(int argc, const char* argv[]){

        // Exit with usage unless at least filename given as argument
        if (argc < 2){
                usage();
        }
	
	string input_file = "", constraints_file = "";
        int i = 1, cleavage = 1, seqlen = 0;
	float sig_total = 0;
        
	// Process command line arguments
        while( i < argc){

                if( argv[i][0] == '-'){
    
                        i++;
                        switch(argv[i-1][1]){

                                // Put some more argument options in here later
				case 'c' : {constraints_file = argv[i]; break;}
				case 'f' : {format = atoi(argv[i]); break;}
				case 'm' : {MINHSC = atoi(argv[i]); break;}
				case 'r' : {MINRHSC = atoi(argv[i]); break;}
				case 's' : {SIGNAL = atoi(argv[i]); break;}
				case 'v' : {VERBOSE = argv[i]; break;}
                                case 'h' : {usage();}
                                default:   {usage();}
                        }
   
                }else{
                        input_file = argv[i];
                }
                i++;
        }
	
	// Parse the input file, fill vectors with SVM scores, get sequence length
        parse_input_file(input_file, seqlen);
	
	// If constraints file is parsed, modify SVM score vectors
	vector<string> constraints_vector;
	for (int c = 0; c < seqlen; c++){
		constraints_vector.push_back("x");
	}

	if (constraints_file != "") parse_constraints_file(constraints_file,constraints_vector,seqlen);

	
	// Process signal peptide, get total signal score and cleavage position
	if (seqlen >= 40 && SIGNAL == 1) process_signal_peptide(sig_total, cleavage);
	
	vector<Topology> topologies;
 	float  maxsc = -BIG, finalsc[50][2], score;
	int maxt, maxnh, h, nhelix = 1; 
	bool inflg = true;
	
	// 1 helix in
    	calcmat(inflg = true, nhelix, seqlen);
    	firstlp(0, inflg, nhelix, seqlen);
    	score = trace_back(0, inflg, nhelix, topologies, seqlen);
    	if (score > maxsc){
    		//memcpy(bestst, structure, sizeof(structure));
    		maxsc = score;
    		maxt = inflg;
    		maxnh = 1;
    	}
    	finalsc[0][inflg] = score;
	printf("Score = %f\n", score);
	cout << endl;

	// 1 helix out
   	calcmat(inflg = false, nhelix, seqlen);
    	firstlp(0, inflg, nhelix, seqlen);
    	score = trace_back(0, inflg, nhelix, topologies, seqlen);
    	if (score > maxsc){
    		//memcpy(bestst, structure, sizeof(structure));
    		maxsc = score;
    		maxt = inflg;
    		maxnh = 1;
    	}
    	finalsc[0][inflg] = score;
	printf("Score = %f\n", score);
	cout << endl;

	// Work out how many helices we can fit in the sequence
	nhelix = max(1,int(min((float)MAXNHEL,(float)seqlen/HELDIV)));

    	calcmat(inflg = true, nhelix, seqlen);
    	for (h=0; h<nhelix-1; h++, inflg = !inflg){
	
    		firstlp(h, inflg, nhelix, seqlen);
    		score = trace_back(h, inflg, nhelix, topologies, seqlen);
		if (score > maxsc){
	    		//memcpy(bestst, structure, sizeof(structure));
	    		maxsc = score;
	    		maxt = inflg;
	    		maxnh = nhelix-h;
		}
		finalsc[nhelix-h-1][inflg] = score;
	        printf("Score = %f\n\n", score);
    	}

    	calcmat(inflg = false, nhelix, seqlen);
    	for (h=0; h<nhelix-1; h++, inflg = !inflg){

    		firstlp(h, inflg, nhelix, seqlen);
    		score = trace_back(h, inflg, nhelix, topologies, seqlen);
		if (score > maxsc){
	   		//memcpy(bestst, structure, sizeof(structure));
	    		maxsc = score;
	    		maxt = inflg;
	    		maxnh = nhelix-h;
		}
		finalsc[nhelix-h-1][inflg] = score;
	    	printf("Score = %f\n\n", score);
    	}

	float best_score = -BIG;
	int best_topology;
	for (unsigned int i = 0; i < topologies.size(); i++){
	
		//topologies[i].Print();
		
		// Score down any topologies that don't work with constraints
		if (constraints_file != ""){
		
			vector<string> compare = topologies[i].CompareToConstraintsVector();
			for (unsigned int c = 0; c < compare.size(); c++){
			
				if (constraints_vector[c] != "x"){
				
					if (constraints_vector[c] != compare[c]){
									
						topologies[i].AlterScore(-100000);
						
						if (topologies[i].Nterm() == "in"){
							finalsc[topologies[i].HelixCount()-1][1] = -100000;
						}else{
							finalsc[topologies[i].HelixCount()-1][0] = -100000;
						}						
					}
				}			
			}		
		}		
	
	 	if (topologies[i].Score() > best_score){

			best_score = topologies[i].Score();
			best_topology = i;
		}
	}

    	printf("Summary of topology analysis:\n");
    	for (h = 1; h <= nhelix; h++){
		if (h == 1){
	   	 	printf("%2d helix   (+) : Score = %g\n", h, finalsc[h-1][1]);
		}else{
	    		printf("%2d helices (+) : Score = %g\n", h, finalsc[h-1][1]);
		}		
		if (h == 1){
	    		printf("%2d helix   (-) : Score = %g\n", h, finalsc[h-1][0]);
		}else{
	    		printf("%2d helices (-) : Score = %g\n", h, finalsc[h-1][0]);
		}
   	}
	cout << endl;
	
	cout << "Results:" << endl;
	if (best_score > -100000){

		// If signal score is above threshold, force the N-terminal to be outside
		if (topologies[best_topology].Nterm() == "in" && sig_total >= SIGNAL_THRESHOLD && SIGNAL == 1){
			topologies[best_topology].AlterNterm("out");
		}

		if (SIGNAL == 1 && sig_total > 0 && topologies[best_topology].Nterm() == "out"){
		
			// Make sure cleavage point isn't after 1st helix
			if ((int)cleavage > (int)topologies[best_topology].FirstHelixBoundary() - MINLLEN){			 
				cleavage = topologies[best_topology].FirstHelixBoundary() - MINLLEN;			 
			}
		
			if (cleavage > MIN_CLEAVAGE){
				if (format == 3){			
					cout << "Signal peptide:\t\to1-" << cleavage << "o" << endl;				
				}else if (format == 2){			
					cout << "Signal peptide:\t\tA.1," << cleavage << endl;				
				}else{			
					cout << "Signal peptide:\t\t1-" << cleavage << endl;
				}					
				cout << "Signal score:\t\t" << sig_total << endl;
				
			}else{
				cout << "Signal peptide:\t\tNot detected." << endl;
				cout << "Signal score:\t\t" << sig_total << endl;
			}			
		}else{		
			cout << "Signal peptide:\t\tNot detected." << endl;
			cout << "Signal score:\t\t" << sig_total << endl;			
		}
		
		topologies[best_topology].Print();
  		cout << endl;
		
	}else{
		cout << "No transmembrane helices predicted!" << endl;
	}
  
	exit(1);

}


