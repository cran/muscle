/* Functions to parse FASTA sequence files.
 *
 * Author: Alex T. Kalinka
 *
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include <R.h>


extern "C" {

using namespace std;


void read_fasta(char **file)
	{

	bool start = TRUE;
	string line, first;

	ifstream infile;
	ofstream outfile;

	infile.open(*file, ios::in | ios::binary);
	outfile.open("temp.rafa", ios::out | ios::binary);

	if(! infile.is_open()){
		Rprintf("\nERROR: %s not found!\n",*file); return;
		}

	while(getline(infile, line)){
		first = line.substr(0,1);
		if(first.compare(">") == 0){
			if(start){
				outfile << line.substr(1,line.length()) << "\t";
				start = FALSE;
			}else{
				outfile << "\n" << line.substr(1,line.length()) << "\t";
				}
		}else{
			outfile << line;
			}
		}

	outfile << endl;

	infile.close();
	outfile.close();

	}


void write_fasta(char **seqs, char **file, int *num, int *len)
	{

	int i, ind;
	bool done;
	ofstream outfile;

	outfile.open(*file, ios::out | ios::binary);

	if(! outfile.is_open()){
		Rprintf("\nERROR: %s not found!\n",*file); return;
		}

	for(i = 0; i < *num; i++){
		if(i % 2 == 0){
			outfile << ">" << seqs[i] << "\n";
		}else{
			string tseq(seqs[i]);
			ind = 0;
			done = FALSE;
			while(!done){
				if((ind+59) <= tseq.length()){
					outfile << tseq.substr(ind,60) << "\n";
				}else{
					outfile << tseq.substr(ind,tseq.length()) << "\n";
					done = TRUE;
					}
				ind = ind + 60;
				}
			tseq.clear();
			}
		}
	outfile << endl;

	outfile.close();

	}
			


}





