/* Functions to parse FASTA sequence files.
 *
 * Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
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

	infile.open(*file, ios::in );
	outfile.open("temp.rafa", ios::out );

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
				outfile << endl << line.substr(1,line.length()) << "\t";
				}
		}else{
			outfile << line;
			}
		}

	outfile << endl;

	infile.close();
	outfile.close();

	}


void write_fasta(char **seqs, char **file, int *num)
	{

	int i, ind;
	bool done;
	ofstream outfile;

	outfile.open(*file, ios::out );

	if(! outfile.is_open()){
		Rprintf("\nERROR: %s not found!\n",*file); return;
		}

	for(i = 0; i < *num; i++){
		if(i % 2 == 0){
			outfile << ">" << seqs[i] << endl;
		}else{
			string tseq(seqs[i]);
			ind = 0;
			done = FALSE;
			while(!done){
				if((ind+59) <= tseq.length()){
					outfile << tseq.substr(ind, 60) << endl;
				}else if(ind <= tseq.length()){
					outfile << tseq.substr(ind) << endl;
					done = TRUE;
				}else{
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





