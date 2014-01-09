#include<cstdio>
#include<cctype>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<iomanip>
#include<string>
#include<vector>
#include<algorithm>
using namespace std;

typedef unsigned int INTEGER;

const int STRLEN = 1005;

INTEGER M;
int k; // k-mer size
vector<string> names;
vector<string> seqs;
vector<INTEGER> effL;

// tid starts from 1
struct ReadType {
  INTEGER tid, pos;

  ReadType(INTEGER tid, INTEGER pos) {
    this->tid = tid;
    this->pos = pos;
  }

  bool operator< (const ReadType& o) const {
    string& a = seqs[tid];
    string& b = seqs[o.tid];
    for (int i = 0; i < k; i++) {
      if (a[pos + i] != b[o.pos + i]) {
	return a[pos + i] < b[o.pos + i];
      }
    }
    return tid < o.tid;
  }

  bool seq_equal(const ReadType& o) const {
    string& a = seqs[tid];
    string& b = seqs[o.tid];
    for (int i = 0; i < k; i++) 
      if (a[pos + i] != b[o.pos + i]) return false;
    return true;
  }
};

vector<ReadType> cands;
vector<double> clusteringInfo; 

string convert(const string& rawseq) {
  int size = (int)rawseq.size();
  string seq = rawseq;
  for (int i = 0; i < size; i++) {
    seq[i] = toupper(rawseq[i]);
    if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') seq[i] = 'N';
  }
  return seq;
}

void loadRef(char* inpF) {
  ifstream fin(inpF);
  string tag, line, rawseq;

  assert(fin.is_open());

  names.clear(); names.push_back("");
  seqs.clear(); seqs.push_back("");
  
  getline(fin, line);
  while ((fin) && (line[0] == '>')) {
    tag = line.substr(1);
    rawseq = "";
    while((getline(fin, line)) && (line[0] != '>')) {
      rawseq += line;
    }
    if (rawseq.size() <= 0) {
      printf("Warning: Fasta entry %s has an empty sequence! It is omitted!\n", tag.c_str());
      continue;
    }
    names.push_back(tag);
    seqs.push_back(convert(rawseq));
  }

  fin.close();

  M = names.size() - 1;

  printf("The reference is loaded.\n");
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    printf("Usage: rsem-for-ebseq-calculate-clustering-info k input_reference_fasta_file output_file\n");
    exit(-1);
  }

  k = atoi(argv[1]);
  loadRef(argv[2]);

  cands.clear();
  effL.assign(M + 1, 0);
  for (INTEGER i = 1; i <= M; i++) {
    effL[i] = seqs[i].length() - k + 1;
    if (effL[i] <= 0) effL[i] = 0; // effL should be non-negative
    for (INTEGER j = 0; j < effL[i]; j++) 
      cands.push_back(ReadType(i, j));
  }
  printf("All possbile %d mers are generated.\n", k);

  sort(cands.begin(), cands.end());
  printf("All %d mers are sorted.\n", k);
 
  size_t p = 0;
  clusteringInfo.assign(M + 1, 0.0);

  for (size_t i = 1; i <= cands.size(); i++)
    if (i == cands.size() || !cands[p].seq_equal(cands[i])) {
      size_t denominator = i - p;
      size_t q = p; 
      for (size_t j = p + 1; j <= i; j++)
	if (j == i || cands[q].tid != cands[j].tid) {
	  size_t numerator = j - q;
	  //double prob = numerator * 1.0 / denominator;
	  //clusteringInfo[cands[q].tid] += (double)numerator * prob * (1.0 - prob);
	  if (numerator < denominator) clusteringInfo[cands[q].tid] += numerator;
	  q = j;
	}
      p = i;
    }

  for (INTEGER i = 1; i <= M; i++) 
    if (effL[i] == 0) clusteringInfo[i] = -1.0;
    else clusteringInfo[i] /= effL[i];

  printf("Clustering information is calculated.\n");


  ofstream fout(argv[3]);
  for (INTEGER i = 1; i <= M; i++) fout<<names[i]<<"\t"<<setprecision(6)<<clusteringInfo[i]<<endl;
  fout.close();

  return 0;
}
