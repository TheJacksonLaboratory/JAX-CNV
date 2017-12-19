#ifndef _DATASTRUCT_H_
#define _DATASTRUCT_H_


struct SReadDepth {
	SReadDepth(const int & i_pos, const unsigned int & i_count) : pos(i_pos), count(i_count){}
	int pos = 0;
	unsigned int count = 0;
	unsigned int n_count = 0; // how many N's in a region of reference genome.
};

struct SHmmStats {
	SHmmStats(){}
	SHmmStats(const std::string & a, const unsigned int b, const unsigned int c, const unsigned int d): chr(a), pos(b), stats(c), length(d){}
	SHmmStats(const unsigned int a, const unsigned int b, const unsigned int c): pos(a), stats(b), length(c){}
	std::string chr;
	unsigned int pos = 0;
	unsigned int stats = 3;
	unsigned int length = 0;
};

struct SHmmStatsHeap {
	SHmmStatsHeap(){}
	SHmmStatsHeap(const SHmmStats & a, const unsigned int & b) : hmm_stats(a), id(b){}
	SHmmStats hmm_stats;
	unsigned int id = 0;
	bool merged = false;
};

#endif
