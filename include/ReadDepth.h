#ifndef _SREADDEPTH_H_
#define _SREADDEPTH_H_


struct SReadDepth {
	SReadDepth(const int & i_pos, const unsigned int & i_count) : pos(i_pos), count(i_count){}
	int pos = 0;
	unsigned int count = 0;
};

#endif
