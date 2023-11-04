/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "bridger.h"
#include "config.h"
#include "util.h"

int entry::print()
{
	printf("entry: length = %d, trace = (%d, %d), stack = (", length, trace1, trace2);
	printv(stack);
	printf(")\n");
	return 0;
}

bool entry_compare(const entry &x, const entry &y)
{
	for(int i = 0; i < x.stack.size() && i < y.stack.size(); i++)
	{
		if(x.stack[i] > y.stack[i]) return true;
		if(x.stack[i] < y.stack[i]) return false;
	}
	if(x.length < y.length) return true;
	else return false;
}

bridger::bridger(bundle_bridge *b)
{
	bd = b;
	max_pnode_length = 50;
}

int bridger::bridge_normal_fragments()
{
	/*
	printf("before bridging ... \n");
	for(int i = 0; i < bd->fragments.size(); i++)
	{
		bd->fragments[i].print(i);
	}
	printf("===\n");
	*/

	set_normal_length();
	int n = bd->fragments.size();

	bridge_overlapped_fragments(bd->fragments);
	filter_paths(bd->fragments);
	int n1 = get_paired_fragments(bd->fragments);

	vector<fcluster> open_fclusters;
	cluster_open_fragments(open_fclusters, bd->fragments);

	// 1st round of briding hard fragments
	build_junction_graph(bd->fragments);
	bridge_hard_fragments_normal(open_fclusters);
	filter_paths(bd->fragments);
	int n2 = get_paired_fragments(bd->fragments);

	// 2nd round of briding hard fragments
	build_junction_graph(bd->fragments);
	bridge_hard_fragments_normal(open_fclusters);
	filter_paths(bd->fragments);
	int n3 = get_paired_fragments(bd->fragments);

	// recluster open fragments
	open_fclusters.clear();
	cluster_open_fragments(open_fclusters, bd->fragments);
	bridge_phased_fragments(open_fclusters);
	filter_paths(bd->fragments);
	int n4 = get_paired_fragments(bd->fragments);

	double r1 = n1 * 100.0 / n;
	double r2 = n2 * 100.0 / n;
	double r3 = n3 * 100.0 / n;
	double r4 = n4 * 100.0 / n;

	vector<int> ct = get_bridged_fragments_type(bd->fragments);	// ct<ct1, ct2, ct3> paired-end, UMI-linked, both
	if(verbose >= 1)
	{
		printf("#normal fragments = %d, #fixed = %d -> %d -> %d -> %d, ratio = %.2lf -> %.2lf -> %.2lf -> %.2lf, #remain = %d, length = (%d, %d, %d), total paired-end = %d, UMI-linked only = %d, intersection: %d, bridged paired-end = %d, UMI-linked only = %d, intersection: %d\n", 
				n, n1, n2, n3, n4, r1, r2, r3, r4, n - n4, length_low, length_median, length_high, ct[3], ct[4], ct[5], ct[0], ct[1], ct[2]);
	}

	/*
	printf("after bridging ... \n");
	for(int i = 0; i < bd->fragments.size(); i++)
	{
		bd->fragments[i].print(i);
	}
	printf("===\n");
	*/

	return 0;
}

int bridger::bridge_circ_fragments()
{
	/*
	printf("before bridging ... \n");
	for(int i = 0; i < bd->fragments.size(); i++)
	{
		bd->fragments[i].print(i);
	}
	printf("===\n");
	*/

	set_normal_length();
	int n = bd->circ_fragments.size();

	bridge_overlapped_fragments(bd->circ_fragments);
	int n1 = get_paired_fragments(bd->circ_fragments);

	vector<fcluster> open_fclusters;
	cluster_open_fragments(open_fclusters, bd->circ_fragments);

	build_junction_graph(bd->fragments);
	bridge_hard_fragments_circ(open_fclusters);

	int n2 = get_paired_fragments(bd->circ_fragments);

	bridge_phased_fragments(open_fclusters);
	int n3 = get_paired_fragments(bd->circ_fragments);
	int n4 = n3;

	pick_bridge_path(bd->circ_fragments);

	double r1 = n1 * 100.0 / n;
	double r2 = n2 * 100.0 / n;
	double r3 = n3 * 100.0 / n;
	double r4 = n4 * 100.0 / n;

	vector<int> ct = get_bridged_fragments_type(bd->circ_fragments);	// ct<ct1, ct2, ct3> paired-end, UMI-linked, both
	if(verbose >= 1)
	{
		printf("#CIRC fragments = %d, #fixed = %d -> %d -> %d -> %d, ratio = %.2lf -> %.2lf -> %.2lf -> %.2lf, #remain = %d, length = (%d, %d, %d), total paired-end = %d, UMI-linked only = %d, intersection: %d, bridged paired-end = %d, UMI-linked only = %d, intersection: %d\n", n, n1, n2, n3, n4, r1, r2, r3, r4, n - n4, length_low, length_median, length_high, ct[3], ct[4], ct[5], ct[0], ct[1], ct[2]);
	}
	/*
	printf("after bridging ... \n");
	for(int i = 0; i < bd->fragments.size(); i++)
	{
		bd->fragments[i].print(i);
	}
	printf("===\n");
	*/

	return 0;
}

int bridger::bridge_clip(int32_t p1, int32_t p2, circular_transcript &circ)
{
	// locate the regions for p1 and p2
	int x1 = -1, x2 = -1; //region ids of boundaries that match p1 and p2
	for(int i=0;i<bd->regions.size();i++)
	{
		region r = bd->regions[i];
		if(r.lpos == p1 && x1 == -1)
		{
			x1 = i;
		}
		if(r.rpos == p2 && x2 == -1)
		{
			x2 = i;
		}
		if(x1 != -1 && x2 != -1) break;
	}


	//printf("Printing x1 and x2: x1 = %d, x2 = %d\n",x1,x2);
	
	
	if(x1 != -1 && x2 != -1 && x1 > x2) return -1;

	if(x1 != -1 && x2 != -1 && x1 == x2) //single exon circrna
	{
		//printf("Entering case equal\n");
		circ.start = bd->regions[x1].lpos;
		circ.end = bd->regions[x1].rpos;
		circ.circ_path.push_back(x1);
		//printf("check equal x1=%d, x2=%d\n",x1,x2);
		circ.circ_path_regions.push_back(bd->regions[x1]);
		circ.merged_regions.push_back(bd->regions[x1]);
	}

	if(x1 != -1 && x2 != -1 && x1 < x2)
	{
		build_junction_graph(bd->fragments);

		vector< vector<entry> > table;
		table.resize(bd->regions.size());
		printf("check dp x1=%d, x2=%d\n",x1,x2);

		dynamic_programming(x1, x2, table);
		vector< vector<int> > pb = trace_back(x2, table);
		if(pb.size() == 0) return -1;

		// just consider the best path
		vector<int> pp = pb[0];

		/*printf("Print path:");
		for(int j=0;j<pp.size();j++)
		{
			printf("%d ",pp[j]);
		}*/

		// then pp is path from x1 to x2 (I belive it includes x1 and x2)
    	// now the circRNA consists the regions in the path pp

		circ.start = bd->regions[x1].lpos;
		circ.end = bd->regions[x2].rpos;
		circ.circ_path.insert(circ.circ_path.end(),pp.begin(),pp.end());
		
		for(int i=0;i<circ.circ_path.size();i++)
		{
			circ.circ_path_regions.push_back(bd->regions[circ.circ_path[i]]);
		}

		join_interval_map jmap;
		for(int k = 0; k < circ.circ_path_regions.size(); k++)
		{
			int32_t p1 = circ.circ_path_regions[k].lpos;
			int32_t p2 = circ.circ_path_regions[k].rpos;
			jmap += make_pair(ROI(p1, p2), 1);
		}

		for(JIMI it = jmap.begin(); it != jmap.end(); it++)
		{
			region r(lower(it->first), upper(it->first), '.', '.');
			circ.merged_regions.push_back(r);
		}
	}
	
	return 0;
}

int bridger::bridge_overlapped_fragments(vector<fragment> &frags)
{
	for(int i = 0; i < frags.size(); i++)
	{
		fragment &fr = frags[i];
		bridge_overlapped_fragment(fr, 0, 0);
	}
	return 0;
}

int bridger::bridge_overlapped_fragment(fragment &fr, int ex1, int ex2)
{
	vector<int> v1 = decode_vlist(fr.h1->vlist);
	vector<int> v2 = decode_vlist(fr.h2->vlist);

	assert(v1.size() > ex1);
	assert(v2.size() > ex2);

	vector<int>::iterator t1 = v1.end() - ex1;
	vector<int>::iterator t2 = v2.begin() + ex2;

	int x1 = v1[v1.size() - 1 - ex1];
	int x2 = v2[ex2];

	
	if(x1 < x2) return 0;

	vector<int>::iterator it = find(t2, v2.end(), x1);
	if(it == v2.end()) return 0;

	/*if(strcmp(fr.h1->qname.c_str(),"simulate:173319") == 0)
	{
		printf("simulate:173319 is in bridge_overlapped_fragments\n");
	}*/

	// calculate the maximum abundence of overlapped vertices
	double max_abd = 0;
	vector<int>::iterator j1, j2;
	for(j1 = t1 - 1, j2 = it; j1 >= v1.begin() && j2 >= t2; j1--, j2--)
	{
		if((*j1) != (*j2)) return 0;
		double w = bd->regions[*j1].ave;
		if(w > max_abd) max_abd = w;
	}

	path p;
	p.ex1 = ex1;
	p.ex2 = ex2;
	p.v.insert(p.v.end(), v1.begin(), t1);
	p.v.insert(p.v.end(), it + 1, v2.end());
	p.length = bd->compute_aligned_length(fr.k1l, fr.k2r, p.v);
	p.v = encode_vlist(p.v);
	if(p.length >= length_low && p.length <= length_high)
	{
		p.type = 3;
		//bd->breads.insert(fr.h1->qname);
	}
	else p.type = 4;

	p.score = max_abd;

	fr.paths.push_back(p);

	return 0;
}

int bridger::bridge_phased_fragments(vector<fcluster> &fclusters)
{
	//vector<fcluster> fclusters;
	//cluster_open_fragments(fclusters);

	for(int k = 0; k < fclusters.size(); k++)
	{
		fcluster &fc = fclusters[k];

		if(fc.v1.size() <= 0) continue;
		if(fc.v2.size() <= 0) continue;

		phase_cluster(fc);

		if(fc.phase.size() <= 0) continue;
		bridge_phased_cluster(fc);
	}
	return 0;
}

int bridger::cluster_open_fragments(vector<fcluster> &fclusters, vector<fragment> &frags)
{
	vector<fragment*> open;
	for(int i = 0; i < frags.size(); i++)
	{
		fragment &fr = frags[i];
		if(fr.paths.size() >= 1) continue;
		//if(fr.paths.size() == 1 && fr.paths[0].type == 1) continue;
		int last1 = fr.h1->vlist[fr.h1->vlist.size() - 2] + fr.h1->vlist.back() - 1;
		int last2 = fr.h2->vlist[fr.h2->vlist.size() - 2] + fr.h2->vlist.back() - 1;
		if(last1 >= last2) continue;
		//if(fr.h1->vlist.back() >= fr.h2->vlist.front()) continue;
		open.push_back(&(frags[i]));
	}

	if(open.size() == 0) return 0;

	//sort(open.begin(), open.end(), compare_fragment_v3);
	sort(open.begin(), open.end(), compare_fragment_v3_flank);

	fcluster fc;
	vector<int> vv1;
	vector<int> vv2;

	int32_t flank1 = 0 - max_clustering_flank;
	int32_t flank2 = 0 - max_clustering_flank;
	for(int k = 0; k < open.size(); k++)
	{
		fragment *fr = open[k];
		int32_t f1 = fr->k1l + fr->k2l;
		int32_t f2 = fr->k1r + fr->k2r;
		int diff = (int)(fabs(f1 - flank1) + fabs(f2 - flank2));
		flank1 = f1;
		flank2 = f2;
		if(fr->h1->vlist == vv1 && fr->h2->vlist == vv2 && diff <= max_clustering_flank)
		{
			//printf("flank1 = %d, flank2 = %d, f1 = %d, f2 = %d, diff = %d\n", flank1, flank2, f1, f2, diff);
			fc.fset.push_back(fr);
		}
		else
		{
			if(fc.fset.size() >= 1) fclusters.push_back(fc);
			vv1 = fr->h1->vlist;
			vv2 = fr->h2->vlist;
			fc.clear();
			fc.type = 0;
			fc.fset.push_back(fr);
			fc.v1 = decode_vlist(vv1);
			fc.v2 = decode_vlist(vv2);
		}
	}
	if(fc.fset.size() >= 1) fclusters.push_back(fc);
	return 0;
}

int bridger::build_junction_graph(vector<fragment> &frags)
{
	pnodes.clear();
	build_path_nodes(2, frags);
	add_consecutive_path_nodes();

	int n = bd->regions.size();
	jsetx.clear();
	jsety.clear();
	jsetx.resize(n);
	jsety.resize(n);

	for(int i = 0; i < pnodes.size(); i++)
	{
		vector<int> &v = pnodes[i].v;
		assert(v.size() <= 2);
		if(v.size() <= 1) continue;
		int w = (int)(pnodes[i].score);
		int x = v[0];
		int y = v[1];

		assert(jsetx[x].find(y) == jsetx[x].end());
		assert(jsety[y].find(x) == jsety[y].end());
		jsetx[x].insert(PI(y, w));
		jsety[y].insert(PI(x, w));
	}
	
	return 0;
}

int bridger::add_consecutive_path_nodes()
{
	set<PI> s;
	for(int i = 0; i < pnodes.size(); i++)
	{
		vector<int> &v = pnodes[i].v;
		for(int k = 0; k < v.size() - 1; k++)
		{
			PI p(v[k], v[k + 1]);
			if(s.find(p) == s.end()) s.insert(p);
		}
	}

	for(int i = 0; i < bd->regions.size() - 1; i++)
	{
		if(s.find(PI(i, i + 1)) != s.end()) continue;
		region &r1 = bd->regions[i + 0];
		region &r2 = bd->regions[i + 1];
		if(r1.rpos != r2.lpos) continue;

		path p;
		p.v.push_back(i + 0);
		p.v.push_back(i + 1);
		p.score = 1;
		p.acc = bd->build_accumulate_length(p.v);
		pnodes.push_back(p);
	}
	return 0;
}

int bridger::build_path_nodes(int low, int high, vector<fragment> &frags)
{
	int m = (low + high) / 2;
	build_path_nodes(m, frags);

	//printf("build path nodes: pnodes = %lu, low = %d, high = %d\n", pnodes.size(), low, high);

	if(high - low <= 1) return 0;
	if(pnodes.size() > max_num_path_nodes) return build_path_nodes(low, max_pnode_length, frags);
	if(pnodes.size() <= max_num_path_nodes) return build_path_nodes(max_pnode_length, high, frags);

	return 0;
}

int bridger::build_path_nodes(int max_len, vector<fragment> &frags)
{
	max_pnode_length = max_len;
	map<vector<int>, int> m;
	set<int> hs;
	for(int i = 0; i < frags.size(); i++)
	{
		// TODO, also check length
		fragment &fr = frags[i];
		if(fr.paths.size() == 1 && fr.paths[0].type == 1)
		{
			vector<int> v = decode_vlist(fr.paths[0].v);
			if(v.size() <= 1) continue;
			build_path_nodes(m, v, fr.cnt);		// TODO: consider cnt of fragments
			hs.insert(fr.h1->hid);
			hs.insert(fr.h2->hid);
		}
		/*
		else
		{
			vector<int> v1 = decode_vlist(fr.h1->vlist);
			vector<int> v2 = decode_vlist(fr.h2->vlist);
			build_path_nodes(m, v1, fr.cnt);
			build_path_nodes(m, v2, fr.cnt);
		}
		*/
	}

	for(int i = 0; i < bd->bb.hits.size(); i++)
	{
		// TODO, also check length
		hit &h = bd->bb.hits[i];
		if(hs.find(h.hid) != hs.end()) continue;
		vector<int> v1 = decode_vlist(h.vlist);
		build_path_nodes(m, v1, 1);
	}


	pnodes.clear();
	for(map<vector<int>, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		path p;
		p.v = it->first;
		p.score = it->second;
		p.acc = bd->build_accumulate_length(p.v);
		pnodes.push_back(p);
		//adjust_path_score(p);
	}
	
	return 0;
}

int bridger::adjust_path_score(path &p)
{
	double m = p.score;
	for(int k = 0; k < p.v.size(); k++)
	{
		region &r = bd->regions[p.v[k]];
		if(r.ave < m) m = r.ave;
		//if(r.ltype != LEFT_SPLICE || r.rtype != RIGHT_SPLICE) continue;
		//if(r.ave - r.dev < m) m = r.ave - r.dev;
	}
	p.score = m * 100;
	//if(p.score < 0) p.score = 0;
	return 0;
}

int bridger::build_path_nodes(map<vector<int>, int> &m, const vector<int> &v, int cnt)
{
	cnt = 1;
	if(v.size() <= 0) return 0;
	int n = v.size();
	if(n > max_pnode_length) n = max_pnode_length;

	for(int i = 0; i <= v.size() - n; i++)
	{
		int j = i + n;
		vector<int> s(v.begin() + i, v.begin() + j);
		if(m.find(s) == m.end()) m.insert(pair<vector<int>, int>(s, cnt));
		else m[s] += cnt;
	}

	//printf(" current map contains %lu pnodes\n", m.size());
	return 0;
}

int bridger::phase_cluster(fcluster &fc)
{
	fc.phase.clear();
	map<int, int> xm;
	vector<PI> xp = bd->ref_index[fc.v1.front()];
	for(int j = 0; j < xp.size(); j++)
	{
		int ti = xp[j].first;
		int ki = xp[j].second;
		bool b = true;
		for(int k = 0; k < fc.v1.size(); k++)
		{
			if(ki + k >= bd->ref_phase[ti].size() || fc.v1[k] != bd->ref_phase[ti][ki + k])
			{
				b = false;
				break;
			}
		}
		if(b == false) continue;
		xm.insert(pair<int, int>(ti, ki));
	}

	vector<PI> yp = bd->ref_index[fc.v2.front()];
	for(int j = 0; j < yp.size(); j++)
	{
		int ti = yp[j].first;
		int ki = yp[j].second;
		if(xm.find(ti) == xm.end()) continue;
		if(ki < xm[ti] + fc.v1.size()) continue;

		bool b = true;
		for(int k = 0; k < fc.v2.size(); k++)
		{
			if(ki + k >= bd->ref_phase[ti].size() || fc.v2[k] != bd->ref_phase[ti][ki + k])
			{
				b = false;
				break;
			}
		}
		if(b == false) continue;

		vector<int> vv(bd->ref_phase[ti].begin() + xm[ti], bd->ref_phase[ti].begin() + ki + fc.v2.size());
		fc.add_phase(vv);
	}

	return 0;
}

int bridger::bridge_phased_cluster(fcluster &fc)
{
	for(int i = 0; i < fc.fset.size(); i++)
	{
		fragment *fr = fc.fset[i];

		for(int k = 0; k < fc.phase.size(); k++)
		{
			path p;
			p.ex1 = p.ex2 = 0;
			p.v = fc.phase[k];
			p.length = bd->compute_aligned_length(fr->k1l, fr->k2r, p.v);
			p.v = encode_vlist(p.v);
			p.score = 1;
			if(p.length >= length_low && p.length <= length_high) p.type = 1;
			else p.type = 2;
			fr->paths.push_back(p);
		}
	}
	return 0;
}

int bridger::remove_tiny_boundary(vector<fragment> &frags)
{
	for(int i = 0; i < frags.size(); i++)
	{
		fragment &fr = frags[i];

		if(fr.paths.size() == 1 && fr.paths[0].type == 1) continue;

		//if(fr.h1->bridged == true) continue;
		//if(fr.h2->bridged == true) continue;

		vector<int> v1 = decode_vlist(fr.h1->vlist);
		int n1 = v1.size();
		if(n1 >= 2 && v1[n1 - 2] + 1 == v1[n1 - 1])
		{
			int k = v1[n1 - 1];
			int32_t total = bd->regions[k].rpos - bd->regions[k].lpos;
			int32_t flank = fr.h1->rpos - bd->regions[k].lpos;

			if(flank <= flank_tiny_length && 1.0 * flank / total < flank_tiny_ratio)
			{
				vector<int> v(v1.begin(), v1.begin() + n1 - 1);
				assert(v.size() + 1 == v1.size());
				fr.h1->vlist = encode_vlist(v);
				fr.h1->rpos = bd->regions[k].lpos;
			}
		}

		vector<int> v2 = decode_vlist(fr.h2->vlist);
		int n2 = v2.size();
		if(n2 >= 2 && v2[0] + 1 == v2[1])
		{
			int k = v2[0];
			int32_t total = bd->regions[k].rpos - bd->regions[k].lpos;
			int32_t flank = bd->regions[k].rpos - fr.h2->pos;

			if(flank <= flank_tiny_length && 1.0 * flank / total < flank_tiny_ratio)
			{
				vector<int> v(v2.begin() + 1, v2.end());
				assert(v.size() + 1 == v2.size());
				fr.h2->vlist = encode_vlist(v);
				fr.h2->pos = bd->regions[k].rpos;
			}
		}
	}
	return 0;
}

int bridger::build_path_nodes(vector<fragment> &frags)
{
	int low = 10;
	int high = 50;
	build_path_nodes(high, frags);
	if(pnodes.size() > max_num_path_nodes)
	{
		build_path_nodes(low, frags);
		if(pnodes.size() < max_num_path_nodes) build_path_nodes(low, high, frags);
	}
	sort(pnodes.begin(), pnodes.end(), compare_path_vertices);
	return 0;
}

int bridger::bridge_hard_fragments_normal(vector<fcluster> &open)
{
	/*
	if(use_overlap_scoring == true)
	{
		build_path_nodes(frags);
		build_overlap_index();
	}
	*/

	//vector<fcluster> open;
	//cluster_open_fragments(open);
	sort(open.begin(), open.end(), compare_fcluster_v1_v2);

	//print open clusters
	if(verbose >= 1)
	{
		for(int k = 0; k < open.size(); k++)
		{
			open[k].print(k);
		}
	}

	/*
	printf("print pnode bridge...\n");
	for(int k = 0; k < bd->regions.size(); k++) bd->regions[k].print(k);
	for(int k = 0; k < pnodes.size(); k++) pnodes[k].print_bridge(k);
	*/

	vector< set<int> > affected(bd->regions.size());
	vector<int> max_needed(bd->regions.size(), -1);
	for(int k = 0; k < open.size(); k++)
	{
		fcluster &fc = open[k];
		int x1 = fc.v1.back();
		int x2 = fc.v2.front();
		assert(x1 >= 0 && x1 < bd->regions.size());
		assert(x2 >= 0 && x2 < bd->regions.size());
		affected[x1].insert(k);
		if(max_needed[x1] < x2) max_needed[x1] = x2;
	}

	vector< vector<entry> > table;
	table.resize(bd->regions.size());
	for(int k = 0; k < bd->regions.size(); k++)
	{
		if(affected[k].size() <= 0) continue;
		if(max_needed[k] < k) continue;

		dynamic_programming(k, max_needed[k], table);

		// print table
		/*
		   printf("table from vertex %d to %d\n", k, max_needed);
		   for(int j = k; j <= max_needed; j++)
		   {
		   for(int i = 0; i < table[j].size(); i++)
		   {
		   entry &e = table[j][i];
		   printf("vertex %d, solution %d: ", j, i);
		   e.print();
		   }
		   }
		 */

		for(set<int>::iterator it = affected[k].begin(); it != affected[k].end(); it++)
		{
			fcluster &fc = open[*it];

			int j = fc.v2.front();
			assert(k == fc.v1.back());
			assert(j <= max_needed[k]);

			if(j < k) continue;
			if(table[j].size() == 0) continue;

			vector< vector<int> > pb = trace_back(j, table);
			vector< vector<int> > pn;
			vector<int> ps;

			for(int e = 0; e < pb.size(); e++)
			{
				vector<int> px = fc.v1;
				if(pb[e].size() >= 2) px.insert(px.end(), pb[e].begin() + 1, pb[e].end() - 1);
				px.insert(px.end(), fc.v2.begin(), fc.v2.end());
				int s = (int)(min_bridging_score) + 2;
				if(use_overlap_scoring) s = evaluate_bridging_path(px);
				pn.push_back(px);
				ps.push_back(s);
			}

			/* not used
			int best_path = 0;
			for(int e = 1; e < pb.size(); e++)
			{
				if(ps[e] > ps[best_path])
				{
					best_path = e;
				}
				if(ps[e] == ps[best_path] && compare_stack(table[j][e].stack, table[j][best_path].stack) >= 1)
				{
					best_path = e;
				}
			}
			*/

			vector<int> votes;
			votes.resize(pb.size(), 0);
			for(int i = 0; i < fc.fset.size(); i++)
			{
				fragment *fr = fc.fset[i];

				vector<int> best_stack;
				int best_score = -1;
				int best_index = -1;

				for(int e = 0; e < pb.size(); e++)
				{
					int32_t length = bd->compute_aligned_length(fr->k1l, fr->k2r, pn[e]);
					//printf(" fragment %d length = %d using path %d\n", i, p.length, e);

					// note by Mingfu
					// length_low and length_high is set for paired-end reads, which admits a fragment distribution
					// for smart-seq 3, such distribution may not make sense
					// by commenting the two lines below, the best single path (among the n candicates) will be used

					// try using distribution
					// only fragments bridged within reasonable 
					// range are entitled to vote; others don't vote

					// TODO, vote for the path with best insertsize?
					if(length < length_low) continue;
					if(length > length_high) continue;

					if(ps[e] > best_score)
					{
						best_score = ps[e];
						best_stack = table[j][e].stack;
						best_index = e;
					}
					else if(ps[e] == best_score && compare_stack(table[j][e].stack, best_stack) >= 1)
					{
						best_stack = table[j][e].stack;
						best_index = e;
					}
				}
				if(best_index >= 0) votes[best_index]++;
			}

			int be = 0;
			int voted = votes[0];
			for(int i = 1; i < votes.size(); i++)
			{
				voted += votes[i];
				if(votes[i] > votes[be]) be = i;
			}

			// don't do these -- even if no one votes, still keep the path
			/*
			if(votes[be] <= 0) continue;
			if(voted <= 0) continue;
			*/


			/*
			double voting_ratio = 100.0 * voted / fc.fset.size();
			double best_ratio = 100.0 * votes[be] / voted;

			printf("total %lu fragments, %d voted, best = %d, voting-ratio = %.2lf, best-ratio = %.2lf ( ", 
					fc.fset.size(), voted, be, voting_ratio, best_ratio);
			printv(votes);
			printf(")\n");
			*/

			//if(voting_ratio <= 0.49) continue;
			//if(best_ratio < 0.8 && be != best_path) continue;

			/*
			printf("fcluster with %lu fragments, total %lu paths, best = %d, from %d to %d, v1 = (", fc.fset.size(), pb.size(), be, k, j);
			printv(fc.v1);
			printf("), v2 = ( ");
			printv(fc.v2);
			printf(")\n");
			for(int e = 0; e < pb.size(); e++)
			{
				printf(" path %d, votes = %d, score = %d, stack = (", e, votes[e], ps[e]); 
				printv(table[j][e].stack);
				printf("), pb = (");
				printv(pb[e]);
				printf("), pn = (");
				printv(pn[e]);
				printf(")\n");
			}
			*/

			for(int i = 0; i < fc.fset.size(); i++)
			{
				fragment *fr = fc.fset[i];

				path p;
				p.ex1 = p.ex2 = 0;
				p.v = pn[be];
				p.length = bd->compute_aligned_length(fr->k1l, fr->k2r, p.v);
				p.v = encode_vlist(p.v);
				p.score = ps[be];

				if(p.length >= length_low && p.length <= length_high)
				{
					//bd->breads.insert(fr->h1->qname);
					p.type = 3;
				}
				else p.type = 4;

				fr->paths.push_back(p);
				//printf(" fragment %d length = %d using path %d, p.type = %d\n", i, p.length, be, p.type);
			}
		}
	}
	return 0;
}

int bridger::bridge_hard_fragments_circ(vector<fcluster> &open)
{
	/*
	if(use_overlap_scoring == true)
	{
		build_path_nodes(frags);
		build_overlap_index();
	}
	*/

	//vector<fcluster> open;
	//cluster_open_fragments(open);
	sort(open.begin(), open.end(), compare_fcluster_v1_v2);

	//print open clusters
	if(verbose >= 1)
	{
		for(int k = 0; k < open.size(); k++)
		{
			open[k].print(k);
		}
	}

	/*
	printf("print pnode bridge...\n");
	for(int k = 0; k < bd->regions.size(); k++) bd->regions[k].print(k);
	for(int k = 0; k < pnodes.size(); k++) pnodes[k].print_bridge(k);
	*/

	vector< set<int> > affected(bd->regions.size());
	vector<int> max_needed(bd->regions.size(), -1);
	for(int k = 0; k < open.size(); k++)
	{
		fcluster &fc = open[k];
		int x1 = fc.v1.back();
		int x2 = fc.v2.front();
		assert(x1 >= 0 && x1 < bd->regions.size());
		assert(x2 >= 0 && x2 < bd->regions.size());
		affected[x1].insert(k);
		if(max_needed[x1] < x2) max_needed[x1] = x2;
	}

	vector< vector<entry> > table;
	table.resize(bd->regions.size());
	for(int k = 0; k < bd->regions.size(); k++)
	{
		if(affected[k].size() <= 0) continue;
		if(max_needed[k] < k) continue;

		dynamic_programming(k, max_needed[k], table);

		// print table
		/*
		   printf("table from vertex %d to %d\n", k, max_needed);
		   for(int j = k; j <= max_needed; j++)
		   {
		   for(int i = 0; i < table[j].size(); i++)
		   {
		   entry &e = table[j][i];
		   printf("vertex %d, solution %d: ", j, i);
		   e.print();
		   }
		   }
		 */

		for(set<int>::iterator it = affected[k].begin(); it != affected[k].end(); it++)
		{
			fcluster &fc = open[*it];

			int j = fc.v2.front();
			assert(k == fc.v1.back());
			assert(j <= max_needed[k]);

			if(j < k) continue;
			if(table[j].size() == 0) continue;

			vector< vector<int> > pb = trace_back(j, table);
			vector< vector<int> > pn;
			vector<int> ps;

			for(int e = 0; e < pb.size(); e++)
			{
				vector<int> px = fc.v1;
				if(pb[e].size() >= 2) px.insert(px.end(), pb[e].begin() + 1, pb[e].end() - 1);
				px.insert(px.end(), fc.v2.begin(), fc.v2.end());
				//int s = (int)(min_bridging_score) + 2;
				int s = table[j][e].stack.front();
				if(use_overlap_scoring) s = evaluate_bridging_path(px);
				pn.push_back(px);
				ps.push_back(s);
			}

			for(int i = 0; i < fc.fset.size(); i++)
			{
				fragment *fr = fc.fset[i];

				for(int e = 0; e < pb.size(); e++)
				{
					path p;
					p.ex1 = p.ex2 = 0;
					p.v = pn[e];
					p.length = bd->compute_aligned_length(fr->k1l, fr->k2r, p.v);
					p.v = encode_vlist(p.v);
					p.score = ps[e];

					// compare score with fset.size
					//printf("FSET-SCORE: fset %lu, score %.1lf, read %s, lpos = %d/%d, length %d\n",fc.fset.size(), p.score, fr->h1->qname.c_str(), fr->h1->pos, fr->h2->pos, p.length);

					double fset_score = log(1 + fc.fset.size()) - log(1 + p.score);
					if(fset_score > max_fset_score) continue;

					if(p.length >= length_low && p.length <= length_high)
					{
						//bd->breads.insert(fr->h1->qname);
						p.type = 3;
					}
					else p.type = 4;

					fr->paths.push_back(p);
					//printf(" fragment %d length = %d using path %d, p.type = %d\n", i, p.length, be, p.type);
				}
			}
		}
	}
	return 0;
}

int bridger::compare_stack(const vector<int> &x, const vector<int> &y)
{
	if(x.size() != y.size())
	{
		printf("|");
		printv(x);
		printf("|");
		printv(y);
		printf("|\n");
	}
	assert(x.size() == y.size());
	for(int i = 0; i < x.size() - 1; i++) assert(x[i] <= x[i + 1]);
	for(int i = 0; i < y.size() - 1; i++) assert(y[i] <= y[i + 1]);
	for(int i = 0; i < x.size(); i++)
	{
		if(x[i] > y[i]) return +1;
		if(x[i] < y[i]) return -1;
	}
	return 0;
}

vector<int> bridger::update_stack(const vector<int> &v, int s)
{
	vector<int> stack(v.size(), 0);
	for(int i = 0, j = 0; i < v.size() && j < v.size(); i++, j++)
	{
		if(i == j && v[i] > s)
		{
			stack[j] = s;
			j++;
		}
		stack[j] = v[i];
	}
	return stack;
}

int bridger::dynamic_programming(int k1, int k2, vector< vector<entry> > &table)
{
	int n = bd->regions.size();
	assert(k1 >= 0 && k1 < n);
	assert(k2 >= 0 && k2 < n);

	table.clear();
	table.resize(n);

	table[k1].resize(1);
	table[k1][0].stack.assign(dp_stack_size, 999999);
	table[k1][0].length = bd->regions[k1].rpos - bd->regions[k1].lpos;
	table[k1][0].trace1 = -1;
	table[k1][0].trace2 = -1;

	for(int k = k1 + 1; k <= k2; k++)
	{
		vector<entry> v;
		int32_t len = bd->regions[k].rpos - bd->regions[k].lpos;
		for(map<int, int>::reverse_iterator it = jsety[k].rbegin(); it != jsety[k].rend(); it++)
		{
			int j = it->first;
			int w = it->second;
			if(j < k1) continue;
			if(table[j].size() == 0) continue;

			for(int i = 0; i < table[j].size(); i++)
			{
				entry e;
				e.stack = update_stack(table[j][i].stack, w);
				e.length = table[j][i].length + len;
				e.trace1 = j;
				e.trace2 = i;
				v.push_back(e);
			}
		}

		sort(v.begin(), v.end(), entry_compare);
		if(v.size() > dp_solution_size) v.resize(dp_solution_size);
		table[k] = v;
	}
	return 0;
}

int bridger::evaluate_bridging_path(const vector<int> &pb)
{
	int max_score = 0;
	vector<int> table(pnodes.size(), 0);
	for(int k = 0; k < pnodes.size(); k++)
	{
		const vector<int> &v = pnodes[k].v;
		if(v.back() < pb.front()) continue;
		if(v.front() > pb.back()) break;

		int s = (int)(pnodes[k].score);

		PI p;
		int t = determine_overlap(v, pb, p);

		if(t == 1)
		{
			table[k] = s;
			if(v.back() >= pb.back() && s > max_score) max_score = s;
		}
		else if(t == 2)
		{
			if(s > max_score) max_score = s;
		}
		else if(t == 3 || t == 4)
		{
			int max = 0;
			for(map<int, int>::reverse_iterator it = psety[k].rbegin(); it != psety[k].rend(); it++)
			{
				int j = it->first;
				if(table[j] > max) max = table[j];
			}
			table[k] = max < s ? max : s;
			if(t == 3 && table[k] > max_score) max_score = table[k];
		}
	}
	return max_score;
}

int bridger::dynamic_programming(int k1, int k2, vector<int> &trace, vector< vector<int> > &table_cov, vector<int32_t> &table_len)
{
	assert(k1 >= 0 && k1 < pnodes.size());
	assert(k2 >= 0 && k2 < pnodes.size());

	for(int k = 0; k < table_cov.size(); k++) table_cov[k].assign(dp_stack_size, -1);
	table_len.assign(table_len.size(), 888888);
	trace.assign(pnodes.size(), -1);

	table_cov[k1].assign(dp_stack_size, 999999);
	table_len[k1] = pnodes[k1].acc.back();
	trace[k1] = -1;
	for(int k = k1 + 1; k <= k2; k++)
	{
		vector<int> stack;
		stack.assign(dp_stack_size, -1);
		int32_t blen = 888888;
		int back = -1;
		for(map<int, int>::reverse_iterator it = psety[k].rbegin(); it != psety[k].rend(); it++)
		{
			int j = it->first;
			if(j < k1) continue;
			if(table_cov[j][0] < 0) continue;

			int jx = it->second;
			int32_t len = table_len[j] + get_extended_length2(jx, j, k);
			//if(len > length_high) continue;

			int s = (int)(pnodes[j].score);
			if(s < stack[0] && table_cov[j][0] < stack[0]) continue;

			vector<int> v = update_stack(table_cov[j], s);

			int b = compare_stack(stack, v);

			if(b == -1)
			{
				stack = v;
				blen = len;
				back = j;
			}
			else if(b == 0 && len < blen)
			{
				blen = len;
				back = j;
			}

			/*
			double cov1 = table_cov[j] < pnodes[j].score ? table_cov[j] : pnodes[j].score;
			double cov2 = table_cov[j] > pnodes[j].score ? table_cov[j] : pnodes[j].score;
			if(cov1 >= bcov1 + 0.0001)
			{
				bcov1 = cov1;
				bcov2 = cov2;
				blen = len;
				back = j;
			}
			else if(fabs(cov1 - bcov1) < 0.0001 && cov2 >= bcov2 + 0.0001)
			{
				bcov1 = cov1;
				bcov2 = cov2;
				blen = len;
				back = j;
			}
			else if(fabs(cov1 - bcov1) < 0.0001 && fabs(cov2 - bcov2) < 0.0001 && len < blen)
			{
				//printf("cov = %.2f, bcov = %.2f, len = %d, blen = %d, k1 = %d, k2 = %d, s = %d, t = %d, back = %d\n", cov, bcov, len, blen, k1, k2, s, t, back);
				blen = len;
				back = j;
			}
			*/
		}
		table_cov[k] = stack;
		table_len[k] = blen;
		trace[k] = back;
	}
	return 0;
}

int32_t bridger::get_extended_length1(int k2, int p1, int p2)
{
	path &px = pnodes[p1];
	path &py = pnodes[p2];
	int l = k2 + 1;
	int k1 = py.v.size() - l;
	if(k1 == 0) return 0;
	else return px.acc[k1 - 1];
}

int32_t bridger::get_extended_length2(int k1, int p1, int p2)
{
	path &px = pnodes[p1];
	path &py = pnodes[p2];
	int l = px.v.size() - k1;
	int k2 = l - 1;
	return py.acc.back() - py.acc[k2];
}

vector< vector<int> > bridger::trace_back(int k, const vector< vector<entry> > &table)
{
	vector< vector<int> > vv;
	for(int i = 0; i < table[k].size(); i++)
	{
		vector<int> v;
		int p = k;
		int q = i;
		while(true)
		{
			v.push_back(p);
			const entry &e = table[p][q];
			p = e.trace1;
			q = e.trace2;
			if(p < 0) break;
		}
		reverse(v);
		vv.push_back(v);
	}
	return vv;
}

vector<int> bridger::trace_back(int k1, int k2, const vector<int> &trace)
{
	vector<int> v;
	int k = k2;
	while(true)
	{
		v.push_back(k);
		k = trace[k];
		if(k < 0) break;
		assert(k >= k1);
	}
	reverse(v);
	return v;
}

vector<int> bridger::get_bridge(const vector<int> &vv, const vector<int> &v1, const vector<int> &v2)
{
	vector<int> v;
	if(vv.size() == 0) return v;

	v = v1;
	for(int k = 0; k < vv.size(); k++)
	{
		const path &p = pnodes[vv[k]];
		int x = v.back();
		vector<int>::const_iterator it = lower_bound(p.v.begin(), p.v.end(), x);
		assert(it != p.v.end());
		assert((*it) == x);
		v.insert(v.end(), it + 1, p.v.end());
	}

	int x = v.back();
	vector<int>::const_iterator it = lower_bound(v2.begin(), v2.end(), x);
	assert(it != v2.end());
	assert((*it) == x);
	v.insert(v.end(), it + 1, v2.end());
	return v;
}

int bridger::build_overlap_index()
{
	/*
	pset1.resize(pnodes.size());
	pset2.resize(pnodes.size());
	for(int k = 0; k < pnodes.size(); k++) pset1[k].clear();
	for(int k = 0; k < pnodes.size(); k++) pset2[k].clear();
	*/

	psetx.clear();
	psety.clear();
	psetx.resize(pnodes.size());
	psety.resize(pnodes.size());
	for(int k = 0; k < pnodes.size(); k++) psetx[k].clear();
	for(int k = 0; k < pnodes.size(); k++) psety[k].clear();

	int cnt1 = 0;
	for(int i = 0; i < pnodes.size(); i++)
	{
		const vector<int> &vx = pnodes[i].v;
		for(int j = i + 1; j < pnodes.size(); j++)
		{
			const vector<int> &vy = pnodes[j].v;

			PI p;
			int t = determine_overlap(vx, vy, p);
			if(t != 1) continue;
			assert(p.first >= 0 && p.first < vx.size());
			assert(p.second >= 0 && p.second < vy.size());

			//int32_t len2 = pnodes[j].acc.back() - pnodes[j].acc[p.second];
			//int32_t len1 = (p.first == 0) ? 0 : pnodes[i].acc[p.first - 1];

			// pset1 & 2 are full
			//pset1[i].insert(pair<int, int>(j, p.second));
			//pset2[j].insert(pair<int, int>(i, p.first));

			// psetx & y are not full
			// consider only (1,2,3) -> (2,3,4)
			//if(p.first == 0) continue;
			//if(p.second == vy.size() - 1) continue;
			psetx[i].insert(pair<int, int>(j, p.second));
			psety[j].insert(pair<int, int>(i, p.first));
			cnt1++;
		}
	}

	// do not use this filtering when with fixed kmer-size
	/*
	// filtering: (1,2,3,4) -> (2,3,4,5,6), (3,4,5,6), (4,5,6):
	// only keep the one with maximized score
	int cnt3 = cnt2;
	for(int i = 0; i < pnodes.size(); i++)
	{
		map<size_t, int> m;
		const path &pi = pnodes[i];
		for(map<int, int>::iterator x1 = psetx[i].begin(); x1 != psetx[i].end(); x1++)
		{
			int j = x1->first;
			int jx = x1->second;
			if(jx < 0) continue;

			const path &pj = pnodes[j];
			//int32_t l = pj.acc.back() - pj.acc[jx];
			size_t l = hash_range(pj.v.begin() + jx + 1, pj.v.end());

			map<size_t, int>::iterator it = m.find(l);

			if(it == m.end())
			{
				m.insert(pair<size_t, int>(l, j));
			}
			else
			{
				int k = it->second;
				assert(j != k);

				const path &pk = pnodes[k];
				if(pj.score >= pk.score) 
				{
					it->second = j;
					psetx[i][k] = -1;
					psety[k][i] = -1;
					//printf("B: remove edge %d -> %d (j = %d, hash = %lu)\n", i, k, j, l);
					cnt3--;
				}
				else
				{
					psetx[i][j] = -1;
					psety[j][i] = -1;
					//printf("C: remove edge %d -> %d\n", i, j);
					cnt3--;
				}
			}
		}
	}
	*/

	// this is not always true
	/*
	// filtering: (1,2,3), (1,2,3,4), (1,2,3,4,5) -> (3,4,5,6)
	// only keep the one with maximized score
	for(int i = 0; i < pnodes.size(); i++)
	{
		map<int32_t, int> m;
		const path &pi = pnodes[i];
		for(map<int, int>::iterator x1 = psety[i].begin(); x1 != psety[i].end(); x1++)
		{
			int j = x1->first;
			int jx = x1->second;
			if(jx < 0) continue;
			assert(jx >= 1);

			const path &pj = pnodes[j];
			int32_t l = pj.acc[jx - 1];
			map<int32_t, int>::iterator it = m.find(l);

			if(it == m.end())
			{
				m.insert(pair<int32_t, int>(l, j));
			}
			else
			{
				int k = it->second;
				assert(j != k);

				const path &pk = pnodes[k];
				if(pj.score >= pk.score) 
				{
					it->second = j;
					psety[i][k] = -1;
					psetx[k][i] = -1;
					cnt3--;
				}
				else
				{
					psety[i][j] = -1;
					psetx[j][i] = -1;
					cnt3--;
				}
			}
		}
	}
	*/

	// i = (1,2,3), j = (2,3,4), k = (3,4,5)
	// if w(i) <= w(j) AND w(k) <= w(j) => remove (i,k)
	for(int i = 0; i < pnodes.size(); i++)
	{
		const path &pi = pnodes[i];
		for(map<int, int>::iterator x1 = psetx[i].begin(); x1 != psetx[i].end(); x1++)
		{
			int j = x1->first;
			if(x1->second < 0) continue;

			const path &pj = pnodes[j];
			if(pi.score > pj.score) continue;

			for(map<int, int32_t>::iterator x2 = psetx[j].begin(); x2 != psetx[j].end(); x2++)
			{
				int k = x2->first;
				if(x2->second < 0) continue;

				map<int, int32_t>::iterator it = psetx[i].find(k);
				if(it == psetx[i].end()) continue;
				if(it->second < 0) continue;

				assert(psety[k].find(i) != psety[k].end());
				assert(psety[k][i] >= 0);

				const path &pk = pnodes[k];
				if(pk.score > pj.score) continue;

				psetx[i][k] = -1;
				psety[k][i] = -1;
			}
		}
	}

	// copy and shrink psetx and psety
	int cnt2 = 0;
	for(int i = 0; i < pnodes.size(); i++)
	{
		map<int, int> m;
		for(map<int, int>::iterator x1 = psetx[i].begin(); x1 != psetx[i].end(); x1++)
		{
			int j = x1->first;
			int jx = x1->second;
			if(jx < 0) continue;
			m.insert(pair<int, int>(j, jx));
			cnt2++;
		}
		psetx[i] = m;
	}

	int cnt3 = 0;
	for(int i = 0; i < pnodes.size(); i++)
	{
		map<int, int> m;
		for(map<int, int>::iterator x1 = psety[i].begin(); x1 != psety[i].end(); x1++)
		{
			int j = x1->first;
			int jx = x1->second;
			if(jx < 0) continue;
			m.insert(pair<int, int>(j, jx));
			cnt3++;
		}
		psety[i] = m;
	}

	// TODO: test later
	//int cnt7 = cnt1;
	//int cnt8 = cnt1;

	/*
	// simplify pset1
	// i = (1,2,3), j = (2,3,4), k = (3,4,5)
	// (i,j) and (i,k) in pset1
	// (j,k) in psetx
	// if w(i) <= w(j) and w=> remove (i,k) in pset1
	for(int i = 0; i < pnodes.size(); i++)
	{
		const path &pi = pnodes[i];
		for(map<int, int>::iterator x1 = pset1[i].begin(); x1 != pset1[i].end(); x1++)
		{
			int j = x1->first;
			if(x1->second < 0) continue;

			const path &pj = pnodes[j];
			if(pi.score > pj.score) continue;

			for(map<int, int32_t>::iterator x2 = psetx[j].begin(); x2 != psetx[j].end(); x2++)
			{
				int k = x2->first;
				if(x2->second < 0) continue;

				map<int, int32_t>::iterator it = pset1[i].find(k);
				if(it == pset1[i].end()) continue;
				if(it->second < 0) continue;
				pset1[i][k] = -1;
				cnt7--;
			}
		}
	}
	*/
	
	/*
	pset4 = pset2;
	// simplify pset4
	// they are not necessarily symmetric
	// k = (1,2,3), j = (2,3,4), i = (3,4,5)
	// (j,i) and (k,i) in pset4
	// (k,j) in psety
	// if w(i) <= w(j) => remove (k,i) in pset4
	for(int i = 0; i < pnodes.size(); i++)
	{
		const path &pi = pnodes[i];
		for(map<int, int>::iterator x1 = pset4[i].begin(); x1 != pset4[i].end(); x1++)
		{
			int j = x1->first;
			if(x1->second < 0) continue;

			const path &pj = pnodes[j];
			for(map<int, int32_t>::iterator x2 = psety[j].begin(); x2 != psety[j].end(); x2++)
			{
				int k = x2->first;
				if(x2->second < 0) continue;

				const path &pk = pnodes[k];
				if(pk.score > pj.score) continue;

				map<int, int32_t>::iterator it = pset4[i].find(k);
				if(it == pset4[i].end()) continue;
				if(it->second < 0) continue;
				pset4[i][k] = -1;
				cnt8--;
			}
		}
	}
	*/

	printf("build overlap index with %lu nodes, max-pnode-length = %d, and %d -> %d / %d edges\n", 
			pnodes.size(), max_pnode_length, cnt1, cnt2, cnt3);

	// print
	/*
	printf("psetx:\n");
	for(int i = 0; i < pnodes.size(); i++)
	{
		printf("path node %d, score = %.0lf, list = ( ", i, pnodes[i].score);
		printv(pnodes[i].v);
		printf(")\n");
		for(map<int, int>::iterator it = psetx[i].begin(); it != psetx[i].end(); it++)
		{
			printf(" edge to %d, index = %d, score = %.0lf\n", it->first, it->second, pnodes[it->first].score);
		}
	}
	*/

	return 0;
}

vector<int> bridger::get_prefix(const vector<int> &v)
{
	if(max_pnode_length >= v.size()) return v;
	vector<int> s(v.begin(), v.begin() + max_pnode_length);
	return s;
}

vector<int> bridger::get_suffix(const vector<int> &v)
{
	if(max_pnode_length >= v.size()) return v;
	int k = v.size() - max_pnode_length;
	vector<int> s(v.begin() + k, v.end());
	return s;
}

bool bridger::determine_identical(const vector<int> &vx, const vector<int> &vy, int x1, int x2, int y1, int y2)
{
	assert(x1 <= x2);
	assert(y1 <= y2);
	if(x2 - x1 != y2 - y1) return false;
	if(vx[x1] != vy[y1]) return false;
	if(vx[x2] != vy[y2]) return false;
	int xx = (x1 + x2) / 2;
	int yy = (y1 + y2) / 2;
	if(vx[xx] != vy[yy]) return false;
	if(x1 + 1 <= xx - 1 && determine_identical(vx, vy, x1 + 1, xx - 1, y1 + 1, yy - 1) == false) return false;
	if(xx + 1 <= x2 - 1 && determine_identical(vx, vy, xx + 1, x2 - 1, yy + 1, y2 - 1) == false) return false;
	return true;
}

int bridger::determine_overlap(const vector<int> &vx, const vector<int> &vy, PI &p) 
{
	// assume that x and y are sorted
	// x = (x1, x2), y = (y1, y2)
	// if x1 <= y1 && x2 <= y2: type 1
	// if x1 <= y1 && x2 >  y2: type 2
	// if x1 >  y1 && x2 >= y2: type 3
	// if x1 >  y1 && x2 <  y2: type 4

	if(vx.size() == 0) return -1;
	if(vy.size() == 0) return -1;

	int t = determine_overlap1(vx, vy, p);
	if(t == 1 || t == 2) return t;

	PI q;
	t = determine_overlap1(vy, vx, q);
	if(t == 1)
	{
		p = PI(q.second, p.first);
		return 3;
	}

	if(t == 2)
	{
		p = q;
		return 4;
	}
	return -1;
}

int bridger::determine_overlap1(const vector<int> &vx, const vector<int> &vy, PI &p) 
{
	vector<int>::const_iterator ix;
	vector<int>::const_iterator iy;

	ix = lower_bound(vx.begin(), vx.end(), vy[0]);
	if(ix == vx.end()) return -1;

	int kx = ix - vx.begin();
	iy = lower_bound(vy.begin(), vy.end(), vx.back());
	if(iy != vy.end())
	{
		int ky = iy - vy.begin();
		p = PI(kx, ky);
		bool b = determine_identical(vx, vy, kx, vx.size() - 1, 0, ky);
		if(b == true) return 1;
		else return -1;
	}

	iy = lower_bound(vx.begin(), vx.end(), vy.back());
	if(iy != vx.end())
	{
		int ky = iy - vx.begin();
		p = PI(kx, ky);
		bool b = determine_identical(vx, vy, kx, ky, 0, vy.size() - 1);
		if(b == true) return 2;
		else return -1;
	}

	return -1;
}

int bridger::set_thresholds()
{
	return 0;
}

int bridger::set_normal_length()
{
	length_median = insertsize_median;

	length_high = length_median * 3.0;
	if(length_high > insertsize_high) length_high = insertsize_high;

	length_low = insertsize_low * 0.5;
	//length_low = length_median * 0.3;
	//if(length_low < insertsize_low) length_low = insertsize_low;
	return 0;
}

int bridger::set_circ_length()
{
	length_median = insertsize_median;
	length_high = 999999999;
	length_low = 0;

	/*length_median = insertsize_median;
	length_high = 500;
	length_low = 0;*/
	return 0;
}

int bridger::filter_paths(vector<fragment> &frags)
{
	for(int k = 0; k < frags.size(); k++)
	{
		fragment &fr = frags[k];

		if(fr.paths.size() <= 0) continue;
		
		int minp = -1;
		double mind = 9999999999999;

		for(int i = 0; i < fr.paths.size(); i++)
		{
			//if(fr.paths[i].type != 1 && fr.paths[i].type != 3) continue;
			if(fr.paths[i].length < length_low) continue;
			if(fr.paths[i].length > length_high) continue;
			double d = fabs(fr.paths[i].length - length_median);
			if(d >= mind) continue;
			mind = d;
			minp = i;
		}

		//fr.print(789);

		if(minp == -1)
		{
			fr.paths.clear();
			fr.set_bridged(false);
			//fr.h1->bridged = false;
			//fr.h2->bridged = false;
		}
		else
		{
			fr.paths[0] = fr.paths[minp];
			fr.paths.resize(1);
			assert(fr.paths.size() == 1);
			if(fr.paths[0].type == 1) fr.set_bridged(true);
			//fr.h1->bridged = true;
			//fr.h2->bridged = true;
		}
	}
	return 0;
}

int bridger::pick_bridge_path(vector<fragment> &frags)
{
	//printf("1.5*length_high = %lf\n",1.5*length_high);
	for(int k = 0; k < frags.size(); k++)
	{
		fragment &fr = frags[k];

		if(fr.paths.size() <= 0) 
		{
			fr.set_bridged(false);
			continue;
		}

		/*for(int i=0;i<fr.paths.size();i++)
		{
			int32_t len = fr.paths[i].length;
			
			if(strcmp(fr.h1->qname.c_str(),"ST-E00299:245:HKTJJALXX:6:2205:11647:15936") == 0)
			{
				printf("ST-E00299:245:HKTJJALXX:6:2205:11647:15936\n");
				vector<int> path_v = decode_vlist(fr.paths[i].v);
				printv(path_v);
				printf("score = %lf\n",fr.paths[i].score);
				printf("exon len: %d",len);

				if(fr.paths[i].type == 1 || fr.paths[i].type == 2)
				{
					printf(" ref path\n");
				}
				else
				{
					printf(" read path\n");
				}

				for(int j=0;j<path_v.size();j++)
				{
					printf("%d-%d, ",bd->regions[path_v[j]].lpos, bd->regions[path_v[j]].rpos);
				}
				printf("\n");
			}

			if(strcmp(fr.h1->qname.c_str(),"ST-E00299:245:HKTJJALXX:6:2107:25256:70486") == 0)
			{
				printf("ST-E00299:245:HKTJJALXX:6:2107:25256:70486\n");
				vector<int> path_v = decode_vlist(fr.paths[i].v);
				printv(path_v);
				printf("score = %lf\n",fr.paths[i].score);
				printf("exon len: %d",len);

				if(fr.paths[i].type == 1 || fr.paths[i].type == 2)
				{
					printf(" ref path\n");
				}
				else
				{
					printf(" read path\n");
				}

				for(int j=0;j<path_v.size();j++)
				{
					printf("%d-%d, ",bd->regions[path_v[j]].lpos, bd->regions[path_v[j]].rpos);
				}
				printf("\n");
			}

			if(strcmp(fr.h1->qname.c_str(),"simulate:689721") == 0)
			{
				printf("simulate:689721\n");
				vector<int> path_v = decode_vlist(fr.paths[i].v);
				printv(path_v);
				printf("score = %lf\n",fr.paths[i].score);
				printf("exon len: %d",len);

				if(fr.paths[i].type == 1 || fr.paths[i].type == 2)
				{
					printf(" ref path\n");
				}
				else
				{
					printf(" read path\n");
				}

				for(int j=0;j<path_v.size();j++)
				{
					printf("%d-%d, ",bd->regions[path_v[j]].lpos, bd->regions[path_v[j]].rpos);
				}
				printf("\n");
			}
		}*/

		//set regions, merged regions and junctions for paths
		for(int i=0;i<fr.paths.size();i++)
		{
			path *p1 = &fr.paths[i];
			//find juncs of p1

			vector<int> p1_v = decode_vlist(p1->v);

			//convert path vertices to regions
			for(int i=0;i<p1_v.size();i++)
			{
				p1->path_regions.push_back(bd->regions[p1_v[i]]);
			}

			//merge consecutive regions
			join_interval_map jmap;
			for(int i = 0; i < p1->path_regions.size(); i++)
			{
				int32_t left = p1->path_regions[i].lpos;
				int32_t right = p1->path_regions[i].rpos;
				jmap += make_pair(ROI(left, right), 1);
			}

			for(JIMI it = jmap.begin(); it != jmap.end(); it++)
			{
				region r(lower(it->first), upper(it->first), '.', '.');
				p1->merged_regions.push_back(r);
			}
			
			int intron_cnt = 0;
			for(int i=1;i<p1->merged_regions.size();i++)
			{
				if(p1->merged_regions[i].lpos != p1->merged_regions[i-1].rpos)
				{
					p1->junc_regions.push_back(pair<int32_t,int32_t>(p1->merged_regions[i-1].rpos+1, p1->merged_regions[i].lpos-1));
					intron_cnt++;
				}
			}
			p1->exon_count = intron_cnt + 1;

			/*if(strcmp(fr.h1->qname.c_str(),"E00511:127:HJN5NALXX:5:1105:12672:71260") == 0)
			{
				printv(p1_v);
				printf("\n");
				for(int i=0;i<p1->merged_regions.size();i++)
				{
					printf("%d-%d, ",p1->merged_regions[i].lpos, p1->merged_regions[i].rpos);
				}
				printf("\nfinding juncs\n");
				for(int i=0;i<p1->junc_regions.size();i++)
				{
					printf("%d-%d, ",p1->junc_regions[i].first,p1->junc_regions[i].second);
				}
				printf("\n");
			}*/
		}

		//print path information of fragments
		/*printf("\nPrinting path info:\n");
		printf("chrm = %s\n",bd->bb.chrm.c_str());
		printf("fragment read = %s, h1 pos = %d, h2 pos = %d\n",fr.h1->qname.c_str(),fr.h1->pos,fr.h2->pos);
		if((fr.h1->flag & 0x800) >= 1) printf("fr.h1 is supplementary\n");
		if((fr.h2->flag & 0x800) >= 1) printf("fr.h2 is supplementary\n");
		if(fr.h1->is_fake == true) printf("fr.h1 is fake\n");
		if(fr.h2->is_fake == true) printf("fr.h2 is fake\n");

		for(int i=0;i<fr.paths.size();i++)
		{
			path p1 = fr.paths[i];

			printf("path %d\n",i+1);
			printf("path vertices encoded: ");
			//vector<int> path_v = decode_vlist(p1.v);
			printv(p1.v);
			printf("\nscore = %lf,",p1.score);
			printf("length = %d,",p1.length);
			printf("type = %d,",p1.type);
			printf("kind = ");
			if(p1.type == 1 || p1.type == 2)
			{
				printf("ref path\n");
			}
			else
			{
				printf("read path\n");
			}
			
			printf("path regions: \n");
			for(int j=0;j<p1.path_regions.size();j++)
			{
				printf("%d-%d ",p1.path_regions[j].lpos, p1.path_regions[j].rpos);
				if(p1.path_regions[j].gapped == true)
				{
					printf("gapped true\n");
				}
				else
				{
					printf("gapped false\n");
				}
			}
		}
		printf("\n");*/

		vector<path> primary_selected_paths;
		vector<path> remove_list;
		remove_list.clear();
		primary_selected_paths.clear();

		//for read paths, discard path if middle region has gap
		for(int i=0;i<fr.paths.size();i++)
		{
			path p1 = fr.paths[i];

			if(p1.type == 1 || p1.type == 2) continue; //exclude ref paths

			if(p1.path_regions.size() > 1)
			{
				for(int j=1;j<p1.path_regions.size()-1;j++) //exclude first and last exon for gapped checking
				{
					region r = p1.path_regions[j];
					if(r.gapped == true)
					{
						remove_list.push_back(p1);
						break;
					}
				}
			}
		}

		//remove exon retention path, compare all paths pairwise
		for(int i=0;i<fr.paths.size();i++)
		{
			path p1 = fr.paths[i];

			//printf("p1 junc size; %lu\n",p1.junc_regions.size());

			for(int j=0;j<p1.junc_regions.size();j++)
			{
				pair<int32_t,int32_t> junc = p1.junc_regions[j];
				
				for(int s=0;s<fr.paths.size();s++)
				{
					path p2 = fr.paths[s];
					//printf("p2 merged size; %lu\n",p2.merged_regions.size());

					for(int t=0;t<p2.merged_regions.size();t++)
					{
						if(junc.first >= p2.merged_regions[t].lpos && junc.second <= p2.merged_regions[t].rpos)
						{
							remove_list.push_back(p2);
							break;
						}
					}
				}
			}
		}

		map<string,int> remove_map;
		remove_map.clear();

		for(int i=0;i<remove_list.size();i++)
		{
			path p = remove_list[i];
			string hash = "";
			
			for(int j=0;j<p.v.size();j++)
			{
				hash = hash + tostring(p.v[j]) + "|";
			}
			hash = hash + tostring(p.score) + "|" + tostring(p.length) + "|" + tostring(p.type);

			if(remove_map.find(hash) != remove_map.end()) //path present already in map
			{
				remove_map[hash]++;
			}
			else //path not present in map
			{
				remove_map.insert(pair<string,int>(hash,1));
			}
		}

		/*if(remove_map.size() > 0)
		{
			printf("initial remove map:\n");
			map<string,int>::iterator itn;
			for(itn = remove_map.begin(); itn != remove_map.end(); itn++)
			{
				printf("%s,",itn->first.c_str());
			}
			printf("\n");
		}*/

		//vector<path> selected_paths;
		//selected_paths.clear();

		for(int i=0;i<fr.paths.size();i++)
		{
			path p = fr.paths[i];
			string hash = "";
			for(int j=0;j<p.v.size();j++)
			{
				hash = hash + tostring(p.v[j]) + "|";
			}
			hash = hash + tostring(p.score) + "|" + tostring(p.length) + "|" + tostring(p.type);

			if(remove_map.find(hash) == remove_map.end()) //path not present in remove map
			{
				primary_selected_paths.push_back(p);
			}
		}

		//remove paths with score <= min_pathscore(1) if there exists higher score paths
		vector<path> selected_paths;
		selected_paths.clear();
		remove_list.clear();
		remove_map.clear();
		bool has_higher_score = false;

		for(int i=0;i<primary_selected_paths.size();i++)
		{
			path p = primary_selected_paths[i];
			if((p.type == 3 || p.type == 4) && p.score > min_pathscore) has_higher_score = true;
			break;
		}

		if(has_higher_score == true)
		{
			for(int i=0;i<primary_selected_paths.size();i++)
			{
				path p = primary_selected_paths[i];
				if((p.type == 3 || p.type == 4) && p.score <= min_pathscore) remove_list.push_back(p);
			}
		}

		for(int i=0;i<remove_list.size();i++)
		{
			path p = remove_list[i];
			string hash = "";
			
			for(int j=0;j<p.v.size();j++)
			{
				hash = hash + tostring(p.v[j]) + "|";
			}
			hash = hash + tostring(p.score) + "|" + tostring(p.length) + "|" + tostring(p.type);

			if(remove_map.find(hash) != remove_map.end()) //path present already in map
			{
				remove_map[hash]++;
			}
			else //path not present in map
			{
				remove_map.insert(pair<string,int>(hash,1));
			}
		}

		/*if(remove_map.size() > 0)
		{
			printf("final remove map:\n");
			map<string,int>::iterator itn;
			for(itn = remove_map.begin(); itn != remove_map.end(); itn++)
			{
				printf("%s,",itn->first.c_str());
			}
			printf("\n");
		}*/

		for(int i=0;i<primary_selected_paths.size();i++)
		{
			path p = primary_selected_paths[i];
			string hash = "";
			for(int j=0;j<p.v.size();j++)
			{
				hash = hash + tostring(p.v[j]) + "|";
			}
			hash = hash + tostring(p.score) + "|" + tostring(p.length) + "|" + tostring(p.type);

			if(remove_map.find(hash) == remove_map.end()) //path not present in remove map
			{
				selected_paths.push_back(p);
			}
		}

		//check if selected paths is zero
		if(selected_paths.size() == 0)
		{
			//printf("selected path size zero\n");
			//selected_paths = fr.paths;
			fr.set_bridged(false);
			fr.paths.resize(0);
			continue;
		}

		fr.candidate_path_count = selected_paths.size();

		// let A be the set of b-paths whose type is either 1 or 2 -- ref
		// let B be the set of b-paths whose type is either 3 or 4 -- reads

		vector<path> ref_paths;
		vector<path> read_paths;

		for(int i=0;i<selected_paths.size();i++)
		{
			path p = selected_paths[i];
			if(p.type == 1 || p.type == 2)
			{
				ref_paths.push_back(p);
			}
			else if(p.type == 3 || p.type == 4)
			{
				read_paths.push_back(p);
			}
		}

		if(read_paths.size() > 1)
		{
			//printf("ref_paths size = %lu, read_paths size = %lu\n",ref_paths.size(),read_paths.size());
		}
		//printf("fragment paths size: %lu\n",fr.paths.size());

		map<string,pair<path,int>> ref_paths_map;
		map<string,pair<path,int>> read_paths_map;

		for(int i=0;i<ref_paths.size();i++)
		{
			path p = ref_paths[i];
			string hash = "";

			for(int j=0;j<p.v.size();j++)
			{
				hash = hash + tostring(p.v[j]) + "|";
			}

			if(ref_paths_map.find(hash) != ref_paths_map.end()) //path present already in map
			{
				ref_paths_map[hash].second++;
			}
			else //path not present in map
			{
				ref_paths_map.insert(pair<string,pair<path,int>>(hash,pair<path,int>(p,1)));
			}
		}

		for(int i=0;i<read_paths.size();i++)
		{
			path p = read_paths[i];
			string hash = "";

			for(int j=0;j<p.v.size();j++)
			{
				hash = hash + tostring(p.v[j]) + "|";
			}

			if(read_paths_map.find(hash) != read_paths_map.end()) //path present already in map
			{
				read_paths_map[hash].second++;
			}
			else //path not present in map
			{
				read_paths_map.insert(pair<string,pair<path,int>>(hash,pair<path,int>(p,1)));
			}
		}

		//printing selected read/ref paths
		/*map<string, pair<path, int>>::iterator itn;
		for(itn = ref_paths_map.begin(); itn != ref_paths_map.end(); itn++)
		{
			printf("ref_path_key = %s, count = %d\n",itn->first.c_str(),itn->second.second);
		}
		for(itn = read_paths_map.begin(); itn != read_paths_map.end(); itn++)
		{
			printf("read_path_key = %s, count = %d\n",itn->first.c_str(),itn->second.second);
		}*/

		vector<path> intersection;
		map<string, pair<path, int>>::iterator itn1;
		map<string, pair<path, int>>::iterator itn2;

		for(itn1 = read_paths_map.begin(); itn1 != read_paths_map.end(); itn1++)
		{
			string p1 = itn1->first;
			for(itn2 = ref_paths_map.begin(); itn2 != ref_paths_map.end(); itn2++)
			{
				string p2 = itn2->first;
				if(strcmp(p1.c_str(),p2.c_str()) == 0)
				{
					intersection.push_back(itn1->second.first);
				}
			}
		}

		// find overlap betwen A and B
		int max_score = -1000000;
		path best_path;

		if(intersection.size() > 0)
		{	
			// if A overlaps with B
			// then we only consider the intersection
			// and we pick one whose score is maximized among 3/4 types
			for(int i=0;i<intersection.size();i++)
			{
				path p = intersection[i];
				if(p.score > max_score)
				{
					max_score = p.score;
					best_path = p;
				}
			}
		}
		else
		{	
			// if intersection is empty, we give priority to 3/4
			// in this case, pick one whose score is maximized among 3/4 types

			if(read_paths_map.size() > 0)
			{
				map<string, pair<path, int>>::iterator itn;
				for(itn = read_paths_map.begin(); itn != read_paths_map.end(); itn++)
				{
					if(itn->second.first.score > max_score)
					{
						max_score = itn->second.first.score;
						best_path = itn->second.first;
					}
				}
			}
			else
			{
				// if no 3/4 types, pick one randomly from 1/2
				// later on we can take the #counts in reference into account

				best_path = ref_paths_map.begin()->second.first;

			}
		}

		/*if(strcmp(fr.h1->qname.c_str(),"simulate:448707") == 0 && ref_paths_map.size() > 0)
		{
			best_path = ref_paths_map.begin()->second.first;
		}*/

		// printf("best path:\n");
		// printv(best_path.v);
		// printf("\n");

		fr.paths[0] = best_path;
		fr.paths.resize(1);
		assert(fr.paths.size() == 1);
		fr.set_bridged(true);
	}
	return 0;
}

int bridger::get_paired_fragments(vector<fragment> &frags)
{
	int n1 = 0;
	int n2 = 0;
	for(int k = 0; k < frags.size(); k++)
	{
		if(frags[k].paths.size() >= 1) n1++;
		if(frags[k].paths.size() != 1) continue;
		if(frags[k].paths[0].type != 1) continue;
		n2++;
	}
	return n1;
}

vector<int> bridger::get_bridged_fragments_type(vector<fragment> &frags)
{
	vector<int> ct(6, 0);
	
	for(int k = 0; k < frags.size(); k++)
	{
		if(frags[k].type == 0) 
		{
			ct[3]++;
			/*
			printf("paired-end fragments #%d, range = [%d, %d], h1 = %d, h2 = %d\n", k, frags[k].lpos, frags[k].rpos, frags[k].h1->hid, frags[k].h2->hid);
                        vector<int> v1 = decode_vlist(frags[k].h1->vlist);
                        vector<int> v2 = decode_vlist(frags[k].h2->vlist);
                        printf("h1 vlist: ");
                        for(int id_v1 = 0; id_v1 < v1.size(); id_v1++)
                        {
                                printf("%d ", v1[id_v1]);
                        }
                        printf("\nh2 vlist: ");
                        for(int id_v2 = 0; id_v2 < v2.size(); id_v2++)
                        {
                                printf("%d ", v2[id_v2]);
                        }
			printf("\n");
			*/


			if(frags[k].paths.size() == 1) 
			{
				//printf("Success bridge paired-end fragments\n\n");
				ct[0]++;
			}
			//else printf("Fail bridge paired-end fragments\n\n");
		}
		else if(frags[k].type == 1) 
		{
			/*
			printf("UMI-linked fragments #%d, umi = %s, range = [%d, %d], h1 = %d, h2 = %d\n", k, frags[k].h2->umi.c_str(), frags[k].lpos, frags[k].rpos, frags[k].h1->hid, frags[k].h2->hid);
			vector<int> v1 = decode_vlist(frags[k].h1->vlist);
			vector<int> v2 = decode_vlist(frags[k].h2->vlist);
			printf("h1 vlist: ");
			for(int id_v1 = 0; id_v1 < v1.size(); id_v1++)
			{
				printf("%d ", v1[id_v1]);
			}
                        printf("\nh2 vlist: ");
                        for(int id_v2 = 0; id_v2 < v2.size(); id_v2++)
                        {
                                printf("%d ", v2[id_v2]);
                        }
			printf("\n");
			*/

			ct[4]++;
			if(frags[k].paths.size() == 1)
			{
				//printf("Success bridge UMI-linked fragments\n\n");
				ct[1]++;
			}
			//else printf("Fail bridge UMI-linked fragments\n\n");
		}
		else if(frags[k].type == 2)
                {
                        ct[5]++;
                        if(frags[k].paths.size() == 1) ct[2]++;
                }

	}
	return ct;
}

int bridger::print(vector<fragment> &frags)
{
	int n = 0;
	/*
	for(int k = 0; k < fclusters.size(); k++)
	{
		n += fclusters[k].fset.size();
	}
	printf("#fragments = %lu, #open-fragments = %d, #fclusters = %lu\n", frags.size(), n, fclusters.size());
	*/

	for(int k = 0; k < frags.size(); k++)
	{
		if(frags[k].paths.size() == 1) n++;
	}

	int total = frags.size();
	int remain = total - n;
	double ratio = n * 100.0 / total;

	printf("#fragments = %d, #fixed = %d, #remain = %d, ratio = %.1lf, length = (%d, %d, %d)\n", total, n, remain, ratio, length_low, length_median, length_high);

	//for(int k = 0; k < fclusters.size(); k++) fclusters[k].print(k);
	return 0;
}

bool compare_fragment_v1(fragment *f1, fragment *f2)
{
	if(f1->h1->vlist.size() < f2->h1->vlist.size()) return true;
	if(f1->h1->vlist.size() > f2->h1->vlist.size()) return false;

	for(int k = 0; k < f1->h1->vlist.size(); k++)
	{
		if(f1->h1->vlist[k] < f2->h1->vlist[k]) return true;
		if(f1->h1->vlist[k] > f2->h1->vlist[k]) return false;
	}

	return (f1->lpos < f2->lpos);
}

bool compare_fragment_v2(fragment *f1, fragment *f2)
{
	if(f1->h2->vlist.size() < f2->h2->vlist.size()) return true;
	if(f1->h2->vlist.size() > f2->h2->vlist.size()) return false;

	for(int k = 0; k < f1->h2->vlist.size(); k++)
	{
		if(f1->h2->vlist[k] < f2->h2->vlist[k]) return true;
		if(f1->h2->vlist[k] > f2->h2->vlist[k]) return false;
	}

	return (f1->lpos < f2->lpos);
}

bool compare_fragment_v3(fragment *f1, fragment *f2)
{
	if(f1->h1->vlist.size() < f2->h1->vlist.size()) return true;
	if(f1->h1->vlist.size() > f2->h1->vlist.size()) return false;
	if(f1->h2->vlist.size() < f2->h2->vlist.size()) return true;
	if(f1->h2->vlist.size() > f2->h2->vlist.size()) return false;

	for(int k = 0; k < f1->h1->vlist.size(); k++)
	{
		if(f1->h1->vlist[k] < f2->h1->vlist[k]) return true;
		if(f1->h1->vlist[k] > f2->h1->vlist[k]) return false;
	}

	for(int k = 0; k < f1->h2->vlist.size(); k++)
	{
		if(f1->h2->vlist[k] < f2->h2->vlist[k]) return true;
		if(f1->h2->vlist[k] > f2->h2->vlist[k]) return false;
	}

	return (f1->lpos < f2->lpos);
}

bool compare_fragment_v3_flank(fragment *f1, fragment *f2)
{
	if(f1->h1->vlist.size() < f2->h1->vlist.size()) return true;
	if(f1->h1->vlist.size() > f2->h1->vlist.size()) return false;
	if(f1->h2->vlist.size() < f2->h2->vlist.size()) return true;
	if(f1->h2->vlist.size() > f2->h2->vlist.size()) return false;

	for(int k = 0; k < f1->h1->vlist.size(); k++)
	{
		if(f1->h1->vlist[k] < f2->h1->vlist[k]) return true;
		if(f1->h1->vlist[k] > f2->h1->vlist[k]) return false;
	}

	for(int k = 0; k < f1->h2->vlist.size(); k++)
	{
		if(f1->h2->vlist[k] < f2->h2->vlist[k]) return true;
		if(f1->h2->vlist[k] > f2->h2->vlist[k]) return false;
	}

	if(f1->k1l + f1->k2l < f2->k1l + f2->k2l) return true;
	if(f1->k1l + f1->k2l > f2->k1l + f2->k2l) return false;

	return (f1->lpos < f2->lpos);
}

bool compare_fragment_path(fragment *f1, fragment *f2)
{
	assert(f1->paths.size() >= 1);
	assert(f2->paths.size() >= 1);

	int n1 = f1->paths[0].v.size();
	int n2 = f2->paths[0].v.size();

	for(int k = 0; k < n1 && k < n2; k++)
	{
		if(f1->paths[0].v[k] < f2->paths[0].v[k]) return true;
		if(f1->paths[0].v[k] > f2->paths[0].v[k]) return false;
	}

	if(n1 < n2) return true;
	if(n1 > n2) return false;

	return (f1->lpos < f2->lpos);
}

bool check_suffix(const vector<int> &vx, const vector<int> &vy) 
{
	if(vx.size() == 0) return true;
	if(vy.size() == 0) return true;

	double overlap = 0.4;
	if(vx.size() <= vy.size())
	{
		vector<int>::const_iterator it = lower_bound(vx.begin(), vx.end(), vy[0]);
		if(it == vx.end()) return false;
		if(vx.end() - it < overlap * vx.size()) return false;
		for(int kx = vx.size() - 1, ky = (vx.end() - it) - 1; kx >= 0 && ky >= 0; kx--, ky--)
		{
			if(vx[kx] != vy[ky]) return false;
		}
		/*
		for(int k = 0; it != vx.end(); it++, k++)
		{
			assert(k < vy.size());
			if((*it) != vy[k]) return false;
		}
		*/
	}
	else
	{
		vector<int>::const_iterator it = lower_bound(vy.begin(), vy.end(), vx.back());
		if(it == vy.end()) return false;
		if(it - vy.begin() + 1 < overlap * vy.size()) return false;
		for(int kx = vx.size() - (it - vy.begin()) - 1, ky = 0; kx < vx.size() && ky < vy.size(); kx++, ky++)
		{
			if(vx[kx] != vy[ky]) return false;
		}
		/*
		for(int k = vx.size() - 1; k >= 0; it--, k--)
		{
			if((*it) != vx[k]) return false;
			if(it == vy.begin()) break;
		}
		*/
	}
	return true;
}

