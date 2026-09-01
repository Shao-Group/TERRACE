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
#include "generated_pool_gate.h"
#include "generated_selection_scorer.h"
#include <iostream>
#include <algorithm>
#include <set>
#include <climits>
#include <cfloat>
#include <queue>
#include <tuple>
#include <cctype>
#include <sys/stat.h>
#include <dirent.h>
#include <ctime>
#include <cstdlib>

map<int, vector<read_info>> bridger::bundle_chimeric_reads;
map<int, vector<read_info>> bridger::bundle_outward_reads;
map<int, map<string, chimeric_bsj_info>> bridger::bundle_chimeric_bsj;
map<int, map<string, vector<string>>> bridger::bundle_chimeric_support_names;
map<int, map<string, vector<pair<vector<int>, vector<int>>>>> bridger::bundle_chimeric_support_paths;
map<int, map<string, map<pair<int32_t,int32_t>, vector<vector<int>>>>> bridger::bundle_outward_bsj_paths;
map<int, map<string, vector<vector<int>>>> bridger::bundle_chimeric_merged_paths;

const int TOLERANCE_GAP = 0;

// Primary assembled chromosomes (chr1-22, chrX, chrY). Unplaced/decoy contigs
// (GL*, chrUn_*, *_random, *_alt, etc.) and chrM are excluded from circRNA
// output: they are dominated by highly repetitive sequence that produces
// spurious chimeric/outward alignments, and no annotated circRNA ground truth
// is ever placed on them.
static bool is_primary_chromosome(const string &chrm)
{
	static const set<string> primary = []() {
		set<string> s;
		for (int i = 1; i <= 22; i++)
		{
			s.insert("chr" + to_string(i));
		}
		s.insert("chrX");
		s.insert("chrY");
		return s;
	}();

	return primary.count(chrm) > 0;
}

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

	if(v1.size() <= ex1) return 0;
	if(v2.size() <= ex2) return 0;

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
		if(fr.h1->vlist.size() < 2) continue;
		if(fr.h2->vlist.size() < 2) continue;
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

	for (int i = 0; i < pnodes.size(); i++)
	{
		vector<int> &v = pnodes[i].v;
		assert(v.size() <= 2);
		if (v.size() <= 1) continue;
		int w = pnodes[i].score;
		int x = v[0];
		int y = v[1];

		assert(jsetx[x].find(y) == jsetx[x].end());
		assert(jsety[y].find(x) == jsety[y].end());
		jsetx[x].insert(PI(y, w));
		jsety[y].insert(PI(x, w));
	}

    splice_graph_adj.clear();

    for (int i = 0; i < jsetx.size(); i++)
    {
        if (!jsetx[i].empty())
        {
            vector<pair<int, int>> neighbors;
            for (pair<const int, int> const& neighbor : jsetx[i])
                neighbors.push_back(make_pair(neighbor.first, neighbor.second));
            splice_graph_adj[i] = neighbors;
        }
    }

	return 0;
}

vector<int> bridger::get_bsj_vlist_for_reads(int bsj_start_pos, int bsj_end_pos, const vector<read_info>& alignments_for_qname) {
    vector<int> bsj_vlist_calculated;
    int start_vertex_idx = -1;
    int end_vertex_idx = -1;

    for (const read_info& read : alignments_for_qname)
	{
        if (read.pos == bsj_start_pos && !read.vlist.empty())
		{
            if (start_vertex_idx == -1 || read.vlist.front() < start_vertex_idx)
                start_vertex_idx = read.vlist.front();
        }
    }

    for (const read_info& read : alignments_for_qname)
	{
        if (read.rpos == bsj_end_pos && !read.vlist.empty())
		{
            if (end_vertex_idx == -1 || read.vlist.back() > end_vertex_idx)
                end_vertex_idx = read.vlist.back();
        }
    }

    if (start_vertex_idx != -1 && end_vertex_idx != -1)
	{
        bsj_vlist_calculated.push_back(start_vertex_idx);
        bsj_vlist_calculated.push_back(end_vertex_idx);
    }

    return bsj_vlist_calculated;
}

static bool is_vlist_spanning_introns_only(const vector<int>& vlist, const vector<region>& regions)
{
    if (vlist.empty()) return false;
    
	for (int i = 0; i < vlist.size(); i++)
    {
        int v = vlist[i];
		if (v < 0 || v >= regions.size()) return false;
		if (regions[v].ltype != LEFT_SPLICE || regions[v].rtype != RIGHT_SPLICE) return false;
    }
	
    return true;
}

int bridger::collect_bundle_reads(int bundle_index)
{
    vector<read_info> ow_reads;
    for (const fragment& fr : bd->outward_fragments)
	{
        if (fr.h1)
		{
			vector<int> vlist1 = decode_vlist(fr.h1->vlist);
            if (!is_vlist_spanning_introns_only(vlist1, bd->regions))
            {
                bool is_r1 = fr.h1->flag & 0x40;
                bool is_r2 = fr.h1->flag & 0x80;
                bool is_supp = fr.h1->flag & 0x800;
                bool is_prim = !is_supp;
                bool is_rev = fr.h1->flag & 0x10;
                ow_reads.push_back({fr.h1->qname, static_cast<int>(fr.h1->pos), static_cast<int>(fr.h1->rpos), vlist1, is_r1, is_r2, is_prim, is_supp, is_rev, fr.h1});
            }
        }

        if (fr.h2)
		{
			vector<int> vlist2 = decode_vlist(fr.h2->vlist);
			if (!is_vlist_spanning_introns_only(vlist2, bd->regions))
            {
                bool is_r1 = fr.h2->flag & 0x40;
                bool is_r2 = fr.h2->flag & 0x80;
                bool is_supp = fr.h2->flag & 0x800;
                bool is_prim = !is_supp;
                bool is_rev = fr.h2->flag & 0x10;
                ow_reads.push_back({fr.h2->qname, static_cast<int>(fr.h2->pos), static_cast<int>(fr.h2->rpos), vlist2, is_r1, is_r2, is_prim, is_supp, is_rev, fr.h2});
            }
        }
    }

    bundle_outward_reads[bundle_index] = ow_reads;

    map<string, vector<int>> qname_to_hit_indices;
    for (int i = 0; i < bd->bb.hits.size(); i++)
        qname_to_hit_indices[bd->bb.hits[i].qname].push_back(i);

    map<string, vector<int>> qname_to_fake_hit_indices;
    for (int i = 0; i < bd->bb.fake_hits.size(); i++)
        qname_to_fake_hit_indices[bd->bb.fake_hits[i].qname].push_back(i);

    vector<read_info> ch_reads;
    set<string> processed_qnames;
    for (int i = 0; i < bd->bb.hits.size(); i++)
	{
        hit &h = bd->bb.hits[i];
        if (h.suppl != NULL && processed_qnames.find(h.qname) == processed_qnames.end())
		{
            vector<read_info> alignments_for_qname;
            int min_pos = INT_MAX;
            int max_rpos = 0;

            map<string, vector<int>>::iterator hit_idx_it = qname_to_hit_indices.find(h.qname);
            if (hit_idx_it != qname_to_hit_indices.end())
            {
                for (int j : hit_idx_it->second)
				{
                    hit &chim_h = bd->bb.hits[j];
                    vector<int> chim_vlist = decode_vlist(chim_h.vlist);
					if (is_vlist_spanning_introns_only(chim_vlist, bd->regions)) continue;

                    bool is_r1 = chim_h.flag & 0x40;
                    bool is_r2 = chim_h.flag & 0x80;
                    bool is_supp = chim_h.flag & 0x800;
                    bool is_prim = !is_supp;
                    bool is_rev = chim_h.flag & 0x10;

                    read_info current_read_info = {chim_h.qname, static_cast<int>(chim_h.pos), static_cast<int>(chim_h.rpos), chim_vlist, is_r1, is_r2, is_prim, is_supp, is_rev, &chim_h};
                    ch_reads.push_back(current_read_info);
                    alignments_for_qname.push_back(current_read_info);

                    min_pos = min(min_pos, current_read_info.pos);
                    max_rpos = max(max_rpos, current_read_info.rpos);
                }
            }

            map<string, vector<int>>::iterator fake_idx_it = qname_to_fake_hit_indices.find(h.qname);
            for (int j : (fake_idx_it != qname_to_fake_hit_indices.end() ? fake_idx_it->second : vector<int>()))
            {
                hit &fake_h = bd->bb.fake_hits[j];
                vector<int> fake_vlist = decode_vlist(fake_h.vlist);
                if (fake_vlist.empty()) continue;
				if (is_vlist_spanning_introns_only(fake_vlist, bd->regions)) continue;

                bool is_r1 = fake_h.flag & 0x40;
                bool is_r2 = fake_h.flag & 0x80;
                bool is_supp = true;
                bool is_prim = false;
                bool is_rev = fake_h.flag & 0x10;
				
                read_info fake_read_info = {fake_h.qname, static_cast<int>(fake_h.pos), static_cast<int>(fake_h.rpos), fake_vlist, is_r1, is_r2, is_prim, is_supp, is_rev, &fake_h};
                ch_reads.push_back(fake_read_info);
                alignments_for_qname.push_back(fake_read_info);

                min_pos = min(min_pos, fake_read_info.pos);
                max_rpos = max(max_rpos, fake_read_info.rpos);
            }

            processed_qnames.insert(h.qname);

            vector<int> bsj_vlist_calculated = get_bsj_vlist_for_reads(min_pos, max_rpos, alignments_for_qname);
            bundle_chimeric_bsj[bundle_index][h.qname] = {min_pos, max_rpos, bsj_vlist_calculated};
        }
    }

    bundle_chimeric_reads[bundle_index] = ch_reads;

    return 0;
}

void bridger::find_all_paths(int start, int end, vector<int>& current_path, vector<vector<int>>& all_paths, int max_paths)
{
	if (all_paths.size() >= max_paths) return;
	
	current_path.push_back(start);
	
	if (start == end)
		all_paths.push_back(current_path);
	
	else
	{
		map<int, vector<pair<int, int>>>::iterator it = splice_graph_adj.find(start);

		if (it != splice_graph_adj.end())
		{
			for (const pair<int, int>& neighbor_pair : it->second)
			{
				int neighbor = neighbor_pair.first;
				if (neighbor <= end)
					find_all_paths(neighbor, end, current_path, all_paths, max_paths);
			}
		}
	}

	current_path.pop_back();
}

static int find_min_edge_weight_along_path(const vector<int>& path, const map<int, vector<pair<int, int>>>& adj)
{
	int min_ew = INT_MAX;
	for (int j = 0; j + 1 < path.size(); j++)
	{
		int u = path[j];
		int v = path[j + 1];

		map<int, vector<pair<int, int>>>::const_iterator it = adj.find(u);
		if (it != adj.end())
		{
			for (const pair<int, int>& np : it->second)
			{
				if (np.first == v)
				{
					if (np.second < min_ew)
						min_ew = np.second;
					break;
				}
			}
		}
	}

	if (min_ew == INT_MAX)
		return -1;
	else
		return min_ew;
}

// Widest path (maximum-bottleneck route: the path whose weakest edge is as strong as possible, ties broken by fewest hops) from start to end over splice_graph_adj.
vector<int> bridger::widest_path(int start, int end)
{
	if (start == end) return vector<int>{start};

	map<int, int> best;   // vertex -> strongest bottleneck weight reachable from start
	map<int, int> hops;   // vertex -> number of edges on that best path
	map<int, int> prev;
	set<int> finalized;

	best[start] = INT_MAX;
	hops[start] = 0;

	// (bottleneck, -hops, vertex): max-heap pops the strongest bottleneck first, then (on ties) the path with fewest hops.
	priority_queue<tuple<int, int, int>> pq;
	pq.push(make_tuple(best[start], 0, start));

	while (!pq.empty())
	{
		tuple<int, int, int> top = pq.top();
		pq.pop();

		int u = get<2>(top);
		if (finalized.count(u)) continue;
		finalized.insert(u);
		if (u == end) break;
		map<int, vector<pair<int, int>>>::iterator it = splice_graph_adj.find(u);
		if (it == splice_graph_adj.end()) continue;

		for (const pair<int, int>& np : it->second)
		{
			int v = np.first;
			int w = np.second;

			if (finalized.count(v)) continue;

			int cand_best = min(best[u], w);
			int cand_hops = hops[u] + 1;

			map<int, int>::iterator bit = best.find(v);
			bool better = (bit == best.end()) || (cand_best > bit->second) || (cand_best == bit->second && cand_hops < hops[v]);

			if (better)
			{
				best[v] = cand_best;
				hops[v] = cand_hops;
				prev[v] = u;
				pq.push(make_tuple(cand_best, -cand_hops, v));
			}
		}
	}

	if (!best.count(end)) return vector<int>{};

	vector<int> path;
	int v = end;

	while (true)
	{
		path.push_back(v);
		if (v == start) break;
		map<int, int>::iterator pit = prev.find(v);
		if (pit == prev.end()) return vector<int>{};
		v = pit->second;
	}

	reverse(path.begin(), path.end());
	return path;
}

void bridger::write_splice_graph(std::ofstream& fout)
{
	for (int i = 0; i < bd->regions.size(); i++)
	{
		region& r = bd->regions[i];
		fout << "[" << i << "] " << "(" << bd->bb.chrm << ":" << r.lpos << "-" << r.rpos << ") → ";

		map<int, vector<pair<int, int>>>::iterator it = splice_graph_adj.find(i);

		if (it != splice_graph_adj.end() && !it->second.empty())
		{
			const vector<pair<int, int>>& neighbors = it->second;
			for (int j = 0; j < neighbors.size(); j++)
			{
				fout << neighbors[j].first << " (w = " << neighbors[j].second << ")";
				if (j < neighbors.size() - 1)
					fout << ", ";
			}
		}

		else
			fout << "∅";

		fout << "\n";
	}
}

void bridger::write_chimeric_read_paths(int bundle_idx, const std::string& outdir)
{
    std::string chim_filename = outdir + "/chimeric_possible_configs.txt";
    std::ofstream cfout(chim_filename, std::ios_base::app);

    cfout << "=================================================\n";
    cfout << "Bundle " << bundle_idx << "\n";
    cfout << "=================================================\n\n";

    cfout << "Splice Graph\n";
    cfout << "---------------------------------\n";
    write_splice_graph(cfout);
    cfout << "\n";

    map<string, vector<read_info>> grouped_chimeric;
    for (const read_info& read : bundle_chimeric_reads[bundle_idx])
        grouped_chimeric[read.qname].push_back(read);

    cfout << "Chimeric Reads (" << grouped_chimeric.size() << ")\n";
    cfout << "---------------------------------\n";

    auto process_chimeric_variant = [&](vector<read_info> reads, const string& qname, const string& variant_label) -> bool
    {
        if (reads.size() < 2) return false;

        sort(
            reads.begin(), reads.end(), [](const read_info& a, const read_info& b)
            {
                return a.pos < b.pos;
            }
        );

        int n = reads.size();

        cfout << "\t[" << variant_label << "]\n";

        for (const read_info& read : reads)
        {
            string s = "\t\t" + bd->bb.chrm + ":" + to_string(read.pos) + "-" + to_string(read.rpos) + "  (";
            for (int i = 0; i < read.vlist.size(); i++)
                s += to_string(read.vlist[i]) + (i == read.vlist.size() - 1 ? "" : " → ");
            s += ")";

            string lbl;
            if (read.is_read1)
				lbl = "R1";
            else if (read.is_read2)
				lbl = "R2";
            else
				lbl = "Unknown";

            lbl += " ";

            if (read.is_primary)
				lbl += "Primary";
            else if (read.is_supplementary)
				lbl += "Supplementary";
            else
				lbl += "Unknown";

            cfout << s << "  (" << lbl << ")\n";
        }

        bool any_empty = false;
        for (const read_info& r : reads)
		{
			if (r.vlist.empty())
			{
				any_empty = true;
				break;
			}
		}

        if (any_empty)
		{
			cfout << "\n";
			return false;
		}

        vector<vector<CandidatePath>> bridges(n - 1);
        bool bridge_failed = false;
        for (int i = 0; i < n - 1; i++)
        {
            int v_back  = reads[i].vlist.back();
            int v_front = reads[i + 1].vlist.front();
            if (v_back > v_front)
			{
				CandidatePath empty_bridge;
				empty_bridge.path = {};
				empty_bridge.min_edge_weight = INT_MAX;
				empty_bridge.min_vertex_weight = -1.0;
				empty_bridge.length = 0;
				bridges[i] = {empty_bridge};   // BSJ crossing (no bridge needed)
			}
            else
			{
				bridges[i] = find_pareto_paths(v_back, v_front, -1, true);
				if (bridges[i].empty())
				{
					bridge_failed = true;
					break;
				}
			}
        }

        if (bridge_failed)
		{
			cfout << "\t\t[skipped: no connecting path exists for at least one segment gap]\n\n";
			return false;
		}

		// Cartesian product across every gap's non-dominated bridge set (odometer-style index
		// increment; n is typically 2 segments/1 gap, occasionally 3/2 gaps).
		vector<int> idx(bridges.size(), 0);
		vector<vector<int>> merged_paths;
		while (true)
		{
			vector<int> merged = reads[0].vlist;
			for (int i = 0; i < (int)bridges.size(); i++)
			{
				const CandidatePath& b = bridges[i][idx[i]];
				if (b.path.size() >= 2)
					merged.insert(merged.end(), b.path.begin() + 1, b.path.end() - 1);
				merged.insert(merged.end(), reads[i + 1].vlist.begin(), reads[i + 1].vlist.end());
			}

			set<int> seen;
			vector<int> deduped;
			for (int v : merged)
			{
				if (seen.insert(v).second)
					deduped.push_back(v);
			}
			merged_paths.push_back(std::move(deduped));

			int k = (int)bridges.size() - 1;
			while (k >= 0)
			{
				idx[k]++;
				if (idx[k] < (int)bridges[k].size()) break;
				idx[k] = 0;
				k--;
			}
			if (k < 0) break;
		}

        cfout << "\n\t\tSegments: " << n << "\n";
		for (int i = 0; i < (int)bridges.size(); i++)
			cfout << "\t\tGap " << i << " candidates: " << bridges[i].size() << "\n";
        cfout << "\t\tCombinations: " << merged_paths.size() << "\n\n";

		for (const vector<int>& merged : merged_paths)
			bundle_chimeric_merged_paths[bundle_idx][qname].push_back(merged);
			return true;
    };

    for (pair<const string, vector<read_info>>& chim_entry : grouped_chimeric)
    {
        const vector<read_info>& all_reads_for_qname = chim_entry.second;

        // Require at least one supplementary alignment as BSJ evidence.
        bool has_suppl = false;
        bool suppl_is_read1 = false;
        bool suppl_is_read2 = false;
        for (const read_info& r : all_reads_for_qname)
		{
            if (r.is_supplementary)
			{
				has_suppl = true;
				suppl_is_read1 = r.is_read1;
				suppl_is_read2 = r.is_read2;
				break;
			}
        }

        if (!has_suppl || all_reads_for_qname.size() < 2) continue;

        cfout << chim_entry.first << "\n";

        bool ok = process_chimeric_variant(all_reads_for_qname, chim_entry.first, "all reads for qname");

        if (!ok)
        {
            vector<read_info> chimeric_segments;
            for (const read_info& r : all_reads_for_qname)
            {
                if (r.is_read1 == suppl_is_read1 && r.is_read2 == suppl_is_read2)
                    chimeric_segments.push_back(r);
            }

            if (chimeric_segments.size() != all_reads_for_qname.size())
                process_chimeric_variant(chimeric_segments, chim_entry.first, "mate-excluded fallback after failed mate-included bridge");
        }

        cfout << "\n";
    }

    cfout.close();

}

void bridger::write_simplified_chimeric_paths(int bundle_idx, const std::string& outdir)
{
    std::string simplify_chim_filename = outdir + "/simplify_chimeric_full_sequence_paths.txt";
    std::ofstream scfout(simplify_chim_filename, std::ios_base::app);

    scfout << "=================================================\n";
    scfout << "Bundle " << bundle_idx << "\n";
    scfout << "=================================================\n\n";

    map<int, map<string, vector<vector<int>>>>::iterator chim_bundle_it = bundle_chimeric_merged_paths.find(bundle_idx);
    if (chim_bundle_it != bundle_chimeric_merged_paths.end())
    {
        for (const pair<const string, vector<vector<int>>>& qname_entry : chim_bundle_it->second)
        {
            const string& qname = qname_entry.first;
            scfout << qname << "\n";
            scfout << "---------------------------------\n";

            vector<read_info> segs;
            for (const read_info& read : bundle_chimeric_reads[bundle_idx])
            {
                if (read.qname == qname)
                    segs.push_back(read);
            }

            sort(
				segs.begin(), segs.end(), [](const read_info& a, const read_info& b)
				{
					return a.pos < b.pos;
				}
			);

            for (const read_info& seg : segs)
            {
                string lbl;
                if (seg.is_read1)
					lbl = "R1";

                else if (seg.is_read2)
					lbl = "R2";

                else
					lbl = "Unknown";

                lbl += " ";
                
				if (seg.is_primary)
					lbl += "Primary";
                
				else if (seg.is_supplementary)
					lbl += "Supplementary";
                
				else
					lbl += "Unknown";

                scfout << "\t" << lbl << ": (";
                for (int i = 0; i < seg.vlist.size(); i++)
                    scfout << seg.vlist[i] << (i == seg.vlist.size() - 1 ? "" : " → ");

                scfout << ")\n";
            }

            scfout << "\n";

            vector<vector<int>> unique_paths;
            set<vector<int>> seen_paths;
            for (const vector<int>& path : qname_entry.second)
            {
                if (seen_paths.find(path) == seen_paths.end())
                {
                    seen_paths.insert(path);
                    unique_paths.push_back(path);
                }
            }

            scfout << "\tSupporting Paths (" << unique_paths.size() << "):\n";
            for (int i = 0; i < unique_paths.size(); i++)
            {
                scfout << "\t\tP" << i + 1 << " = (";
                for (int j = 0; j < unique_paths[i].size(); j++)
                    scfout << unique_paths[i][j] << (j == unique_paths[i].size() - 1 ? "" : " → ");
                scfout << ")\n";
            }

            scfout << "\n";
        }
    }

    else
        scfout << "No full sequence paths\n";

    scfout << "\n\n";
    scfout.close();
}

/*
Given the 3 segments of a chimeric read (sorted ascending by pos), identify:
   S – the single-alignment read (the mate with exactly one alignment)
   T_first – the T-read alignment hit first when traveling from S
   T_second – the T-read alignment hit second
   start_seg / mid_seg / end_seg – the three segments in left-to-right order for constructing the path
   frag_start / frag_end – actual read boundary positions

Algorithm:
   Forward S: path expressed L→R is S → T_first → T_second
     frag_start = S.pos, frag_end = T_second.rpos
   Backward S:  path expressed L→R is T_second → T_first → S
     frag_start = T_second.pos, frag_end = S.rpos

T_first / T_second are determined by clockwise traversal order from S along the circular path (stored as seg1 → seg2 → seg3 → [BSJ] → seg1):
   S=seg1 fwd → right: seg2 (T_first), seg3 (T_second)
   S=seg2 fwd → right: seg3 (T_first), seg1 (T_second)
   S=seg3 fwd → right: seg1 (T_first), seg2 (T_second)
   S=seg1 bwd → left: seg3 (T_first), seg2 (T_second)
   S=seg2 bwd → left: seg1 (T_first), seg3 (T_second)
   S=seg3 bwd → left: seg2 (T_first), seg1 (T_second)
*/

static bool identify_chimeric_traversal(const read_info& seg1, const read_info& seg2, const read_info& seg3, const read_info*& S, const read_info*& T_first, const read_info*& T_second, const read_info*& start_seg, const read_info*& mid_seg, const read_info*& end_seg, int& frag_start, int& frag_end)
{
	int r1 = 0; 
	int r2 = 0;
	for (const read_info* r : {&seg1, &seg2, &seg3})
	{
		if (r->is_read1)
			r1++;

		else if (r->is_read2)
			r2++;
	}

	S = nullptr;
	if (r1 == 1 && r2 == 2)
	{
		for (const read_info* r : {&seg1, &seg2, &seg3})
		{
			if (r->is_read1)
			{
				S = r;
				break;
			}
		}
	}

	else if (r1 == 2 && r2 == 1)
	{
		for (const read_info* r : {&seg1, &seg2, &seg3})
		{
			if (r->is_read2)
			{
				S = r;
				break;
			}
		}
	}

	if (!S) return false;

	// Collect T alignments in position order (ascending, since seg1 < seg2 < seg3).
	const read_info* T_left = nullptr;
	const read_info* T_right = nullptr;
	for (const read_info* r : {&seg1, &seg2, &seg3})
	{
		if (r == S) continue;

		if (!T_left)
			T_left = r;
		else
			T_right = r;
	}

	bool fwd = !S->is_reverse;

	if (S == &seg1 && fwd)
	{
		T_first = T_left;
		T_second = T_right;
	}

	else if (S == &seg2 && fwd)
	{
		T_first = T_right;
		T_second = T_left;
	}

	else if (S == &seg3 && fwd)
	{
		T_first = T_left;
		T_second = T_right;
	}

	else if (S == &seg1 && !fwd)
	{
		T_first = T_right;
		T_second = T_left;
	}

	else if (S == &seg2 && !fwd) 
	{ 
		T_first = T_left;
		T_second = T_right; 
	}
	
	else {                         
		T_first = T_right;
		T_second = T_left;  
	}

	if (fwd)
	{
		start_seg = S;
		mid_seg = T_first;
		end_seg = T_second;
		frag_start = S->pos;
		frag_end = T_second->rpos;
	}

	else
	{
		start_seg = T_second;
		mid_seg = T_first;
		end_seg = S;
		frag_start = T_second->pos;
		frag_end = S->rpos;
	}

	return true;
}

static vector<int> extract_chimeric_subpath(const vector<int>& path, int start_v, int end_v)
{
	int i_start = -1;
	int i_end = -1;

	for (int j = 0; j < path.size(); j++)
	{
		if (path[j] == start_v && i_start < 0)
			i_start = j;

		if (path[j] == end_v && i_end < 0)
			i_end = j;
	}

	if (i_start < 0 || i_end < 0) return {};

	vector<int> result;

	if (i_start == i_end)
	{
		result.insert(result.end(), path.begin() + i_start, path.end());
		result.insert(result.end(), path.begin(), path.begin() + i_start);
		result.push_back(path[i_start]);
	}

	// Does not cross the BSJ
	else if (i_start < i_end)
		result.assign(path.begin() + i_start, path.begin() + i_end + 1);

	// Wraps around and crosses the BSJ
	else
	{
		result.insert(result.end(), path.begin() + i_start, path.end());
		result.insert(result.end(), path.begin(), path.begin() + i_end + 1);
	}

	return result;
}

static int compute_chimeric_insert_size(const vector<int>& path, int frag_start, int frag_end, const vector<region>& regions, vector<int>* parts = nullptr)
{
    int len = 0;
    for (int i = 0; i < path.size(); i++)
    {
        int v = path[i];
		int lpos;
		int rpos;

		if (i == 0)
			lpos = frag_start;

		else
			lpos = regions[v].lpos;

		if (i == path.size() - 1)
			rpos = frag_end;

		else
			rpos = regions[v].rpos;

        int contrib = rpos - lpos;
        len += contrib;

        if (parts)
			parts->push_back(contrib);
    }

    return len;
}

// The "single" read S (the mate with exactly one alignment) sets the traversal direction.
// Starting from S's entry end, we travel along the circular path: we cross T_first (the first T alignment encountered) and continue to T_second (the second), stopping at T_second's far end.
// The total exon length of this traversal is the insert size.
void bridger::write_chimeric_insert_sizes(int bundle_idx, const std::string& outdir)
{
    map<string, vector<read_info>> grouped;
    for (const read_info& read : bundle_chimeric_reads[bundle_idx])
        grouped[read.qname].push_back(read);

    bool has_output = false;
    for (const pair<const string, vector<read_info>>& e : grouped)
    {
        for (const read_info& r : e.second)
        {
            if (r.is_supplementary)
			{
				has_output = true;
				break;
			}
        }

        if (has_output) break;
    }

    if (!has_output) return;

    std::string filename = outdir + "/chimeric_insert_sizes.txt";
    std::ofstream fout(filename, std::ios_base::app);

    fout << "=================================================\n";
    fout << "Bundle " << bundle_idx << "\n";
    fout << "=================================================\n\n";

    for (pair<const string, vector<read_info>>& chim_entry : grouped)
    {
        vector<read_info>& reads = chim_entry.second;

        bool has_suppl = false;
        for (const read_info& r : reads)
		{
            if (r.is_supplementary)
			{
				has_suppl = true;
				break;
			}
		}

        if (!has_suppl) continue;

        sort(
            reads.begin(), reads.end(), [](const read_info& a, const read_info& b)
            {
                return a.pos < b.pos;
            }
        );

        int n = reads.size();

        fout << "\t" << chim_entry.first << "\n";

        for (const read_info& read : reads)
        {
            string lbl;
            if (read.is_read1)
				lbl = "R1";
            else if (read.is_read2)
				lbl = "R2";
            else
				lbl = "Unknown";

            lbl += " ";

            if (read.is_primary)
				lbl += "Primary";
            else if (read.is_supplementary)
				lbl += "Supplementary";
            else
				lbl += "Unknown";

            fout << "\t\t" << bd->bb.chrm << ":" << read.pos << "-" << read.rpos << "  (";
            for (int i = 0; i < read.vlist.size(); i++)
            {
                fout << read.vlist[i];
                if (i < read.vlist.size() - 1)
					fout << " → ";
            }
            fout << ")  (" << lbl << ")\n";
        }

        bool any_empty = false;
        for (const read_info& r : reads)
		{
            if (r.vlist.empty())
			{
				any_empty = true;
				break;
			}
		}

        if (any_empty)
		{
			fout << "\n";
			continue;
        }

        vector<int> merged;
        merged.insert(merged.end(), reads[0].vlist.begin(), reads[0].vlist.end());

        for (int i = 0; i < n - 1; i++)
        {
            int v_back = reads[i].vlist.back();
            int v_front = reads[i + 1].vlist.front();

            vector<int> bridge;
            if (v_back <= v_front)
                bridge = widest_path(v_back, v_front);

            vector<int> bridge_internal;

            if (bridge.size() >= 2)
                bridge_internal.assign(bridge.begin() + 1, bridge.end() - 1);
			
            merged.insert(merged.end(), bridge_internal.begin(), bridge_internal.end());
            merged.insert(merged.end(), reads[i + 1].vlist.begin(), reads[i + 1].vlist.end());
        }

		set<int> seen; vector<int> deduped;
		for (int v : merged) 
		{
			if (seen.insert(v).second)
				deduped.push_back(v);
		}
		merged = std::move(deduped);

        int32_t frag_start = reads[0].pos;
        int32_t frag_end = reads[n - 1].rpos;

        vector<int> parts;
        int ins = compute_chimeric_insert_size(merged, frag_start, frag_end, bd->regions, &parts);

        fout << "\n";

        if (ins < 0)
            fout << "\t\t(could not compute insert size)\n\n";
        
		else
        {
            fout << "\t\tPath: (";
            for (int i = 0; i < merged.size(); i++)
                fout << merged[i] << (i < merged.size() - 1 ? " → " : "");
            fout << ")\n";

            fout << "\t\tLength: ";

            for (int i = 0; i < parts.size(); i++)
            {
                if (i > 0)
					fout << " + ";
                fout << parts[i];
            }

            if (parts.size() > 1)
				fout << " = " << ins;
			
            fout << "\n\n";
        }
    }

    fout.close();
}

void bridger::write_simplified_outward_paths(int bundle_idx, const std::string& outdir)
{
	std::string grouped_filename = outdir + "/simplify_outward_full_sequence_paths.txt";
    std::ofstream gout(grouped_filename, std::ios_base::app);

    gout << "=================================================\n";
    gout << "Bundle " << bundle_idx << "\n";
    gout << "=================================================\n\n";

    map<int, map<string, map<pair<int32_t,int32_t>, vector<vector<int>>>>>::iterator bundle_it = bundle_outward_bsj_paths.find(bundle_idx);

    if (bundle_it != bundle_outward_bsj_paths.end())
    {
        for (const pair<const string, map<pair<int32_t,int32_t>, vector<vector<int>>>>& qname_entry : bundle_it->second)
        {
            const string& qname = qname_entry.first;
            gout << qname << "\n";
            gout << "---------------------------------\n";

            for (const read_info& read : bundle_outward_reads[bundle_idx])
            {
                if (read.qname != qname) continue;

                string label = read.is_read1 ? "R1" : (read.is_read2 ? "R2" : "?");

                gout << "\t" << label << ": (";
                for (int i = 0; i < read.vlist.size(); i++)
                    gout << read.vlist[i] << (i == read.vlist.size() - 1 ? "" : " → ");

                gout << ")\n";
            }

			gout << "\n";

            for (const pair<const pair<int32_t,int32_t>, vector<vector<int>>>& pos_entry : qname_entry.second)
            {
                int32_t L_pos = pos_entry.first.first;
                int32_t R_pos = pos_entry.first.second;
                const vector<vector<int>>& paths = pos_entry.second;

                gout << "\t(L = " << bd->bb.chrm << ":" << L_pos << ", R = " << bd->bb.chrm << ":" << R_pos << ")\n";
                gout << "\tSupporting Paths (" << paths.size() << "):\n";

                for (int i = 0; i < paths.size(); i++)
                {
                    gout << "\t\tP" << i + 1 << " =  (";
                    for (int j = 0; j < paths[i].size(); j++)
                        gout << paths[i][j] << (j == paths[i].size() - 1 ? "" : " → ");

                    gout << ")\n";
                }
            }
            gout << "\n";
        }
    }

    else
        gout << "No Outward BSJ Candidates\n";

    gout << "\n\n";
    gout.close();


}

map<string, vector<vector<int>>> bridger::build_chimeric_bg(int bundle_idx, map<string, vector<string>>& rep_members)
{
    rep_members.clear();
    map<string, vector<vector<int>>> chim_bg;
	map<int, map<string, vector<vector<int>>>>::iterator it = bundle_chimeric_merged_paths.find(bundle_idx);
	if (it != bundle_chimeric_merged_paths.end())
		chim_bg = it->second;

	map<string, vector<read_info>> chim_segs_map;
	for (const read_info& r : bundle_chimeric_reads[bundle_idx])
	{
		if (chim_bg.count(r.qname))
			chim_segs_map[r.qname].push_back(r);
	}

	for (pair<const string, vector<read_info>>& e : chim_segs_map)
	{
		sort(
			e.second.begin(), e.second.end(), [](const read_info& a, const read_info& b)
			{
				return a.pos < b.pos;
			}
		);
	}

	typedef vector<tuple<int32_t, int32_t, vector<int>>> SegKey;
	map<SegKey, string> key_to_rep;
	map<string, string> chim_rep_map;

	for (const pair<const string, vector<vector<int>>>& qname_entry : chim_bg)
	{
		const string& qname = qname_entry.first;
		map<string, vector<read_info>>::iterator segs_it = chim_segs_map.find(qname);
		if (segs_it == chim_segs_map.end())
		{
			chim_rep_map[qname] = qname;
			continue;
		}
		
		const vector<read_info>& segs = segs_it->second;
		SegKey key;
		for (const read_info& seg : segs)
			key.emplace_back(seg.pos, seg.rpos, seg.vlist);

		map<SegKey, string>::iterator kit = key_to_rep.find(key);
		if (kit != key_to_rep.end())
			chim_rep_map[qname] = kit->second;
		
		else
		{
			chim_rep_map[qname] = qname;
			key_to_rep[key] = qname;
		}
	}

	map<string, vector<vector<int>>> consolidated_chim_bg;
	for (const pair<const string, vector<vector<int>>>& qname_entry : chim_bg)
	{
		const string& rep = chim_rep_map[qname_entry.first];
		consolidated_chim_bg[rep].insert(consolidated_chim_bg[rep].end(), qname_entry.second.begin(), qname_entry.second.end());
		rep_members[rep].push_back(qname_entry.first);
	}

	chim_bg = consolidated_chim_bg;

    return chim_bg;
}

map<string, vector<vector<int>>> bridger::build_outward_bg(int bundle_idx, map<string, vector<string>>& rep_members)
{
    // Consolidate read pairs whose individual reads share identical vlists and pos/rpos so duplicates don't inflate left-side node count.
    rep_members.clear();
    map<string, vector<vector<int>>> outward_bg;
    {
        map<int, map<string, map<pair<int32_t,int32_t>, vector<vector<int>>>>>::iterator bundle_it_bg = bundle_outward_bsj_paths.find(bundle_idx);
        if (bundle_it_bg != bundle_outward_bsj_paths.end())
        {
            struct OutwardKey
            {
                vector<int> left_vlist;
                vector<int> right_vlist;
                int32_t left_pos, left_rpos;
                int32_t right_pos, right_rpos;
            };

            vector<pair<string, OutwardKey>> representatives;
            map<string, string> outward_rep_map; // qname -> representative qname

            for (const pair<const string, map<pair<int32_t,int32_t>, vector<vector<int>>>>& qname_entry : bundle_it_bg->second)
            {
                const string& qname = qname_entry.first;

                const read_info* left_rd  = nullptr;
                const read_info* right_rd = nullptr;
                for (const read_info& r : bundle_outward_reads[bundle_idx])
                {
                    if (r.qname != qname) continue;

                    if (!left_rd || r.pos < left_rd->pos)
						left_rd  = &r;

                    if (!right_rd || r.rpos > right_rd->rpos)
						right_rd = &r;
                }

                if (!left_rd || !right_rd || left_rd == right_rd)
                {
                    outward_rep_map[qname] = qname;
                    continue;
                }

                bool matched = false;
                for (const pair<string, OutwardKey>& rep : representatives)
                {
                    const OutwardKey& key = rep.second;
                    if (left_rd->vlist  == key.left_vlist && right_rd->vlist == key.right_vlist &&
                        abs(left_rd->pos - key.left_pos) <= TOLERANCE_GAP &&
                        abs(left_rd->rpos - key.left_rpos) <= TOLERANCE_GAP &&
                        abs(right_rd->pos - key.right_pos) <= TOLERANCE_GAP &&
                        abs(right_rd->rpos - key.right_rpos) <= TOLERANCE_GAP)
                    {
                        outward_rep_map[qname] = rep.first;
                        matched = true;
                        break;
                    }
                }

                if (!matched)
                {
                    outward_rep_map[qname] = qname;
                    OutwardKey key;
                    key.left_vlist = left_rd->vlist;
                    key.right_vlist = right_rd->vlist;
                    key.left_pos = left_rd->pos;
                    key.left_rpos = left_rd->rpos;
                    key.right_pos = right_rd->pos;
                    key.right_rpos = right_rd->rpos;
                    representatives.push_back({qname, key});
                }
            }

            for (const pair<const string, map<pair<int32_t,int32_t>, vector<vector<int>>>>& qname_entry : bundle_it_bg->second)
            {
                const string& rep = outward_rep_map[qname_entry.first];
                for (const pair<const pair<int32_t,int32_t>, vector<vector<int>>>& pos_entry : qname_entry.second)
				{
                    for (const vector<int>& path : pos_entry.second)
                        outward_bg[rep].push_back(path);
                }
                rep_members[rep].push_back(qname_entry.first);
            }
        }
    }

    return outward_bg;
}

static int32_t bridge_interior_len(bundle_bridge *bd, const vector<int>& path, int v1, int v2)
{
	int i1 = -1;
	int i2 = -1;
	for (int j = 0; j < path.size(); j++)
	{
		if (path[j] == v1 && i1 < 0)
			i1 = j;

		if (path[j] == v2 && i2 < 0)
			i2 = j;
	}

	if (i1 < 0 || i2 <= i1 + 1) return 0;

	vector<int> interior(path.begin() + i1 + 1, path.begin() + i2);
	return bd->compute_aligned_length(0, 0, interior);
}

static map<int, int> assign_ranks(const vector<int>& candidates, const vector<double>& values, bool higher_is_better, double invalid_val)
{
	vector<pair<double, int>> valid_items, invalid_items;
	for (int idx : candidates)
	{
		double v = values[idx];
		(v == invalid_val ? invalid_items : valid_items).push_back({v, idx});
	}

	if (higher_is_better)
	{
		sort(
			valid_items.begin(), valid_items.end(), [](const pair<double,int>& a, const pair<double,int>& b)
			{
				return a.first > b.first;
			}
		);
	}

	else
	{
		sort(
			valid_items.begin(), valid_items.end(), [](const pair<double,int>& a, const pair<double,int>& b)
			{
				return a.first < b.first;
			}
		);
	}

	map<int, int> ranks;
	int r = 1;
	for (int i = 0; i < valid_items.size(); i++)
	{
		if (i > 0 && valid_items[i].first != valid_items[i - 1].first)
			r = i + 1;
		ranks[valid_items[i].second] = r;
	}

	int worst = valid_items.size() + 1;

	for (const pair<double, int>& item : invalid_items)
		ranks[item.second] = worst;

	return ranks;
}

// Rank-sum tie-break: unchanged Borda count over 5 features, restricted to whichever
// candidates are still tied after the learned selection_score cut below.
static pair<int, double> pick_best_ranked_borda(const vector<int>& candidates, const vector<set<string>>& path_chim_reads, const vector<set<string>>& path_outward_reads, const set<string>& explained, const vector<int>& min_edge_weight, const vector<int32_t>& insert_sizes, const vector<double>& min_vertex_weight, int32_t length_median)
{
	int n = min_edge_weight.size();
	vector<double> chim_vals(n, 0.0), out_vals(n, 0.0), ew_vals(n, -1.0), dev_vals(n, -1.0), vw_vals(n, -1.0);

	for (int idx : candidates)
	{
		int cc = 0;
		for (const string& r : path_chim_reads[idx])
		{
			if (!explained.count(r))
				cc++;
		}

		chim_vals[idx] = (double)cc;

		int oc = 0;
		for (const string& r : path_outward_reads[idx])
		{
			if (!explained.count(r))
				oc++;
		}

		out_vals[idx] = (double)oc;

		ew_vals[idx] = (double)min_edge_weight[idx];

		dev_vals[idx] = (insert_sizes[idx] >= 0 && length_median > 0)
		    ? (double)abs(insert_sizes[idx] - (int32_t)length_median)
		    : -1.0;

		vw_vals[idx] = min_vertex_weight[idx];
	}

	map<int, int> r_chim = assign_ranks(candidates, chim_vals, true, -1.0);
	map<int, int> r_out = assign_ranks(candidates, out_vals, true, -1.0);
	map<int, int> r_ew = assign_ranks(candidates, ew_vals, true, -1.0);
	map<int, int> r_dev = assign_ranks(candidates, dev_vals, false, -1.0);
	map<int, int> r_vw = assign_ranks(candidates, vw_vals, true, -1.0);

	int best_idx = -1;
	int best_rank_sum = INT_MAX;

	for (int idx : candidates)
	{
		int rank_sum = r_chim[idx] + r_out[idx] + r_ew[idx] + r_dev[idx] + r_vw[idx];
		if (best_idx < 0 || rank_sum < best_rank_sum)
		{
			best_rank_sum = rank_sum;
			best_idx = idx;
		}
	}

	return {best_idx, (double)best_rank_sum};
}

static pair<int, double> pick_best_ranked(const vector<int>& candidates, const vector<set<string>>& path_chim_reads, const vector<set<string>>& path_outward_reads, const set<string>& explained, const vector<int>& min_edge_weight, const vector<int32_t>& insert_sizes, const vector<double>& min_vertex_weight, int32_t length_median, const vector<int>& n_vertices_feat, const vector<int>& exon_count_feat, const vector<int>& weak_edge_count, const vector<int>& circ_start, const vector<int>& circ_end, const vector<int>& cross_tissue_repro_count, const vector<int>& cross_tissue_repro_full_count, const vector<double>& circ_ratio, const vector<int>& chimeric_support_weighted, const vector<int>& cross_tissue_chim_repro_count, const vector<int>& fake_supple_count, const vector<int>& supple_len)
{
	if (candidates.empty()) return {-1, 0.0};

	vector<int> effective_candidates;
	{
		auto score_of = [&](int i) {
			int cc = 0;
			for (const string& r : path_chim_reads[i])
				if (!explained.count(r)) cc++;
			int oc = 0;
			for (const string& r : path_outward_reads[i])
				if (!explained.count(r)) oc++;
			bool ew_valid = min_edge_weight[i] >= 0;
			bool vw_valid = min_vertex_weight[i] >= 0.0;
			int32_t dev = (insert_sizes[i] >= 0 && length_median > 0)
			    ? abs(insert_sizes[i] - (int32_t)length_median) : 0;
			bool dev_valid = (insert_sizes[i] >= 0 && length_median > 0);
			return selection_score(
				n_vertices_feat[i], exon_count_feat[i], cc, oc,
				ew_valid ? min_edge_weight[i] : 0, ew_valid,
				vw_valid ? min_vertex_weight[i] : 0.0, vw_valid,
				dev, dev_valid, weak_edge_count[i], cross_tissue_repro_count[i], cross_tissue_repro_full_count[i],
				circ_ratio[i], chimeric_support_weighted[i], cross_tissue_chim_repro_count[i], fake_supple_count[i], supple_len[i]);
		};

		map<pair<int,int>, int> best_by_boundary;
		for (int idx : candidates)
		{
			pair<int,int> key = {circ_start[idx], circ_end[idx]};
			map<pair<int,int>, int>::iterator it = best_by_boundary.find(key);
			if (it == best_by_boundary.end())
			{
				best_by_boundary[key] = idx;
				continue;
			}

			int incumbent = it->second;
			if (min_edge_weight[idx] > min_edge_weight[incumbent])
			{
				best_by_boundary[key] = idx;
			}
			else if (min_edge_weight[idx] == min_edge_weight[incumbent])
			{
				double score_idx = score_of(idx);
				double score_incumbent = score_of(incumbent);
				if (score_idx > score_incumbent)
				{
					best_by_boundary[key] = idx;
				}
				else if (score_idx == score_incumbent)
				{
					pair<int, double> tie_break = pick_best_ranked_borda(
						{incumbent, idx}, path_chim_reads, path_outward_reads, explained,
						min_edge_weight, insert_sizes, min_vertex_weight, length_median);
					best_by_boundary[key] = tie_break.first;
				}
			}
		}
		for (const pair<const pair<int,int>, int>& kv : best_by_boundary)
			effective_candidates.push_back(kv.second);
	}
	const vector<int>& scored_candidates = effective_candidates;

	double best_score = -1.0;
	for (int idx : scored_candidates)
	{
		int cc = 0;
		for (const string& r : path_chim_reads[idx])
			if (!explained.count(r)) cc++;

		int oc = 0;
		for (const string& r : path_outward_reads[idx])
			if (!explained.count(r)) oc++;

		bool ew_valid = min_edge_weight[idx] >= 0;
		bool vw_valid = min_vertex_weight[idx] >= 0.0;
		int32_t dev = (insert_sizes[idx] >= 0 && length_median > 0)
		    ? abs(insert_sizes[idx] - (int32_t)length_median) : 0;
		bool dev_valid = (insert_sizes[idx] >= 0 && length_median > 0);

		double score = selection_score(
			n_vertices_feat[idx], exon_count_feat[idx], cc, oc,
			ew_valid ? min_edge_weight[idx] : 0, ew_valid,
			vw_valid ? min_vertex_weight[idx] : 0.0, vw_valid,
			dev, dev_valid, weak_edge_count[idx], cross_tissue_repro_count[idx], cross_tissue_repro_full_count[idx],
			circ_ratio[idx], chimeric_support_weighted[idx], cross_tissue_chim_repro_count[idx], fake_supple_count[idx], supple_len[idx]);

		if (score > best_score)
			best_score = score;
	}

	vector<int> tied;
	for (int idx : scored_candidates)
	{
		int cc = 0;
		for (const string& r : path_chim_reads[idx])
			if (!explained.count(r)) cc++;

		int oc = 0;
		for (const string& r : path_outward_reads[idx])
			if (!explained.count(r)) oc++;

		bool ew_valid = min_edge_weight[idx] >= 0;
		bool vw_valid = min_vertex_weight[idx] >= 0.0;
		int32_t dev = (insert_sizes[idx] >= 0 && length_median > 0)
		    ? abs(insert_sizes[idx] - (int32_t)length_median) : 0;
		bool dev_valid = (insert_sizes[idx] >= 0 && length_median > 0);

		double score = selection_score(
			n_vertices_feat[idx], exon_count_feat[idx], cc, oc,
			ew_valid ? min_edge_weight[idx] : 0, ew_valid,
			vw_valid ? min_vertex_weight[idx] : 0.0, vw_valid,
			dev, dev_valid, weak_edge_count[idx], cross_tissue_repro_count[idx], cross_tissue_repro_full_count[idx],
			circ_ratio[idx], chimeric_support_weighted[idx], cross_tissue_chim_repro_count[idx], fake_supple_count[idx], supple_len[idx]);

		if (score == best_score)
			tied.push_back(idx);
	}

	if (tied.size() == 1)
		return {tied[0], best_score};

	pair<int, double> tie_break = pick_best_ranked_borda(tied, path_chim_reads, path_outward_reads, explained, min_edge_weight, insert_sizes, min_vertex_weight, length_median);
	return {tie_break.first, best_score};
}

static bool is_shifted_duplicate_bsj(bundle_bridge* bd, int s1, int e1, int s2, int e2)
{
	int shift = s2 - s1;
	if (shift == 0 || (e2 - e1) != shift || abs(shift) > read_length) return false;

	string seq_a, seq_b;
	if (shift > 0)
	{
		seq_a = bd->get_fasta_seq(s1, s1 + shift - 1);
		seq_b = bd->get_fasta_seq(e1, e1 + shift - 1);
	}

	else
	{
		int k = -shift;
		seq_a = bd->get_fasta_seq(e1 - k, e1 - 1);
		seq_b = bd->get_fasta_seq(s1 - k, s1 - 1);
	}

	if (seq_a.empty() || seq_a.size() != seq_b.size()) return false;

	for (char& c : seq_a)
		c = toupper(c);

	for (char& c : seq_b)
		c = toupper(c);

	return seq_a == seq_b;
}

static int uf_find(map<int,int>& uf_parent, int x)
{
	while (uf_parent[x] != x)
	{
		uf_parent[x] = uf_parent[uf_parent[x]];
		x = uf_parent[x];
	}

	return x;
}

static bool is_subset(const set<string>& sub, const set<string>& sup)
{
	for (const string& x : sub)
		if (!sup.count(x)) return false;

		return true;
}

static vector<pair<int32_t,int32_t>> merge_path_to_exons(bundle_bridge* bd, const vector<int>& vpath)
{
	vector<pair<int32_t,int32_t>> exons;
	for (int v : vpath)
	{
		int32_t lp = bd->regions[v].lpos;
		int32_t rp = bd->regions[v].rpos;

		if (!exons.empty() && lp <= exons.back().second)
			exons.back().second = max(exons.back().second, rp);
		else
			exons.push_back({lp, rp});
	}
	return exons;
}

struct SiblingReproCounts
{
	map<tuple<string, int32_t, int32_t>, int> boundary;
	map<tuple<string, int32_t, int32_t, string, string>, int> full_structure;
	map<tuple<string, int32_t, int32_t>, int> chim_boundary;
};

static const SiblingReproCounts& get_sibling_tissue_repro_counts(const string& outdir)
{
	static bool loaded = false;
	static SiblingReproCounts counts;

	if (loaded) return counts;
	loaded = true;

	return counts;

	char own_resolved[4096];
	string own_pool_path = outdir + "/circrna_candidate_pool.tsv";
	bool have_own = (realpath(own_pool_path.c_str(), own_resolved) != NULL);
	string own_resolved_str = have_own ? string(own_resolved) : string();

	DIR* dir = opendir("..");
	if (dir == NULL) return counts;

	struct dirent* entry;
	while ((entry = readdir(dir)) != NULL)
	{
		string name = entry->d_name;
		if (name == "." || name == "..") continue;

		string sibling_pool_path = string("../") + name + "/results/circrna_candidate_pool.tsv";

		char sibling_resolved[4096];
		if (realpath(sibling_pool_path.c_str(), sibling_resolved) == NULL) continue;
		if (have_own && own_resolved_str == string(sibling_resolved)) continue;

		std::ifstream fin(sibling_pool_path.c_str());
		if (!fin.is_open()) continue;

		set<tuple<string, int32_t, int32_t>> this_tissue_boundaries;
		set<tuple<string, int32_t, int32_t, string, string>> this_tissue_full;
		set<tuple<string, int32_t, int32_t>> this_tissue_chim_boundaries;
		string line;
		while (std::getline(fin, line))
		{
			if (line.empty() || line[0] == '#') continue;

			size_t p1 = line.find('\t');
			if (p1 == string::npos) continue;
			size_t p2 = line.find('\t', p1 + 1);
			if (p2 == string::npos) continue;
			size_t p3 = line.find('\t', p2 + 1);
			if (p3 == string::npos) continue;
			size_t p4 = line.find('\t', p3 + 1);
			if (p4 == string::npos) continue;
			size_t p5 = line.find('\t', p4 + 1);
			if (p5 == string::npos) continue;
			size_t p6 = line.find('\t', p5 + 1);
			if (p6 == string::npos) continue;
			size_t p7 = line.find('\t', p6 + 1);
			if (p7 == string::npos) continue;
			size_t p8 = line.find('\t', p7 + 1);
			if (p8 == string::npos) continue;
			size_t p9 = line.find('\t', p8 + 1);
			if (p9 == string::npos) continue;
			size_t p10 = line.find('\t', p9 + 1);
			if (p10 == string::npos) continue;
			size_t p11 = line.find('\t', p10 + 1);
			if (p11 == string::npos) continue;

			string chrm = line.substr(p2 + 1, p3 - p2 - 1);
			string cs_str = line.substr(p4 + 1, p5 - p4 - 1);
			string ce_str = line.substr(p5 + 1, p6 - p5 - 1);
			// field6 (p6..p7) is exon_count, skipped; field7/8 are exon_starts/exon_ends.
			string exon_starts_str = line.substr(p7 + 1, p8 - p7 - 1);
			string exon_ends_str = line.substr(p8 + 1, p9 - p8 - 1);
			// field9 (p9..p10) is n_vertices, skipped; field10 is chimeric_support_reads.
			string chim_support_str = line.substr(p10 + 1, p11 - p10 - 1);

			char* endptr1 = NULL;
			char* endptr2 = NULL;
			long cs = strtol(cs_str.c_str(), &endptr1, 10);
			long ce = strtol(ce_str.c_str(), &endptr2, 10);
			if (endptr1 == cs_str.c_str() || endptr2 == ce_str.c_str()) continue;

			this_tissue_boundaries.insert(make_tuple(chrm, (int32_t)cs, (int32_t)ce));
			this_tissue_full.insert(make_tuple(chrm, (int32_t)cs, (int32_t)ce, exon_starts_str, exon_ends_str));

			char* endptr3 = NULL;
			long chim_support = strtol(chim_support_str.c_str(), &endptr3, 10);
			if (endptr3 != chim_support_str.c_str() && chim_support > 0)
				this_tissue_chim_boundaries.insert(make_tuple(chrm, (int32_t)cs, (int32_t)ce));
		}

		for (const tuple<string, int32_t, int32_t>& key : this_tissue_boundaries)
			counts.boundary[key]++;
		for (const tuple<string, int32_t, int32_t, string, string>& key : this_tissue_full)
			counts.full_structure[key]++;
		for (const tuple<string, int32_t, int32_t>& key : this_tissue_chim_boundaries)
			counts.chim_boundary[key]++;
	}
	closedir(dir);

	return counts;
}

static bool is_canonical_bsj_boundary(bundle_bridge* bd, int32_t circ_start, int32_t circ_end, char strand, char* resolved_strand_out = nullptr, bool treat_missing_as_canonical = false)
{
	string acceptor_seq = bd->get_fasta_seq(circ_start - 2, circ_start - 1);
	string donor_seq = bd->get_fasta_seq(circ_end, circ_end + 1);
	if (acceptor_seq.size() != 2 || donor_seq.size() != 2) return treat_missing_as_canonical;

	for (char &c : acceptor_seq) c = toupper(c);
	for (char &c : donor_seq) c = toupper(c);

	bool matches_plus = (acceptor_seq == "AG" && donor_seq == "GT");
	bool matches_minus = (acceptor_seq == "AC" && donor_seq == "CT");

	if (strand == '+') return matches_plus;
	if (strand == '-') return matches_minus;

	if (resolved_strand_out)
	{
		if (matches_plus) *resolved_strand_out = '+';
		else if (matches_minus) *resolved_strand_out = '-';
	}
	return matches_plus || matches_minus;
}

static int count_noncanonical_internal_junctions(bundle_bridge* bd, const vector<pair<int32_t, int32_t>>& exons, char strand)
{
	if (exons.size() < 2 || (strand != '+' && strand != '-')) return 0;

	int count = 0;
	for (size_t k = 0; k + 1 < exons.size(); k++)
	{
		int32_t up_end = exons[k].second;
		int32_t down_start = exons[k + 1].first;
		if (!is_canonical_bsj_boundary(bd, down_start, up_end, strand, nullptr, true))
			count++;
	}
	return count;
}

static int inverted_repeat_kmer_matches(bundle_bridge* bd, int32_t circ_start, int32_t circ_end)
{
	const int32_t FLANK = max_multi_exon_length;
	const int K = 16;

	string up = bd->get_fasta_seq(circ_start - 1 - FLANK, circ_start - 2);
	string down = bd->get_fasta_seq(circ_end, circ_end + FLANK - 1);
	if ((int)up.size() < K || (int)down.size() < K) return 0;

	for (char &c : up) c = toupper(c);
	for (char &c : down) c = toupper(c);

	string down_rc(down.rbegin(), down.rend());
	for (char &c : down_rc)
	{
		if (c == 'A') c = 'T';
		else if (c == 'T') c = 'A';
		else if (c == 'C') c = 'G';
		else if (c == 'G') c = 'C';
	}

	set<string> up_kmers;
	for (int i = 0; i + K <= (int)up.size(); i++)
		up_kmers.insert(up.substr(i, K));

	set<string> matched;
	for (int i = 0; i + K <= (int)down_rc.size(); i++)
	{
		string km = down_rc.substr(i, K);
		if (up_kmers.count(km)) matched.insert(km);
	}

	return (int)matched.size();
}

static bool find_softclip_confirmed_outward_bsj(bundle_bridge* bd, const hit* left_h, const hit* right_h,
	const map<int32_t, int>& lpos_to_v, const map<int32_t, int>& rpos_to_v,
	int32_t& out_bsj_start, int32_t& out_bsj_end)
{
	// Case 1: left mate's leading soft-clip -- its own alignment start (left_h->pos) is the BSJ
	// start; the clipped bases are hypothesized to belong at a junction to the right of the pair.
	if (left_h && !left_h->cigar_vector.empty() && left_h->cigar_vector.front().first == 'S'
		&& !left_h->soft_left_clip_seqs.empty())
	{
		int32_t len = (int32_t)left_h->soft_left_clip_seqs[0].size();
		if (len >= min_soft_clip_len)
		{
			for (const junction& junc : bd->junctions)
			{
				if (junc.junc_type != 1) continue;
				if (right_h && junc.lpos < right_h->rpos) continue;
				if (rpos_to_v.find(junc.lpos) == rpos_to_v.end()) continue;

				int32_t circ_start = left_h->pos;
				int32_t circ_end = junc.lpos;
				if (!is_canonical_bsj_boundary(bd, circ_start, circ_end, bd->bb.strand)) continue;

				string region_seq = bd->get_fasta_seq(circ_end - len, circ_end - 1);
				if ((int32_t)region_seq.size() != len) continue;

				int mapped_len = max(1, (int)left_h->qlen - len);
				int allowed = max(1, (int)((left_h->nm * (int64_t)len + mapped_len - 1) / mapped_len));

				if (bd->get_edit_distance(left_h->soft_left_clip_seqs[0], region_seq) <= allowed)
				{
					out_bsj_start = circ_start;
					out_bsj_end = circ_end;
					return true;
				}
			}
		}
	}

	// Case 2: right mate's trailing soft-clip -- symmetric; right_h->rpos is the BSJ end, clipped
	// bases hypothesized to belong at a junction to the left of the pair.
	if (right_h && !right_h->cigar_vector.empty() && right_h->cigar_vector.back().first == 'S'
		&& !right_h->soft_right_clip_seqs.empty())
	{
		int32_t len = (int32_t)right_h->soft_right_clip_seqs[0].size();
		if (len >= min_soft_clip_len)
		{
			for (const junction& junc : bd->junctions)
			{
				if (junc.junc_type != 1) continue;
				if (left_h && junc.rpos > left_h->pos) continue;
				if (lpos_to_v.find(junc.rpos) == lpos_to_v.end()) continue;

				int32_t circ_start = junc.rpos;
				int32_t circ_end = right_h->rpos;
				if (!is_canonical_bsj_boundary(bd, circ_start, circ_end, bd->bb.strand)) continue;

				string region_seq = bd->get_fasta_seq(circ_start - 1, circ_start + len - 2);
				if ((int32_t)region_seq.size() != len) continue;

				int mapped_len = max(1, (int)right_h->qlen - len);
				int allowed = max(1, (int)((right_h->nm * (int64_t)len + mapped_len - 1) / mapped_len));

				if (bd->get_edit_distance(right_h->soft_right_clip_seqs[0], region_seq) <= allowed)
				{
					out_bsj_start = circ_start;
					out_bsj_end = circ_end;
					return true;
				}
			}
		}
	}

	return false;
}

static bool dominates(const CandidatePath& a, const CandidatePath& b)
{
	bool ge = a.min_edge_weight >= b.min_edge_weight && a.min_vertex_weight >= b.min_vertex_weight && a.length <= b.length;
	bool strictly_better = a.min_edge_weight > b.min_edge_weight || a.min_vertex_weight > b.min_vertex_weight || a.length < b.length;
	return ge && strictly_better;
}

// Inserts cand into frontier unless an existing entry already dominates it; on insertion, drops
// any existing entries cand itself dominates.
static void insert_if_nondominated(vector<CandidatePath>& frontier, const CandidatePath& cand)
{
	for (const CandidatePath& existing : frontier)
		if (dominates(existing, cand)) return;

	vector<CandidatePath> kept;
	for (CandidatePath& existing : frontier)
		if (!dominates(cand, existing))
			kept.push_back(existing);

	kept.push_back(cand);
	frontier = kept;
}

vector<CandidatePath> bridger::find_pareto_paths(int start, int end, int32_t len_cap, bool require_canonical_junctions)
{
	map<int, vector<CandidatePath>> frontier;

	double start_vw = bd->regions[start].ave;
	CandidatePath seed;
	seed.path = {start};
	seed.min_edge_weight = INT_MAX;
	seed.min_vertex_weight = (start_vw > 0.0) ? start_vw : -1.0;
	seed.length = bd->regions[start].rpos - bd->regions[start].lpos;
	frontier[start].push_back(seed);

	for (int v = start; v <= end; v++)
	{
		map<int, vector<CandidatePath>>::iterator fit = frontier.find(v);
		if (fit == frontier.end()) continue;

		map<int, vector<pair<int, int>>>::iterator adj_it = splice_graph_adj.find(v);
		if (adj_it == splice_graph_adj.end()) continue;

		for (const CandidatePath& state : fit->second)
		{
			for (const pair<int, int>& np : adj_it->second)
			{
				int w = np.first;
				int edge_w = np.second;
				if (w > end) continue;

				CandidatePath ext = state;
				ext.path.push_back(w);
				ext.min_edge_weight = min(state.min_edge_weight, edge_w);

				double w_vw = bd->regions[w].ave;
				if (w_vw > 0.0)
					ext.min_vertex_weight = (state.min_vertex_weight < 0.0) ? w_vw : min(state.min_vertex_weight, w_vw);

				ext.length = state.length + (bd->regions[w].rpos - bd->regions[w].lpos);

				if (len_cap >= 0 && ext.length > len_cap) continue;

				if (require_canonical_junctions && bd->regions[w].lpos > bd->regions[v].rpos)
				{
					if (!is_canonical_bsj_boundary(bd, bd->regions[w].lpos, bd->regions[v].rpos, bd->bb.strand, nullptr, true))
						continue;
				}

				insert_if_nondominated(frontier[w], ext);
			}
		}
	}

	map<int, vector<CandidatePath>>::iterator eit = frontier.find(end);
	if (eit == frontier.end()) return {};
	return eit->second;
}

void bridger::build_outward_read_paths(int bundle_idx, const std::string& outdir)
{
	std::string filename = outdir + "/outward_possible_configs.txt";
	std::ofstream fout(filename, std::ios_base::app);

	fout << "=================================================\n";
	fout << "Bundle " << bundle_idx << "\n";
	fout << "=================================================\n\n";

	fout << "Splice Graph\n";
	fout << "---------------------------------\n";
	write_splice_graph(fout);
	fout << "\n";

	fout << "Outward Reads (" << bundle_outward_reads[bundle_idx].size() << ")\n";
	fout << "---------------------------------\n";

	const int32_t OUTWARD_SEGMENT_LEN_CAP = 2 * length_median;

	map<string, vector<read_info>> grouped_outward;
	for (const read_info& read : bundle_outward_reads[bundle_idx])
		grouped_outward[read.qname].push_back(read);

	map<int32_t, int> rpos_to_v, lpos_to_v;
	for (int i = 0; i < bd->regions.size(); i++)
	{
		rpos_to_v[bd->regions[i].rpos] = i;
		lpos_to_v[bd->regions[i].lpos] = i;
	}

	for (const std::pair<const string, vector<read_info>>& pr : grouped_outward)
	{
		fout << pr.first << "\n";

		for (const read_info& read : pr.second)
		{
			string label = read.is_read1 ? "(R1)" : (read.is_read2 ? "(R2)" : " -");
			fout << "\t" << bd->bb.chrm << ":" << read.pos << "-" << read.rpos << "  " << label << "\n";
		}

		if (pr.second.size() != 2) { fout << "\n"; continue; }

		const read_info* left_read = (pr.second[0].pos <= pr.second[1].pos) ? &pr.second[0] : &pr.second[1];
		const read_info* right_read = (pr.second[0].rpos >= pr.second[1].rpos) ? &pr.second[0] : &pr.second[1];

		if (left_read->vlist.empty() || right_read->vlist.empty()) { fout << "\n"; continue; }

		int v_left_first = left_read->vlist.front();
		int v_left_last = left_read->vlist.back();
		int v_right_first = right_read->vlist.front();
		int v_right_last = right_read->vlist.back();

		int32_t leftmost_pos = left_read->pos;
		int32_t rightmost_rpos = right_read->rpos;

		// Prospective soft-clip short-circuit.
		int32_t sc_start, sc_end;
		if (find_softclip_confirmed_outward_bsj(bd, left_read->src, right_read->src, lpos_to_v, rpos_to_v, sc_start, sc_end))
		{
			vector<CandidatePath> mid_candidates = find_pareto_paths(v_left_last, v_right_first, -1, true);
			for (const CandidatePath& mid : mid_candidates)
			{
				vector<int> mid_internal;
				if (mid.path.size() >= 2)
					mid_internal.assign(mid.path.begin() + 1, mid.path.end() - 1);

				vector<int> merged_path;
				map<int32_t, int>::iterator lit = lpos_to_v.find(sc_start);
				if (lit != lpos_to_v.end() && lit->second != v_left_first)
					merged_path.push_back(lit->second);

				merged_path.insert(merged_path.end(), left_read->vlist.begin(), left_read->vlist.end());
				merged_path.insert(merged_path.end(), mid_internal.begin(), mid_internal.end());
				merged_path.insert(merged_path.end(), right_read->vlist.begin(), right_read->vlist.end());

				map<int32_t, int>::iterator rit = rpos_to_v.find(sc_end);
				if (rit != rpos_to_v.end() && rit->second != v_right_last)
					merged_path.push_back(rit->second);

				merged_path.erase(unique(merged_path.begin(), merged_path.end()), merged_path.end());

				bundle_outward_bsj_paths[bundle_idx][pr.first][{sc_start, sc_end}].push_back(merged_path);
			}

			if (!mid_candidates.empty())
			{
				fout << "\t[soft-clip confirmed BSJ: L = " << bd->bb.chrm << ":" << sc_start
					<< ", R = " << bd->bb.chrm << ":" << sc_end << "]\n\n";
			}

			continue;
		}

		map<int, vector<CandidatePath>> L_candidates, R_candidates;
		set<int32_t> seen_L_pos, seen_R_pos;

		for (const junction& junc : bd->junctions)
		{
			// skip BSJ junctions
			if (junc.junc_type != 1) continue;

			// L candidate: junc.rpos as L; v_L = vertex whose lpos == junction.rpos; path v_L -> v_left_first
			if (junc.rpos <= leftmost_pos && seen_L_pos.find(junc.rpos) == seen_L_pos.end())
			{
				map<int, int>::iterator it_L = lpos_to_v.find(junc.rpos);
				if (it_L != lpos_to_v.end())
				{
					seen_L_pos.insert(junc.rpos);
					int v_L = it_L->second;
					vector<CandidatePath> raw = find_pareto_paths(v_L, v_left_first, OUTWARD_SEGMENT_LEN_CAP, true);
					if (!raw.empty())
						L_candidates[v_L] = raw;
				}
			}

			// R candidate: junc.lpos as R; v_R = vertex whose rpos == junction.lpos; path v_right_last -> v_R
			if (junc.lpos >= rightmost_rpos && seen_R_pos.find(junc.lpos) == seen_R_pos.end())
			{
				map<int, int>::iterator it_R = rpos_to_v.find(junc.lpos);
				if (it_R != rpos_to_v.end())
				{
					seen_R_pos.insert(junc.lpos);
					int v_R = it_R->second;
					vector<CandidatePath> raw = find_pareto_paths(v_right_last, v_R, OUTWARD_SEGMENT_LEN_CAP, true);
					if (!raw.empty())
						R_candidates[v_R] = raw;
				}
			}
		}

		int n_L_total = 0, n_R_total = 0;
		for (const pair<const int, vector<CandidatePath>>& lp : L_candidates) n_L_total += lp.second.size();
		for (const pair<const int, vector<CandidatePath>>& rp : R_candidates) n_R_total += rp.second.size();

		fout << "\n\tNumber of Left Splice Candidates: " << n_L_total << " (" << L_candidates.size() << " vertices)\n";
		fout << "\tNumber of Right Splice Candidates: " << n_R_total << " (" << R_candidates.size() << " vertices)\n";

		if (L_candidates.empty() || R_candidates.empty())
		{
			fout << "\n";
			continue;
		}

		// Mid segment: connects the two mates' own alignments directly; not subject to the 2M
		// cap (unlike L and R) -- full Pareto set retained, admissible only through edges whose
		// splice junctions (real intron-skips) carry the canonical donor/acceptor motif.
		vector<CandidatePath> mid_candidates = find_pareto_paths(v_left_last, v_right_first, -1, true);
		if (mid_candidates.empty())
		{
			fout << "\tNo Middle Path Found\n\n";
			continue;
		}

		fout << "\tNumber of Middle Path Candidates: " << mid_candidates.size() << "\n";

		int combos_accepted = 0;
		int combos_total = 0;
		for (const pair<const int, vector<CandidatePath>>& lp : L_candidates)
		{
			for (const CandidatePath& l : lp.second)
			{
				for (const pair<const int, vector<CandidatePath>>& rp : R_candidates)
				{
					for (const CandidatePath& r : rp.second)
					{
						for (const CandidatePath& mid : mid_candidates)
						{
							combos_total++;

							vector<int> mid_internal;
							if (mid.path.size() >= 2)
								mid_internal.assign(mid.path.begin() + 1, mid.path.end() - 1);

							vector<int> merged_path = l.path;
							for (int k = 1; k < left_read->vlist.size(); k++)
								merged_path.push_back(left_read->vlist[k]);

							merged_path.insert(merged_path.end(), mid_internal.begin(), mid_internal.end());
							for (int k = 0; k < (int)right_read->vlist.size() - 1; k++)
								merged_path.push_back(right_read->vlist[k]);

							merged_path.insert(merged_path.end(), r.path.begin(), r.path.end());
							merged_path.erase(unique(merged_path.begin(), merged_path.end()), merged_path.end());

							int32_t full_seq_length_o = bd->compute_aligned_length(0, 0, merged_path);
							int32_t mid_internal_len = bd->compute_aligned_length(0, 0, mid_internal);
							int32_t outward_insert_size = full_seq_length_o - mid_internal_len;

							if (outward_insert_size > 0 && outward_insert_size <= OUTWARD_SEGMENT_LEN_CAP)
							{
								int32_t l_junc_pos = bd->regions[lp.first].lpos;
								int32_t r_junc_pos = bd->regions[rp.first].rpos;
								bundle_outward_bsj_paths[bundle_idx][pr.first][{l_junc_pos, r_junc_pos}].push_back(merged_path);
								combos_accepted++;
							}
						}
					}
				}
			}
		}

		fout << "\tCombinations Accepted: " << combos_accepted << " / " << combos_total << "\n\n";
	}

	fout.close();
}

void bridger::write_bipartite_graph_file(int bundle_idx, const std::string& outdir, const map<string, vector<vector<int>>>& chim_bg, const map<string, vector<vector<int>>>& outward_bg, const map<string, vector<vector<int>>>& bg, const map<string, vector<string>>& chim_rep_members, const map<string, vector<string>>& outward_rep_members)
{
    std::string bg_filename = outdir + "/bipartite_graph.txt";
    std::ofstream bgout(bg_filename, std::ios_base::app);

    bgout << "=================================================\n";
    bgout << "Bundle " << bundle_idx << "\n";
    bgout << "=================================================\n\n";

    if (bg.empty())
        bgout << "No Bipartite Graph\n\n\n";

	else
	{
		map<string, string> qname_to_chim_rep;
		for (const pair<const string, vector<string>>& rep_entry : chim_rep_members)
			for (const string& member : rep_entry.second)
				qname_to_chim_rep[member] = rep_entry.first;

		map<string, string> qname_to_outward_rep;
		for (const pair<const string, vector<string>>& rep_entry : outward_rep_members)
			for (const string& member : rep_entry.second)
				qname_to_outward_rep[member] = rep_entry.first;

		double bundle_expr_ceiling = 0.0;
		for (const region& r : bd->regions)
			if (r.ave > bundle_expr_ceiling)
				bundle_expr_ceiling = r.ave;

		map<vector<int>, int> path_to_idx;
		vector<vector<int>> idx_to_path;
		for (const pair<const string, vector<vector<int>>>& qname_entry : bg)
		{
			for (const vector<int>& path : qname_entry.second)
			{
				if (path_to_idx.find(path) == path_to_idx.end())
				{
					path_to_idx[path] = idx_to_path.size();
					idx_to_path.push_back(path);
				}
			}
		}

		vector<int> pre_gate_cc_size(idx_to_path.size(), 1);
		{
			map<int,int> pg_uf_parent;
			for (int i = 0; i < (int)idx_to_path.size(); i++)
				pg_uf_parent[i] = i;

			for (const pair<const string, vector<vector<int>>>& qname_entry : bg)
			{
				set<int> active_paths;
				for (const vector<int>& path : qname_entry.second)
				{
					map<vector<int>, int>::iterator it = path_to_idx.find(path);
					if (it != path_to_idx.end())
						active_paths.insert(it->second);
				}
				if (active_paths.size() <= 1) continue;

				int first = uf_find(pg_uf_parent, *active_paths.begin());
				for (int idx : active_paths)
				{
					int r = uf_find(pg_uf_parent, idx);
					if (r != first)
						pg_uf_parent[r] = first;
				}
			}

			map<int, int> pg_comp_count;
			for (int i = 0; i < (int)idx_to_path.size(); i++)
				pg_comp_count[uf_find(pg_uf_parent, i)]++;

			for (int i = 0; i < (int)idx_to_path.size(); i++)
				pre_gate_cc_size[i] = pg_comp_count[uf_find(pg_uf_parent, i)];
		}

		vector<int> outward_support(idx_to_path.size(), 0);
		vector<int> distinct_chimeric_support(idx_to_path.size(), 0);
		vector<int> distinct_outward_support(idx_to_path.size(), 0);
		for (const pair<const string, vector<vector<int>>>& qname_entry : outward_bg)
		{
			map<string, vector<string>>::const_iterator mit = outward_rep_members.find(qname_entry.first);
			int n_members = (mit != outward_rep_members.end()) ? (int)mit->second.size() : 1;

			set<int> seen_for_qname;
			for (const vector<int>& path : qname_entry.second)
			{
				map<vector<int>, int>::iterator it = path_to_idx.find(path);
				if (it == path_to_idx.end()) continue;

				if (seen_for_qname.insert(it->second).second)
				{
					outward_support[it->second] += n_members;
					distinct_outward_support[it->second] += 1;
				}
			}
		}

		vector<int> chimeric_support(idx_to_path.size(), 0);
		for (const pair<const string, vector<vector<int>>>& qname_entry : chim_bg)
		{
			map<string, vector<string>>::const_iterator mit = chim_rep_members.find(qname_entry.first);
			int n_members = (mit != chim_rep_members.end()) ? (int)mit->second.size() : 1;

			set<int> seen_for_qname;
			for (const vector<int>& path : qname_entry.second)
			{
				map<vector<int>, int>::iterator it = path_to_idx.find(path);
				if (it == path_to_idx.end()) continue;

				if (seen_for_qname.insert(it->second).second)
				{
					chimeric_support[it->second] += n_members;
					distinct_chimeric_support[it->second] += 1;
				}
			}
		}

		vector<double> min_vertex_weight(idx_to_path.size(), -1.0);
		vector<int> zero_coverage_vertices(idx_to_path.size(), 0);
		vector<int> min_edge_weight(idx_to_path.size(), -1);
		vector<int> fragment_lengths(idx_to_path.size(), -1);
		vector<int> weak_edge_count(idx_to_path.size(), 0);
		vector<int> n_vertices_feat(idx_to_path.size(), 0);
		vector<int> exon_count_feat(idx_to_path.size(), 0);

		for (int i = 0; i < idx_to_path.size(); i++)
		{
			const vector<int>& p = idx_to_path[i];
			if (p.empty()) continue;

			n_vertices_feat[i] = (int)p.size();
			exon_count_feat[i] = (int)merge_path_to_exons(bd, p).size();

			// Skip regions with zero coverage
			double min_vw = DBL_MAX;
			int zero_cov = 0;
			for (int v : p)
			{
				if (v >= 0 && v < bd->regions.size())
				{
					double w = bd->regions[v].ave;

					if (w == 0.0)
					{
						zero_cov++;
						continue;
					}

					if (w < min_vw)
						min_vw = w;
				}
			}

			min_vertex_weight[i] = (min_vw == DBL_MAX) ? -1.0 : min_vw;
			zero_coverage_vertices[i] = zero_cov;

			int min_ew = INT_MAX;
			int n_weak = 0;
			for (int j = 0; j + 1 < p.size(); j++)
			{
				int u = p[j];
				int v = p[j + 1];
				map<int, vector<pair<int, int>>>::iterator it_adj = splice_graph_adj.find(u);
				if (it_adj != splice_graph_adj.end())
				{
					for (const pair<int, int>& np : it_adj->second)
					{
						if (np.first == v)
						{
							if (np.second < min_ew)
								min_ew = np.second;
							if (np.second < min_junction_count)
								n_weak++;

							break;
						}
					}
				}
			}

			min_edge_weight[i] = (min_ew == INT_MAX) ? -1 : min_ew;
			weak_edge_count[i] = n_weak;
			fragment_lengths[i] = bd->compute_aligned_length(0, 0, p);
		}

		map<string, set<int>> chim_path_indices;
		for (const pair<const string, vector<vector<int>>>& e : chim_bg)
		{
			for (const vector<int>& cp : e.second)
			{
				map<vector<int>, int>::iterator it = path_to_idx.find(cp);
				if (it != path_to_idx.end())
					chim_path_indices[e.first].insert(it->second);
			}
		}

		map<int, vector<pair<string, bool>>> path_supporters;
		for (const pair<const string, vector<vector<int>>>& qname_entry : bg)
		{
			const string& qname = qname_entry.first;

			map<string, vector<string>>::const_iterator cmit = chim_rep_members.find(qname);
			map<string, vector<string>>::const_iterator omit = outward_rep_members.find(qname);
			const vector<string>* members =
				(cmit != chim_rep_members.end()) ? &cmit->second :
				(omit != outward_rep_members.end()) ? &omit->second : nullptr;

			set<int> seen;
			for (const vector<int>& q_path : qname_entry.second)
			{
				map<vector<int>, int>::iterator it = path_to_idx.find(q_path);
				if (it == path_to_idx.end()) continue;

				int idx = it->second;

				if (!seen.insert(idx).second) continue;

				bool is_chim_path = chim_path_indices.count(qname) && chim_path_indices.at(qname).count(idx);

				if (members)
				{
					for (const string& member_qname : *members)
						path_supporters[idx].push_back(make_pair(member_qname, is_chim_path));
				}
				else
				{
					path_supporters[idx].push_back(make_pair(qname, is_chim_path));
				}
			}
		}

		map<string, vector<read_info>> qname_segs;
		for (const read_info& r : bundle_chimeric_reads[bundle_idx])
			qname_segs[r.qname].push_back(r);

		for (pair<const string, vector<read_info>>& e : qname_segs)
		{
			sort(
				e.second.begin(), e.second.end(), [](const read_info& a, const read_info& b)
				{
					return a.pos < b.pos;
				}
			);
		}

		map<string, pair<const read_info*, const read_info*>> qname_outward;
		for (const read_info& r : bundle_outward_reads[bundle_idx])
		{
			pair<const read_info*, const read_info*>& lr = qname_outward[r.qname];
			if (!lr.first || r.pos < lr.first->pos)
				lr.first  = &r;

			if (!lr.second || r.rpos > lr.second->rpos)
				lr.second = &r;
		}

		vector<int> insert_sizes(idx_to_path.size(), -1);
		for (int i = 0; i < idx_to_path.size(); i++)
		{
			const vector<int>& p = idx_to_path[i];
			if (p.empty()) continue;

			map<int, vector<pair<string, bool>>>::iterator supp_it = path_supporters.find(i);
			if (supp_it == path_supporters.end()) continue;

			int full_len = fragment_lengths[i];

			int insert_sum = 0;
			int insert_count = 0;

			for (const pair<string, bool>& supporter : supp_it->second)
			{
				const string& qname = supporter.first;
				const bool is_chim_path = supporter.second;
				int insert_val = -1;

				if (is_chim_path)
				{
					map<string, vector<read_info>>::iterator segs_it = qname_segs.find(qname);
					if (segs_it == qname_segs.end() || segs_it->second.size() < 2) continue;

					const vector<read_info>& segs = segs_it->second;

					bool any_vlist_empty = false;
					for (const read_info& seg : segs)
					{
						if (seg.vlist.empty())
						{
							any_vlist_empty = true;
							break;
						}
					}

					if (any_vlist_empty) continue;

					if (segs.size() == 2)
						insert_val = (segs[0].rpos - segs[0].pos) + (segs[1].rpos - segs[1].pos);

					else if (segs.size() == 3)
					{
						const read_info& seg1 = segs[0];
						const read_info& seg2 = segs[1];
						const read_info& seg3 = segs[2];

						const read_info *S_ptr;
						const read_info *T_first_ptr;
						const read_info *T_second_ptr;

						const read_info *start_ptr;
						const read_info *mid_ptr;
						const read_info *end_ptr;

						int c_frag_start;
						int c_frag_end;

						if (!identify_chimeric_traversal(seg1, seg2, seg3, S_ptr, T_first_ptr, T_second_ptr, start_ptr, mid_ptr, end_ptr, c_frag_start, c_frag_end)) continue;

						int start_v = start_ptr->vlist.front();
						int end_v = end_ptr->vlist.back();

						if (start_v == end_v)
							insert_val = std::abs(c_frag_end - c_frag_start);

						else
						{
							vector<int> sub_path = extract_chimeric_subpath(p, start_v, end_v);
							if (sub_path.empty()) continue;
							insert_val = compute_chimeric_insert_size(sub_path, c_frag_start, c_frag_end, bd->regions);
						}
					}

					else
					{
						int c_frag_start = segs.front().pos;
						int c_frag_end = segs.back().rpos;

						int start_v = segs.front().vlist.front();
						int end_v = segs.back().vlist.back();

						if (start_v == end_v)
							insert_val = std::abs(c_frag_end - c_frag_start);
						
						else
						{
							vector<int> sub_path = extract_chimeric_subpath(p, start_v, end_v);
							if (sub_path.empty()) continue;
							insert_val = compute_chimeric_insert_size(sub_path, c_frag_start, c_frag_end, bd->regions);
						}
					}
				}

				else
				{
					map<string, pair<const read_info*, const read_info*>>::iterator ow_it = qname_outward.find(qname);
					if (ow_it == qname_outward.end()) continue;

					const read_info* left_rd = ow_it->second.first;
					const read_info* right_rd = ow_it->second.second;

					if (!left_rd || !right_rd || left_rd == right_rd) continue;
					if (left_rd->vlist.empty() || right_rd->vlist.empty()) continue;

					int middle = bridge_interior_len(bd, p, left_rd->vlist.back(), right_rd->vlist.front());
					insert_val = full_len - middle;
				}

				if (insert_val >= 0 && insert_val <= 2 * length_median)
				{
					insert_sum += insert_val;
					insert_count++;
				}
			}

			if (insert_count > 0)
				insert_sizes[i] = static_cast<int>(insert_sum / insert_count);
		}

		vector<int> kept_indices;
		for (int i = 0; i < idx_to_path.size(); i++)
		{
			if (!idx_to_path[i].empty())
				kept_indices.push_back(i);
		}

		vector<set<string>> path_chim_reads(idx_to_path.size());
		vector<set<string>> path_outward_reads(idx_to_path.size());
		for (const pair<const int, vector<pair<string, bool>>>& kv : path_supporters)
		{
			int idx = kv.first;
			for (const pair<string, bool>& sup : kv.second)
			{
				if (sup.second)
					path_chim_reads[idx].insert(sup.first);

				else
					path_outward_reads[idx].insert(sup.first);
			}
		}

		vector<bool> pruned(idx_to_path.size(), false);

		vector<int> circ_start(idx_to_path.size(), 0), circ_end(idx_to_path.size(), 0);
		for (int i : kept_indices)
		{
			const vector<int>& vp = idx_to_path[i];
			circ_start[i] = bd->regions[vp.front()].lpos;
			circ_end[i] = bd->regions[vp.back()].rpos;
		}

		vector<int> cross_tissue_repro_count(idx_to_path.size(), 0);
		vector<int> cross_tissue_repro_full_count(idx_to_path.size(), 0);
		vector<int> cross_tissue_chim_repro_count(idx_to_path.size(), 0);
		{
			const SiblingReproCounts& sibling_counts = get_sibling_tissue_repro_counts(outdir);
			for (int i : kept_indices)
			{
				tuple<string, int32_t, int32_t> key(bd->bb.chrm, circ_start[i] + 1, circ_end[i]);
				map<tuple<string, int32_t, int32_t>, int>::const_iterator it = sibling_counts.boundary.find(key);
				if (it != sibling_counts.boundary.end())
					cross_tissue_repro_count[i] = it->second;

				vector<pair<int32_t, int32_t>> exons = merge_path_to_exons(bd, idx_to_path[i]);
				if (exons.empty()) continue;

				string starts_str, ends_str;
				for (int ei = 0; ei < (int)exons.size(); ei++)
				{
					starts_str += std::to_string(exons[ei].first + 1);
					if (ei + 1 != (int)exons.size()) starts_str += ",";
				}
				for (int ei = 0; ei < (int)exons.size(); ei++)
				{
					ends_str += std::to_string(exons[ei].second);
					if (ei + 1 != (int)exons.size()) ends_str += ",";
				}

				tuple<string, int32_t, int32_t, string, string> full_key(bd->bb.chrm, circ_start[i] + 1, circ_end[i], starts_str, ends_str);
				map<tuple<string, int32_t, int32_t, string, string>, int>::const_iterator fit = sibling_counts.full_structure.find(full_key);
				if (fit != sibling_counts.full_structure.end())
					cross_tissue_repro_full_count[i] = fit->second;

				map<tuple<string, int32_t, int32_t>, int>::const_iterator cit = sibling_counts.chim_boundary.find(key);
				if (cit != sibling_counts.chim_boundary.end())
					cross_tissue_chim_repro_count[i] = cit->second;
			}
		}

		vector<double> circ_ratio(idx_to_path.size(), 0.0);
		for (int i : kept_indices)
			if (bundle_expr_ceiling > 0.0)
				circ_ratio[i] = (double)path_chim_reads[i].size() / bundle_expr_ceiling;

		// inverted_repeat_matches: dump-only (see inverted_repeat_kmer_matches's own comment
		// and the plan doc for the offline validation).
		vector<int> inverted_repeat_matches(idx_to_path.size(), 0);
		for (int i : kept_indices)
			inverted_repeat_matches[i] = inverted_repeat_kmer_matches(bd, circ_start[i], circ_end[i]);

		// internal_junction_noncanonical_count: dump-only (see count_noncanonical_internal_junctions's
		// own comment and the plan doc for the offline validation).
		vector<int> internal_junction_noncanonical_count(idx_to_path.size(), 0);
		for (int i : kept_indices)
		{
			vector<pair<int32_t, int32_t>> exons = merge_path_to_exons(bd, idx_to_path[i]);
			internal_junction_noncanonical_count[i] = count_noncanonical_internal_junctions(bd, exons, bd->bb.strand);
		}

		vector<int> anchor_length(idx_to_path.size(), 0);
		vector<int> fake_supple_count(idx_to_path.size(), 0);
		vector<int> supple_len(idx_to_path.size(), 0);
		{
			map<string, const hit*> qname_to_hit;
			for (const read_info& r : bundle_chimeric_reads[bundle_idx])
				if (r.src != NULL && r.src->suppl != NULL)
					qname_to_hit[r.qname] = r.src;

			for (int i : kept_indices)
			{
				int best = 0;
				int fake_count = 0;
				int best_supple_len = 0;
				for (const string& qn : path_chim_reads[i])
				{
					map<string, const hit*>::iterator it = qname_to_hit.find(qn);
					if (it == qname_to_hit.end()) continue;

					const hit* h = it->second;
					int h_anchor = max(h->left_cigar_len, h->right_cigar_len);
					int s_anchor = max(h->suppl->left_cigar_len, h->suppl->right_cigar_len);
					int read_anchor = min(h_anchor, s_anchor);

					if (read_anchor > best)
						best = read_anchor;

					if (h->is_fake || h->suppl->is_fake)
						fake_count++;

					int32_t this_supple_len = 0;
					for (size_t ci = 0; ci < h->suppl->cigar_vector.size(); ci++)
						if (h->suppl->cigar_vector[ci].first == 'M')
							this_supple_len += h->suppl->cigar_vector[ci].second;
					if (this_supple_len > best_supple_len)
						best_supple_len = this_supple_len;
				}
				anchor_length[i] = best;
				fake_supple_count[i] = fake_count;
				supple_len[i] = best_supple_len;
			}
		}

		map<int,int> uf_parent;
		for (int i : kept_indices)
			uf_parent[i] = i;

		vector<int> sorted_by_start = kept_indices;
		sort(sorted_by_start.begin(), sorted_by_start.end(), [&](int x, int y) {
			return circ_start[x] < circ_start[y];
		});

		for (int ai = 0; ai < (int)sorted_by_start.size(); ai++)
		{
			int a = sorted_by_start[ai];
			for (int bi = ai + 1; bi < (int)sorted_by_start.size(); bi++)
			{
				int b = sorted_by_start[bi];
				if (circ_start[b] - circ_start[a] > read_length) break;
				if (uf_find(uf_parent, a) == uf_find(uf_parent, b)) continue;

				if (is_shifted_duplicate_bsj(bd, circ_start[a], circ_end[a], circ_start[b], circ_end[b]))
					uf_parent[uf_find(uf_parent, b)] = uf_find(uf_parent, a);
			}
		}

		map<int, vector<int>> groups;
		for (int i : kept_indices)
			groups[uf_find(uf_parent, i)].push_back(i);

		for (const pair<const int, vector<int>>& g : groups)
		{
			if (g.second.size() < 2) continue;

			int rep = pick_best_ranked(g.second, path_chim_reads, path_outward_reads, set<string>(), min_edge_weight, insert_sizes, min_vertex_weight, length_median, n_vertices_feat, exon_count_feat, weak_edge_count, circ_start, circ_end, cross_tissue_repro_count, cross_tissue_repro_full_count, circ_ratio, chimeric_support, cross_tissue_chim_repro_count, fake_supple_count, supple_len).first;

			for (int i : g.second)
			{
				if (i == rep) continue;

				path_chim_reads[rep].insert(path_chim_reads[i].begin(), path_chim_reads[i].end());
				path_outward_reads[rep].insert(path_outward_reads[i].begin(), path_outward_reads[i].end());
				pruned[i] = true;
			}
		}

		// Dominance pruning: discard path a if there exists path b such that all hold:
		//   (1) chim_reads(a) ⊆ chim_reads(b) AND outward_reads(a) ⊆ outward_reads(b)
		//   (2) min_edge_weight[a] < min_edge_weight[b]
		//   (3) |insert_sizes[a] - length_median| > |insert_sizes[b] - length_median|
		// Conditions checked cheapest-first (2, 3, then 1) to terminate early.
		for (int ai = 0; ai < kept_indices.size(); ai++)
		{
			int a = kept_indices[ai];
			if (pruned[a]) continue;

			for (int bi = 0; bi < kept_indices.size(); bi++)
			{
				if (ai == bi) continue;

				int b = kept_indices[bi];

				// Condition 2: both edge weights valid, a's strictly lower
				if (min_edge_weight[a] < 0 || min_edge_weight[b] < 0) continue;

				if (min_edge_weight[a] >= min_edge_weight[b]) continue;

				// Condition 3: a's insert size farther from median than b's
				if (insert_sizes[a] < 0 || insert_sizes[b] < 0) continue;

				int dist_a = abs(insert_sizes[a] - (int)length_median);
				int dist_b = abs(insert_sizes[b] - (int)length_median);

				if (dist_a <= dist_b) continue;

				// Condition 1: a's read sets are subsets of b's
				if (! is_subset(path_chim_reads[a], path_chim_reads[b])) continue;

				if (! is_subset(path_outward_reads[a], path_outward_reads[b])) continue;

				pruned[a] = true;
				break;
			}
		}

		vector<int> surviving;
		surviving.reserve(kept_indices.size());
		for (int i : kept_indices)
		{
			if (!pruned[i])
				surviving.push_back(i);
		}

		kept_indices = std::move(surviving);

		vector<int> gated;
		gated.reserve(kept_indices.size());

		{
			std::string pool_filename = outdir + "/circrna_candidate_pool.tsv";
			std::ofstream pool_out(pool_filename, std::ios_base::app);
			pool_out << "#bundle_idx\tpath_idx\tchrm\tstrand\tcirc_start\tcirc_end\t"
			         << "exon_count\texon_starts\texon_ends\tn_vertices\t"
			         << "chimeric_support_reads\toutward_support_reads\t"
			         << "chimeric_support_weighted\toutward_support_weighted\t"
			         << "min_edge_weight\tmin_edge_weight_valid\t"
			         << "min_vertex_weight\tmin_vertex_weight_valid\t"
			         << "zero_coverage_vertices\t"
			         << "insert_size\tinsert_size_deviation\tinsert_size_valid\t"
			         << "length_median\tfragment_length\t"
			         << "total_exon_len\tmax_exon_len\tmin_exon_len\tavg_exon_len\t"
			         << "weak_edge_count\tanchor_length\tcross_tissue_repro_count\tcross_tissue_repro_full_count\tbundle_expr_ceiling\tcirc_ratio\tpre_gate_cc_size\tcross_tissue_chim_repro_count\tinverted_repeat_matches\tinternal_junction_noncanonical_count\tfake_supple_count\tsupple_len\tdistinct_chimeric_support\tdistinct_outward_support\n";

			for (int i : kept_indices)
			{
				vector<pair<int32_t, int32_t>> exons = merge_path_to_exons(bd, idx_to_path[i]);
				if (exons.empty()) continue;

				int32_t total_exon_len = 0, max_exon_len = 0;
				int32_t min_exon_len = exons[0].second - exons[0].first;
				for (const pair<int32_t, int32_t>& ex : exons)
				{
					int32_t elen = ex.second - ex.first;
					total_exon_len += elen;
					if (elen > max_exon_len) max_exon_len = elen;
					if (elen < min_exon_len) min_exon_len = elen;
				}
				double avg_exon_len = (double)total_exon_len / exons.size();

				pool_out << bundle_idx << "\t" << i << "\t" << bd->bb.chrm << "\t" << bd->bb.strand << "\t"
				         << (circ_start[i] + 1) << "\t" << circ_end[i] << "\t" << exons.size() << "\t";

				for (int ei = 0; ei < (int)exons.size(); ei++)
					pool_out << (exons[ei].first + 1) << (ei + 1 == (int)exons.size() ? "" : ",");
				pool_out << "\t";

				for (int ei = 0; ei < (int)exons.size(); ei++)
					pool_out << exons[ei].second << (ei + 1 == (int)exons.size() ? "" : ",");
				pool_out << "\t";

				pool_out << idx_to_path[i].size() << "\t"
				         << path_chim_reads[i].size() << "\t" << path_outward_reads[i].size() << "\t"
				         << chimeric_support[i] << "\t" << outward_support[i] << "\t";

				if (min_edge_weight[i] >= 0)
					pool_out << min_edge_weight[i] << "\t1\t";
				else
					pool_out << "NA\t0\t";

				if (min_vertex_weight[i] >= 0.0)
					pool_out << min_vertex_weight[i] << "\t1\t";
				else
					pool_out << "NA\t0\t";

				pool_out << zero_coverage_vertices[i] << "\t";

				if (insert_sizes[i] >= 0)
					pool_out << insert_sizes[i] << "\t" << abs(insert_sizes[i] - (int32_t)length_median) << "\t1\t";
				else
					pool_out << "NA\tNA\t0\t";

				pool_out << length_median << "\t" << fragment_lengths[i] << "\t"
				         << total_exon_len << "\t" << max_exon_len << "\t" << min_exon_len << "\t" << avg_exon_len << "\t"
				         << weak_edge_count[i] << "\t" << anchor_length[i] << "\t" << cross_tissue_repro_count[i] << "\t" << cross_tissue_repro_full_count[i] << "\t" << bundle_expr_ceiling << "\t" << circ_ratio[i] << "\t" << pre_gate_cc_size[i] << "\t" << cross_tissue_chim_repro_count[i] << "\t" << inverted_repeat_matches[i] << "\t" << internal_junction_noncanonical_count[i] << "\t" << fake_supple_count[i] << "\t" << supple_len[i] << "\t" << distinct_chimeric_support[i] << "\t" << distinct_outward_support[i] << "\n";

				bool mew_valid = min_edge_weight[i] >= 0;
				bool mvw_valid = min_vertex_weight[i] >= 0.0;
				bool isd_valid = insert_sizes[i] >= 0;

				bool gate_pass = passes_pool_gate(
					mew_valid ? min_edge_weight[i] : 0, mew_valid,
					mvw_valid ? min_vertex_weight[i] : 0.0, mvw_valid,
					isd_valid ? (int)abs(insert_sizes[i] - (int32_t)length_median) : 0, isd_valid,
					(int)path_chim_reads[i].size(), (int)path_outward_reads[i].size(),
					zero_coverage_vertices[i], (int)idx_to_path[i].size(), (int)exons.size(),
					total_exon_len, max_exon_len, min_exon_len, avg_exon_len,
					fragment_lengths[i], weak_edge_count[i], cross_tissue_repro_count[i], cross_tissue_repro_full_count[i],
					circ_ratio[i], chimeric_support[i], cross_tissue_chim_repro_count[i], fake_supple_count[i], supple_len[i],
					anchor_length[i], distinct_chimeric_support[i]);

				if (weak_edge_count[i] > 5)
					gate_pass = false;

				if (exons.size() > 21)
					gate_pass = false;

				bool is_single_exon_chimeric_only = (exons.size() == 1)
					&& !path_chim_reads[i].empty() && path_outward_reads[i].empty();
				if (is_single_exon_chimeric_only
					&& (int)path_chim_reads[i].size() < 3
					&& !(mvw_valid && min_vertex_weight[i] >= 5.0))
					gate_pass = false;

				if (gate_pass)
					gated.push_back(i);
			}
		}

		kept_indices = std::move(gated);

		if (kept_indices.empty())
		{
			bgout << "No Bipartite Graph\n\n\n";
			bgout.close();
			return;
		}

		map<int, int> old_to_new;
		for (int ni = 0; ni < kept_indices.size(); ni++)
			old_to_new[kept_indices[ni]] = ni;

		int global_max_ew = 0;
		for (int i : kept_indices)
		{
			if (min_edge_weight[i] > global_max_ew)
				global_max_ew = min_edge_weight[i];
		}

		bgout << "Unique Paths (" << kept_indices.size() << ")\n";
		bgout << "---------------------------------\n";
		for (int ni = 0; ni < kept_indices.size(); ni++)
		{
			int i = kept_indices[ni];
			bgout << "P" << ni + 1 << " = (";
			for (int j = 0; j < idx_to_path[i].size(); j++)
				bgout << idx_to_path[i][j] << (j == idx_to_path[i].size() - 1 ? "" : " → ");

			bgout << ")\n";
			bgout << "\tSupporting (Merged) Outward Reads = " << outward_support[i] << "\n";
			bgout << "\tSupporting (Merged) Chimeric Reads = " << chimeric_support[i] << "\n";
			
			bgout << "\tMinimum Vertex Weight = ";
			
			if (min_vertex_weight[i] < 0.0)
				bgout << "N/A";

			else
				bgout << min_vertex_weight[i];
			
			bgout << "\n";
			bgout << "\tZero-Coverage Vertices = " << zero_coverage_vertices[i] << "\n";
			bgout << "\tMinimum Edge Weight = ";
			
			if (min_edge_weight[i] < 0)
				bgout << "N/A";

			else
				bgout << min_edge_weight[i];
			
			bgout << "\n";
			bgout << "\tInsert Size = " << insert_sizes[i] << "\n";
			bgout << "\tFull-Sequence Length = " << fragment_lengths[i] << "\n";
		}

		bgout << "\n";

		map<string, set<int>> qname_to_active_paths;
		for (const pair<const string, vector<vector<int>>>& qname_entry : bg)
		{
			const string& qname = qname_entry.first;
			for (const vector<int>& path : qname_entry.second)
			{
				map<vector<int>, int>::iterator pit = path_to_idx.find(path);
				if (pit == path_to_idx.end()) continue;
				int idx = pit->second;
				if (!old_to_new.count(idx)) continue;
				qname_to_active_paths[qname].insert(idx);
			}
		}

		uf_parent.clear();
		for (int i : kept_indices)
			uf_parent[i] = i;

		for (const pair<const string, set<int>>& qp : qname_to_active_paths)
		{
			const set<int>& path_set = qp.second;

			if (path_set.size() <= 1) continue;

			int first = uf_find(uf_parent, *path_set.begin());
			for (int idx : path_set)
			{
				int r = uf_find(uf_parent, idx);
				if (r != first)
					uf_parent[r] = first;
			}
		}

		map<int, vector<int>> comp_paths;
		map<int, vector<string>> comp_qnames;

		for (int i : kept_indices)
			comp_paths[uf_find(uf_parent, i)].push_back(i);

		for (const pair<const string, set<int>>& qp : qname_to_active_paths)
		{
			if (qp.second.empty()) continue;
			comp_qnames[uf_find(uf_parent, *qp.second.begin())].push_back(qp.first);
		}

		bgout << "CCs (" << comp_paths.size() << ")\n\n";
		int comp_num = 1;
		for (pair<const int, vector<int>>& comp_entry : comp_paths)
		{
			int root = comp_entry.first;
			vector<int>& comp_path_indices = comp_entry.second;

			sort(
				comp_path_indices.begin(), comp_path_indices.end(), [&](int a, int b)
				{
					return old_to_new[a] < old_to_new[b];
				}
			);

			const vector<string>& comp_qname_list = comp_qnames.count(root) ? comp_qnames.at(root) : vector<string>{};
			int n_comp_reads = comp_qname_list.size();

			bgout << "CC " << comp_num
					<< " (" << comp_path_indices.size() << " path"
					<< (comp_path_indices.size() == 1 ? "" : "s") << ", "
					<< n_comp_reads << " read pair"
					<< (n_comp_reads == 1 ? "" : "s") << ")\n";
			bgout << "---------------------------------\n";

			bgout << "\tPaths: ";
			for (int ci = 0; ci < comp_path_indices.size(); ci++)
			{
				bgout << "P" << old_to_new[comp_path_indices[ci]] + 1;
				if (ci + 1 < comp_path_indices.size())
					bgout << ", ";
			}
			
			bgout << "\n\n";

			bgout << "\tBipartite Graph (" << n_comp_reads << " read pair" << (n_comp_reads == 1 ? "" : "s") << ")\n";

			bgout << "\t---------------------------------\n";

			for (const string& qname : comp_qname_list)
			{
				const set<int>& active_paths = qname_to_active_paths.at(qname);
				vector<int> surviving_new_indices;
				for (int idx : active_paths)
					surviving_new_indices.push_back(old_to_new[idx]);

				sort(surviving_new_indices.begin(), surviving_new_indices.end());

				string bracket;

				bool is_chimeric = bundle_chimeric_merged_paths.count(bundle_idx) && bundle_chimeric_merged_paths.at(bundle_idx).count(qname);

				if (is_chimeric)
				{
					vector<read_info> chim_segs;
					for (const read_info& read : bundle_chimeric_reads[bundle_idx])
					{
						if (read.qname == qname)
							chim_segs.push_back(read);
					}

					sort(
						chim_segs.begin(), chim_segs.end(), [](const read_info& a, const read_info& b)
						{
							return a.pos < b.pos;
						}
					);

					for (const read_info& seg : chim_segs)
					{
						string lbl;
						if (seg.is_read1)
							lbl = "R1";

						else if (seg.is_read2)
							lbl = "R2";

						else
							lbl = "Unknown";

						lbl += " ";

						if (seg.is_primary)
							lbl += "Primary";

						else if (seg.is_supplementary)
							lbl += "Supplementary";

						string loc = bd->bb.chrm + ":" + to_string(seg.pos) + "-" + to_string(seg.rpos);

						if (!bracket.empty())
							bracket += ", ";

						bracket += lbl + " " + loc;
					}
				}

				else
				{
					string r1_str, r2_str;
					for (const read_info& read : bundle_outward_reads[bundle_idx])
					{
						if (read.qname != qname) continue;

						string loc = bd->bb.chrm + ":" + to_string(read.pos) + "-" + to_string(read.rpos);
						
						if (read.is_read1)
							r1_str = "R1 " + loc;

						else if (read.is_read2)
							r2_str = "R2 " + loc;
					}

					if (!r1_str.empty())
						bracket += r1_str;

					if (!r2_str.empty())
					{
						if (!bracket.empty())
							bracket += ", ";

						bracket += r2_str;
					}
				}

				bgout << "\t" << qname;

				if (!bracket.empty())
					bgout << " [" << bracket << "]";

				bgout << "\n";
				for (int ni : surviving_new_indices)
					bgout << "\t\t- P" << ni + 1 << "\n";
			}

			bgout << "\n";
			comp_num++;
		}

		std::string circrna_filename = outdir + "/circrna_paths.txt";
		std::ofstream circ_out(circrna_filename, std::ios_base::app);

		static int gtf_circ_id = 0;
		std::string gtf_filename = outdir + "/terrace_2.0.gtf";
		std::ofstream gtf_out(gtf_filename, std::ios_base::app);

		std::string diag_filename = outdir + "/circrna_diagnostics.tsv";
		std::ofstream diag_out(diag_filename, std::ios_base::app);
		diag_out << "#bundle_idx\tcc_id\titeration\tpath_id\tchrm\tstrand\t"
		         << "circ_start\tcirc_end\tn_vertices\tresidual_chimeric\tresidual_outward\t"
		         << "min_edge_weight\tmin_vertex_weight\tinsert_size\tinsert_size_deviation\tgate_result\n";

		std::string feat_filename = outdir + "/circrna_features.tsv";
		std::ofstream feat_out(feat_filename, std::ios_base::app);
		feat_out << "#tid\tchrm\tstrand\tcirc_start\tcirc_end\tn_vertices\t"
		         << "residual_chimeric\tresidual_outward\tmin_edge_weight\t"
		         << "min_vertex_weight\tinsert_size\tinsert_size_deviation\n";

		circ_out << "=================================================\n";
		circ_out << "Bundle " << bundle_idx << "\n";
		circ_out << "=================================================\n\n";

		int total_selected = 0;
		int comp_num_circ = 1;

		struct BundleSelectedEntry
		{
			int idx;
			int comp_num_circ;
			int iteration;
			set<string> residual_chim;
			set<string> residual_out;
			char strand;
			double score;
		};
		vector<BundleSelectedEntry> bundle_wide_selected;

		for (pair<const int, vector<int>>& comp_entry : comp_paths)
		{
			const vector<int>& comp_path_indices = comp_entry.second;

			set<int> available(comp_path_indices.begin(), comp_path_indices.end());
			set<string> explained;
			const set<string> orig_explained;

			vector<pair<int, double>> selected;
			// Reads at time of selection (chimeric, outward)
			vector<pair<set<string>, set<string>>> selected_reads;
			vector<char> selected_strand;

			int iteration = 0;

			while (!available.empty())
			{
				// Stop once every available path has all its reads already explained.
				bool any_unexplained = false;
				for (int idx : available)
				{
					for (const string& r : path_chim_reads[idx])
					{
						if (!explained.count(r))
						{
							any_unexplained = true;
							break;
						}
					}

					if (any_unexplained) break;

					for (const string& r : path_outward_reads[idx])
					{
						if (!explained.count(r))
						{
							any_unexplained = true;
							break;
						}
					}
					
					if (any_unexplained) break;
				}

				if (!any_unexplained) break;

				vector<int> avail_vec;
				for (int idx : available)
				{
					bool has_unexplained = false;
					for (const string& r : path_chim_reads[idx])
					{
						if (!explained.count(r))
						{
							has_unexplained = true;
							break;
						}
					}

					if (!has_unexplained)
					{
						for (const string& r : path_outward_reads[idx])
						{
							if (!explained.count(r))
							{
								has_unexplained = true;
								break;
							}
						}
					}

					if (has_unexplained)
						avail_vec.push_back(idx);
				}

				if (avail_vec.empty()) break;

				pair<int, double> best = pick_best_ranked(avail_vec, path_chim_reads, path_outward_reads, orig_explained, min_edge_weight, insert_sizes, min_vertex_weight, length_median, n_vertices_feat, exon_count_feat, weak_edge_count, circ_start, circ_end, cross_tissue_repro_count, cross_tissue_repro_full_count, circ_ratio, chimeric_support, cross_tissue_chim_repro_count, fake_supple_count, supple_len);
				int best_path = best.first;
				double best_rank_sum = best.second;

				{
					std::string sel_filename = outdir + "/circrna_selection_pool.tsv";
					std::ofstream sel_out(sel_filename, std::ios_base::app);
					sel_out << "#bundle_idx\tcc_id\titeration\tpath_id\tis_winner\tchrm\tstrand\t"
					        << "circ_start\tcirc_end\texon_starts\texon_ends\tn_vertices\t"
					        << "unexplained_chim_reads\tunexplained_outward_reads\t"
					        << "min_edge_weight\tmin_vertex_weight\tinsert_size\tinsert_size_deviation\t"
					        << "weak_edge_count\tcross_tissue_repro_count\tcross_tissue_repro_full_count\tbundle_expr_ceiling\tcirc_ratio\tpre_gate_cc_size\tcross_tissue_chim_repro_count\tinverted_repeat_matches\tinternal_junction_noncanonical_count\tchimeric_support_weighted\toutward_support_weighted\tfake_supple_count\tsupple_len\n";

					for (int idx : avail_vec)
					{
						vector<pair<int32_t, int32_t>> sel_exons = merge_path_to_exons(bd, idx_to_path[idx]);
						if (sel_exons.empty()) continue;

						int sel_cc = 0;
						for (const string& r : path_chim_reads[idx])
							if (!explained.count(r)) sel_cc++;

						int sel_oc = 0;
						for (const string& r : path_outward_reads[idx])
							if (!explained.count(r)) sel_oc++;

						sel_out << bundle_idx << "\t" << comp_num_circ << "\t" << iteration + 1 << "\t" << idx << "\t"
						        << (idx == best_path ? 1 : 0) << "\t" << bd->bb.chrm << "\t" << bd->bb.strand << "\t"
						        << (circ_start[idx] + 1) << "\t" << circ_end[idx] << "\t";

						for (int ei = 0; ei < (int)sel_exons.size(); ei++)
							sel_out << (sel_exons[ei].first + 1) << (ei + 1 == (int)sel_exons.size() ? "" : ",");
						sel_out << "\t";

						for (int ei = 0; ei < (int)sel_exons.size(); ei++)
							sel_out << sel_exons[ei].second << (ei + 1 == (int)sel_exons.size() ? "" : ",");
						sel_out << "\t";

						sel_out << idx_to_path[idx].size() << "\t" << sel_cc << "\t" << sel_oc << "\t";

						if (min_edge_weight[idx] >= 0) sel_out << min_edge_weight[idx]; else sel_out << "NA";
						sel_out << "\t";

						if (min_vertex_weight[idx] >= 0.0) sel_out << min_vertex_weight[idx]; else sel_out << "NA";
						sel_out << "\t";

						if (insert_sizes[idx] >= 0)
							sel_out << insert_sizes[idx] << "\t" << abs(insert_sizes[idx] - (int32_t)length_median);
						else
							sel_out << "NA\tNA";

						sel_out << "\t" << weak_edge_count[idx] << "\t" << cross_tissue_repro_count[idx] << "\t" << cross_tissue_repro_full_count[idx] << "\t" << bundle_expr_ceiling << "\t" << circ_ratio[idx] << "\t" << pre_gate_cc_size[idx] << "\t" << cross_tissue_chim_repro_count[idx] << "\t" << inverted_repeat_matches[idx] << "\t" << internal_junction_noncanonical_count[idx] << "\t" << chimeric_support[idx] << "\t" << outward_support[idx] << "\t" << fake_supple_count[idx] << "\t" << supple_len[idx] << "\n";
					}
				}

				if (best_path < 0) break;

				set<string> residual_chim;
				set<string> residual_out;

				for (const string& r : path_chim_reads[best_path])
				{
					if (!explained.count(r))
						residual_chim.insert(r);
				}

				for (const string& r : path_outward_reads[best_path])
				{
					if (!explained.count(r))
						residual_out.insert(r);
				}

				const vector<int>& best_vpath = idx_to_path[best_path];
				int32_t circ_start = bd->regions[best_vpath.front()].lpos;
				int32_t circ_end = bd->regions[best_vpath.back()].rpos;

				size_t total_support = residual_chim.size() + residual_out.size();
				bool fails_evidence_gate = (total_support < 2);
				if (fails_evidence_gate && cross_tissue_repro_full_count[best_path] >= 2)
					fails_evidence_gate = false;
				bool has_chimeric_bsj_support = !path_chim_reads[best_path].empty();
				if (!fails_evidence_gate && !has_chimeric_bsj_support && min_edge_weight[best_path] >= 0 && min_edge_weight[best_path] < min_junction_count)
					fails_evidence_gate = true;

				vector<pair<int32_t, int32_t>> best_exons = merge_path_to_exons(bd, best_vpath);

				bool fails_splice_signal_gate = false;
				char resolved_strand = bd->bb.strand;
				if (!fails_evidence_gate)
				{
					fails_splice_signal_gate = !is_canonical_bsj_boundary(bd, circ_start, circ_end, bd->bb.strand, &resolved_strand, true);

					char junction_strand = (bd->bb.strand == '.') ? resolved_strand : bd->bb.strand;
					if (!fails_splice_signal_gate && best_exons.size() > 1
						&& (junction_strand == '+' || junction_strand == '-'))
					{
						for (size_t k = 0; k + 1 < best_exons.size() && !fails_splice_signal_gate; k++)
						{
							int32_t up_end = best_exons[k].second;
							int32_t down_start = best_exons[k + 1].first;

							if (!is_canonical_bsj_boundary(bd, down_start, up_end, junction_strand, nullptr, true))
								fails_splice_signal_gate = true;
						}
					}
				}

				bool fails_contig_gate = !is_primary_chromosome(bd->bb.chrm);

				bool fails_structural_gate = false;

				bool structural_rescue = cross_tissue_repro_full_count[best_path] >= 1;

				if ((int)best_exons.size() > max_circ_vsize)
				{
					if (!structural_rescue)
						fails_structural_gate = true;
				}
				else if (best_exons.size() == 1)
				{
					if (best_exons[0].second - best_exons[0].first > max_single_exon_length)
						fails_structural_gate = true;
				}
				else if (!structural_rescue)
				{
					for (const pair<int32_t, int32_t>& ex : best_exons)
					{
						if (ex.second - ex.first > max_multi_exon_length)
						{
							fails_structural_gate = true;
							break;
						}
					}
				}

				iteration++;

				{
					string gate_result = "none";
					if (fails_evidence_gate) gate_result = "evidence";
					else if (fails_splice_signal_gate) gate_result = "splice_signal";
					else if (fails_contig_gate) gate_result = "contig";
					else if (fails_structural_gate) gate_result = "structural";

					diag_out << bundle_idx << "\t" << comp_num_circ << "\t" << iteration << "\t"
					         << best_path << "\t" << bd->bb.chrm << "\t" << bd->bb.strand << "\t"
					         << circ_start << "\t" << circ_end << "\t" << best_vpath.size() << "\t"
					         << residual_chim.size() << "\t" << residual_out.size() << "\t";

					if (min_edge_weight[best_path] >= 0) diag_out << min_edge_weight[best_path]; else diag_out << "NA";
					diag_out << "\t";

					if (min_vertex_weight[best_path] >= 0.0) diag_out << min_vertex_weight[best_path]; else diag_out << "NA";
					diag_out << "\t";

					if (insert_sizes[best_path] >= 0)
					{
						diag_out << insert_sizes[best_path] << "\t" << abs(insert_sizes[best_path] - (int32_t)length_median);
					}
					else
					{
						diag_out << "NA\tNA";
					}

					diag_out << "\t" << gate_result << "\n";
				}

				if (!fails_evidence_gate && !fails_splice_signal_gate
					&& !fails_contig_gate && !fails_structural_gate)
				{
					selected.push_back({best_path, best_rank_sum});
					selected_reads.push_back({residual_chim, residual_out});
					selected_strand.push_back(resolved_strand);

					// Mark all reads supporting the chosen path as explained
					for (const string& r : path_chim_reads[best_path])
						explained.insert(r);

					for (const string& r : path_outward_reads[best_path])
						explained.insert(r);
				}

				available.erase(best_path);
			}

			circ_out << "CC " << comp_num_circ
			         << " (" << selected.size() << " Selected Path"
			         << (selected.size() == 1 ? ")" : "s)") << "\n";
			circ_out << "---------------------------------\n";

			for (int rank = 0; rank < selected.size(); rank++)
			{
				int idx = selected[rank].first;
				double score = selected[rank].second;
				int ni = old_to_new[idx];

				circ_out << "  Rank " << rank + 1 << ": P" << ni + 1 << " = (";
				for (int j = 0; j < idx_to_path[idx].size(); j++)
				{
					circ_out << idx_to_path[idx][j];
					if (j + 1 < idx_to_path[idx].size())
						circ_out << " → ";
				}

				circ_out << ")\n";

				circ_out << "    Rank Sum          = " << score << "\n";
				const set<string>& residual_cr = selected_reads[rank].first;
				const set<string>& residual_or = selected_reads[rank].second;
				circ_out << "    Chimeric Support  = " << residual_cr.size()  << "\n";
				circ_out << "    Outward Support   = " << residual_or.size()  << "\n";

				circ_out << "    Min Edge Weight   = ";

				if (min_edge_weight[idx] >= 0)
					circ_out << min_edge_weight[idx];

				else
					circ_out << "N/A";

				circ_out << "\n";

				circ_out << "    Insert Size       = " << insert_sizes[idx] << " (|deviation| = ";

				if (insert_sizes[idx] >= 0)
					circ_out << abs(insert_sizes[idx] - length_median);

				else
					circ_out << "N/A";

				circ_out << ")\n";

				circ_out << "    Full-Seq Length   = " << fragment_lengths[idx] << "\n";

				circ_out << "    Chimeric Reads (" << residual_cr.size() << "): [";

				bool first = true;
				for (const string& r : residual_cr)
				{
					if (!first)
						circ_out << ", ";

					circ_out << r;
					first = false;
				}

				circ_out << "]\n";

				circ_out << "    Outward Reads  (" << residual_or.size() << "): [";

				bool first_or = true;
				for (const string& r : residual_or)
				{
					if (!first_or)
						circ_out << ", ";

					circ_out << r;
					first_or = false;
				}
				circ_out << "]\n\n";

				if (!idx_to_path[idx].empty())
					bundle_wide_selected.push_back({idx, comp_num_circ, rank + 1, residual_cr, residual_or, selected_strand[rank], score});
			}

			total_selected += selected.size();
			comp_num_circ++;
		}

		vector<bool> keep_bundle_wide(bundle_wide_selected.size(), true);
		vector<string> bundle_wide_drop_reason(bundle_wide_selected.size(), "");
		{
			map<pair<int,int>, vector<int>> groups;
			for (int i = 0; i < (int)bundle_wide_selected.size(); i++)
			{
				int idx = bundle_wide_selected[i].idx;
				groups[{circ_start[idx], circ_end[idx]}].push_back(i);
			}

			for (const pair<const pair<int,int>, vector<int>>& g : groups)
			{
				if (g.second.size() <= 1) continue;

				int best_i = g.second[0];
				for (int i : g.second)
					if (min_edge_weight[bundle_wide_selected[i].idx] > min_edge_weight[bundle_wide_selected[best_i].idx])
						best_i = i;

				for (int i : g.second)
				{
					if (i == best_i) continue;
					if (min_edge_weight[bundle_wide_selected[i].idx] < min_junction_count
						&& cross_tissue_repro_full_count[bundle_wide_selected[i].idx] < 2)
					{
						keep_bundle_wide[i] = false;
						bundle_wide_drop_reason[i] = "cross_cc_duplicate";
					}
				}
			}
		}

		for (int i = 0; i < (int)bundle_wide_selected.size(); i++)
		{
			const BundleSelectedEntry& e = bundle_wide_selected[i];
			int idx = e.idx;
			const set<string>& residual_cr = e.residual_chim;
			const set<string>& residual_or = e.residual_out;

			if (!keep_bundle_wide[i])
			{
				diag_out << bundle_idx << "\t" << e.comp_num_circ << "\t" << e.iteration << "\t"
				         << idx << "\t" << bd->bb.chrm << "\t" << bd->bb.strand << "\t"
				         << circ_start[idx] << "\t" << circ_end[idx] << "\t" << idx_to_path[idx].size() << "\t"
				         << residual_cr.size() << "\t" << residual_or.size() << "\t";

				if (min_edge_weight[idx] >= 0) diag_out << min_edge_weight[idx]; else diag_out << "NA";
				diag_out << "\t";
				if (min_vertex_weight[idx] >= 0.0) diag_out << min_vertex_weight[idx]; else diag_out << "NA";
				diag_out << "\t";
				if (insert_sizes[idx] >= 0)
					diag_out << insert_sizes[idx] << "\t" << abs(insert_sizes[idx] - (int32_t)length_median);
				else
					diag_out << "NA\tNA";
				diag_out << "\t" << bundle_wide_drop_reason[i] << "\n";
				continue;
			}

			const vector<int>& vpath = idx_to_path[idx];

			// Merge adjacent/overlapping regions into exon intervals
			vector<pair<int32_t, int32_t>> exons = merge_path_to_exons(bd, vpath);

			int32_t out_circ_start = exons.front().first;
			int32_t out_circ_end = exons.back().second;
			int support = residual_cr.size() + residual_or.size();

			string tid;
			bool first = true;

			for (const string& r : residual_cr)
			{
				if (!first)
					tid += "|";
				tid += r;
				first = false;
			}

			for (const string& r : residual_or)
			{
				if (!first)
					tid += "|";
				tid += r;
				first = false;
			}

			gtf_out << bd->bb.chrm << "\tTERRACE_2.0\tcircRNA\t"
					<< (out_circ_start + 1) << "\t" << out_circ_end << "\t"
					<< support << "\t" << e.strand << "\t.\t"
					<< "gene_id \"gene\"; "
					<< "transcript_id \"" << tid << "\"; "
					<< "cov \"" << support << "\";\n";

			for (int ei = 0; ei < exons.size(); ei++)
			{
				gtf_out << bd->bb.chrm << "\tTERRACE_2.0\texon\t"
						<< (exons[ei].first + 1) << "\t" << exons[ei].second << "\t"
						<< support << "\t" << e.strand << "\t.\t"
						<< "gene_id \"gene\"; "
						<< "transcript_id \"" << tid << "\"; "
						<< "exon_number \"" << (ei + 1) << "\";\n";
			}

			feat_out << tid << "\t" << bd->bb.chrm << "\t" << e.strand << "\t"
			         << (out_circ_start + 1) << "\t" << out_circ_end << "\t" << vpath.size() << "\t"
			         << residual_cr.size() << "\t" << residual_or.size() << "\t";

			if (min_edge_weight[idx] >= 0)
				feat_out << min_edge_weight[idx];
			else
				feat_out << "NA";

			feat_out << "\t";

			if (min_vertex_weight[idx] >= 0.0)
				feat_out << min_vertex_weight[idx];
			else
				feat_out << "NA";

			feat_out << "\t";

			if (insert_sizes[idx] >= 0)
				feat_out << insert_sizes[idx] << "\t" << abs(insert_sizes[idx] - (int32_t)length_median);
			else
				feat_out << "NA\tNA";

			feat_out << "\n";
		}

		circ_out << "Total Selected circRNA Paths: " << total_selected << "\n\n\n";
		circ_out.close();
		gtf_out.close();
		diag_out.close();
		feat_out.close();
	}

    bgout << "\n\n";
    bgout.close();
}

int bridger::print_splice_graph()
{
	build_junction_graph(bd->fragments);
	static int bundle_counter = 0;
	int current_bundle_index = bundle_counter++;

	collect_bundle_reads(current_bundle_index);

	// NO USEFUL INFORMATION IF NEITHER CHIMERIC READS NOR OUTWARD READS IN BUNDLE
	if (bundle_chimeric_reads[current_bundle_index].empty() && bundle_outward_reads[current_bundle_index].empty()) return 0;

	static std::string outdir;
	static bool dir_created = false;

	if (!dir_created)
	{
		time_t now = time(nullptr);
		char buf[32];
		strftime(buf, sizeof(buf), "results_%Y%m%d_%H%M%S", localtime(&now));
		outdir = buf;
		mkdir(outdir.c_str(), 0755);

		dir_created = true;
	}

	write_chimeric_read_paths(current_bundle_index, outdir);
	write_simplified_chimeric_paths(current_bundle_index, outdir);
	write_chimeric_insert_sizes(current_bundle_index, outdir);
	build_outward_read_paths(current_bundle_index, outdir);
	// DISABLED (Round 15/21 perf fix): this diagnostic dump
	// (simplify_outward_full_sequence_paths.txt) is unbounded and reached
	// 29.4GB on liver (vs 36MB on kidney), confirmed unused by any evaluate/
	// script. Re-enable only for targeted debugging, never for a full run.
	// write_simplified_outward_paths(current_bundle_index, outdir);

	map<string, vector<string>> chim_rep_members;
	map<string, vector<string>> outward_rep_members;
	map<string, vector<vector<int>>> chim_bg = build_chimeric_bg(current_bundle_index, chim_rep_members);
	map<string, vector<vector<int>>> outward_bg = build_outward_bg(current_bundle_index, outward_rep_members);

	map<string, vector<vector<int>>> bg;
	for (const pair<const string, vector<vector<int>>>& e : chim_bg)
		bg[e.first].insert(bg[e.first].end(), e.second.begin(), e.second.end());

	for (const pair<const string, vector<vector<int>>>& e : outward_bg)
		bg[e.first].insert(bg[e.first].end(), e.second.begin(), e.second.end());

	write_bipartite_graph_file(current_bundle_index, outdir, chim_bg, outward_bg, bg, chim_rep_members, outward_rep_members);

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
	int all_count = 0;
	int only_ref_count = 0;
	int single_ref_count = 0;
	int multi_ref_count = 0;

	//printf("1.5*length_high = %lf\n",1.5*length_high);
	for(int k = 0; k < frags.size(); k++)
	{
		all_count++;

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

		//not use discard path with score < min_pathscore
		/*selected_paths.clear();
		for(int i=0;i<primary_selected_paths.size();i++)
		{
			selected_paths.push_back(primary_selected_paths[i]);
		}*/

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
				only_ref_count++;

				if(ref_paths_map.size() == 1) single_ref_count++;
				else if(ref_paths_map.size() > 1) multi_ref_count++;

				//discard this frag if comes to only ref
				// fr.set_bridged(false);
				// fr.paths.resize(0);
				// continue;

				//discard this frag if comes to only ref and ref size > 1
				if(ref_paths_map.size() > 1)
				{
					fr.set_bridged(false);
					fr.paths.resize(0);
					continue;
				}
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

	bd->total_frag_count = all_count;
	bd->only_ref_path_frag_count = only_ref_count;
	bd->single_ref_chosen_count = single_ref_count;
	bd->multi_ref_chosen_count = multi_ref_count;
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