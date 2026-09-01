/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BRIDGER_H__
#define __BRIDGER_H__

#include "bundle_bridge.h"
#include "fcluster.h"
#include <map>
#include <vector>
#include <utility>
#include <string>
#include <fstream>
#include <queue> // Added for BFS
#include <tuple> // Added for storing support info

using namespace std;

struct read_info {
    string qname;
    int pos;
    int rpos;
    vector<int> vlist;
    bool is_read1;
    bool is_read2;
    bool is_primary;
    bool is_supplementary;
    bool is_reverse;    // true if aligned to the − strand (SAM flag 0x10)
    const hit* src;     // back-pointer to the underlying alignment record
};

struct chimeric_bsj_info {
    int32_t bsj_start_pos;
    int32_t bsj_end_pos;
    vector<int> bsj_vlist;
};

// A candidate path from a Pareto (multi-criteria) search over the splice graph -- see
// bridger::find_pareto_paths in bridger.cc.
struct CandidatePath {
    vector<int> path;
    int min_edge_weight;
    double min_vertex_weight;   // -1.0 sentinel if every vertex on the path had zero coverage
    int32_t length;
};

class entry
{
public:
	vector<int> stack;
	int32_t length;
	int trace1;
	int trace2;

public:
	int print();
};

bool entry_compare(const entry &x, const entry &y);

class bridger
{
public:
	bridger(bundle_bridge *b);

public:
	bundle_bridge *bd;				// parent bundle
	vector<path> pnodes;			// path nodes (not used)
	vector< map<int, int> > jsetx;	// junction graph (out) 
	vector< map<int, int> > jsety;	// junction graph (in)
	map<int, vector<pair<int, int>>> splice_graph_adj; // splice graph adjacency list
	vector< map<int, int> > psetx;	// path graph (out) (not used)
	vector< map<int, int> > psety;	// path graph (in) (not used)
	int max_pnode_length;			// kmer size
	int32_t length_median;
	int32_t length_low;	//DISTRBN OF FRGAMNET length 0
	int32_t length_high; //DISTRBN OF FRGAMNET length 10000

    static map<int, vector<read_info>> bundle_chimeric_reads;
    static map<int, vector<read_info>> bundle_outward_reads;
    static map<int, map<string, chimeric_bsj_info>> bundle_chimeric_bsj;
    static map<int, map<string, vector<string>>> bundle_chimeric_support_names;
    static map<int, map<string, vector<pair<vector<int>, vector<int>>>>> bundle_chimeric_support_paths;
    static map<int, map<string, map<pair<int32_t,int32_t>, vector<vector<int>>>>> bundle_outward_bsj_paths;
	static map<int, map<string, vector<vector<int>>>> bundle_chimeric_merged_paths;


public:
	int bridge_normal_fragments();
	int bridge_circ_fragments();
	int bridge_clip(int32_t p1, int32_t p2, circular_transcript &circ);
	int pick_bridge_path(vector<fragment> &frags);
	int print(vector<fragment> &frags);

public:
	int bridge_overlapped_fragments(vector<fragment> &frags);
	int bridge_overlapped_fragment(fragment &fr, int ex1, int ex2);

	int bridge_phased_fragments(vector<fcluster> &fclusters);
	int phase_cluster(fcluster &fc);
	int bridge_phased_cluster(fcluster &fc);

	int remove_tiny_boundary(vector<fragment> &frags);

	int build_junction_graph(vector<fragment> &frags);
	int print_splice_graph();
	void write_splice_graph(std::ofstream& fout);
	void write_chimeric_read_paths(int bundle_idx, const std::string& outdir);
	void write_simplified_chimeric_paths(int bundle_idx, const std::string& outdir);
	void write_chimeric_insert_sizes(int bundle_idx, const std::string& outdir);
	void build_outward_read_paths(int bundle_idx, const std::string& outdir);
	void write_simplified_outward_paths(int bundle_idx, const std::string& outdir);
	map<string, vector<vector<int>>> build_chimeric_bg(int bundle_idx, map<string, vector<string>>& rep_members);
	map<string, vector<vector<int>>> build_outward_bg(int bundle_idx, map<string, vector<string>>& rep_members);
	void write_bipartite_graph_file(int bundle_idx, const std::string& outdir,
	    const map<string, vector<vector<int>>>& chim_bg,
	    const map<string, vector<vector<int>>>& outward_bg,
	    const map<string, vector<vector<int>>>& bg,
	    const map<string, vector<string>>& chim_rep_members,
	    const map<string, vector<string>>& outward_rep_members);
    int collect_bundle_reads(int bundle_index);
    vector<int> get_bsj_vlist_for_reads(int bsj_start_pos, int bsj_end_pos, const vector<read_info>& alignments_for_qname);
	int bridge_hard_fragments_normal(vector<fcluster> &open);
	int bridge_hard_fragments_circ(vector<fcluster> &open);
	int bridge_hard_fragments_outward(vector<fcluster> &open);
	int dynamic_programming(int k1, int k2, vector< vector<entry> > &table);
	vector< vector<int> > trace_back(int k, const vector< vector<entry> > &table);
	int evaluate_bridging_path(const vector<int> &pb);
	int determine_overlap(const vector<int> &vx, const vector<int> &vy, PI &p);
	int determine_overlap1(const vector<int> &vx, const vector<int> &vy, PI &p);
	bool determine_identical(const vector<int> &vx, const vector<int> &vy, int x1, int x2, int y1, int y2);

	int build_overlap_index();
	int dynamic_programming(int k1, int k2, vector<int> &trace, vector< vector<int> > &table_cov, vector<int32_t> &table_len);
	int compare_stack(const vector<int> &x, const vector<int> &y);
	vector<int> update_stack(const vector<int> &v, int s);

	vector<int> trace_back(int k1, int k2, const vector<int> &trace);
	vector<int> get_bridge(const vector<int> &vv, const vector<int> &v1, const vector<int> &v2);
	int32_t get_extended_length1(int k2, int p1, int p2);
	int32_t get_extended_length2(int k1, int p1, int p2);
	vector<int> get_suffix(const vector<int> &v);
	vector<int> get_prefix(const vector<int> &v);

	int cluster_open_fragments(vector<fcluster> &fclusters, vector<fragment> &frags);
	int build_path_nodes(vector<fragment> &frags);
	int build_path_nodes(int max_len, vector<fragment> &frags);
	int build_path_nodes(int low, int high, vector<fragment> &fragments);
	int build_path_nodes(map<vector<int>, int> &m, const vector<int> &v, int cnt);
	int add_consecutive_path_nodes();
	int adjust_path_score(path &p);

	int set_thresholds();
	int set_normal_length();
	int set_circ_length();
	int filter_paths(vector<fragment> &frags);
	int get_paired_fragments(vector<fragment> &frags);
	vector<int> get_bridged_fragments_type(vector<fragment> &frags);
private:
	void find_all_paths(int start, int end, vector<int>& current_path, vector<vector<int>>& all_paths, int max_paths = 150);
	vector<CandidatePath> find_pareto_paths(int start, int end, int32_t len_cap, bool require_canonical_junctions = false);
	vector<int> widest_path(int start, int end);
};

bool compare_fragment_v1(fragment *f1, fragment *f2);
bool compare_fragment_v2(fragment *f1, fragment *f2);
bool compare_fragment_v3(fragment *f1, fragment *f2);
bool compare_fragment_v3_flank(fragment *f1, fragment *f2);
bool compare_fragment_path(fragment *f1, fragment *f2);
bool check_suffix(const vector<int> &vx, const vector<int> &vy);

#endif