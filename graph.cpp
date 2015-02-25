// graph.cpp

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdint>
#include <bitset>
#include <algorithm>
#include <chrono>
#include <random>
#include <intrin.h>

typedef std::uint64_t index_t;

class graph
{
public:
	graph(index_t n) : adj_(n) { }

	graph& operator=(const graph&) = delete;

	~graph() { }

	void add_edge(index_t u, index_t v)
	{
		assert(u >= 0 &&
			v >= 0 &&
			u != v &&
			u < adj_.size() &&
			v < adj_.size());

		adj_[u] |= (1 << v);
		adj_[v] |= (1 << u);
	}

	index_t get_degree(index_t u) const
	{
		return __popcnt64(adj_[u]);
	}

	index_t num_vertices() const
	{
		return adj_.size();
	}

	std::vector<index_t> adj_;
};

// TODO: not suitable for generation for large graphs.
graph random_graph(index_t n, double p)
{
	assert(n >= 2 && n <= 64);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::bernoulli_distribution d(p);

	graph g(n);
	
	for (index_t i = 0; i < n; ++i)
	{
		for (index_t j = i + 1; j < n; ++j)
		{
			if (d(gen))
			{
				g.add_edge(i, j);
			}
		}
	}

	return g;
}

void complement(graph& g)
{
	for (index_t i = 0; i < g.adj_.size(); ++i)
	{
		g.adj_[i] = ~g.adj_[i];
		g.adj_[i] ^= (1 << i);
	}
}

bool is_adjacent(const graph& g, int u, int v)
{
	assert(u >= 0 && v >= 0 && u < g.adj_.size() && v < g.adj_.size());
	return g.adj_[u] & (1 << v);
}

template <typename OutputIterator>
void bron_kerbosch_recursive(const graph& g, index_t r, index_t p, index_t x, OutputIterator out)
{
	if (p == 0 && x == 0)
	{
		*out = r;
	}

	for (unsigned long v; p != 0; p &= ~(1 << v))
	{
		_BitScanForward64(&v, p);

		bron_kerbosch_recursive(g, r | (1 << v), p & g.adj_[v], x & g.adj_[v], out);
		p &= ~(1 << v);
		x |= (1 << v);
	}
}

template <typename OutputIterator>
void bron_kerbosch(const graph& g, OutputIterator out)
{
	index_t p = (std::numeric_limits<index_t>::max() >>
		(std::numeric_limits<index_t>::digits - g.num_vertices()));

	bron_kerbosch_recursive(g, 0, p, 0, out);
}

index_t independent_domination_number(const graph& g)
{
	graph h(g);
	complement(h);

	std::vector<index_t> doms;
	bron_kerbosch(h, std::back_inserter(doms));

	return *std::min_element(doms.cbegin(), doms.cend(), [](index_t x, index_t y) { return __popcnt64(x) < __popcnt64(y); });
}
