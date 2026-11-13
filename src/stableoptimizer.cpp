// This file is part of meshoptimizer library; see meshoptimizer.h for version/license details
#include "meshoptimizer.h"

#include <assert.h>
#include <float.h>
#include <string.h>

// This work is based on:
// Fabian Giesen. Decoding Morton codes. 2009
namespace meshopt
{

// "Insert" two 0 bits after each of the 20 low bits of x
inline unsigned long long part1By2(unsigned long long x)
{
	x &= 0x000fffffull;                          // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- jihg fedc ba98 7654 3210
	x = (x ^ (x << 32)) & 0x000f00000000ffffull; // x = ---- ---- ---- jihg ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
	x = (x ^ (x << 16)) & 0x000f0000ff0000ffull; // x = ---- ---- ---- jihg ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- 7654 3210
	x = (x ^ (x << 8)) & 0x000f00f00f00f00full;  // x = ---- ---- ---- jihg ---- ---- fedc ---- ---- ba98 ---- ---- 7654 ---- ---- 3210
	x = (x ^ (x << 4)) & 0x00c30c30c30c30c3ull;  // x = ---- ---- ji-- --hg ---- fe-- --dc --b- --98 --7- -6-- 5--4 --3- -2-- 1--0
	x = (x ^ (x << 2)) & 0x0249249249249249ull;  // x = ---- --j- -i-- h--g ---- fe-- --dc --b- --a- -9-- 8--7 -6-- 5--4 --3- -2-- 1--0
	return x;
}

struct StableVertexLess
{
	const meshopt_Stream* streams;
	size_t stream_count;
	const unsigned int* representative; // source vertex for each unique vertex

	bool operator()(unsigned int lhs, unsigned int rhs) const
	{
		for (size_t i = 0; i < stream_count; ++i)
		{
			const meshopt_Stream& s = streams[i];
			const unsigned char* data = static_cast<const unsigned char*>(s.data);

			int cmp = memcmp(data + representative[lhs] * s.stride, data + representative[rhs] * s.stride, s.size);

			if (cmp != 0)
				return cmp < 0;
		}

		return false;
	}
};

struct StableTriangleLess
{
	const unsigned long long* keys;
	const unsigned int* indices; // canonical vertex ids with normalized rotation

	bool operator()(unsigned int lhs, unsigned int rhs) const
	{
		if (keys[lhs] != keys[rhs])
			return keys[lhs] < keys[rhs];

		const unsigned int* a = indices + lhs * 3;
		const unsigned int* b = indices + rhs * 3;

		if (a[0] != b[0])
			return a[0] < b[0];
		if (a[1] != b[1])
			return a[1] < b[1];

		return a[2] < b[2];
	}
};

// compare a triangle rotation (x0, x1, x2) with the current best rotation t lexicographically
static bool rotationLess(unsigned int x0, unsigned int x1, unsigned int x2, const unsigned int* t)
{
	return x0 < t[0] || (x0 == t[0] && (x1 < t[1] || (x1 == t[1] && x2 < t[2])));
}

// stable bottom-up merge sort; the result only depends on the set of elements, not their initial order
template <typename T, typename Less>
static void stableSort(T* data, T* scratch, size_t count, const Less& less)
{
	if (count < 2)
		return;

	T* src = data;
	T* dst = scratch;

	for (size_t width = 1; width < count; width *= 2)
	{
		for (size_t i = 0; i < count; i += 2 * width)
		{
			size_t mid = i + width < count ? i + width : count;
			size_t end = i + 2 * width < count ? i + 2 * width : count;

			size_t l = i, r = mid;

			for (size_t o = i; o < end; ++o)
				dst[o] = (l < mid && (r >= end || !less(src[r], src[l]))) ? src[l++] : src[r++];
		}

		T* swap = src;
		src = dst;
		dst = swap;
	}

	if (src != data)
		memcpy(data, src, count * sizeof(T));
}

} // namespace meshopt

size_t meshopt_optimizeMeshStable(unsigned int* destination, unsigned int* vertex_remap, const unsigned int* indices, size_t index_count, const float* vertex_positions, size_t vertex_count, size_t vertex_positions_stride, const meshopt_Stream* streams, size_t stream_count, float max_cache_cost_ratio)
{
	using namespace meshopt;

	assert(index_count % 3 == 0);
	assert(vertex_positions_stride >= 12 && vertex_positions_stride % sizeof(float) == 0);
	assert(stream_count > 0 && stream_count <= 16);
	assert(max_cache_cost_ratio >= 0.f);

	// empty meshes consistently produce an empty result
	if (index_count == 0 || vertex_count == 0)
		return 0;

	meshopt_Allocator allocator;

	// weld vertices that are binary equivalent across all streams; vertices that are not referenced by the index buffer are excluded
	unsigned int* remap = allocator.allocate<unsigned int>(vertex_count);
	size_t unique_count = meshopt_generateVertexRemapMulti(remap, indices, index_count, vertex_count, streams, stream_count);
	assert(unique_count > 0);

	// representative source vertex for each unique vertex; all vertices in a weld group are binary equivalent
	unsigned int* representative = allocator.allocate<unsigned int>(unique_count);
	memset(representative, -1, unique_count * sizeof(unsigned int));

	for (size_t i = 0; i < vertex_count; ++i)
		if (remap[i] != ~0u && representative[remap[i]] == ~0u)
			representative[remap[i]] = unsigned(i);

	// canonical vertex order: sort unique vertices by stream contents, which makes the order independent of input vertex numbering
	unsigned int* vertex_order = allocator.allocate<unsigned int>(unique_count);
	unsigned int* vertex_scratch = allocator.allocate<unsigned int>(unique_count);

	for (size_t i = 0; i < unique_count; ++i)
		vertex_order[i] = unsigned(i);

	meshopt::StableVertexLess vertex_less = {streams, stream_count, representative};
	meshopt::stableSort(vertex_order, vertex_scratch, unique_count, vertex_less);

	// canonical vertex id for each unique vertex
	unsigned int* vertex_rank = allocator.allocate<unsigned int>(unique_count);

	for (size_t i = 0; i < unique_count; ++i)
		vertex_rank[vertex_order[i]] = unsigned(i);

	// canonical index buffer; each triangle is rotated to the lexicographically smallest rotation, which preserves winding
	unsigned int* canonical = allocator.allocate<unsigned int>(index_count);

	for (size_t i = 0; i < index_count; i += 3)
	{
		unsigned int a = vertex_rank[remap[indices ? indices[i + 0] : unsigned(i + 0)]];
		unsigned int b = vertex_rank[remap[indices ? indices[i + 1] : unsigned(i + 1)]];
		unsigned int c = vertex_rank[remap[indices ? indices[i + 2] : unsigned(i + 2)]];

		unsigned int t[3] = {a, b, c};

		// rotate to the lexicographically smallest rotation so that the result does not depend on the input vertex order
		if (meshopt::rotationLess(b, c, a, t))
		{
			t[0] = b, t[1] = c, t[2] = a;
		}

		if (meshopt::rotationLess(c, a, b, t))
		{
			t[0] = c, t[1] = a, t[2] = b;
		}

		canonical[i + 0] = t[0];
		canonical[i + 1] = t[1];
		canonical[i + 2] = t[2];
	}

	// mesh bounds over referenced vertices; used to quantize triangle centroids into Morton keys
	size_t vertex_stride_float = vertex_positions_stride / sizeof(float);

	float minv[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
	float maxv[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

	for (size_t i = 0; i < unique_count; ++i)
	{
		const float* v = vertex_positions + representative[i] * vertex_stride_float;

		for (int j = 0; j < 3; ++j)
		{
			minv[j] = minv[j] > v[j] ? v[j] : minv[j];
			maxv[j] = maxv[j] < v[j] ? v[j] : maxv[j];
		}
	}

	float extent = 0.f;

	extent = (maxv[0] - minv[0]) < extent ? extent : (maxv[0] - minv[0]);
	extent = (maxv[1] - minv[1]) < extent ? extent : (maxv[1] - minv[1]);
	extent = (maxv[2] - minv[2]) < extent ? extent : (maxv[2] - minv[2]);

	// rescale each axis to 16 bits to get 48-bit Morton keys
	float scale = extent == 0 ? 0.f : 65535.f / extent;

	size_t face_count = index_count / 3;

	unsigned long long* triangle_keys = allocator.allocate<unsigned long long>(face_count);

	for (size_t i = 0; i < face_count; ++i)
	{
		const float* p0 = vertex_positions + representative[vertex_order[canonical[i * 3 + 0]]] * vertex_stride_float;
		const float* p1 = vertex_positions + representative[vertex_order[canonical[i * 3 + 1]]] * vertex_stride_float;
		const float* p2 = vertex_positions + representative[vertex_order[canonical[i * 3 + 2]]] * vertex_stride_float;

		float centroid[3] = {
		    (p0[0] + p1[0] + p2[0]) / 3.f,
		    (p0[1] + p1[1] + p2[1]) / 3.f,
		    (p0[2] + p1[2] + p2[2]) / 3.f,
		};

		int x = int((centroid[0] - minv[0]) * scale + 0.5f);
		int y = int((centroid[1] - minv[1]) * scale + 0.5f);
		int z = int((centroid[2] - minv[2]) * scale + 0.5f);

		triangle_keys[i] = meshopt::part1By2(x) | (meshopt::part1By2(y) << 1) | (meshopt::part1By2(z) << 2);
	}

	// canonical triangle order: sort by spatial key, breaking ties between identical keys with the canonical triangle vertices
	unsigned int* triangle_order = allocator.allocate<unsigned int>(face_count);
	unsigned int* triangle_scratch = allocator.allocate<unsigned int>(face_count);

	for (size_t i = 0; i < face_count; ++i)
		triangle_order[i] = unsigned(i);

	meshopt::StableTriangleLess triangle_less = {triangle_keys, canonical};
	meshopt::stableSort(triangle_order, triangle_scratch, face_count, triangle_less);

	// renumber vertices in first-use order and build the resulting index buffer
	unsigned int* vertex_renumber = allocator.allocate<unsigned int>(unique_count);
	memset(vertex_renumber, -1, unique_count * sizeof(unsigned int));

	unsigned int* result = allocator.allocate<unsigned int>(index_count);
	unsigned int next_vertex = 0;

	for (size_t i = 0; i < index_count; ++i)
	{
		unsigned int v = canonical[triangle_order[i / 3] * 3 + (i % 3)];

		if (vertex_renumber[v] == ~0u)
			vertex_renumber[v] = next_vertex++;

		result[i] = vertex_renumber[v];
	}

	assert(next_vertex == unique_count);

	// map output vertices to source vertices
	unsigned int* output_remap = allocator.allocate<unsigned int>(unique_count);

	for (size_t i = 0; i < unique_count; ++i)
		output_remap[vertex_renumber[i]] = representative[vertex_order[i]];

	// compare the cache cost of the result against the existing vertex cache optimizer for the same mesh
	unsigned int* reference = allocator.allocate<unsigned int>(index_count);
	meshopt_optimizeVertexCache(reference, result, index_count, unique_count);

	meshopt_VertexCacheStatistics reference_stats = meshopt_analyzeVertexCache(reference, index_count, unique_count, 16, 0, 0);
	meshopt_VertexCacheStatistics result_stats = meshopt_analyzeVertexCache(result, index_count, unique_count, 16, 0, 0);

	// the budget can not be satisfied; leave destination/vertex_remap unmodified
	if (reference_stats.acmr > 0 && result_stats.acmr > reference_stats.acmr * max_cache_cost_ratio)
		return 0;

	memcpy(destination, result, index_count * sizeof(unsigned int));
	memcpy(vertex_remap, output_remap, unique_count * sizeof(unsigned int));

	return unique_count;
}
