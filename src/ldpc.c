#include "ldpc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int index;
    float residual;
} residual_entry_t;

static int cmp_residual_desc(const void *a, const void *b) {
    const residual_entry_t *ra = (const residual_entry_t *)a;
    const residual_entry_t *rb = (const residual_entry_t *)b;
    if (ra->residual < rb->residual) return 1;
    if (ra->residual > rb->residual) return -1;
    return 0;
}

static uint32_t xorshift32(uint32_t *state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

static float uniform01(uint32_t *seed) {
    return (xorshift32(seed) >> 8) * (1.0f / 16777216.0f);
}

static float gaussian01(uint32_t *seed) {
    const float u1 = fmaxf(uniform01(seed), 1e-7f);
    const float u2 = uniform01(seed);
    return sqrtf(-2.0f * logf(u1)) * cosf(6.28318530718f * u2);
}

int ldpc_generate_regular(ldpc_matrix_t *h, int n, int m, int dv, uint32_t seed) {
    (void)seed;
    if (!h || n <= 0 || m <= 0 || dv <= 1 || (n * dv) % m != 0 || (n % m) != 0) {
        return -1;
    }

    memset(h, 0, sizeof(*h));
    h->n = n;
    h->m = m;
    h->edges = n * dv;
    const int dc = h->edges / m;

    h->cn_start = (int *)calloc((size_t)m + 1, sizeof(int));
    h->vn_of_edge = (int *)calloc((size_t)h->edges, sizeof(int));
    h->vn_degree = (int *)calloc((size_t)n, sizeof(int));
    h->vn_edges_start = (int *)calloc((size_t)n + 1, sizeof(int));
    h->vn_edges = (int *)calloc((size_t)h->edges, sizeof(int));
    if (!h->cn_start || !h->vn_of_edge || !h->vn_degree || !h->vn_edges_start || !h->vn_edges) {
        ldpc_free_matrix(h);
        return -1;
    }

    int *used = (int *)calloc((size_t)m, sizeof(int));
    uint8_t *adj = (uint8_t *)calloc((size_t)n * (size_t)m, sizeof(uint8_t));
    int *cursor = (int *)calloc((size_t)((n > m) ? n : m), sizeof(int));
    if (!used || !adj || !cursor) {
        free(used);
        free(adj);
        free(cursor);
        ldpc_free_matrix(h);
        return -1;
    }

    const int step = 5;
    const int jump = m / dv;
    for (int v = 0; v < n; ++v) {
        for (int k = 0; k < dv; ++k) {
            int c = (v * step + k * jump) % m;
            if (adj[(size_t)v * (size_t)m + (size_t)c]) {
                free(used); free(adj); free(cursor); ldpc_free_matrix(h); return -1;
            }
            adj[(size_t)v * (size_t)m + (size_t)c] = 1;
            used[c]++;
            h->vn_degree[v]++;
        }
    }

    for (int c = 0; c < m; ++c) {
        if (used[c] != dc) {
            free(used); free(adj); free(cursor); ldpc_free_matrix(h); return -1;
        }
    }

    h->cn_start[0] = 0;
    for (int c = 0; c < m; ++c) {
        h->cn_start[c + 1] = h->cn_start[c] + used[c];
        cursor[c] = h->cn_start[c];
    }

    for (int v = 0; v < n; ++v) {
        for (int c = 0; c < m; ++c) {
            if (adj[(size_t)v * (size_t)m + (size_t)c]) {
                int e = cursor[c]++;
                h->vn_of_edge[e] = v;
            }
        }
    }

    h->vn_edges_start[0] = 0;
    for (int v = 0; v < n; ++v) {
        h->vn_edges_start[v + 1] = h->vn_edges_start[v] + h->vn_degree[v];
        cursor[v] = h->vn_edges_start[v];
    }

    for (int e = 0; e < h->edges; ++e) {
        int v = h->vn_of_edge[e];
        h->vn_edges[cursor[v]++] = e;
    }

    free(used);
    free(adj);
    free(cursor);
    return 0;
}

void ldpc_free_matrix(ldpc_matrix_t *h) {
    if (!h) {
        return;
    }
    free(h->cn_start);
    free(h->vn_of_edge);
    free(h->vn_degree);
    free(h->vn_edges_start);
    free(h->vn_edges);
    memset(h, 0, sizeof(*h));
}

void ldpc_awgn_bpsk_llr(const uint8_t *bits, int n, float ebn0_db, float rate, uint32_t *seed, float *llr) {
    const float ebn0 = powf(10.0f, ebn0_db / 10.0f);
    const float sigma2 = 1.0f / (2.0f * rate * ebn0);
    const float sigma = sqrtf(sigma2);
    for (int i = 0; i < n; ++i) {
        float x = bits[i] ? -1.0f : 1.0f;
        float y = x + sigma * gaussian01(seed);
        llr[i] = 2.0f * y / sigma2;
    }
}

int ldpc_check_syndrome(const ldpc_matrix_t *h, const uint8_t *bits) {
    for (int c = 0; c < h->m; ++c) {
        int p = 0;
        for (int e = h->cn_start[c]; e < h->cn_start[c + 1]; ++e) {
            p ^= bits[h->vn_of_edge[e]];
        }
        if (p) {
            return 0;
        }
    }
    return 1;
}

static void check_node_update(const ldpc_matrix_t *h, int c, float alpha, float beta, float *v2c, float *c2v) {
    float min1 = 1e30f;
    float min2 = 1e30f;
    int min_idx = -1;
    int sign = 0;
    for (int e = h->cn_start[c]; e < h->cn_start[c + 1]; ++e) {
        float val = v2c[e];
        float a = fabsf(val);
        if (a < min1) {
            min2 = min1;
            min1 = a;
            min_idx = e;
        } else if (a < min2) {
            min2 = a;
        }
        sign ^= (val < 0.0f);
    }

    for (int e = h->cn_start[c]; e < h->cn_start[c + 1]; ++e) {
        float mag = (e == min_idx) ? min2 : min1;
        mag = fmaxf(alpha * mag - beta, 0.0f);
        int out_sign = sign ^ (v2c[e] < 0.0f);
        c2v[e] = out_sign ? -mag : mag;
    }
}

static float check_residual(const ldpc_matrix_t *h, int c, float alpha, float beta, const float *v2c, const float *c2v) {
    float min1 = 1e30f;
    float min2 = 1e30f;
    int min_idx = -1;
    int sign = 0;
    for (int e = h->cn_start[c]; e < h->cn_start[c + 1]; ++e) {
        const float val = v2c[e];
        const float a = fabsf(val);
        if (a < min1) {
            min2 = min1;
            min1 = a;
            min_idx = e;
        } else if (a < min2) {
            min2 = a;
        }
        sign ^= (val < 0.0f);
    }

    float res = 0.0f;
    for (int e = h->cn_start[c]; e < h->cn_start[c + 1]; ++e) {
        float mag = (e == min_idx) ? min2 : min1;
        mag = fmaxf(alpha * mag - beta, 0.0f);
        const int out_sign = sign ^ (v2c[e] < 0.0f);
        const float next = out_sign ? -mag : mag;
        const float delta = fabsf(next - c2v[e]);
        if (delta > res) {
            res = delta;
        }
    }
    return res;
}

static void refresh_variable_node(
    const ldpc_matrix_t *h,
    int v,
    const float *channel_llr,
    float damp,
    float *app,
    float *v2c,
    const float *c2v) {
    float sum = channel_llr[v];
    for (int i = h->vn_edges_start[v]; i < h->vn_edges_start[v + 1]; ++i) {
        sum += c2v[h->vn_edges[i]];
    }
    app[v] = damp * app[v] + (1.0f - damp) * sum;
    for (int i = h->vn_edges_start[v]; i < h->vn_edges_start[v + 1]; ++i) {
        const int e = h->vn_edges[i];
        const float extrinsic = app[v] - c2v[e];
        v2c[e] = damp * v2c[e] + (1.0f - damp) * extrinsic;
    }
}

static ldpc_decode_result_t decode_conventional(
    const ldpc_matrix_t *h,
    const float *channel_llr,
    const ldpc_decoder_params_t *p,
    uint8_t *hard_bits) {
    float *v2c = (float *)calloc((size_t)h->edges, sizeof(float));
    float *c2v = (float *)calloc((size_t)h->edges, sizeof(float));
    float *app = (float *)calloc((size_t)h->n, sizeof(float));
    if (!v2c || !c2v || !app) {
        free(v2c);
        free(c2v);
        free(app);
        return (ldpc_decode_result_t){0, 0};
    }

    for (int e = 0; e < h->edges; ++e) {
        v2c[e] = channel_llr[h->vn_of_edge[e]];
    }

    ldpc_decode_result_t r = {p->max_iters, 0};
    for (int it = 0; it < p->max_iters; ++it) {
        for (int c = 0; c < h->m; ++c) {
            check_node_update(h, c, p->alpha, p->beta, v2c, c2v);
        }

        for (int v = 0; v < h->n; ++v) {
            float sum = channel_llr[v];
            for (int i = h->vn_edges_start[v]; i < h->vn_edges_start[v + 1]; ++i) {
                sum += c2v[h->vn_edges[i]];
            }
            app[v] = sum;
            hard_bits[v] = (sum < 0.0f) ? 1U : 0U;
        }

        if (ldpc_check_syndrome(h, hard_bits)) {
            r.iterations_used = it + 1;
            r.success = 1;
            break;
        }

        for (int v = 0; v < h->n; ++v) {
            for (int i = h->vn_edges_start[v]; i < h->vn_edges_start[v + 1]; ++i) {
                int e = h->vn_edges[i];
                v2c[e] = app[v] - c2v[e];
            }
        }
    }

    free(v2c);
    free(c2v);
    free(app);
    return r;
}

static ldpc_decode_result_t decode_rmas1(
    const ldpc_matrix_t *h,
    const float *channel_llr,
    const ldpc_decoder_params_t *p,
    uint8_t *hard_bits) {
    float *v2c = (float *)calloc((size_t)h->edges, sizeof(float));
    float *c2v = (float *)calloc((size_t)h->edges, sizeof(float));
    float *app = (float *)calloc((size_t)h->n, sizeof(float));
    if (!v2c || !c2v || !app) {
        free(v2c);
        free(c2v);
        free(app);
        return (ldpc_decode_result_t){0, 0};
    }

    for (int v = 0; v < h->n; ++v) {
        app[v] = channel_llr[v];
        for (int i = h->vn_edges_start[v]; i < h->vn_edges_start[v + 1]; ++i) {
            v2c[h->vn_edges[i]] = app[v];
        }
    }

    const int group = p->group_size < 1 ? 1 : p->group_size;
    const float damp = fminf(fmaxf(p->damping, 0.0f), 1.0f);
    residual_entry_t *groups = (residual_entry_t *)calloc((size_t)((h->m + group - 1) / group), sizeof(residual_entry_t));
    if (!groups) {
        free(v2c);
        free(c2v);
        free(app);
        return (ldpc_decode_result_t){0, 0};
    }

    const int num_groups = (h->m + group - 1) / group;
    ldpc_decode_result_t r = {p->max_iters, 0};

    for (int it = 0; it < p->max_iters; ++it) {
        for (int g = 0; g < num_groups; ++g) {
            int g0 = g * group;
            int g1 = g0 + group;
            if (g1 > h->m) {
                g1 = h->m;
            }
            float g_res = 0.0f;
            for (int c = g0; c < g1; ++c) {
                float c_res = check_residual(h, c, p->alpha, p->beta, v2c, c2v);
                if (c_res > g_res) {
                    g_res = c_res;
                }
            }
            groups[g].index = g;
            groups[g].residual = g_res;
        }

        qsort(groups, (size_t)num_groups, sizeof(groups[0]), cmp_residual_desc);

        for (int gi = 0; gi < num_groups; ++gi) {
            int g0 = groups[gi].index * group;
            int g1 = g0 + group;
            if (g1 > h->m) {
                g1 = h->m;
            }
            for (int c = g0; c < g1; ++c) {
                check_node_update(h, c, p->alpha, p->beta, v2c, c2v);
            }

            for (int v = 0; v < h->n; ++v) {
                refresh_variable_node(h, v, channel_llr, damp, app, v2c, c2v);
            }
        }

        for (int v = 0; v < h->n; ++v) {
            hard_bits[v] = (app[v] < 0.0f) ? 1U : 0U;
        }

        if (ldpc_check_syndrome(h, hard_bits)) {
            r.iterations_used = it + 1;
            r.success = 1;
            break;
        }
    }

    free(v2c);
    free(c2v);
    free(app);
    free(groups);
    return r;
}

static ldpc_decode_result_t decode_rmas2(
    const ldpc_matrix_t *h,
    const float *channel_llr,
    const ldpc_decoder_params_t *p,
    uint8_t *hard_bits) {
    float *v2c = (float *)calloc((size_t)h->edges, sizeof(float));
    float *c2v = (float *)calloc((size_t)h->edges, sizeof(float));
    float *app = (float *)calloc((size_t)h->n, sizeof(float));
    if (!v2c || !c2v || !app) {
        free(v2c);
        free(c2v);
        free(app);
        return (ldpc_decode_result_t){0, 0};
    }

    for (int v = 0; v < h->n; ++v) {
        app[v] = channel_llr[v];
        for (int i = h->vn_edges_start[v]; i < h->vn_edges_start[v + 1]; ++i) {
            v2c[h->vn_edges[i]] = app[v];
        }
    }

    const float damp = fminf(fmaxf(p->damping, 0.0f), 1.0f);
    int *edge_to_cn = (int *)calloc((size_t)h->edges, sizeof(int));
    residual_entry_t *checks = (residual_entry_t *)calloc((size_t)h->m, sizeof(residual_entry_t));
    if (!edge_to_cn || !checks) {
        free(v2c);
        free(c2v);
        free(app);
        free(edge_to_cn);
        free(checks);
        return (ldpc_decode_result_t){0, 0};
    }
    for (int c = 0; c < h->m; ++c) {
        checks[c].index = c;
        for (int e = h->cn_start[c]; e < h->cn_start[c + 1]; ++e) {
            edge_to_cn[e] = c;
        }
    }

    ldpc_decode_result_t r = {p->max_iters, 0};

    for (int it = 0; it < p->max_iters; ++it) {
        for (int c = 0; c < h->m; ++c) {
            checks[c].residual = check_residual(h, c, p->alpha, p->beta, v2c, c2v);
        }
        qsort(checks, (size_t)h->m, sizeof(checks[0]), cmp_residual_desc);

        for (int i = 0; i < h->m; ++i) {
            const int c = checks[i].index;
            check_node_update(h, c, p->alpha, p->beta, v2c, c2v);

            for (int e = h->cn_start[c]; e < h->cn_start[c + 1]; ++e) {
                int v = h->vn_of_edge[e];
                refresh_variable_node(h, v, channel_llr, damp, app, v2c, c2v);
                for (int j = h->vn_edges_start[v]; j < h->vn_edges_start[v + 1]; ++j) {
                    int neigh = edge_to_cn[h->vn_edges[j]];
                    checks[neigh].residual = check_residual(h, neigh, p->alpha, p->beta, v2c, c2v);
                }
            }
        }

        for (int v = 0; v < h->n; ++v) {
            hard_bits[v] = (app[v] < 0.0f) ? 1U : 0U;
        }

        if (ldpc_check_syndrome(h, hard_bits)) {
            r.iterations_used = it + 1;
            r.success = 1;
            break;
        }
    }

    free(v2c);
    free(c2v);
    free(app);
    free(edge_to_cn);
    free(checks);
    return r;
}


static ldpc_decode_result_t decode_as(
    const ldpc_matrix_t *h,
    const float *channel_llr,
    const ldpc_decoder_params_t *p,
    uint8_t *hard_bits) {
    float *v2c = (float *)calloc((size_t)h->edges, sizeof(float));
    float *c2v = (float *)calloc((size_t)h->edges, sizeof(float));
    float *app = (float *)calloc((size_t)h->n, sizeof(float));
    if (!v2c || !c2v || !app) {
        free(v2c);
        free(c2v);
        free(app);
        return (ldpc_decode_result_t){0, 0};
    }

    for (int e = 0; e < h->edges; ++e) {
        v2c[e] = channel_llr[h->vn_of_edge[e]];
    }

    const float alpha0 = p->alpha;
    const float alpha1 = p->alpha_final;
    ldpc_decode_result_t r = {p->max_iters, 0};
    for (int it = 0; it < p->max_iters; ++it) {
        const float t = (p->max_iters <= 1) ? 1.0f : (float)it / (float)(p->max_iters - 1);
        const float alpha_it = alpha0 + (alpha1 - alpha0) * t;

        for (int c = 0; c < h->m; ++c) {
            check_node_update(h, c, alpha_it, p->beta, v2c, c2v);
        }

        for (int v = 0; v < h->n; ++v) {
            float sum = channel_llr[v];
            for (int i = h->vn_edges_start[v]; i < h->vn_edges_start[v + 1]; ++i) {
                sum += c2v[h->vn_edges[i]];
            }
            app[v] = sum;
            hard_bits[v] = (sum < 0.0f) ? 1U : 0U;
        }

        if (ldpc_check_syndrome(h, hard_bits)) {
            r.iterations_used = it + 1;
            r.success = 1;
            break;
        }

        for (int v = 0; v < h->n; ++v) {
            for (int i = h->vn_edges_start[v]; i < h->vn_edges_start[v + 1]; ++i) {
                int e = h->vn_edges[i];
                v2c[e] = app[v] - c2v[e];
            }
        }
    }

    free(v2c);
    free(c2v);
    free(app);
    return r;
}



ldpc_decode_result_t ldpc_decode(
    const ldpc_matrix_t *h,
    const float *channel_llr,
    ldpc_algorithm_t algorithm,
    const ldpc_decoder_params_t *params,
    uint8_t *hard_bits) {
    if (algorithm == LDPC_ALG_RMAS1) {
        return decode_rmas1(h, channel_llr, params, hard_bits);
    }
    if (algorithm == LDPC_ALG_RMAS2) {
        return decode_rmas2(h, channel_llr, params, hard_bits);
    }
    if (algorithm == LDPC_ALG_AS) {
        return decode_as(h, channel_llr, params, hard_bits);
    }
    return decode_conventional(h, channel_llr, params, hard_bits);
}
