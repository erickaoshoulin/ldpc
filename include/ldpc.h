#ifndef LDPC_H
#define LDPC_H

#include <stdint.h>

typedef enum {
    LDPC_ALG_CONVENTIONAL = 0,
    LDPC_ALG_RMAS1 = 1,
    LDPC_ALG_RMAS2 = 2,
} ldpc_algorithm_t;

typedef struct {
    int n;
    int m;
    int edges;
    int *cn_start;
    int *vn_of_edge;
    int *vn_degree;
    int *vn_edges_start;
    int *vn_edges;
} ldpc_matrix_t;

typedef struct {
    int max_iters;
    float alpha;
    float beta;
    int group_size;
    float damping;
} ldpc_decoder_params_t;

typedef struct {
    int iterations_used;
    int success;
} ldpc_decode_result_t;

int ldpc_generate_regular(ldpc_matrix_t *h, int n, int m, int dv, uint32_t seed);
void ldpc_free_matrix(ldpc_matrix_t *h);

void ldpc_awgn_bpsk_llr(const uint8_t *bits, int n, float ebn0_db, float rate, uint32_t *seed, float *llr);

ldpc_decode_result_t ldpc_decode(
    const ldpc_matrix_t *h,
    const float *channel_llr,
    ldpc_algorithm_t algorithm,
    const ldpc_decoder_params_t *params,
    uint8_t *hard_bits);

int ldpc_check_syndrome(const ldpc_matrix_t *h, const uint8_t *bits);

#endif
