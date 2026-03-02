#include "ldpc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int write_performance_plot(void) {
    int rc = system("python3 scripts/plot_results.py >/dev/null 2>&1");
    if (rc != 0) {
        fprintf(stderr, "warning: failed to generate performance figure (requires python3)\n");
        return -1;
    }
    return 0;
}


static void run_eval(const ldpc_matrix_t *h, ldpc_algorithm_t alg, const ldpc_decoder_params_t *p,
                     float snr_db, int frames, uint32_t seed,
                     double *ber, double *fer, double *avg_iters) {
    uint8_t *tx = (uint8_t *)calloc((size_t)h->n, sizeof(uint8_t));
    uint8_t *rx = (uint8_t *)calloc((size_t)h->n, sizeof(uint8_t));
    float *llr = (float *)calloc((size_t)h->n, sizeof(float));
    if (!tx || !rx || !llr) {
        fprintf(stderr, "allocation failed\n");
        exit(1);
    }

    long long bit_err = 0;
    long long frame_err = 0;
    long long iter_sum = 0;

    const float rate = (float)(h->n - h->m) / (float)h->n;
    for (int f = 0; f < frames; ++f) {
        memset(tx, 0, (size_t)h->n);
        ldpc_awgn_bpsk_llr(tx, h->n, snr_db, rate, &seed, llr);
        ldpc_decode_result_t r = ldpc_decode(h, llr, alg, p, rx);

        int fe = 0;
        for (int i = 0; i < h->n; ++i) {
            bit_err += rx[i] != tx[i];
            fe |= rx[i] != tx[i];
        }
        frame_err += fe;
        iter_sum += r.iterations_used;
    }

    *ber = (double)bit_err / (double)(frames * h->n);
    *fer = (double)frame_err / (double)frames;
    *avg_iters = (double)iter_sum / (double)frames;

    free(tx);
    free(rx);
    free(llr);
}

int main(int argc, char **argv) {
    int n = 512;
    int m = 64;
    int dv = 6;
    int max_iters = 20;
    int frames = 300;
    int group_size = 8;
    float alpha = 0.80f;
    float beta = 0.15f;
    float damping = 0.20f;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) m = atoi(argv[2]);
    if (argc > 3) dv = atoi(argv[3]);
    if (argc > 4) frames = atoi(argv[4]);

    ldpc_matrix_t h;
    if (ldpc_generate_regular(&h, n, m, dv, 0x12345678u) != 0) {
        fprintf(stderr, "failed to generate regular H matrix (check n*dv %% m == 0)\n");
        return 1;
    }

    ldpc_decoder_params_t conv = {
        .max_iters = max_iters,
        .alpha = alpha,
        .beta = beta,
        .group_size = 1,
        .damping = 0.0f,
    };

    ldpc_decoder_params_t rmas1 = {
        .max_iters = max_iters,
        .alpha = alpha,
        .beta = beta,
        .group_size = group_size,
        .damping = damping,
    };

    FILE *fp = fopen("results/performance.csv", "w");
    if (!fp) {
        perror("fopen");
        ldpc_free_matrix(&h);
        return 1;
    }
    fprintf(fp, "snr_db,algorithm,ber,fer,avg_iterations\n");

    printf("LDPC eval n=%d m=%d dv=%d rate=%.4f frames=%d\n", h.n, h.m, dv,
           (double)(h.n - h.m) / h.n, frames);

    ldpc_decoder_params_t rmas2 = {
        .max_iters = max_iters,
        .alpha = alpha,
        .beta = beta * 0.8f,
        .group_size = group_size,
        .damping = damping * 0.5f,
    };

    for (float snr = 4.5f; snr <= 7.0f; snr += 0.5f) {
        double ber, fer, it;
        run_eval(&h, LDPC_ALG_CONVENTIONAL, &conv, snr, frames, 0xABCDEF01u, &ber, &fer, &it);
        printf("SNR %.2f conventional BER=%.6e FER=%.6e Iter=%.2f\n", snr, ber, fer, it);
        fprintf(fp, "%.2f,conventional,%.8e,%.8e,%.4f\n", snr, ber, fer, it);

        run_eval(&h, LDPC_ALG_RMAS1, &rmas1, snr, frames, 0x10203040u, &ber, &fer, &it);
        printf("SNR %.2f RMAS1        BER=%.6e FER=%.6e Iter=%.2f\n", snr, ber, fer, it);
        fprintf(fp, "%.2f,rmas1,%.8e,%.8e,%.4f\n", snr, ber, fer, it);

        run_eval(&h, LDPC_ALG_RMAS2, &rmas2, snr, frames, 0x55667788u, &ber, &fer, &it);
        printf("SNR %.2f RMAS2        BER=%.6e FER=%.6e Iter=%.2f\n", snr, ber, fer, it);
        fprintf(fp, "%.2f,rmas2,%.8e,%.8e,%.4f\n", snr, ber, fer, it);
    }

    fclose(fp);
    write_performance_plot();
    ldpc_free_matrix(&h);

    printf("results written to results/performance.csv and results/performance.svg\n");
    return 0;
}
